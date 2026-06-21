#!/usr/bin/env python3
"""HYP-2701 / T936: true-wide survival middle-mass gate for LRC14.

THM-556 gives U4 = p0 + p5 + 5*p6 for the six inner sectors.  Therefore
HYP-2695's true-wide floor gate

    U4(E) <= floor_k = (k-6)/7

is equivalent to the survival-currency inequality

    p1+p2+p3+p4 - 4*p6 >= (13-k)/7.

The exact cap gate has the same left side and replaces the right side by
1-cap_k.  This scout keeps all arithmetic as Fractions, extends the S60 exact
boxes, and records which true-wide row families make the currency tight.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import comb
import sys

from lrc14_wide_branch_ridge_codex_s47 import (
    CAP,
    Row,
    additive_energy,
    longest_run,
    missed_distribution,
    primitive,
    squarefree_profile,
    sumset_excess,
)

sys.stdout.reconfigure(line_buffering=True)


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.9f})"


def floor_cap(k: int) -> Fraction:
    return Fraction(k - 6, 7)


def bonferroni4(dist: tuple[Fraction, ...]) -> Fraction:
    total = Fraction(1)
    sign = -1
    for r in range(1, 5):
        sr = sum(Fraction(comb(t, r)) * dist[t] for t in range(r, 7))
        total += sign * sr
        sign *= -1
    return total


def survival(dist: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(sum(dist[t] for t in range(b, 7)) for b in range(7))


def far_count(row: Row) -> int:
    return sum(1 for x in row if x > 14)


def bounded_core(row: Row) -> Row:
    return tuple(x for x in row if x <= 14)


def gap_signature(row: Row) -> tuple[int, ...]:
    return tuple(b - a for a, b in zip(row, row[1:]))


@dataclass(frozen=True)
class Audit:
    label: str
    row: Row
    dist: tuple[Fraction, ...]
    u4: Fraction
    p0: Fraction
    p6: Fraction
    middle_mass: Fraction
    currency: Fraction
    floor_required: Fraction
    cap_required: Fraction
    floor_slack: Fraction
    cap_slack: Fraction
    tail45: Fraction
    surv: tuple[Fraction, ...]

    @property
    def k(self) -> int:
        return len(self.row)

    @property
    def floor(self) -> Fraction:
        return floor_cap(self.k)

    @property
    def cap(self) -> Fraction:
        return CAP[self.k]

    @property
    def floor_debt(self) -> Fraction:
        return max(Fraction(0), -self.floor_slack)


def make_audit(label: str, row: Row) -> Audit:
    dist = missed_distribution(row)
    u4 = bonferroni4(dist)
    p0 = dist[0]
    p6 = dist[6]
    middle_mass = sum(dist[1:5])
    currency = middle_mass - 4 * p6
    floor_required = Fraction(13 - len(row), 7)
    cap_required = 1 - CAP[len(row)]
    floor_slack = currency - floor_required
    cap_slack = currency - cap_required
    tail45 = dist[5] + 5 * p6
    audit = Audit(
        label=label,
        row=row,
        dist=dist,
        u4=u4,
        p0=p0,
        p6=p6,
        middle_mass=middle_mass,
        currency=currency,
        floor_required=floor_required,
        cap_required=cap_required,
        floor_slack=floor_slack,
        cap_slack=cap_slack,
        tail45=tail45,
        surv=survival(dist),
    )
    assert u4 == p0 + tail45
    assert audit.floor - u4 == floor_slack
    assert audit.cap - u4 == cap_slack
    return audit


def structural_line(row: Row) -> str:
    gaps = gap_signature(row)
    return (
        f"second={row[-2]} span={row[-1]} far_count={far_count(row)} "
        f"core={bounded_core(row)} gaps={gaps} maxgap={max(gaps) if gaps else 0} "
        f"run={longest_run(row)} exc={sumset_excess(row)} "
        f"energy={additive_energy(row)} sqfree={squarefree_profile(row)}"
    )


def print_audit(prefix: str, audit: Audit) -> None:
    g = audit.surv
    print(f"{prefix} E={audit.row}")
    print(
        f"      U4={fmt(audit.u4)} floor_slack={fmt(audit.floor_slack)} "
        f"cap_slack={fmt(audit.cap_slack)}"
    )
    print(
        f"      currency=M-4p6={fmt(audit.currency)} required_floor={fmt(audit.floor_required)} "
        f"required_cap={fmt(audit.cap_required)}"
    )
    print(
        f"      p0={fmt(audit.p0)} middle={fmt(audit.middle_mass)} "
        f"p5={fmt(audit.dist[5])} p6={fmt(audit.p6)} tail45={fmt(audit.tail45)}"
    )
    print(
        f"      G1={fmt(g[1])} G5={fmt(g[5])} G6={fmt(g[6])} "
        f"dist={[str(x) for x in audit.dist]}"
    )
    print(f"      {structural_line(audit.row)}")


def push_min(bucket: list[Audit], audit: Audit, keep: int, key) -> None:
    bucket.append(audit)
    bucket.sort(key=key)
    del bucket[keep:]


def scan_box(k: int, bound: int, keep: int = 8) -> dict[str, object]:
    counts: Counter[str] = Counter()
    tight_floor: list[Audit] = []
    tight_cap: list[Audit] = []
    max_p6: list[Audit] = []
    floor_fail: list[Audit] = []
    cap_fail: list[Audit] = []
    by_far: dict[int, list[Audit]] = defaultdict(list)
    by_core_size: dict[int, list[Audit]] = defaultdict(list)
    far_counts: Counter[int] = Counter()

    for combo in combinations(range(1, bound + 1), k - 1):
        row = (0,) + combo
        counts["raw"] += 1
        if row[-1] <= 14 or not primitive(row):
            continue
        counts["primitive_span_gt14"] += 1
        if row[-2] <= 14:
            counts["boundary"] += 1
            continue
        counts["truewide"] += 1
        audit = make_audit(f"k{k}-B{bound}", row)
        push_min(tight_floor, audit, keep, lambda a: (a.floor_slack, a.cap_slack, a.row))
        push_min(tight_cap, audit, keep, lambda a: (a.cap_slack, a.floor_slack, a.row))
        push_min(max_p6, audit, keep, lambda a: (-a.p6, a.floor_slack, a.row))
        fc = far_count(row)
        cs = len(bounded_core(row))
        far_counts[fc] += 1
        push_min(by_far[fc], audit, keep, lambda a: (a.floor_slack, a.cap_slack, a.row))
        push_min(by_core_size[cs], audit, keep, lambda a: (a.floor_slack, a.cap_slack, a.row))
        if audit.floor_slack < 0:
            floor_fail.append(audit)
        if audit.cap_slack < 0:
            cap_fail.append(audit)

    return {
        "counts": counts,
        "tight_floor": tight_floor,
        "tight_cap": tight_cap,
        "max_p6": max_p6,
        "floor_fail": floor_fail,
        "cap_fail": cap_fail,
        "by_far": by_far,
        "by_core_size": by_core_size,
        "far_counts": far_counts,
    }


def named_rows() -> list[tuple[str, Row]]:
    return [
        ("k8_floor_exception", (0, 3, 6, 9, 12, 14, 15, 18)),
        ("k9_floor_leader", (0, 4, 6, 8, 10, 12, 14, 15, 16)),
        ("k10_floor_leader", (0, 2, 4, 6, 8, 10, 12, 14, 15, 16)),
        ("k11_s60_leader", (0, 2, 4, 6, 8, 9, 10, 12, 14, 15, 16)),
        ("k12_s60_leader", (0, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16)),
        ("boundary_ap_like", (0, 2, 4, 6, 8, 10, 12, 14, 15)),
        ("three_cluster", (0, 1, 2, 30, 31, 32, 60, 61, 62)),
        ("dyadic_block", (0, 1, 2, 4, 8, 12, 16, 20, 24)),
        ("split_blocks", (0, 1, 2, 3, 50, 51, 52, 300, 301)),
    ]


def tournament_analysis(audits: list[Audit]) -> None:
    """Tournament on row-family obligations rather than runners/arcs."""

    dedup: dict[Row, Audit] = {}
    for audit in audits:
        if audit.row[-2] <= 14:
            continue
        if audit.row not in dedup or audit.floor_slack < dedup[audit.row].floor_slack:
            dedup[audit.row] = audit
    verts = sorted(dedup.values(), key=lambda a: (a.floor_slack, a.cap_slack, a.row))[:14]
    wins: set[tuple[int, int]] = set()
    scores: Counter[int] = Counter()
    for i, a in enumerate(verts):
        for j, b in enumerate(verts):
            if i >= j:
                continue
            winner, loser = (i, j) if (a.floor_slack, a.cap_slack) <= (b.floor_slack, b.cap_slack) else (j, i)
            wins.add((winner, loser))
            scores[winner] += 1
            scores.setdefault(loser, scores[loser])
    cycles = 0
    for i, j, h in combinations(range(len(verts)), 3):
        if ((i, j) in wins and (j, h) in wins and (h, i) in wins) or (
            (i, h) in wins and (h, j) in wins and (j, i) in wins
        ):
            cycles += 1
    path = sorted(range(len(verts)), key=lambda idx: (-scores[idx], verts[idx].floor_slack, verts[idx].row))
    print("Tournament Analysis: true-wide survival-currency obligations")
    print("  vertices: exact true-wide row/family obligations, not runners or arcs")
    print("  pairwise observable: smaller floor slack wins")
    print("  switch/gauge: orient toward the row spending more cap-floor currency")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian risk path:")
    for rank, idx in enumerate(path, 1):
        audit = verts[idx]
        print(
            f"    {rank:2d}. k={audit.k} E={audit.row} "
            f"far={far_count(audit.row)} floor_slack={audit.floor_slack} "
            f"cap_slack={audit.cap_slack} currency={audit.currency}"
        )


def main() -> None:
    print("HYP-2701 / T936 -- true-wide survival middle-mass gate")
    print("Exact arithmetic: Fractions from sector-wall integration.")
    print("Identity checked per row: floor-U4 = (p1+p2+p3+p4-4p6) - (13-k)/7.")
    print()
    print("CAP/FLOOR CURRENCY")
    for k in sorted(CAP):
        print(
            f"  k={k}: floor_required={(13-k)}/7={fmt(Fraction(13-k, 7))} "
            f"cap_required=1-cap={fmt(1-CAP[k])} dividend={fmt(CAP[k]-floor_cap(k))}"
        )

    print()
    print("NAMED ROWS")
    all_tight: list[Audit] = []
    for label, row in named_rows():
        if row[-2] <= 14:
            tag = "boundary"
        else:
            tag = "truewide"
        audit = make_audit(f"{label}-{tag}", row)
        all_tight.append(audit)
        print_audit(f"  {label:<20} [{tag}]", audit)

    print()
    print("TRUE-WIDE BOX SCANS")
    scan_plan = ((8, 18), (9, 18), (10, 17), (11, 17), (12, 17))
    all_floor_fail: list[Audit] = []
    all_cap_fail: list[Audit] = []
    for k, bound in scan_plan:
        print("=" * 100)
        print(f"k={k}, bound={bound}")
        result = scan_box(k, bound)
        counts: Counter[str] = result["counts"]  # type: ignore[assignment]
        floor_fail: list[Audit] = result["floor_fail"]  # type: ignore[assignment]
        cap_fail: list[Audit] = result["cap_fail"]  # type: ignore[assignment]
        by_far: dict[int, list[Audit]] = result["by_far"]  # type: ignore[assignment]
        far_counts: Counter[int] = result["far_counts"]  # type: ignore[assignment]
        all_floor_fail.extend(floor_fail)
        all_cap_fail.extend(cap_fail)
        print(
            f"  raw={counts['raw']} primitive_span_gt14={counts['primitive_span_gt14']} "
            f"truewide={counts['truewide']} boundary_skipped={counts['boundary']}"
        )
        print(f"  floor_failures={len(floor_fail)} cap_failures={len(cap_fail)}")
        print("  tightest true-wide floor slacks:")
        for idx, audit in enumerate(result["tight_floor"][:5], 1):  # type: ignore[index]
            all_tight.append(audit)
            print_audit(f"    {idx:2d}.", audit)
        print("  far-count ledger (minimum floor slack by number of >14 speeds):")
        for fc in sorted(by_far):
            best = by_far[fc][0]
            print(
                f"    far_count={fc}: cases={far_counts[fc]} "
                f"best_slack={fmt(best.floor_slack)} best={best.row}"
            )
        print("  largest p6 rows:")
        for idx, audit in enumerate(result["max_p6"][:3], 1):  # type: ignore[index]
            print_audit(f"    p6-{idx}.", audit)

    print("=" * 100)
    print("SYNTHESIS")
    if all_cap_fail:
        print(f"  CAP FAILURES FOUND: {len(all_cap_fail)}")
    else:
        print("  No true-wide cap failures in the audited boxes.")
    k_ge9_fail = [a for a in all_floor_fail if a.k >= 9]
    if k_ge9_fail:
        print(f"  k>=9 FLOOR FAILURES FOUND: {len(k_ge9_fail)}")
    else:
        print("  No k>=9 true-wide floor failures in the audited boxes.")
    k8_fail = [a for a in all_floor_fail if a.k == 8]
    print(f"  k=8 true-wide floor failures in audited boxes: {len(k8_fail)}")
    if k8_fail:
        worst = max(k8_fail, key=lambda a: (a.floor_debt, -a.cap_slack, a.row))
        print_audit("    worst k=8 floor-debt row", worst)
    print()
    print("  Structural reading:")
    print("    * The risky rows are barely true-wide: the tight leaders have exactly two speeds >14.")
    print("    * Rows with three or more far speeds are much looser in the scanned boxes.")
    print("    * The next proof obligation is therefore a two-far survival-currency lemma,")
    print("      followed by a separate, easier r>=3 far-count margin lemma.")
    print()
    tournament_analysis(all_tight)


if __name__ == "__main__":
    main()
