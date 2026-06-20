#!/usr/bin/env python3
"""HYP-2695 / T930: true-wide cap-floor gate for LRC14.

Refine HYP-2693 by splitting the exact cap into two currencies:

    floor_k    = (k - 6) / 7
    dividend_k = cap_k - floor_k

THM-535 proves cap_k >= floor_k.  The stronger target tested here is

    second_largest(E) > 14, k >= 9  =>  U4(E) <= floor_k,

where THM-556 gives the exact six-sector collapse U4=p0+p5+5p6.
The scan also records the k=8 cap-dividend exception.
"""

from __future__ import annotations

from collections import Counter
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


def bonferroni_upper(dist: tuple[Fraction, ...], level: int = 4) -> Fraction:
    total = Fraction(1)
    sign = -1
    for r in range(1, level + 1):
        sr = sum(Fraction(comb(t, r)) * dist[t] for t in range(r, len(dist)))
        total += sign * sr
        sign *= -1
    return total


def tail45(dist: tuple[Fraction, ...]) -> Fraction:
    return dist[5] + 5 * dist[6]


@dataclass(frozen=True)
class Audit:
    label: str
    row: Row
    dist: tuple[Fraction, ...]

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
    def p0(self) -> Fraction:
        return self.dist[0]

    @property
    def tail(self) -> Fraction:
        return tail45(self.dist)

    @property
    def u4(self) -> Fraction:
        return bonferroni_upper(self.dist, 4)

    @property
    def floor_slack(self) -> Fraction:
        return self.floor - self.u4

    @property
    def cap_slack(self) -> Fraction:
        return self.cap - self.u4

    @property
    def floor_debt(self) -> Fraction:
        return max(Fraction(0), -self.floor_slack)


def assert_tail_identity(audit: Audit) -> None:
    assert audit.u4 == audit.p0 + audit.tail


def moments(dist: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(
        sum(Fraction(comb(t, r)) * dist[t] for t in range(r, len(dist)))
        for r in range(7)
    )


def structural(row: Row) -> str:
    gaps = tuple(b - a for a, b in zip(row, row[1:]))
    return (
        f"second={row[-2]} span={row[-1]} maxgap={max(gaps) if gaps else 0} "
        f"run={longest_run(row)} exc={sumset_excess(row)} "
        f"energy={additive_energy(row)} sqfree={squarefree_profile(row)}"
    )


def print_audit(prefix: str, audit: Audit) -> None:
    ms = moments(audit.dist)
    print(
        f"{prefix} E={audit.row} p0={fmt(audit.p0)} tail45={fmt(audit.tail)} "
        f"U4={fmt(audit.u4)}"
    )
    print(
        f"      floor={fmt(audit.floor)} floor-U4={fmt(audit.floor_slack)} "
        f"cap={fmt(audit.cap)} dividend={fmt(audit.cap-audit.floor)} "
        f"cap-U4={fmt(audit.cap_slack)}"
    )
    print(f"      {structural(audit.row)}")
    print(f"      S1..S4={[str(x) for x in ms[1:5]]} dist={[str(x) for x in audit.dist]}")


def push_top(bucket: list[Audit], audit: Audit, keep: int, key) -> None:
    bucket.append(audit)
    bucket.sort(key=key)
    del bucket[keep:]


def scan_box(k: int, bound: int, keep: int = 8) -> dict[str, object]:
    counts = Counter()
    true_floor: list[Audit] = []
    true_cap: list[Audit] = []
    true_tail: list[Audit] = []
    boundary_floor: list[Audit] = []
    floor_v_true: list[Audit] = []
    cap_v_true: list[Audit] = []
    floor_v_boundary: list[Audit] = []
    cap_v_boundary: list[Audit] = []

    for combo in combinations(range(1, bound + 1), k - 1):
        row = (0,) + combo
        counts["raw"] += 1
        if row[-1] <= 14 or not primitive(row):
            continue
        counts["primitive_span_gt14"] += 1
        stratum = "truewide" if row[-2] > 14 else "boundary"
        counts[stratum] += 1
        audit = Audit(f"k{k}-B{bound}-{stratum}", row, missed_distribution(row))
        assert_tail_identity(audit)
        if stratum == "truewide":
            push_top(true_floor, audit, keep, key=lambda a: (a.floor_slack, a.cap_slack, a.row))
            push_top(true_cap, audit, keep, key=lambda a: (a.cap_slack, a.floor_slack, a.row))
            push_top(true_tail, audit, keep, key=lambda a: (-a.tail, a.floor_slack, a.row))
            if audit.floor_slack < 0:
                floor_v_true.append(audit)
            if audit.cap_slack < 0:
                cap_v_true.append(audit)
        else:
            push_top(boundary_floor, audit, keep, key=lambda a: (a.floor_slack, a.cap_slack, a.row))
            if audit.floor_slack < 0:
                floor_v_boundary.append(audit)
            if audit.cap_slack < 0:
                cap_v_boundary.append(audit)

    return {
        "counts": counts,
        "true_floor": true_floor,
        "true_cap": true_cap,
        "true_tail": true_tail,
        "boundary_floor": boundary_floor,
        "floor_v_true": floor_v_true,
        "cap_v_true": cap_v_true,
        "floor_v_boundary": floor_v_boundary,
        "cap_v_boundary": cap_v_boundary,
    }


def named_rows() -> list[tuple[str, Row]]:
    return [
        ("AP9_floor_reference", tuple(range(9))),
        ("doubled_AP_boundary", (0, 2, 4, 6, 8, 10, 12, 14, 16)),
        ("boundary_leader", (0, 2, 4, 6, 8, 10, 12, 14, 15)),
        ("k8_tight_truewide", (0, 3, 6, 9, 12, 14, 15, 18)),
        ("k9_direct_truewide", (0, 4, 6, 8, 10, 12, 14, 15, 16)),
        ("k10_direct_truewide", (0, 2, 4, 6, 8, 10, 12, 14, 15, 16)),
        ("three_cluster_truewide", (0, 1, 2, 30, 31, 32, 60, 61, 62)),
        ("dyadic_block", (0, 1, 2, 4, 8, 12, 16, 20, 24)),
    ]


def tournament_analysis(audits: list[Audit], limit: int = 12) -> None:
    """Tournament on proof-obligation vertices, not runners or arcs."""

    def risk(a: Audit) -> tuple[Fraction, Fraction, Fraction]:
        return (a.floor_debt, -a.floor_slack, -a.cap_slack)

    dedup: dict[tuple[int, Row], Audit] = {}
    for audit in audits:
        key = (audit.k, audit.row)
        if key not in dedup or risk(audit) > risk(dedup[key]):
            dedup[key] = audit
    verts = sorted(dedup.values(), key=lambda a: (risk(a), tuple(-x for x in a.row)), reverse=True)[:limit]
    wins: set[tuple[int, int]] = set()
    score = Counter()
    for i, a in enumerate(verts):
        for j, b in enumerate(verts):
            if i >= j:
                continue
            winner, loser = (i, j) if (risk(a), tuple(-x for x in a.row)) >= (risk(b), tuple(-x for x in b.row)) else (j, i)
            wins.add((winner, loser))
            score[winner] += 1
            score.setdefault(loser, score[loser])
    cycles = 0
    for i, j, k in combinations(range(len(verts)), 3):
        if ((i, j) in wins and (j, k) in wins and (k, i) in wins) or (
            (i, k) in wins and (k, j) in wins and (j, i) in wins
        ):
            cycles += 1
    path = sorted(range(len(verts)), key=lambda idx: (-score[idx], -risk(verts[idx])[0], verts[idx].row))
    print("Tournament Analysis: cap/floor proof-obligation quotient")
    print("  vertices: row/stratum obligations, not runners/arcs")
    print("  pairwise observable: larger floor debt and smaller floor/cap slack wins")
    print("  switch/gauge: orient toward the row consuming more cap currency")
    print(f"  score_hist={dict(sorted(Counter(score.values()).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian risk path:")
    for rank, idx in enumerate(path, 1):
        audit = verts[idx]
        print(
            f"    {rank:2d}. {audit.label} E={audit.row} "
            f"debt={audit.floor_debt} floor_slack={audit.floor_slack} cap_slack={audit.cap_slack}"
        )


def main() -> None:
    print("HYP-2695 / T930 -- LRC14 true-wide cap-floor gate")
    print("Exact arithmetic: Fractions over the sector-wall common denominator.")
    print("Predicate tested: second-largest > 14, k>=9 => U4 <= floor_k=(k-6)/7.")
    print("U4 is THM-556's six-sector collapse U4=p0+p5+5p6.\n")

    print("CAP CURRENCY TABLE")
    for k in sorted(CAP):
        fl = floor_cap(k)
        print(f"  k={k}: floor={fmt(fl)} cap={fmt(CAP[k])} dividend={fmt(CAP[k]-fl)}")
    print()

    print("NAMED ROWS")
    all_audits: list[Audit] = []
    for label, row in named_rows():
        audit = Audit(label, row, missed_distribution(row))
        assert_tail_identity(audit)
        all_audits.append(audit)
        print_audit(f"  {label:<24}", audit)
    print()

    all_floor_v: list[tuple[int, int, Audit]] = []
    all_cap_v: list[tuple[int, int, Audit]] = []
    for k, bound in ((8, 18), (9, 18), (10, 16), (11, 16), (12, 16)):
        print("=" * 96)
        print(f"BOX k={k}, bound={bound}")
        result = scan_box(k, bound)
        counts: Counter = result["counts"]  # type: ignore[assignment]
        fv: list[Audit] = result["floor_v_true"]  # type: ignore[assignment]
        cv: list[Audit] = result["cap_v_true"]  # type: ignore[assignment]
        bfv: list[Audit] = result["floor_v_boundary"]  # type: ignore[assignment]
        bcv: list[Audit] = result["cap_v_boundary"]  # type: ignore[assignment]
        all_floor_v.extend((k, bound, a) for a in fv)
        all_cap_v.extend((k, bound, a) for a in cv)
        print(
            f"  raw={counts['raw']} primitive_span_gt14={counts['primitive_span_gt14']} "
            f"truewide={counts['truewide']} boundary={counts['boundary']}"
        )
        print(
            f"  truewide floor violations={len(fv)} cap violations={len(cv)}; "
            f"boundary floor violations={len(bfv)} cap violations={len(bcv)}"
        )
        print("  tightest truewide floor slacks:")
        for idx, audit in enumerate(result["true_floor"][:5], 1):  # type: ignore[index]
            all_audits.append(audit)
            print_audit(f"    {idx:2d}.", audit)
        print("  tightest boundary floor slacks:")
        for idx, audit in enumerate(result["boundary_floor"][:3], 1):  # type: ignore[index]
            print_audit(f"    {idx:2d}.", audit)
        print()

    print("=" * 96)
    print("SYNTHESIS")
    if all_cap_v:
        print(f"  TRUEWIDE CAP FAILURES FOUND: {len(all_cap_v)}")
    else:
        print("  No truewide U4>cap violations in the audited exact boxes.")
    k_ge9_floor_v = [item for item in all_floor_v if item[0] >= 9]
    if k_ge9_floor_v:
        print(f"  TRUEWIDE k>=9 FLOOR FAILURES FOUND: {len(k_ge9_floor_v)}")
        for k, bound, audit in k_ge9_floor_v[:10]:
            print_audit(f"    k={k} B={bound}", audit)
    else:
        print("  No truewide k>=9 U4>floor violations in the audited exact boxes.")
    k8_floor_v = [audit for k, _, audit in all_floor_v if k == 8]
    if k8_floor_v:
        worst = max(k8_floor_v, key=lambda a: (a.floor_debt, -a.cap_slack, a.row))
        print(f"  k=8 truewide floor failures exist: {len(k8_floor_v)} rows in B18.")
        print_audit("    worst k=8 floor debt", worst)
        print("  Reading: the k=8 exact cap dividend, not the subadditive floor, is needed here.")
    print()
    tournament_analysis(all_audits)
    print()
    print("PROOF-ROUTE READING")
    print("  1. THM-535 proves cap_k >= floor_k.  A true-wide proof of U4<=floor_k")
    print("     automatically closes the HYP-2693 cap gate without exact cap minimizers.")
    print("  2. The audited exception is k=8: true-wide rows can exceed the floor, but")
    print("     their floor debt is covered by the exact dividend cap_8-floor_8.")
    print("  3. Boundary/AP rows remain outside this gate and route to HYP-2691 templates.")


if __name__ == "__main__":
    main()
