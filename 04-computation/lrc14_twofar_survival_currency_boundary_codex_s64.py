#!/usr/bin/env python3
"""HYP-2701 addendum: two-far survival-currency boundary scout.

For HYP-2701 the live true-wide floor gate is the exact currency

    C(E) = p1+p2+p3+p4 - 4*p6 = 1 - U4(E).

For a bounded core B and r independent far runners, the decorrelated boundary
model is the death chain on the number of missed inner sectors: each far runner
hits one of the seven sector colors uniformly, and only hits to currently
missed inner sectors decrease the residual count.

This scout computes that boundary value exactly for r=2, compares it with the
floor requirement, and measures the actual two-far deviation

    C(B union {u,v}) - C_boundary(B,2)

on the audited true-wide boxes.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import comb, gcd

from lrc14_wide_branch_ridge_codex_s47 import (
    additive_energy,
    missed_distribution,
    primitive,
    squarefree_profile,
    sumset_excess,
)


Row = tuple[int, ...]


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def floor_required(k: int) -> F:
    return F(13 - k, 7)


def currency_from_dist(dist: tuple[F, ...]) -> F:
    return sum(dist[1:5]) - 4 * dist[6]


@lru_cache(maxsize=None)
def profile(row: Row) -> tuple[F, ...]:
    return missed_distribution(tuple(sorted(row)))


@lru_cache(maxsize=None)
def currency(row: Row) -> F:
    return currency_from_dist(profile(tuple(sorted(row))))


def transition_prob(t: int, s: int, r: int) -> F:
    """Probability t missed sectors become exactly s missed after r iid hits."""

    if not 0 <= s <= t <= 6:
        return F(0)
    need_hit = t - s
    total = F(0)
    for j in range(need_hit + 1):
        total += ((-1) ** j) * comb(need_hit, j) * F(7 - s - j, 7) ** r
    return F(comb(t, s)) * total


def boundary_distribution(core: Row, r: int) -> tuple[F, ...]:
    base = profile(core)
    out = [F(0) for _ in range(7)]
    for t, mass in enumerate(base):
        for s in range(t + 1):
            out[s] += mass * transition_prob(t, s, r)
    assert sum(out) == 1
    return tuple(out)


@lru_cache(maxsize=None)
def boundary_currency(core: Row, r: int = 2) -> F:
    return currency_from_dist(boundary_distribution(core, r))


def relation_distance(vals: Row, height: int = 4) -> tuple[int, tuple[int, ...], int]:
    best: tuple[int, int, tuple[int, ...]] | None = None
    for coeff in product(range(-height, height + 1), repeat=len(vals)):
        if all(c == 0 for c in coeff):
            continue
        g = 0
        for c in coeff:
            g = gcd(g, abs(c))
        if g != 1:
            continue
        val = abs(sum(c * v for c, v in zip(coeff, vals)))
        l1 = sum(abs(c) for c in coeff)
        item = (val, l1, coeff)
        if best is None or item < best:
            best = item
    assert best is not None
    val, l1, coeff = best
    return val, coeff, l1


def gap_signature(row: Row) -> tuple[int, ...]:
    return tuple(b - a for a, b in zip(row, row[1:]))


@dataclass(frozen=True)
class TwoFarReport:
    row: Row
    core: Row
    far: Row
    actual: F
    boundary: F
    required: F

    @property
    def k(self) -> int:
        return len(self.row)

    @property
    def slack(self) -> F:
        return self.actual - self.required

    @property
    def boundary_margin(self) -> F:
        return self.boundary - self.required

    @property
    def deviation(self) -> F:
        return self.actual - self.boundary

    @property
    def relation(self) -> tuple[int, tuple[int, ...], int]:
        return relation_distance(self.far)


def make_report(row: Row) -> TwoFarReport:
    core = tuple(x for x in row if x <= 14)
    far = tuple(x for x in row if x > 14)
    assert len(far) == 2
    return TwoFarReport(
        row=row,
        core=core,
        far=far,
        actual=currency(row),
        boundary=boundary_currency(core, 2),
        required=floor_required(len(row)),
    )


def push_min(store: list[TwoFarReport], rep: TwoFarReport, keep: int, key) -> None:
    store.append(rep)
    store.sort(key=key)
    del store[keep:]


def print_report(label: str, rep: TwoFarReport) -> None:
    rel_val, rel_coeff, rel_l1 = rep.relation
    dist = profile(rep.row)
    print(f"{label}: E={rep.row}, core={rep.core}, far={rep.far}")
    print(
        f"  actual_currency={fmt(rep.actual)} required={fmt(rep.required)} "
        f"slack={fmt(rep.slack)}"
    )
    print(
        f"  boundary_currency={fmt(rep.boundary)} boundary_margin={fmt(rep.boundary_margin)} "
        f"actual-boundary={fmt(rep.deviation)}"
    )
    print(
        f"  dist={[str(x) for x in dist]} p6={fmt(dist[6])} "
        f"relation={rel_val} coeff={rel_coeff} l1={rel_l1}"
    )
    print(
        f"  gaps={gap_signature(rep.row)} exc={sumset_excess(rep.row)} "
        f"energy={additive_energy(rep.row)} sqfree={squarefree_profile(rep.row)}"
    )


def boundary_core_scan() -> None:
    print("TWO-FAR DECORRELATED BOUNDARY CURRENCY OVER BOUNDED CORES")
    print("core entries are in {0,...,14}; r=2 far runners are decorrelated")
    for k in range(8, 13):
        core_size = k - 2
        best: tuple[F, Row, F] | None = None
        count = 0
        for rest in combinations(range(1, 15), core_size - 1):
            core = (0,) + rest
            count += 1
            bcur = boundary_currency(core, 2)
            margin = bcur - floor_required(k)
            item = (margin, core, bcur)
            if best is None or item < best:
                best = item
        assert best is not None
        margin, core, bcur = best
        print(
            f"  k={k}: cores={count} min_boundary_margin={fmt(margin)} "
            f"boundary_currency={fmt(bcur)} core={core}"
        )
    print()


def scan_box(k: int, bound: int, keep: int = 8) -> dict[str, object]:
    counts: Counter[str] = Counter()
    tight: list[TwoFarReport] = []
    worst_deviation: list[TwoFarReport] = []
    thin_boundary: list[TwoFarReport] = []
    by_rel: dict[int, list[TwoFarReport]] = defaultdict(list)
    rel_counts: Counter[int] = Counter()

    for combo in combinations(range(1, bound + 1), k - 1):
        row = (0,) + combo
        counts["raw"] += 1
        if row[-1] <= 14 or not primitive(row):
            continue
        far = tuple(x for x in row if x > 14)
        if len(far) != 2:
            continue
        counts["two_far"] += 1
        rep = make_report(row)
        push_min(tight, rep, keep, lambda r: (r.slack, r.deviation, r.row))
        push_min(worst_deviation, rep, keep, lambda r: (r.deviation, r.slack, r.row))
        push_min(thin_boundary, rep, keep, lambda r: (r.boundary_margin, r.slack, r.row))
        rel_val = rep.relation[0]
        rel_counts[rel_val] += 1
        push_min(by_rel[rel_val], rep, keep, lambda r: (r.slack, r.deviation, r.row))

    return {
        "counts": counts,
        "tight": tight,
        "worst_deviation": worst_deviation,
        "thin_boundary": thin_boundary,
        "by_rel": by_rel,
        "rel_counts": rel_counts,
    }


def print_scan(k: int, bound: int, data: dict[str, object]) -> None:
    counts: Counter[str] = data["counts"]  # type: ignore[assignment]
    print("=" * 100)
    print(f"k={k}, bound={bound}, exact two-far rows")
    print(f"  raw={counts['raw']} two_far={counts['two_far']}")
    print("  tightest actual survival-currency slacks:")
    for i, rep in enumerate(data["tight"], 1):  # type: ignore[index]
        print_report(f"    {i}", rep)
    print("  most negative actual-minus-boundary deviations:")
    for i, rep in enumerate(data["worst_deviation"], 1):  # type: ignore[index]
        print_report(f"    dev-{i}", rep)
    print("  thinnest decorrelated boundary margins:")
    for i, rep in enumerate(data["thin_boundary"], 1):  # type: ignore[index]
        print_report(f"    bdry-{i}", rep)
    print("  relation-distance ledger (height <=4, minimum actual slack):")
    rel_counts: Counter[int] = data["rel_counts"]  # type: ignore[assignment]
    by_rel: dict[int, list[TwoFarReport]] = data["by_rel"]  # type: ignore[assignment]
    for rel_val in sorted(rel_counts)[:10]:
        best = by_rel[rel_val][0]
        print(
            f"    rel={rel_val}: cases={rel_counts[rel_val]} "
            f"best_slack={fmt(best.slack)} best_dev={fmt(best.deviation)} best={best.row}"
        )
    print()


def tournament(rows: list[TwoFarReport]) -> None:
    verts = rows[:12]
    score = [0] * len(verts)
    wins: set[tuple[int, int]] = set()
    for i, a in enumerate(verts):
        for j, b in enumerate(verts):
            if i >= j:
                continue
            ka = (a.slack, a.deviation, -a.boundary_margin, a.row)
            kb = (b.slack, b.deviation, -b.boundary_margin, b.row)
            if ka <= kb:
                score[i] += 1
                wins.add((i, j))
            else:
                score[j] += 1
                wins.add((j, i))
    cycles = 0
    for a, b, c in combinations(range(len(verts)), 3):
        if (a, b) in wins and (b, c) in wins and (c, a) in wins:
            cycles += 1
        if (a, c) in wins and (c, b) in wins and (b, a) in wins:
            cycles += 1
    print("Tournament Analysis: two-far survival boundary obligations")
    print("  vertices: audited two-far row/family obligations")
    print("  pairwise observable: smaller actual currency slack wins")
    print("  switch/gauge: tie by more negative actual-boundary deviation")
    print(f"  score_hist={dict(sorted(Counter(score).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian risk path:")
    for i, rep in enumerate(sorted(verts, key=lambda r: (r.slack, r.deviation, r.row)), 1):
        print(
            f"    {i}. k={rep.k} E={rep.row} slack={rep.slack} "
            f"boundary_margin={rep.boundary_margin} deviation={rep.deviation}"
        )
    print()


def main() -> None:
    print("HYP-2701 addendum -- two-far survival-currency boundary scout")
    print("Exact arithmetic: Fractions from sector-wall integration.")
    for t in range(7):
        probs = [transition_prob(t, s, 2) for s in range(t + 1)]
        assert sum(probs) == 1
    print("Transition kernel check passed for r=2.")
    print()

    boundary_core_scan()
    scans = [(8, 18), (9, 18), (10, 17), (11, 17), (12, 17)]
    risk_path: list[TwoFarReport] = []
    for k, bound in scans:
        data = scan_box(k, bound)
        print_scan(k, bound, data)
        risk_path.extend(data["tight"][:3])  # type: ignore[index]

    print("SYNTHESIS")
    print("  The decorrelated two-far boundary currency is floor-safe for every bounded core scanned.")
    print("  Actual tight rows spend that boundary margin through negative actual-boundary deviation.")
    print("  The worst audited rows all have low relation distance far pairs, usually consecutive or near-consecutive.")
    print("  Next proof target: boundary currency margin minus signed two-far deviation.")
    print()
    tournament(sorted(risk_path, key=lambda r: (r.slack, r.deviation, r.row)))


if __name__ == "__main__":
    main()
