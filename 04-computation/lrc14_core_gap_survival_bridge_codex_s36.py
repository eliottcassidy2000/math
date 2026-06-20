#!/usr/bin/env python3
"""HYP-2658/T903: fixed-observer core-gap survival bridge scout.

This script connects HYP-2651's exact positive-core gap atlas with HYP-2644's
far-element plateau idea, KPS HYP-2653/HYP-2655 decorrelation work, and the
THM-541/542/543/544 near-collar layer.  For a positive core C, write

    G_C = {t in [0,1): ||c t|| > 1/14 for every c in C}.

If a new speed w is genuinely far from a bounded base B, the expected
fixed-observer limit is

    meas(G_{B union {w}}) -> (6/7) meas(G_B),

because an independent point survives the forbidden central sector with
probability 6/7.  The scout tests that lower-bound route and audits the first
finite resonance that beats the old B<=19 second value from HYP-2651.  After
THM-543, that resonance is no longer a conjectural exception; it is a
component-ledger bridge explaining the proved one-replacement AP-tail layer.
THM-544 then closes the two-replacement AP-tail layer, while KPS HYP-2655
pushes the genuinely-wide side from a small uniform constant to a joint
plateau/Delta recursion.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache, reduce
from heapq import heappop, heappush
from itertools import combinations
from math import gcd


TARGET_DENOM = 14
CORE_SIZE = 12
COLLAR = Fraction(7, 858)
OLD_B19_SECOND = Fraction(426, 35035)


def merge(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    for a, b in sorted(intervals):
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            if b > out[-1][1]:
                out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def measure(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((b - a for a, b in intervals), Fraction(0))


@lru_cache(maxsize=None)
def danger_arcs(v: int, denom: int = TARGET_DENOM) -> tuple[tuple[Fraction, Fraction], ...]:
    arcs: list[tuple[Fraction, Fraction]] = []
    d = denom * v
    for a in range(v):
        lo_num = denom * a - 1
        hi_num = denom * a + 1
        if lo_num < 0:
            arcs.append((Fraction(0), Fraction(hi_num, d)))
            arcs.append((Fraction(d + lo_num, d), Fraction(1)))
        elif hi_num > d:
            arcs.append((Fraction(lo_num, d), Fraction(1)))
            arcs.append((Fraction(0), Fraction(hi_num - d, d)))
        else:
            arcs.append((Fraction(lo_num, d), Fraction(hi_num, d)))
    return tuple(arcs)


def safe_components(core: tuple[int, ...], denom: int = TARGET_DENOM) -> tuple[tuple[Fraction, Fraction], ...]:
    danger = merge([arc for v in core for arc in danger_arcs(v, denom)])
    safe: list[tuple[Fraction, Fraction]] = []
    prev = Fraction(0)
    for lo, hi in danger:
        if lo > prev:
            safe.append((prev, lo))
        if hi > prev:
            prev = hi
    if prev < 1:
        safe.append((prev, Fraction(1)))
    return tuple(safe)


@lru_cache(maxsize=None)
def safe_measure(core: tuple[int, ...], denom: int = TARGET_DENOM) -> Fraction:
    return measure(list(safe_components(core, denom)))


def primitive(core: tuple[int, ...]) -> bool:
    return reduce(gcd, core, 0) == 1


def sumset_excess(core: tuple[int, ...]) -> int:
    sums = {a + b for a in core for b in core}
    return len(sums) - (2 * len(core) - 1)


def holes(core: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(set(range(min(core), max(core) + 1)) - set(core)))


def fmt(q: Fraction) -> str:
    return f"{q} = {float(q):.9f}"


@dataclass(frozen=True)
class GapRow:
    core: tuple[int, ...]
    safe: Fraction
    components: int
    min_component: Fraction
    excess: int
    holes: tuple[int, ...]

    @classmethod
    def from_core(cls, core: tuple[int, ...]) -> "GapRow":
        comps = safe_components(core)
        return cls(
            core=core,
            safe=measure(list(comps)),
            components=len(comps),
            min_component=min((b - a for a, b in comps), default=Fraction(0)),
            excess=sumset_excess(core),
            holes=holes(core),
        )


def top_push(heap: list[tuple[Fraction, tuple[int, ...], GapRow]], row: GapRow, keep: int) -> None:
    item = (-row.safe, row.core, row)
    heappush(heap, item)
    if len(heap) > keep:
        heappop(heap)


def smallest_rows(size: int, bmax: int, keep: int) -> tuple[int, list[GapRow], GapRow]:
    heap: list[tuple[Fraction, tuple[int, ...], GapRow]] = []
    best: GapRow | None = None
    count = 0
    for core in combinations(range(1, bmax + 1), size):
        if not primitive(core):
            continue
        count += 1
        row = GapRow.from_core(core)
        top_push(heap, row, keep)
        if best is None or row.safe < best.safe:
            best = row
    if best is None:
        raise ValueError("empty primitive scan")
    return count, sorted((item[2] for item in heap), key=lambda r: (r.safe, r.core)), best


def print_multi_far_ledger(max_base_size: int, keep_cap: int) -> None:
    print("multi-far independent-survival ledger")
    print("For an r-core base, m=12-r independent far speeds predict (6/7)^m * meas(G_base).")
    print("The bmax column is the bounded atlas used for the base minimum in this scout.")
    print(" r | bmax | primitive rows | min G_base | argmin | predicted 12-core floor | margin over collar")
    for r in range(1, max_base_size + 1):
        bmax = min(keep_cap, max(12, r + 8))
        count, _, best = smallest_rows(r, bmax, keep=1)
        predicted = best.safe * (Fraction(6, 7) ** (CORE_SIZE - r))
        print(
            f"{r:2d} | {bmax:4d} | {count:14d} | {str(best.safe):>14} | "
            f"{best.core} | {str(predicted):>18} | {fmt(predicted - COLLAR)}"
        )
    print()


def print_append_resonance(base_bmax: int, keep: int, wmax: int) -> list[tuple[GapRow, Fraction, int, GapRow, Fraction, int]]:
    print(f"one-far append resonance scan: top {keep} 11-core bases in [1,{base_bmax}], append w<= {wmax}")
    count, bases, _ = smallest_rows(11, base_bmax, keep)
    print(f"primitive 11-core rows scanned: {count}")
    print("rank | base gap | base | limit 6/7 | best append gap | w | margin over collar | max w*err")
    out: list[tuple[GapRow, Fraction, int, GapRow, Fraction, int]] = []
    for i, base in enumerate(bases, 1):
        limit = base.safe * Fraction(6, 7)
        best_append: GapRow | None = None
        best_w = -1
        max_werr = Fraction(0)
        max_werr_w = -1
        for w in range(max(base.core) + 1, wmax + 1):
            core = tuple(sorted(base.core + (w,)))
            if not primitive(core):
                continue
            row = GapRow.from_core(core)
            if best_append is None or row.safe < best_append.safe:
                best_append = row
                best_w = w
            werr = abs(row.safe - limit) * w
            if werr > max_werr:
                max_werr = werr
                max_werr_w = w
        if best_append is None:
            continue
        out.append((base, limit, best_w, best_append, max_werr, max_werr_w))
        print(
            f"{i:4d} | {str(base.safe):>10} | {base.core} | {str(limit):>10} | "
            f"{str(best_append.safe):>15} | {best_w:3d} | "
            f"{fmt(best_append.safe - COLLAR):>24} | {str(max_werr):>12} at {max_werr_w}"
        )
    print()
    return out


def endpoint_owners(core: tuple[int, ...], t: Fraction) -> tuple[int, ...]:
    owners: list[int] = []
    for v in core:
        n = TARGET_DENOM * v * t
        if n.denominator == 1 and n.numerator % TARGET_DENOM in (1, TARGET_DENOM - 1):
            owners.append(v)
    return tuple(owners)


def print_component_ledger(core: tuple[int, ...], label: str) -> None:
    print(label)
    row = GapRow.from_core(core)
    print(
        f"core={core}; safe={fmt(row.safe)}; components={row.components}; "
        f"excess={row.excess}; holes={row.holes}"
    )
    for lo, hi in safe_components(core):
        print(
            f"  [{lo}, {hi}] len={hi - lo} mid={(lo + hi) / 2} "
            f"owners={endpoint_owners(core, lo)}->{endpoint_owners(core, hi)}"
        )
    print()


def print_double_grafts() -> None:
    base = set(range(1, 14))
    base.remove(6)
    collar = tuple(sorted(base))
    collar_measure = safe_measure(collar)
    print("drop-6 collar double-graft test: replace h by 2h")
    print(f"collar={collar}; safe={collar_measure}")
    print(" h | graft core | safe | delta over collar")
    for h in sorted(base):
        graft = tuple(sorted((base - {h}) | {2 * h}))
        if len(graft) != CORE_SIZE:
            continue
        m = safe_measure(graft)
        print(f"{h:2d} | {graft} | {str(m):>12} | {fmt(m - collar_measure)}")
    print()


def scan_twelve(bmax: int, keep: int) -> list[GapRow]:
    print(f"exact positive 12-core extension through B={bmax}")
    heap: list[tuple[Fraction, tuple[int, ...], GapRow]] = []
    best: GapRow | None = None
    count = 0
    cumulative: list[tuple[int, int, GapRow]] = []
    for B in range(CORE_SIZE, bmax + 1):
        rows_this_B: list[tuple[int, ...]] = []
        if B == CORE_SIZE:
            rows_this_B.append(tuple(range(1, CORE_SIZE + 1)))
        else:
            for core in combinations(range(1, B), CORE_SIZE - 1):
                rows_this_B.append(tuple(core + (B,)))
        for core in rows_this_B:
            if not primitive(core):
                continue
            count += 1
            row = GapRow.from_core(core)
            top_push(heap, row, keep)
            if best is None or row.safe < best.safe:
                best = row
        if best is None:
            raise ValueError("empty scan")
        cumulative.append((B, count, best))
    top = sorted((item[2] for item in heap), key=lambda r: (r.safe, r.core))

    print("  B | primitive rows | current min | argmin")
    for B, nrows, row in cumulative:
        print(f"{B:3d} | {nrows:14d} | {str(row.safe):>11} | {row.core}")
    print()
    print("near-collar rows in the extended 12-core atlas")
    print("rank | safe | margin over collar | core | exc | comps | mincomp | holes")
    for i, row in enumerate(top, 1):
        print(
            f"{i:4d} | {str(row.safe):>15} | {str(row.safe - COLLAR):>18} | "
            f"{row.core} | {row.excess:3d} | {row.components:5d} | "
            f"{str(row.min_component):>8} | {row.holes}"
        )
    if len(top) > 1:
        print()
        print("second-value correction")
        print(f"  old HYP-2651 B<=19 second value: {fmt(OLD_B19_SECOND)}")
        print(f"  extended B<={bmax} second value: {fmt(top[1].safe)}")
        print(f"  extended second minus collar: {fmt(top[1].safe - COLLAR)}")
    print()
    return top


def print_tournament_analysis() -> None:
    print("Tournament Analysis")
    print("Vertices are proof quotients, not runners:")
    vertices = [
        "proved_one_replacement_tail",
        "proved_two_replacement_tail",
        "exact_core_gap_components",
        "drop6_collar_ledger",
        "owner_bubble_ledger",
        "far_survival_multiplier",
        "joint_plateau_delta_recursion",
        "finite_resonance_discrepancy",
        "sumset_excess",
        "raw_far_speed",
    ]
    print("Pairwise observable: which quotient preserves the lower-bound predicate and explains the sub-426/35035 rows.")
    print("Switch/gauge: keep fixed-observer positivity; judge tails after THM-543/544 and the 6/7 survival quotient.")
    print("Hamiltonian path:")
    print("  " + " > ".join(vertices))
    print("Fingerprint: transitive proof-priority tournament; score histogram {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}; directed 3-cycles: 0")
    print()
    print("Assumption challenge")
    print(
        "Alternate vertices considered: runners, far speed w, missing holes, safe components, "
        "component endpoint owners, sumset excess, state words, and proof obligations.  The "
        "winner is not raw w or raw excess: it is the proved AP-tail layers plus the "
        "endpoint-owner component ledger, because it shows the B=20 near-collar row is the old "
        "drop-6 collar plus two owner-addressed bubbles from the 10->20 graft."
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger-max-base-size", type=int, default=11)
    parser.add_argument("--ledger-bmax-cap", type=int, default=20)
    parser.add_argument("--base-bmax", type=int, default=19)
    parser.add_argument("--base-keep", type=int, default=12)
    parser.add_argument("--append-wmax", type=int, default=240)
    parser.add_argument("--scan12-bmax", type=int, default=21)
    parser.add_argument("--scan12-keep", type=int, default=12)
    args = parser.parse_args()

    print("HYP-2658/T903 LRC14 core-gap survival bridge scout")
    print(f"target gap: 1/{TARGET_DENOM}; known collar: {fmt(COLLAR)}")
    print("Integrates HYP-2651 exact core gaps, KPS HYP-2653/HYP-2655 decorrelation, HYP-2654/THM-543/544 near-collar structure, HYP-2648 state words, and HYP-2652 endpoint geometry.")
    print("Caveat: this is exact finite evidence plus a proof route; LRC(14) is not proved here.")
    print()

    print_multi_far_ledger(args.ledger_max_base_size, args.ledger_bmax_cap)
    append_rows = print_append_resonance(args.base_bmax, args.base_keep, args.append_wmax)
    print_double_grafts()

    collar_core = (1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13)
    graft_core = (1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 20)
    old_second_core = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13)
    print_component_ledger(collar_core, "component ledger: drop-6 collar")
    print_component_ledger(graft_core, "component ledger: B=20 near-collar 10->20 graft")
    print_component_ledger(old_second_core, "component ledger: old B<=19 second row")

    top = scan_twelve(args.scan12_bmax, args.scan12_keep)
    if append_rows:
        best_append = min((item[3] for item in append_rows), key=lambda row: row.safe)
        print("cross-check")
        print(f"  best appended row among scanned 11-bases: {best_append.core} with {fmt(best_append.safe)}")
        print(f"  best extended 12-core row after collar: {top[1].core} with {fmt(top[1].safe)}")
        print(f"  match? {best_append.core == top[1].core and best_append.safe == top[1].safe}")
        print()

    print("proof reading")
    print(
        "  The independent far-speed ledger has large margins: in the tested base atlases, even "
        "the closest r=11 one-far limit is 313/11319, about 0.01949 above the collar.  The "
        "danger is therefore not the far limit itself but finite resonance before decorrelation."
    )
    print(
        "  The first resonance is structured: replacing 10 by 20 in the drop-6 collar leaves the "
        "four old collar components intact and adds exactly two symmetric bubbles of length 1/1960. "
        "THM-543 proves this 10->20 row is the unique one-replacement AP-tail exception below 426/35035."
    )
    print(
        "  Updated target after THM-541/542/543/544: keep the proved collar, mouth-retention, "
        "one-replacement, and two-replacement layers as the near-collar base, then route "
        "genuinely far rows through fixed-observer 6/7 survival plus the KPS HYP-2655 joint "
        "plateau/Delta recursion."
    )
    print()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
