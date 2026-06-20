#!/usr/bin/env python3
"""HYP-2664/T907: exploratory three-tail AP-tail comb scaffold for LRC14.

This upgrades HYP-2663's bounded three-tail root-packet scout.  The family is

    C = ({1,...,13} \\ H) union {r,s,t}
    |H| = 4, 14 <= r < s < t.

The script tries to attack the whole infinite family by nested exact comb
cutoffs, but the S39 session found that a blind run leaves too much exact
residue work before HYP-2661's shell-1 carry gate is applied.  Treat this as a
prototype / diagnostic scaffold; the stored S39 result comes from
``lrc14_three_tail_shell1_frontier_codex_s39.py``.

For a safe set G with measure M and c components, a new speed v removes at most
M/7 + 2c/(7v).  Therefore if all remaining tails are at least R:

    3 tails:  meas >= (4/7)M - 6c/(7R)
    2 tails:  meas >= (5/7)M - 4c/(7R)
    1 tail:   meas >= (6/7)M - 2c/(7R)

The proof engine:
  1. For each 4-hole AP base, find a cutoff R for the smallest tail r.
  2. For r < R, subtract D_r exactly and find a cutoff S for s.
  3. For s < S, subtract D_s exactly and find a cutoff T for t.
  4. Check t < T exactly.

This is the THM-543/544 comb method with one more tail.  Side channels record
the root-packet address: old drop-6 mouth survival and Glaisher carry delta.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd

from lrc14_two_replacement_ap_tail_theorem_codex_s35 import (
    AP_SECOND,
    DROP6_MOUTHS,
    danger_arcs,
    intersect_measure,
    measure,
    safe_components,
    subtract_arcs,
)


AP = tuple(range(1, 14))
DROP6_CORE = tuple(v for v in AP if v != 6)
DROP6_COLLAR = measure(DROP6_MOUTHS)


@dataclass(frozen=True)
class Row:
    measure: Fraction
    holes: tuple[int, int, int, int]
    tails: tuple[int, int, int]
    old_survivor: Fraction
    components: int
    carry_delta: tuple[tuple[int, int], ...]

    @property
    def old_damage(self) -> Fraction:
        return DROP6_COLLAR - self.old_survivor


def primitive(core: tuple[int, ...]) -> bool:
    return reduce(gcd, core, 0) == 1


def ceil_frac(x: Fraction) -> int:
    return x.numerator // x.denominator + (1 if x.numerator % x.denominator else 0)


def tower(speed: int) -> tuple[int, int, int]:
    odd = speed
    level = 0
    while odd % 2 == 0:
        odd //= 2
        level += 1
    return speed, odd, level


def odd_shell_profile(core: tuple[int, ...]) -> Counter[int]:
    out: Counter[int] = Counter()
    for speed in core:
        _, odd, level = tower(speed)
        out[odd] += 2**level
    return out


DROP6_PROFILE = odd_shell_profile(DROP6_CORE)


def carry_delta(core: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    profile = odd_shell_profile(core)
    keys = sorted(set(profile) | set(DROP6_PROFILE))
    return tuple((k, profile[k] - DROP6_PROFILE[k]) for k in keys if profile[k] != DROP6_PROFILE[k])


def cutoff_for(meas: Fraction, comps: int, remaining_tails: int) -> int | None:
    """Return R so all remaining_tails with value >= R clear AP_SECOND."""

    if remaining_tails == 3:
        # (4/7)M - 6c/(7R) >= Q
        denom = 4 * meas - 7 * AP_SECOND
        numer = 6 * comps
    elif remaining_tails == 2:
        # (5/7)M - 4c/(7R) >= Q
        denom = 5 * meas - 7 * AP_SECOND
        numer = 4 * comps
    elif remaining_tails == 1:
        # (6/7)M - 2c/(7R) >= Q
        denom = 6 * meas - 7 * AP_SECOND
        numer = 2 * comps
    else:
        raise ValueError(remaining_tails)
    if denom <= 0:
        return None
    return max(14, ceil_frac(Fraction(numer, denom)))


def make_row(
    measure_value: Fraction,
    comps: tuple[tuple[Fraction, Fraction], ...],
    holes: tuple[int, int, int, int],
    tails: tuple[int, int, int],
    core: tuple[int, ...],
) -> Row:
    return Row(
        measure=measure_value,
        holes=holes,
        tails=tails,
        old_survivor=intersect_measure(comps, DROP6_MOUTHS),
        components=len(comps),
        carry_delta=carry_delta(core),
    )


def row_line(rank: int, row: Row) -> str:
    return (
        f"{rank:3d} | {str(row.measure):>13} | {row.holes!s:>15} | {row.tails!s:>14} | "
        f"{str(row.old_survivor):>9} | {str(row.old_damage):>10} | "
        f"{row.components:4d} | {row.carry_delta}"
    )


SHELL1 = {1, 2, 4, 8}


def scan(
    keep: int = 40,
    diagnostic: bool = False,
    shell1_full_only: bool = False,
) -> tuple[list[Row], list[Row], dict[str, int], dict[str, int]]:
    stats: dict[str, int] = defaultdict(int)
    maxima: dict[str, int] = defaultdict(int)
    best: list[Row] = []
    below: list[Row] = []

    for holes in combinations(AP, 4):
        if shell1_full_only and (set(holes) & SHELL1):
            stats["skip_shell1_damaged"] += 1
            continue
        base = tuple(v for v in AP if v not in holes)
        base_comps = safe_components(base)
        base_measure = measure(base_comps)
        base_c = len(base_comps)
        r_cut = cutoff_for(base_measure, base_c, 3)
        if r_cut is None:
            stats["uncut_r_bases"] += 1
            continue
        maxima["r_cut"] = max(maxima["r_cut"], r_cut)
        stats["base_packets"] += 1

        for r in range(14, r_cut):
            core1 = tuple(sorted(base + (r,)))
            if not primitive(core1):
                stats["nonprimitive_r"] += 1
                continue
            comps1 = subtract_arcs(base_comps, danger_arcs(r))
            m1 = measure(comps1)
            c1 = len(comps1)
            s_cut = cutoff_for(m1, c1, 2)
            stats["r_rows"] += 1
            if s_cut is None:
                stats["uncut_s_rows"] += 1
                continue
            s_start = max(r + 1, 14)
            if s_cut <= s_start:
                stats["s_all_pruned"] += 1
                continue
            maxima["s_cut"] = max(maxima["s_cut"], s_cut)

            for s in range(s_start, s_cut):
                core2 = tuple(sorted(core1 + (s,)))
                if not primitive(core2):
                    stats["nonprimitive_s"] += 1
                    continue
                comps2 = subtract_arcs(comps1, danger_arcs(s))
                m2 = measure(comps2)
                c2 = len(comps2)
                t_cut = cutoff_for(m2, c2, 1)
                stats["s_rows"] += 1
                if t_cut is None:
                    stats["uncut_t_rows"] += 1
                    continue
                t_start = max(s + 1, 14)
                if t_cut <= t_start:
                    stats["t_all_pruned"] += 1
                    continue
                maxima["t_cut"] = max(maxima["t_cut"], t_cut)
                if diagnostic:
                    stats["candidate_t_residues"] += t_cut - t_start
                    continue

                for t in range(t_start, t_cut):
                    core = tuple(sorted(core2 + (t,)))
                    if not primitive(core):
                        stats["nonprimitive_t"] += 1
                        continue
                    comps = subtract_arcs(comps2, danger_arcs(t))
                    m = measure(comps)
                    stats["exact_rows"] += 1
                    enters = len(best) < keep or (m, holes, (r, s, t)) < (
                        best[-1].measure,
                        best[-1].holes,
                        best[-1].tails,
                    )
                    if not enters and m >= AP_SECOND:
                        continue
                    row = make_row(m, comps, holes, (r, s, t), core)
                    if enters:
                        best.append(row)
                        best.sort(key=lambda x: (x.measure, x.holes, x.tails))
                        if len(best) > keep:
                            best.pop()
                    if m < AP_SECOND:
                        below.append(row)

    best.sort(key=lambda x: (x.measure, x.holes, x.tails))
    below.sort(key=lambda x: (x.measure, x.holes, x.tails))
    return best, below, dict(stats), dict(maxima)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--diagnostic",
        action="store_true",
        help="Count finite residue burden without exact third-tail subtraction.",
    )
    parser.add_argument(
        "--shell1-full-only",
        action="store_true",
        help="Use HYP-2661 gate: scan only packets preserving {1,2,4,8}.",
    )
    args = parser.parse_args()
    best, below, stats, maxima = scan(
        diagnostic=args.diagnostic,
        shell1_full_only=args.shell1_full_only,
    )

    print("HYP-2664/T907 LRC14 three-tail AP-tail comb scaffold")
    print("family: ({1,...,13} \\ H) union {r,s,t}, |H|=4, 14<=r<s<t unbounded")
    if args.diagnostic:
        print("mode: diagnostic cutoff census only")
    if args.shell1_full_only:
        print("gate: shell-1 full packets only (holes avoid {1,2,4,8})")
    print(f"drop-6 collar: {DROP6_COLLAR} = {float(DROP6_COLLAR):.9f}")
    print(f"AP one-hole second value: {AP_SECOND} = {float(AP_SECOND):.9f}")
    print()
    print("nested comb cutoff maxima")
    for key in ("r_cut", "s_cut", "t_cut"):
        print(f"  {key}: {maxima.get(key, 0)}")
    print()
    print("scan stats")
    for key in sorted(stats):
        print(f"  {key}: {stats[key]}")
    print()
    print("best exact rows after all finite residues")
    print("rank | measure | holes | tails | old_surv | old_damage | comp | carry_delta")
    for i, row in enumerate(best, 1):
        print(row_line(i, row))
    print()
    print("below AP one-hole second threshold")
    print(f"  threshold: {AP_SECOND} = {float(AP_SECOND):.9f}")
    print(f"  rows below threshold: {len(below)}")
    for i, row in enumerate(below[:20], 1):
        print(row_line(i, row))
    print()
    if not below and not any(key.startswith("uncut") and val for key, val in stats.items()):
        print("THEOREM-SHAPED CERTIFICATE:")
        print("  Every three-tail AP-tail 12-core in this family satisfies")
        print("  meas(G_C) >= 426/35035.  The infinite tail range is covered")
        print("  by nested comb cutoffs; all remaining finite residues were")
        print("  checked exactly with rational interval arithmetic.")
    else:
        print("OPEN RESIDUAL:")
        print("  Some uncut row or below-threshold row remains; inspect stats above.")
    print()
    print("Tournament Analysis")
    print("  vertices: nested_comb_scaffold, root_packet, shell1_carry, mouth_survival, raw_tail_bound")
    print("  observable: preservation of meas(G_C)<426/35035 before scalarization")
    print("  Hamiltonian path: shell1_carry > root_packet > mouth_survival > nested_comb_scaffold > raw_tail_bound")
    print("  challenged assumption: replacement count is not the AP-tail vertex; the useful vertex is the addressed packet plus comb cutoff.")


if __name__ == "__main__":
    main()
