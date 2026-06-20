#!/usr/bin/env python3
"""HYP-2663/T906: AP-tail root-packet scout for LRC14.

This script uses the old FKN/root-packet lesson as a near-AP LRC invariant:
perturbations of the AP collar should be read as packets of addressed roots
(holes plus tail insertions), not as uniform Hamming moves.

It stress-tests the first layer beyond THM-544:

    ({1,...,13} \ holes) union {r,s,t}
    |holes|=4, 14 <= r < s < t <= B.

The engine reuses the exact sorted-arc interval arithmetic from the THM-544
two-replacement certificate and adds two side channels:

* Glaisher odd-shell carry delta relative to the drop-6 collar.
* survival/damage of the four drop-6 mouth intervals.

The output is a packet atlas, not a proof.  It is meant to identify the next
mouth-retention lemma with enough exact data to avoid scalar-only guessing.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce

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
class PacketRow:
    measure: Fraction
    holes: tuple[int, int, int, int]
    tails: tuple[int, int, int]
    old_survivor: Fraction
    components: int
    core: tuple[int, ...]
    carry_delta: tuple[tuple[int, int], ...]
    tail_towers: tuple[tuple[int, int, int], ...]

    @property
    def old_damage(self) -> Fraction:
        return DROP6_COLLAR - self.old_survivor

    @property
    def new_mass(self) -> Fraction:
        return self.measure - self.old_survivor

    @property
    def keeps_drop6_hole(self) -> bool:
        return 6 in self.holes

    @property
    def keeps_full_mouth(self) -> bool:
        return self.old_survivor == DROP6_COLLAR


def primitive(core: tuple[int, ...]) -> bool:
    return reduce(gcd, core, 0) == 1


def base_core(holes: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(v for v in AP if v not in holes)


def tower(speed: int) -> tuple[int, int, int]:
    level = 0
    odd = speed
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


def fmt(q: Fraction) -> str:
    return f"{q} = {float(q):.9f}"


def lower_after_two_tails(meas: Fraction, components: int, s: int) -> Fraction:
    """Lower bound after adding two tails both at least s."""

    return Fraction(5, 7) * meas - Fraction(4 * components, 7 * s)


def lower_after_one_tail(meas: Fraction, components: int, t: int) -> Fraction:
    return Fraction(6, 7) * meas - Fraction(2 * components, 7 * t)


def scan(bound: int, keep: int) -> tuple[list[PacketRow], list[PacketRow], dict[str, int]]:
    stats = defaultdict(int)
    best: list[PacketRow] = []
    below: list[PacketRow] = []

    def make_row(
        m: Fraction,
        comps: tuple[tuple[Fraction, Fraction], ...],
        holes: tuple[int, int, int, int],
        tails: tuple[int, int, int],
        core: tuple[int, ...],
    ) -> PacketRow:
        return PacketRow(
            measure=m,
            holes=holes,
            tails=tails,
            old_survivor=intersect_measure(comps, DROP6_MOUTHS),
            components=len(comps),
            core=core,
            carry_delta=carry_delta(core),
            tail_towers=tuple(tower(v) for v in tails),
        )

    for holes in combinations(AP, 4):
        base = base_core(holes)
        base_comps = safe_components(base)
        base_measure = measure(base_comps)
        base_components = len(base_comps)
        for r in range(14, bound + 1):
            core1 = tuple(sorted(base + (r,)))
            if not primitive(core1):
                stats["nonprimitive_after_r"] += 1
                continue
            comps1 = subtract_arcs(base_comps, danger_arcs(r))
            m1 = measure(comps1)
            c1 = len(comps1)
            stats["r_rows"] += 1
            for s in range(r + 1, bound + 1):
                if lower_after_two_tails(m1, c1, s) >= AP_SECOND:
                    stats["two_tail_prune"] += bound - s + 1
                    break
                core2 = tuple(sorted(core1 + (s,)))
                if not primitive(core2):
                    stats["nonprimitive_after_s"] += 1
                    continue
                comps2 = subtract_arcs(comps1, danger_arcs(s))
                m2 = measure(comps2)
                c2 = len(comps2)
                stats["s_rows"] += 1
                for t in range(s + 1, bound + 1):
                    if lower_after_one_tail(m2, c2, t) >= AP_SECOND:
                        stats["one_tail_prune"] += bound - t + 1
                        break
                    core = tuple(sorted(core2 + (t,)))
                    if not primitive(core):
                        stats["nonprimitive_after_t"] += 1
                        continue
                    comps = subtract_arcs(comps2, danger_arcs(t))
                    m = measure(comps)
                    stats["exact_rows"] += 1
                    enters_best = len(best) < keep or (m, holes, (r, s, t)) < (
                        best[-1].measure,
                        best[-1].holes,
                        best[-1].tails,
                    )
                    if not enters_best and m >= AP_SECOND:
                        continue
                    row = make_row(m, comps, holes, (r, s, t), core)
                    if enters_best:
                        best.append(row)
                        best.sort(key=lambda x: (x.measure, x.holes, x.tails))
                        if len(best) > keep:
                            best.pop()
                    if m < AP_SECOND:
                        below.append(row)

        # Record base rows for sanity; these variables also document the comb
        # scale even when every tail is pruned.
        if base_measure and base_components:
            stats["base_packets"] += 1

    below.sort(key=lambda x: (x.measure, x.holes, x.tails))
    return best, below, dict(stats)


def row_line(rank: int, row: PacketRow) -> str:
    return (
        f"{rank:4d} | {str(row.measure):>13} | {row.holes!s:>15} | {row.tails!s:>14} | "
        f"{str(row.old_survivor):>9} | {str(row.old_damage):>10} | "
        f"{row.components:4d} | {row.carry_delta!s:>28} | {row.tail_towers}"
    )


def print_rows(title: str, rows: list[PacketRow], limit: int) -> None:
    print(title)
    print(
        "rank | measure | holes | tails | old_surv | old_damage | comp | "
        "carry_delta | tail_towers"
    )
    for i, row in enumerate(rows[:limit], 1):
        print(row_line(i, row))
    if len(rows) > limit:
        print(f"  ... {len(rows) - limit} more rows omitted")
    print()


def summarize_below(rows: list[PacketRow]) -> None:
    print("below AP one-hole second threshold")
    print(f"  threshold: {fmt(AP_SECOND)}")
    print(f"  rows below threshold: {len(rows)}")
    print(f"  keep drop-6 hole: {sum(r.keeps_drop6_hole for r in rows)}/{len(rows)}")
    print(f"  full old mouth retained: {sum(r.keeps_full_mouth for r in rows)}/{len(rows)}")
    print(f"  old mouth damaged: {sum(r.old_damage > 0 for r in rows)}/{len(rows)}")
    print()
    by_delta: Counter[tuple[tuple[int, int], ...]] = Counter(r.carry_delta for r in rows)
    print("  carry delta packets among below-threshold rows:")
    for delta, count in by_delta.most_common(12):
        print(f"    {count:4d} x {delta}")
    if len(by_delta) > 12:
        print(f"    ... {len(by_delta) - 12} more carry packets")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=34)
    parser.add_argument("--keep", type=int, default=30)
    parser.add_argument("--below-limit", type=int, default=80)
    args = parser.parse_args()

    print("HYP-2663/T906 LRC14 AP-tail root-packet scout")
    print("family: ({1,...,13} \\ 4 holes) union 3 tails")
    print(f"tail bound: 14 <= r<s<t <= {args.bound}")
    print(f"drop-6 collar: {fmt(DROP6_COLLAR)}")
    print(f"AP one-hole second value: {fmt(AP_SECOND)}")
    print(f"drop-6 odd-shell profile: {dict(sorted(DROP6_PROFILE.items()))}")
    print()

    best, below, stats = scan(args.bound, args.keep)
    print("scan stats")
    for key in sorted(stats):
        print(f"  {key}: {stats[key]}")
    print()
    print_rows("best exact rows in scanned packet", best, args.keep)
    summarize_below(below)
    print_rows("below-threshold packet rows", below, args.below_limit)

    print("invariant reading")
    print("1. This is the AP-tail analogue of FKN/root-packet structure: the")
    print("   relevant perturbation is the addressed packet (holes, tails, carry,")
    print("   mouth survival), not the number of replacements.")
    print("2. The proof target is a mouth-damage theorem: once a packet damages the")
    print("   four drop-6 mouths or opens genuinely new odd-shell carry, it should")
    print("   pay the AP-second threshold before scalar measure is used.")
    print("3. Any below-threshold packet retained here should be treated as a finite")
    print("   root-packet template feeding HYP-2654/HYP-2659/HYP-2660, while the")
    print("   pruned tail side routes to HYP-2655/HYP-2658.")
    print()
    print("Tournament Analysis")
    print("  vertices: root_packet, mouth_survival, odd_shell_carry, comb_tail_bound, raw_speed")
    print("  observable: preservation of the predicate meas(G_C)<426/35035")
    print("  Hamiltonian path: root_packet > mouth_survival > odd_shell_carry > comb_tail_bound > raw_speed")
    print("  challenged assumption: AP-tail vertices are not replacement counts; they are addressed root packets.")


if __name__ == "__main__":
    main()
