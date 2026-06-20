#!/usr/bin/env python3
"""HYP-2664/T907: shell-1 gate for the three-tail AP-tail frontier.

HYP-2663 found no below-second rows for the three-tail AP-tail family through
tail bound 35, but the naive unbounded comb proof leaves a large finite residue
set.  HYP-2661 says any sub-second AP-tail row must preserve the dyadic shell-1
tower {1,2,4,8}.  This script measures how that new gate acts on the three-tail
comb frontier before expensive exact tail enumeration.

It computes the first-tail comb cutoff for every 4-hole AP base:

    C = ({1,...,13} \\ H) union {r,s,t}, |H|=4.

If all three tails are at least R, then

    meas(G_C) >= (4/7) meas(G_base) - 6*c_base/(7R).

Thus R is the first crude "residue burden" of the three-tail theorem.  The
point is not that R alone proves the family, but that HYP-2661 removes the worst
root packets before this residue search begins.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations

from lrc14_three_tail_ap_tail_comb_codex_s39 import (
    AP,
    AP_SECOND,
    SHELL1,
    cutoff_for,
    safe_components,
    measure,
)


@dataclass(frozen=True)
class BaseFrontier:
    r_cut: int
    holes: tuple[int, int, int, int]
    measure: Fraction
    components: int

    @property
    def shell1_damaged(self) -> bool:
        return bool(set(self.holes) & SHELL1)

    @property
    def missing_shell1(self) -> tuple[int, ...]:
        return tuple(sorted(set(self.holes) & SHELL1))

    @property
    def residue_load_upper(self) -> int:
        """Crude count of tail triples below the first cutoff."""

        n = max(0, self.r_cut - 14)
        return n * (n - 1) * (n - 2) // 6


def quantiles(values: list[int]) -> dict[str, int]:
    values = sorted(values)
    if not values:
        return {}
    def q(p: float) -> int:
        return values[min(len(values) - 1, int(round(p * (len(values) - 1))))]
    return {
        "min": values[0],
        "q25": q(0.25),
        "median": q(0.50),
        "q75": q(0.75),
        "q90": q(0.90),
        "q99": q(0.99),
        "max": values[-1],
    }


def build_frontier() -> list[BaseFrontier]:
    out: list[BaseFrontier] = []
    for holes in combinations(AP, 4):
        base = tuple(v for v in AP if v not in holes)
        comps = safe_components(base)
        meas = measure(comps)
        r_cut = cutoff_for(meas, len(comps), 3)
        if r_cut is None:
            raise RuntimeError(f"uncut base {holes}")
        out.append(BaseFrontier(r_cut, holes, meas, len(comps)))
    out.sort(key=lambda row: (row.r_cut, row.residue_load_upper), reverse=True)
    return out


def row_line(rank: int, row: BaseFrontier) -> str:
    gate = "damaged" if row.shell1_damaged else "full"
    return (
        f"{rank:3d} | {row.r_cut:5d} | {row.residue_load_upper:10d} | "
        f"{row.holes!s:>15} | {gate:7s} | {row.missing_shell1!s:>9} | "
        f"{str(row.measure):>13} | {row.components:3d}"
    )


def main() -> None:
    rows = build_frontier()
    damaged = [r for r in rows if r.shell1_damaged]
    full = [r for r in rows if not r.shell1_damaged]
    top40 = rows[:40]
    missing_counter = Counter(r.missing_shell1 for r in rows)

    print("HYP-2664/T907 three-tail AP-tail shell-1 frontier")
    print("family: ({1,...,13} \\ H) union {r,s,t}, |H|=4")
    print(f"AP one-hole second value: {AP_SECOND} = {float(AP_SECOND):.9f}")
    print("first-tail cutoff R: all tails >=R are certified by the crude 3-tail comb bound")
    print()
    print("frontier counts")
    print(f"  total 4-hole bases: {len(rows)}")
    print(f"  shell-1 damaged bases: {len(damaged)}")
    print(f"  shell-1 full bases: {len(full)}")
    print(f"  top-40 damaged by HYP-2661 gate: {sum(r.shell1_damaged for r in top40)}/40")
    print(f"  top-100 damaged by HYP-2661 gate: {sum(r.shell1_damaged for r in rows[:100])}/100")
    print()
    print("r_cut quantiles")
    print(f"  all:      {quantiles([r.r_cut for r in rows])}")
    print(f"  damaged:  {quantiles([r.r_cut for r in damaged])}")
    print(f"  shell1+:  {quantiles([r.r_cut for r in full])}")
    print()
    print("missing shell-1 packet counts")
    for miss, count in missing_counter.most_common():
        label = "full" if not miss else str(miss)
        print(f"  {label:12s}: {count}")
    print()
    print("top crude comb bases")
    print("rank | r_cut | load_upper | holes | shell1 | missing | base_measure | comp")
    for i, row in enumerate(rows[:30], 1):
        print(row_line(i, row))
    print()
    print("top shell-1-full bases after applying HYP-2661")
    print("rank | r_cut | load_upper | holes | shell1 | missing | base_measure | comp")
    for i, row in enumerate(full[:30], 1):
        print(row_line(i, row))
    print()
    print("invariant reading")
    print("1. The naive three-tail comb frontier is dominated by packets that delete")
    print("   a shell-1 tower bit. HYP-2661 should remove those before exact tail")
    print("   enumeration; the broad comb was doing avoidable work.")
    print("2. After the shell-1 gate, the worst first-tail cutoff drops from the")
    print("   global frontier to the shell-1-full frontier. The remaining packets")
    print("   should be attacked by mouth survival / endpoint-owner addresses, not")
    print("   by replacement count alone.")
    print("3. This gives a proof-order recipe: tower conservation first, root packet")
    print("   second, nested comb finite residues third.")
    print()
    print("Tournament Analysis")
    print("  vertices: shell1_gate, root_packet, mouth_owner, nested_comb, raw_replacement_count")
    print("  observable: finite residue burden for proving meas(G_C)>=426/35035")
    print("  Hamiltonian path: shell1_gate > root_packet > mouth_owner > nested_comb > raw_replacement_count")
    print("  challenged assumption: the first AP-tail quotient is not the number of tails; it is whether the shell-1 carrier survives.")


if __name__ == "__main__":
    main()
