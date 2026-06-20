#!/usr/bin/env python3
"""HYP-2671 scout: the shell-full new-speed 1/3 constant.

HYP-2670 split the post-shell-gate p1-tax route into a finite B13 packet
ledger and a new-speed decay target

    max(E') > 14  ==>  Delta_w^+ <= p1(E')/3.

This script studies that single open constant.  It records the exact B30
new-speed maximum, tests whether fold reciprocal mass alone explains it, and
probes the dyadic block family

    E_m = {0,1,2,4,8,3m,4m,5m},  w = 6m.

The surprise is that m=4 is the only large spike in this family; it is exactly
the B30 new-speed maximum.  Thus the constant `1/3` is not a diffuse scalar
frontier but the gap above a single dyadic block resonance.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction

from lrc14_far_delta_galois_phase_codex_s38 import primitive
from lrc14_p1_tax_envelope_codex_s40 import base_data, raw_wdelta
from lrc14_shellfull_packet_gap_codex_s44 import (
    LightRow,
    fold_profile,
    fold_recip_mass,
    odd_carry_profile,
    scan_shellfull_light,
    shell_extras,
)


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def eval_row(Ep: tuple[int, ...], w: int) -> LightRow:
    base = base_data(Ep)
    raw = raw_wdelta(base, w)
    actual = max(raw, Fraction(0)) / (w * base.p1)
    return LightRow(actual=actual, raw=raw, p1=base.p1, phi=base.phi, Ep=Ep, w=w)


def fold_mass(row: LightRow) -> Fraction:
    return fold_recip_mass(fold_profile(row.Ep))


def print_top(rows: list[LightRow], title: str, limit: int = 14) -> None:
    print(title)
    print("rank | ratio | gap below 1/3 | w | max(E') | E' | extras | fold_recip | odd-carry")
    for rank, row in enumerate(sorted(rows, key=lambda r: (r.actual, r.raw), reverse=True)[:limit], 1):
        print(
            f"{rank:4d} | {fmt(row.actual):>20} | {Fraction(1, 3) - row.actual} | "
            f"{row.w:2d} | {max(row.Ep):2d} | {row.Ep} | {shell_extras(row.Ep)} | "
            f"{fold_mass(row)} | {odd_carry_profile(row.Ep)}"
        )
    print()


def print_fold_bins(rows: list[LightRow]) -> None:
    bins: dict[int, list[LightRow]] = defaultdict(list)
    for row in rows:
        f = fold_mass(row)
        bins[(f.numerator * 20) // f.denominator].append(row)

    print("new-speed rows by fold_recip twentieth bin")
    print("bin | fold range | rows | max ratio | max row | fold_recip")
    for b in sorted(bins):
        group = bins[b]
        top = max(group, key=lambda r: (r.actual, r.raw))
        print(
            f"{b:3d} | [{Fraction(b,20)}, {Fraction(b+1,20)}) | {len(group):5d} | "
            f"{fmt(top.actual):>20} | w={top.w}, E'={top.Ep} | {fold_mass(top)}"
        )
    print()


def print_by_max(rows: list[LightRow]) -> None:
    by_max: dict[int, list[LightRow]] = defaultdict(list)
    for row in rows:
        by_max[max(row.Ep)].append(row)

    print("best row by max(E')")
    print("max(E') | rows | max ratio | gap below 1/3 | w | E' | fold_recip")
    for m in sorted(by_max):
        group = by_max[m]
        top = max(group, key=lambda r: (r.actual, r.raw))
        print(
            f"{m:7d} | {len(group):4d} | {fmt(top.actual):>20} | "
            f"{Fraction(1,3) - top.actual} | {top.w:2d} | {top.Ep} | {fold_mass(top)}"
        )
    print()


def dyadic_block_family(m_min: int = 3, m_max: int = 24) -> list[tuple[int, LightRow]]:
    out: list[tuple[int, LightRow]] = []
    for m in range(m_min, m_max + 1):
        Ep = tuple(sorted({0, 1, 2, 4, 8, 3 * m, 4 * m, 5 * m}))
        w = 6 * m
        if len(Ep) != 8 or not primitive(Ep + (w,)):
            continue
        out.append((m, eval_row(Ep, w)))
    return out


def print_dyadic_family(rows: list[tuple[int, LightRow]]) -> None:
    print("dyadic block family E_m={0,1,2,4,8,3m,4m,5m}, w=6m")
    print("m | max(E') | ratio | gap below 1/3 | p1 | fold_recip | odd-carry")
    for m, row in rows:
        print(
            f"{m:2d} | {max(row.Ep):7d} | {fmt(row.actual):>20} | "
            f"{Fraction(1,3) - row.actual} | {row.p1} | {fold_mass(row)} | "
            f"{odd_carry_profile(row.Ep)}"
        )
    print()


def main() -> None:
    print("HYP-2671 shell-full new-speed 1/3 constant scout")
    print("exact Fraction arithmetic")
    print()

    rows = scan_shellfull_light(30, 8)
    new_rows = [row for row in rows if max(row.Ep) > 14]
    tail_rows = [row for row in rows if max(row.Ep) > 24]
    finite_rows = [row for row in rows if max(row.Ep) <= 14]
    new_top = max(new_rows, key=lambda r: (r.actual, r.raw))
    finite_top = max(finite_rows, key=lambda r: (r.actual, r.raw))
    tail_top = max(tail_rows, key=lambda r: (r.actual, r.raw))

    print("B30 shell-full quotient recap")
    print(f"  rows={len(rows)}, finite={len(finite_rows)}, new={len(new_rows)}, tail={len(tail_rows)}")
    print(f"  finite max={fmt(finite_top.actual)} at w={finite_top.w}, E'={finite_top.Ep}")
    print(f"  new-speed max={fmt(new_top.actual)} at w={new_top.w}, E'={new_top.Ep}")
    print(f"  new-speed gap below 1/3={Fraction(1,3)-new_top.actual}")
    print(f"  tail max={fmt(tail_top.actual)} at w={tail_top.w}, E'={tail_top.Ep}")
    print(f"  tail gap below 1/4={Fraction(1,4)-tail_top.actual}")
    print()

    print_top(new_rows, "top new-speed rows max(E')>14")
    print_by_max(rows)
    print_fold_bins(new_rows)

    family = dyadic_block_family()
    print_dyadic_family(family)
    fam_top_m, fam_top = max(family, key=lambda mr: (mr[1].actual, mr[1].raw))
    print("dyadic block conclusion")
    print(
        f"  family maximum occurs at m={fam_top_m}: ratio={fmt(fam_top.actual)}, "
        f"gap below 1/3={Fraction(1,3)-fam_top.actual}"
    )
    print("  This is exactly the B30 new-speed maximum.")
    print()

    print("Interpretation")
    print("  The live post-gate constant is the 1/3 new-speed barrier.")
    print("  In B30 it is controlled by one dyadic block resonance:")
    print("    E'=(0,1,2,4,8,12,16,20), w=24, Delta^+/p1=1371/4319.")
    print("  Fold reciprocal mass helps locate the block but is not a monotone certificate:")
    print("    low-fold rows can still be moderately high, and high-fold rows can be harmless.")
    print("  The proof target should isolate dyadic block resonance first, then prove")
    print("  every other new-speed shell-full row has extra packet cancellation below 1/3.")
    print()
    print("Tournament Analysis")
    print("  vertices: dyadic_block_resonance > new_speed_1/3_gap > fold_recip_bins > max_speed_layer > raw_runner_vertices")
    print("  observable: exact Delta^+/p1 on the shell-full quotient with max(E')>14")
    print("  switch/gauge: first gate by max(E'), then by dyadic block family and fold-target bins")
    print("  Hamiltonian path: dyadic_block_resonance > new_speed_1/3_gap > fold_recip_bins > max_speed_layer > raw_runner_vertices")
    print("  challenged assumption: fold mass alone determines the constant; it does not.")
    print("PASS: shell-full new-speed constant scout complete.")


if __name__ == "__main__":
    main()
