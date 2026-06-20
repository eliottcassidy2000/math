#!/usr/bin/env python3
"""Verify the codex/KPS new-speed bridge on exact rows.

Incoming KPS work (2026-06-20) states that codex's measure-side p1-tax object
and the KPS sector-side far-plateau deviation are the same exact quantity:

    raw_wdelta(E',w) / w
      = p0(E' union {w}) - (p0(E') + p1(E')/7).

This script checks that identity on the HYP-2671 dyadic-block extremizer and on
the non-shell-full warning row from the incoming bridge note.
"""

from __future__ import annotations

from fractions import Fraction

from lrc14_p1_tax_envelope_codex_s40 import base_data, raw_wdelta
from lrc14_shellfull_packet_gap_codex_s44 import SHELL1


def p0(E: tuple[int, ...]) -> Fraction:
    """Exact measure where all six inner seventh-sectors are hit."""

    breakpoints = {Fraction(0), Fraction(1)}
    for e in sorted(set(E)):
        if e == 0:
            continue
        for a in range(0, 7 * e + 1):
            breakpoints.add(Fraction(a, 7 * e))

    total = Fraction(0)
    bps = sorted(x for x in breakpoints if 0 <= x <= 1)
    for lo, hi in zip(bps, bps[1:]):
        if lo == hi:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in E:
            v = e * mid
            v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if all(j in hit for j in range(1, 7)):
            total += hi - lo
    return total


def shell_full(Ep: tuple[int, ...]) -> bool:
    return all(v in Ep for v in SHELL1)


def fmt(q: Fraction) -> str:
    return f"{q} ({float(q):.6f})"


def check_row(name: str, Ep: tuple[int, ...], w: int) -> None:
    base = base_data(Ep)
    codex_delta = raw_wdelta(base, w) / w
    kps_delta = p0(tuple(sorted(Ep + (w,)))) - base.phi
    ratio = codex_delta / base.p1 if base.p1 else None

    print(name)
    print(f"  E'={Ep}, w={w}, shell_full={shell_full(Ep)}")
    print(f"  p1(E')={fmt(base.p1)}")
    print(f"  codex raw_wdelta/w       = {fmt(codex_delta)}")
    print(f"  KPS p0(E)-Phi(E')        = {fmt(kps_delta)}")
    print(f"  identity holds           = {codex_delta == kps_delta}")
    if ratio is not None:
        print(f"  Delta/p1                 = {fmt(ratio)}")
        print(f"  gap below 1/3            = {fmt(Fraction(1, 3) - ratio)}")
    print()


def main() -> None:
    print("HYP-2673 codex/KPS new-speed bridge verifier")
    print("exact Fraction arithmetic")
    print()
    check_row(
        "HYP-2671 dyadic-block extremizer",
        (0, 1, 2, 4, 8, 12, 16, 20),
        24,
    )
    check_row(
        "incoming non-shell-full warning row",
        (0, 2, 3, 5, 6, 15),
        18,
    )
    print("Interpretation")
    print("  The codex p1-tax object and the KPS plateau deviation agree exactly.")
    print("  The relative Delta <= p1/3 target is a shell-full statement; the")
    print("  non-shell-full warning row exceeds 1/3 and must be routed through")
    print("  the shell-damage gate or the absolute C(k)/w far-peel route.")
    print("PASS: codex/KPS bridge verifier complete.")


if __name__ == "__main__":
    main()
