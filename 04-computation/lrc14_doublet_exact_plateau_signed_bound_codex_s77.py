#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Actual-size doublet signed-tail audit against the corrected frozen plateau.

Incoming HYP-2797 identifies the genuine-wide doublet family

    E_N(M) = {0,1,...,N-3} union {M,M+1}.

The companion signed-bound script uses a shifted label convention in one table.
This follow-up keeps the actual row size ``N`` fixed, imports the corrected
two-far frozen plateau from HYP-2799, and measures the exact rational signed
tail

    M * (p0(E_N(M)) - D7({0,...,N-3}, 1)).

The point is not to prove the asymptotic bound by computation.  It gives the
right constants and finite-window targets for a proof: once the displayed
positive signed constant is made rigorous for all M>=15, the actual-size
doublet tail is closed against cap_N, and often against Q(N-1), with only a
small finite prefix left.
"""
from __future__ import annotations

import argparse
import functools
import sys
from fractions import Fraction as F

sys.path.insert(0, "04-computation")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
print = functools.partial(print, flush=True)

from lrc14_genuinewide_corrected_freeze_tail_codex_s77 import (  # noqa: E402
    CAP,
    QVAL,
    dblock7_exact,
    p0_fast,
)


def actual_doublet(size: int, M: int) -> tuple[int, ...]:
    return tuple(list(range(size - 2)) + [M, M + 1])


def strict_tail_cutoff(constant: F, room: F) -> int | None:
    """Smallest integer M0 such that constant/M < room for all M>=M0."""
    if room <= 0:
        return None
    if constant <= 0:
        return 1
    return constant // room + 1


def scan_size(size: int, lo: int, hi: int) -> dict[str, object]:
    base = tuple(range(size - 2))
    plateau = dblock7_exact(base, 1)
    cap_room = CAP[size] - plateau
    q_room = QVAL[size] - plateau

    best_p0_M = lo
    best_p0 = F(-1)
    max_pos_M = lo
    max_pos = F(0)
    min_scaled_M = lo
    min_scaled = F(0)
    max_abs_M = lo
    max_abs = F(0)

    for M in range(lo, hi + 1):
        val = p0_fast(actual_doublet(size, M))
        scaled = M * (val - plateau)
        if val > best_p0:
            best_p0 = val
            best_p0_M = M
        if scaled > max_pos:
            max_pos = scaled
            max_pos_M = M
        if scaled < min_scaled:
            min_scaled = scaled
            min_scaled_M = M
        if abs(scaled) > max_abs:
            max_abs = abs(scaled)
            max_abs_M = M

    return {
        "size": size,
        "base": base,
        "plateau": plateau,
        "cap_room": cap_room,
        "q_room": q_room,
        "best_p0_M": best_p0_M,
        "best_p0": best_p0,
        "cap_gap_at_best": CAP[size] - best_p0,
        "q_gap_at_best": QVAL[size] - best_p0,
        "max_pos_M": max_pos_M,
        "max_pos": max_pos,
        "min_scaled_M": min_scaled_M,
        "min_scaled": min_scaled,
        "max_abs_M": max_abs_M,
        "max_abs": max_abs,
        "cap_tail_from": strict_tail_cutoff(max_pos, cap_room),
        "q_tail_from": strict_tail_cutoff(max_pos, q_room),
        "cap_closes_from_lo": max_pos < lo * cap_room,
        "q_closes_from_lo": max_pos < lo * q_room,
    }


def fmt_fraction(x: object) -> str:
    if isinstance(x, F):
        return str(x)
    return str(x)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lo", type=int, default=15)
    parser.add_argument("--hi", type=int, default=600)
    args = parser.parse_args()

    print("=" * 96)
    print("LRC14 actual-size doublet exact plateau signed-tail audit (codex-S77)")
    print("=" * 96)
    print(f"Family: E_N(M)={{0,...,N-3}} union {{M,M+1}}, scanned {args.lo}<=M<={args.hi}.")
    print("Plateau: corrected HYP-2799 D7({0,...,N-3}, gap=1).")
    print()
    print("Exact scan maxima:")
    print(
        f"{'N':>3} {'D7':>18} {'M_p0':>5} {'max p0':>18} "
        f"{'cap_N-p0':>14} {'Q(N-1)-p0':>14}"
    )
    rows = [scan_size(size, args.lo, args.hi) for size in range(8, 13)]
    for row in rows:
        print(
            f"{row['size']:>3} {fmt_fraction(row['plateau']):>18} "
            f"{row['best_p0_M']:>5} {fmt_fraction(row['best_p0']):>18} "
            f"{fmt_fraction(row['cap_gap_at_best']):>14} "
            f"{fmt_fraction(row['q_gap_at_best']):>14}"
        )
    print()
    print("Signed-tail constants C+ = max M*(p0-D7), Cabs = max |M*(p0-D7)|:")
    print(
        f"{'N':>3} {'M_C+':>5} {'C+':>18} {'M_min':>5} {'min scaled':>18} "
        f"{'M_abs':>5} {'Cabs':>18}"
    )
    for row in rows:
        print(
            f"{row['size']:>3} {row['max_pos_M']:>5} {fmt_fraction(row['max_pos']):>18} "
            f"{row['min_scaled_M']:>5} {fmt_fraction(row['min_scaled']):>18} "
            f"{row['max_abs_M']:>5} {fmt_fraction(row['max_abs']):>18}"
        )
    print()
    print("Tail closure if the scanned C+ is proved for all M>=lo:")
    print(
        f"{'N':>3} {'cap_N-D7':>18} {'Q(N-1)-D7':>18} "
        f"{'cap from M>=':>12} {'Q from M>=':>11} {'cap@lo':>7} {'Q@lo':>5}"
    )
    for row in rows:
        print(
            f"{row['size']:>3} {fmt_fraction(row['cap_room']):>18} "
            f"{fmt_fraction(row['q_room']):>18} "
            f"{str(row['cap_tail_from']):>12} {str(row['q_tail_from']):>11} "
            f"{str(row['cap_closes_from_lo']):>7} {str(row['q_closes_from_lo']):>5}"
        )
    print()
    print("Interpretation:")
    print("  * The actual-size cap_N margins are wide; C+ closes cap_N from M=15 for N=8..12.")
    print("  * Q(N-1) is much tighter.  N=10 has the live razor margin, but the signed tail")
    print("    still suggests only a tiny finite prefix once a uniform C+ bound is proved.")
    print("  * This is a cap-branch tool for genuine-wide LRC14; it does not replace the")
    print("    corrected frozen-gap extremality and finite low-f obligations in HYP-2799.")


if __name__ == "__main__":
    main()
