#!/usr/bin/env python3
"""AP-prefix insertion template for the HYP-2691 transfer DP.

The residual-budget audit showed the highest pressure in the named row bank at
AP-like prefixes.  This script isolates the exact template

    d*{0,1,...,m-1}  ->  add d*m

and records the one-missed-sector transfer residual and component-budget
pressure.  Integer dilation preserves the measured transition data, while the
run count and inserted speed both scale by d, so the pressure is unchanged.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd

from lrc14_state_transfer_dp_codex_s58 import state_profile, transfer_step
from lrc14_transfer_residual_budget_codex_s58 import one_missed_run_count
from lrc14_wide_branch_ridge_codex_s47 import Row


def template_row(m: int, d: int) -> Row:
    return tuple(d * i for i in range(m))


def template_stats(m: int, d: int) -> tuple[Fraction, Fraction, Fraction, int, Fraction, Fraction]:
    prefix = template_row(m, d)
    added = d * m
    step = transfer_step(prefix, added)
    p1 = state_profile(prefix).p_by_count[1]
    decor = p1 / 7
    residual = step.delta_p0 - decor
    V = one_missed_run_count(prefix)
    bound = Fraction(6 * V, 49 * added) if added else Fraction(0)
    pressure = abs(residual) / bound if bound else Fraction(0)
    return p1, step.delta_p0, residual, V, bound, pressure


def main() -> None:
    print("HYP-2691 AP-prefix transfer template")
    print("template: d*{0,...,m-1} -> add d*m")
    print("pressure = |delta-p1/7| / ((6/49)*V/(d*m))")
    print()

    print("primitive scale d=1")
    print(" m p1 delta residual V bound pressure scaled_abs")
    leaders: list[tuple[Fraction, int, Fraction]] = []
    for m in range(1, 14):
        p1, delta, residual, V, bound, pressure = template_stats(m, 1)
        leaders.append((pressure, m, abs(residual) * m))
        print(
            f"{m:2d} {p1} {delta} {residual} {V:2d} {bound} "
            f"{pressure} {abs(residual) * m}"
        )

    print()
    print("pressure leaders")
    for pressure, m, scaled in sorted(leaders, reverse=True)[:6]:
        print(f"  m={m:2d} pressure={pressure} ({float(pressure):.9f}) scaled_abs={scaled}")

    print()
    print("dilation check for the live peak m=6")
    base = template_stats(6, 1)
    for d in range(1, 15):
        p1, delta, residual, V, bound, pressure = template_stats(6, d)
        relation = "same" if (p1, delta, residual, pressure) == (base[0], base[1], base[2], base[5]) else "changed"
        print(
            f"  d={d:2d} gcd(d,7)={gcd(d,7)} p1={p1} delta={delta} residual={residual} "
            f"V={V:3d} bound={bound} pressure={pressure} {relation}"
        )

    print()
    print("SYNTHESIS")
    print("  1. Among AP-prefix appends through m=13, the exact pressure peak is m=6:")
    print("     residual=27/245, V=13, bound=13/49, pressure=27/65.")
    print("  2. Integer dilation preserves p1, delta, residual, and pressure,")
    print("     while V and the added speed both scale by d.")
    print("  3. This explains why AP9, boundary doubled-AP, and related row-bank steps")
    print("     share the same top pressure.  It is a finite AP-prefix template, not a")
    print("     true-wide obstruction.")


if __name__ == "__main__":
    main()
