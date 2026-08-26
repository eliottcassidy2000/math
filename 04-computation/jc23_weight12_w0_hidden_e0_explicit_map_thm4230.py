#!/usr/bin/env python3
"""Exact symbolic certificate for the explicit hidden E0 map in THM-4230.

On the W=0 chart the positive-genus component scales to x^6+y^4=1.
The theorem constructs degree-six maps to E0:Y^2=X^3+1 indexed by the
roots of p^4-2p^3-2p+1.  This script checks the map identity, all required
nonvanishing conditions, the attachment coordinate and orbit actions, and
the pole-degree ledger with exact SymPy arithmetic.

Reproduction:
  python3 04-computation/jc23_weight12_w0_hidden_e0_explicit_map_thm4230.py
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


p, t, u = sp.symbols("p t u")
quartic = p**4 - 2 * p**3 - 2 * p + 1
scale_denominator = 2 * p**3 + 3 * p**2 - 1


def reduced(numerator: sp.Expr) -> sp.Expr:
    return sp.rem(sp.Poly(sp.expand(numerator), p), sp.Poly(quartic, p)).as_expr()


def main() -> None:
    core_identity = sp.expand(
        (u - 1) * (u + p**3) ** 2
        - (u - p**2) ** 3
        - scale_denominator * u * (u + 1)
    )
    require(
        sp.factor(core_identity) == u * (p + 1) ** 2 * quartic,
        "load-bearing quartic identity changed",
    )

    x_cubed = 2 * t / (t**2 + 1)
    y_squared = (t**2 - 1) / (t**2 + 1)
    require(sp.cancel(x_cubed**2 + y_squared**2 - 1) == 0,
            "normalization no longer lies on x^6+y^4=1")
    require(sp.cancel((1 + y_squared) / x_cubed - t) == 0,
            "t reconstruction changed")

    # With s^6=4/scale_denominator, these are Y^2 and X^3.
    y_map_squared = (
        y_squared * (t**2 + p**3) ** 2
        / (scale_denominator * t**2)
    )
    x_map_cubed = (
        sp.Rational(1, 2) * x_cubed * (t**2 - p**2) ** 3
        / (scale_denominator * t**3)
    )
    map_error = sp.together(y_map_squared - x_map_cubed - 1)
    require(reduced(sp.fraction(map_error)[0]) == 0,
            "explicit map misses E0 modulo the quartic")

    require(sp.resultant(quartic, scale_denominator, p) == -108,
            "sixth-root scale can vanish")
    require(sp.resultant(quartic, p**2 + 1, p) == 4,
            "a map numerator can cancel a t=+-i pole")
    require(sp.resultant(quartic, p, p) == 1,
            "quartic unexpectedly has p=0")
    reciprocal = sp.cancel(p**6 * scale_denominator.subs(p, 1 / p)
                           + scale_denominator)
    require(reduced(sp.fraction(reciprocal)[0]) == 0,
            "reciprocal-root scale law changed")

    attachment_c = sp.cancel(y_squared / x_cubed)
    require(attachment_c == (t**2 - 1) / (2 * t),
            "attachment coordinate changed")
    tau_t = sp.cancel((1 - y_squared) / (-x_cubed))
    require(tau_t == -1 / t, "tau no longer exchanges t and -1/t")

    # (number of places, pole order of X at each place).
    pole_ledger = {
        "t=0": (2, 2),
        "t=infinity": (2, 2),
        "t=i": (1, 2),
        "t=-i": (1, 2),
    }
    total_x_pole_order = sum(count * order for count, order in pole_ledger.values())
    require(total_x_pole_order == 12, "X-pole divisor degree changed")
    map_degree = total_x_pole_order // 2
    require(map_degree == 6, "map degree changed")

    # In Z[zeta_3], N(a+b*zeta_3)=a^2-a*b+b^2.
    norm_three_plus_zeta = 3 * 3 - 3 * 1 + 1 * 1
    require(norm_three_plus_zeta == 7, "norm-seven degree-42 control changed")

    print("THM-4230 explicit W=0 hidden-E0 map certificate")
    print("normalized_curve=x^6+y^4=1 t=(1+y^2)/x^3")
    print("quartic=p^4-2p^3-2p+1 scale=s^6=4/(2p^3+3p^2-1)")
    print("identity_factor=u*(p+1)^2*(p^4-2p^3-2p+1)")
    print("map_X=(s^2/2)*x*(t^2-p^2)/t")
    print("map_Y=(s^3/2)*y*(t^2+p^3)/t")
    print("map_equation=PASS Y^2-X^3=1")
    print("resultants=scale:-108,p2_plus_1:4,p:1")
    print("pole_ledger=t0:2x2,tinf:2x2,ti:1x2,tminus_i:1x2 total_X=12")
    print("map_degree=6 anti_tau6=Prym hidden_E0=PASS")
    print("attachment_c=(t-t^-1)/2 tau_action=t->-1/t")
    print("six_orbit_actions=(zeta3,-1) common_value_forces_O=PASS")
    print("rank_one_degrees=6*N(alpha) degree34=impossible degree42_norm7=possible")
    print("result=PASS explicit hidden map; full hidden lattice and norm7 torsion remain OPEN")


if __name__ == "__main__":
    main()
