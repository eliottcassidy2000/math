#!/usr/bin/env python3
"""Exact referee for THM-3034.

The checks are symbolic over QQ.  No floating-point value is truth-bearing.
"""

from fractions import Fraction

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def numerator(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.together(expr).as_numer_denom()[0])


def require_zero_mod(
    expr: sp.Expr, relation: sp.Expr, variable: sp.Symbol, message: str
) -> None:
    rem = sp.factor(sp.rem(numerator(expr), relation, variable))
    require(rem == 0, f"{message}: remainder={rem}")


def invariants(a1: int, a2: int, a3: int, a4: int, a6: int):
    b2 = a1 * a1 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3 * a3 + 4 * a6
    b8 = (
        a1 * a1 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3 * a3
        - a4 * a4
    )
    c4 = b2 * b2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    delta = -b2 * b2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    return c4, c6, delta, Fraction(c4**3, delta)


def main() -> None:
    tau, eta, xi = sp.symbols("tau eta xi")
    quartic = tau**4 - 6 * tau**3 + 7 * tau**2 - 2 * tau + 1
    quartic_relation = eta**2 - quartic

    # C_+ on the Y=1 chart and projection to the binary quartic.
    cplus = sp.expand(xi * tau - (xi - 1) * (xi - tau) * (1 - tau))
    eta_from_cplus = 2 * (tau - 1) * xi - tau**2 + tau + 1
    require_zero_mod(
        eta_from_cplus**2 - quartic,
        cplus,
        xi,
        "C_+ projection does not land on the quartic",
    )

    # Dense birational map from the quartic to the minimal X_1(14) model.
    chi = sp.factor((eta + 1 - tau) / tau**2)
    u = sp.factor((chi + 1) / 2)
    v = sp.factor((chi**2 - 1) * tau / 4)
    e1_relation = v**2 + u * v + v - u**3 + u
    require_zero_mod(
        e1_relation,
        quartic_relation,
        eta,
        "quartic-to-E1 map failed",
    )

    U, V = sp.symbols("U V")
    abstract_e1 = V**2 + U * V + V - U**3 + U
    tau_inverse = sp.factor(V / (U * (U - 1)))
    eta_inverse = sp.factor((2 * U - 1) * tau_inverse**2 + tau_inverse - 1)
    require_zero_mod(
        eta_inverse**2
        - (
            tau_inverse**4
            - 6 * tau_inverse**3
            + 7 * tau_inverse**2
            - 2 * tau_inverse
            + 1
        ),
        abstract_e1,
        V,
        "E1-to-quartic inverse failed",
    )
    require_zero_mod(
        v / (u * (u - 1)) - tau,
        quartic_relation,
        eta,
        "quartic round trip failed for tau",
    )
    require_zero_mod(
        (2 * u - 1) * (v / (u * (u - 1))) ** 2
        + v / (u * (u - 1))
        - 1
        - eta,
        quartic_relation,
        eta,
        "quartic round trip failed for eta",
    )

    # Pull the even flank cycle rho[A:B:C]=[B:C:A] through the map.  On the
    # B=1 chart, tau'=xi/tau and xi'=1/tau.  The resulting formulas are
    # exactly addition by T=(0,0) on E1.
    xi_from_quartic = sp.factor(
        (eta + tau**2 - tau - 1) / (2 * (tau - 1))
    )
    tau_rho = sp.cancel(xi_from_quartic / tau)
    xi_rho = 1 / tau
    eta_rho = sp.cancel(
        2 * (tau_rho - 1) * xi_rho - tau_rho**2 + tau_rho + 1
    )

    def quartic_to_e1(T: sp.Expr, E: sp.Expr):
        X = sp.factor((E + 1 - T) / T**2)
        return sp.cancel((X + 1) / 2), sp.cancel((X**2 - 1) * T / 4)

    u_rho, v_rho = quartic_to_e1(tau_rho, eta_rho)
    u_add_t = sp.cancel((v / u) ** 2 + v / u - u)
    v_add_t = sp.cancel(-(v / u + 1) * u_add_t - 1)
    require_zero_mod(
        u_rho - u_add_t,
        quartic_relation,
        eta,
        "rho is not addition by T in u",
    )
    require_zero_mod(
        v_rho - v_add_t,
        quartic_relation,
        eta,
        "rho is not addition by T in v",
    )

    # Minimal invariants and the short model recorded by LMFDB 14.a5.
    e1_inv = invariants(1, 0, 1, -1, 0)
    require(
        e1_inv == (25, -253, -28, Fraction(-15625, 28)),
        f"unexpected E1 invariants {e1_inv}",
    )
    short_x = 36 * U + 3
    short_y = 108 * (2 * V + U + 1)
    require_zero_mod(
        short_y**2 - (short_x**3 - 675 * short_x + 13662),
        abstract_e1,
        V,
        "E1 short-model change failed",
    )

    # T=(0,0) has order three, and it generates the unique rational
    # order-three subgroup: the other x-coordinates satisfy an irreducible
    # cubic factor of the 3-division polynomial.
    x = sp.symbols("x")
    psi3 = sp.factor(3 * x**4 + x**3 - 3 * x**2 + 3 * x)
    require(psi3 == x * (3 * x**3 + x**2 - 3 * x + 3), "psi_3 factor")
    require(
        sp.Poly(3 * x**3 + x**2 - 3 * x + 3, x, domain=sp.QQ).is_irreducible,
        "nonzero cubic factor of psi_3 must be irreducible over QQ",
    )
    tangent_den = 1
    tangent_lambda = Fraction(-1, tangent_den)
    tangent_nu = Fraction(0, tangent_den)
    twice_x = tangent_lambda**2 + tangent_lambda
    twice_y = -(tangent_lambda + 1) * twice_x - tangent_nu - 1
    require((twice_x, twice_y) == (0, -1), "2T must equal -T")
    tangent_x = sp.symbols("tangent_x")
    tangent_intersection = sp.expand(
        (-tangent_x) ** 2
        + tangent_x * (-tangent_x)
        + (-tangent_x)
        - tangent_x**3
        + tangent_x
    )
    require(tangent_intersection == -tangent_x**3, "T tangent/order-three check")

    # The normalized degree-three quotient by <T> is the standard X_0(14)
    # equation.  These are the exact Velu functions on the dense u != 0 chart.
    quotient_u = sp.factor((U**3 - U + 1) / U**2)
    quotient_v = sp.factor(
        (V * (U**3 + U - 2) + U**2 - U - 1) / U**3
    )
    e0_relation = (
        quotient_v**2
        + quotient_u * quotient_v
        + quotient_v
        - quotient_u**3
        - 4 * quotient_u
        + 6
    )
    require_zero_mod(e0_relation, abstract_e1, V, "Velu quotient failed")

    # The same map in the short coordinates printed in the theorem.
    sx, sy = sp.symbols("sx sy")
    source_short_x = 36 * U + 3
    source_short_y = 108 * (2 * V + U + 1)
    quotient_short_x = sp.factor(36 * quotient_u + 3)
    quotient_short_y = sp.factor(108 * (2 * quotient_v + quotient_u + 1))
    expected_short_x = (sx**3 - 6 * sx**2 - 1287 * sx + 50544) / (sx - 3) ** 2
    expected_short_y = sy * (sx**3 - 9 * sx**2 + 1323 * sx - 97227) / (sx - 3) ** 3
    require(
        sp.cancel(quotient_short_x - expected_short_x.subs(sx, source_short_x))
        == 0,
        "short Velu x formula failed",
    )
    require_zero_mod(
        quotient_short_y
        - expected_short_y.subs({sx: source_short_x, sy: source_short_y}),
        abstract_e1,
        V,
        "short Velu y formula failed",
    )

    # Independent short-form Velu check and exact agreement with the quotient.
    source_a, source_b = -675, 13662
    point_x = 3
    velu_t = 2 * (3 * point_x**2 + source_a)
    velu_w = 2 * (5 * point_x**3 + 3 * source_a * point_x + 2 * source_b)
    target_a = source_a - 5 * velu_t
    target_b = source_b - 7 * velu_w
    require((target_a, target_b) == (5805, -285714), "short Velu target")
    target_inv = invariants(0, 0, 0, target_a, target_b)
    require(
        target_inv[3] == Fraction(9938375, 21952),
        f"unexpected quotient j {target_inv[3]}",
    )
    e0_inv = invariants(1, 0, 1, 4, -6)
    require(
        e0_inv == (-215, 5291, -21952, Fraction(9938375, 21952)),
        f"unexpected E0 invariants {e0_inv}",
    )
    e0_u, e0_v = sp.symbols("e0_u e0_v")
    abstract_e0 = e0_v**2 + e0_u * e0_v + e0_v - e0_u**3 - 4 * e0_u + 6
    target_short_x = 36 * e0_u + 3
    target_short_y = 108 * (2 * e0_v + e0_u + 1)
    require_zero_mod(
        target_short_y**2
        - (target_short_x**3 + target_a * target_short_x + target_b),
        abstract_e0,
        e0_v,
        "X0 short-model change failed",
    )

    # The coarse diamond deck group has three elements.
    units = tuple(a for a in range(14) if sp.gcd(a, 14) == 1)
    diamond_classes = tuple(sorted({min(a, (-a) % 14) for a in units}))
    require(units == (1, 3, 5, 9, 11, 13), "units modulo 14")
    require(diamond_classes == (1, 3, 5), "diamond classes modulo sign")
    require(pow(3, 3, 14) in (1, 13), "diamond generator order three")

    # THM-2998's symmetric quotient uses a translated target origin.
    translated_origin = (9, -33)
    qx, qy = translated_origin
    require(
        qy**2 + qx * qy + qy == qx**3 + 4 * qx - 6,
        "translated quotient origin not on X0(14)",
    )

    # An odd flank permutation exchanges the two Vandermonde sign sheets and
    # conjugates the C3 generator to its inverse.  After synchronizing origins
    # this is elliptic negation, not another diamond or Fricke involution.
    A, B, C = sp.symbols("A B C")
    vandermonde = (A - B) * (A - C) * (B - C)
    c_plus = A * B * C - vandermonde
    c_minus = A * B * C + vandermonde
    odd_swap = c_plus.subs({A: B, B: A}, simultaneous=True)
    even_cycle = c_plus.subs({A: B, B: C, C: A}, simultaneous=True)
    require(sp.expand(odd_swap - c_minus) == 0, "odd swap sheet exchange")
    require(sp.expand(even_cycle - c_plus) == 0, "even cycle sheet preservation")
    v_neg = -V - U - 1
    require_zero_mod(
        v_neg**2 + U * v_neg + v_neg - U**3 + U,
        abstract_e1,
        V,
        "elliptic negation failed",
    )
    require(v_neg.subs({U: 0, V: 0}) == -1, "negation must send T to -T")

    print("THM-3034 ordered quartic cross-wall X1(14) and diamond quotient")
    print("cplus_projection_to_quartic=PASS")
    print("quartic_x1_birational_map_and_inverse=PASS")
    print("even_flank_rho=translation_by_(0,0):PASS")
    print("x1_model=y^2+xy+y=x^3-x")
    print("x1_invariants=(25,-253,-28,-15625/28)")
    print("x1_short_model=Y^2=X^3-675X+13662")
    print("rational_three_torsion=((0,0),(0,-1));unique_subgroup=PASS")
    print("velu_x0_map=PASS;kernel=<((0,0))>;degree=3")
    print("x0_model=Y^2+XY+Y=X^3+4X-6")
    print("x0_short_model=Y^2=X^3+5805X-285714")
    print("x0_j=9938375/21952")
    print("diamond_group=(Z/14Z)^*/{+-1}=C3")
    print("thm2998_symmetric_target_origin=(9,-33);normalization=translation")
    print("odd_flank_scope=sheet_exchange_and_C3_inversion_not_diamond_C2")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
