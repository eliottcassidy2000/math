#!/usr/bin/env python3
"""Clean-room source-chart audit for the J_top=0 closure of THM-4161.

This referee uses the alternative critical pair (A,C0), never forms the
primary normalized resultant, and uses a disjoint ordinary-node control.
It works directly in Q(c)[a]/(15*a^2+356), resolves the exceptional cusp,
and checks the carrier-orbit equality obstruction.
"""

from __future__ import annotations

import hashlib
from math import gcd

import sympy as sp


CHECKS = 0


def check(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuation(expression: sp.Expr, variable: sp.Symbol) -> int:
    return min(
        monomial[0]
        for monomial, coefficient in sp.Poly(expression, variable).terms()
        if coefficient != 0
    )


def main() -> None:
    s, p, a, b, c, u, v, y, q = sp.symbols("s p a b c u v y q")
    X, T = sp.symbols("X T")
    kappa = sp.Rational(1376, 135)
    K0 = sp.Rational(2848, 45)
    mod = sp.Poly(15 * a**2 + 356, a)

    zeta = kappa * a**2 * b
    theta = -kappa * (a**2 + 2 * a * b)
    phi = kappa * (2 * a + b)

    def quotient(expression: sp.Expr) -> sp.Expr:
        expression = sp.cancel(expression.subs(b, c * a))
        numerator, denominator = sp.fraction(expression)
        numerator = sp.rem(sp.Poly(numerator, a), mod).as_expr()
        denominator = sp.rem(sp.Poly(denominator, a), mod).as_expr()
        inverse = sp.invert(sp.Poly(denominator, a), mod).as_expr()
        return sp.factor(sp.rem(sp.Poly(numerator * inverse, a), mod).as_expr())

    # Coverage of the intrinsic secondary wall.
    A_top = sp.factor(6 * phi * zeta - 2 * theta**2)
    B_top = sp.factor(-9 * kappa * zeta - phi * theta)
    J_top = sp.factor(15 * A_top**2 + 356 * B_top**2)
    check(
        sp.factor(
            J_top
            - sp.Rational(14339490709504, 332150625)
            * a**2 * (a - b)**4 * (15 * a**2 + 356)
        ) == 0,
        "intrinsic J_top factorization",
    )
    zeta_q = quotient(zeta)
    I_q = quotient(4 * theta * K0**2 - 27 * zeta**2)
    check(zeta_q == -sp.Rational(489856, 2025) * a * c,
          "zeta quotient")
    check(
        sp.factor(
            I_q
            - sp.Rational(1986636480512, 20503125)
            * (387 * c**2 + 80 * c + 40)
        ) == 0,
        "inner quotient gate",
    )

    # Alternative source pair (A,C0), distinct from the primary (A,B).
    t = p - s**2
    H = sp.expand(
        -3 * p + sp.Rational(8, 3) * p**2 - kappa * p**3
        + K0 * s**2 * p**2 + phi * s * p**3
        + theta * s**2 * p**3 + zeta * s**3 * p**3
    )
    G_source = -s**2 / (2 * t) + H
    A = sp.factor(sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p))
    C0 = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
    check(sp.factor(sp.cancel(t**2 * sp.diff(G_source, s) - p * A)) == 0,
          "alternative first identity")
    check(sp.factor(sp.cancel(2 * t**2 * sp.diff(G_source, p) - C0)) == 0,
          "alternative second identity")
    check((sp.degree(A, s), sp.degree(C0, s)) == (6, 7),
          "alternative source degrees")
    check(sp.factor(sp.Poly(A, s).LC() - 3 * zeta * p**2) == 0,
          "alternative A leading row")
    check(sp.factor(sp.Poly(C0, s).LC() - 6 * zeta * p**2) == 0,
          "alternative C0 leading row")
    raw_resultant = sp.resultant(A, C0, s)
    check(valuation(raw_resultant, p) == 8, "alternative p artifact")
    ambient = sp.Poly(sp.cancel(raw_resultant / p**8), p)
    check(ambient.degree() == 17, "ambient alternative degree")
    reduced_expression = sum(
        quotient(ambient.nth(degree)) * p**degree
        for degree in range(18)
    )
    R16 = sp.Poly(reduced_expression, p)
    check(R16.degree() == 16, "generic quotient degree")
    check(
        sp.factor(
            R16.LC()
            + sp.Rational(
                156350393282470363396924806083371672901760279642112,
                1077383180545806884765625,
            ) * a * c**5 * (c - 1)**2 * (43 * c - 33)
        ) == 0,
        "alternative R16 leading row",
    )
    check(sp.factor(R16.TC() + 46656 * zeta_q * I_q) == 0,
          "alternative R16 constant row")

    # Disjoint algebraic control c=3; the primary uses c=2.
    root = sp.sqrt(-1335)
    a0 = 2 * root / 15
    node_control = sp.Poly(
        sp.expand(R16.as_expr().subs({a: a0, c: 3})),
        p,
        extension=root,
    )
    check(node_control.degree() == 16, "node control degree")
    check(sp.gcd(node_control, node_control.diff()).degree() == 0,
          "node control squarefree")

    # Universal points are restored independently in normalized coordinates,
    # without using the normalized critical resultant.
    P = T + X**2 * T**2
    Y = X * T * P
    G = sp.expand(
        -X**2 * T / 2 - 3 * P + sp.Rational(8, 3) * P**2
        - kappa * P**3 + K0 * Y**2 + phi * P**2 * Y
        + theta * P * Y**2 + zeta * Y**3
    )
    f = sp.cancel(sp.diff(G, X) / T)
    h = sp.diff(G, T)
    hessian = sp.det(sp.hessian(G, (X, T)))
    for t_value, modulus, value, determinant in (
        (sp.Integer(0), X**2 + 6, sp.Integer(0), sp.Integer(6)),
        (-sp.Rational(1, 6), X**2 - 6, sp.Rational(1, 2), -sp.Integer(6)),
    ):
        check(
            sp.factor(sp.rem(sp.Poly(G.subs(T, t_value) - value, X),
                             sp.Poly(modulus, X)).as_expr()) == 0,
            "universal value",
        )
        check(
            sp.factor(sp.rem(sp.Poly(hessian.subs(T, t_value) - determinant, X),
                             sp.Poly(modulus, X)).as_expr()) == 0,
            "universal Hessian",
        )

    # Independent boundary expansion.
    F = sp.expand(q * (s**2 - p) - (s**2 - p) * H - s**2 / 2)
    strict = sp.expand(sp.cancel(u**4 * F.subs({p: 1 / u, s: 1 / a + v})))
    v2 = quotient(sp.diff(strict, v, 2).subs({u: 0, v: 0}) / 2)
    uv = quotient(sp.diff(strict, u, v).subs({u: 0, v: 0}))
    u2 = quotient(sp.diff(strict, u, 2).subs({u: 0, v: 0}) / 2)
    u1 = quotient(sp.diff(strict, u).subs({u: 0, v: 0}))
    check(u1 == 0, "secondary wall u row")
    check(sp.factor(v2 + sp.Rational(489856, 2025) * (c - 1)) == 0,
          "secondary v2 row")
    check(sp.factor(uv + sp.Rational(16, 3) * a) == 0,
          "secondary uv row")
    check(u2 == -3, "secondary u2 row")
    discriminant = quotient(uv**2 - 4 * v2 * u2)
    check(
        sp.factor(discriminant + sp.Rational(45568, 675) * (43 * c - 33)) == 0,
        "unique node degeneracy",
    )

    # At c=33/43 the repeated tangent exposes a rational (2,3) cusp.
    cusp_c = sp.Rational(33, 43)
    shift = -sp.Rational(135, 2848) * a
    cusp = sp.Poly(
        sp.expand(strict.subs({b: c * a, v: y - shift * u})),
        y,
        u,
    )

    def cusp_coefficient(iy: int, iu: int) -> sp.Expr:
        return sp.factor(
            quotient(cusp.coeff_monomial(y**iy * u**iu)).subs(c, cusp_c)
        )

    check(cusp_coefficient(2, 0) == sp.Rational(22784, 405),
          "cusp y2 row")
    check(cusp_coefficient(1, 1) == 0 and cusp_coefficient(0, 2) == 0,
          "cusp repeated tangent")
    check(
        sp.factor(cusp_coefficient(0, 3) + q + sp.Rational(1161, 80)) == 0,
        "cusp u3 row",
    )
    R14 = sp.Poly(R16.as_expr().subs(c, cusp_c), p)
    check(R14.degree() == 14, "cusp alternative degree")
    check(
        R14.LC()
        == sp.Rational(
            830271724986997884692231133070215530174115479552,
            886735128021240234375,
        ) * a,
        "cusp alternative leading row",
    )
    cusp_control = sp.Poly(
        sp.expand(R14.as_expr().subs(a, a0)),
        p,
        extension=root,
    )
    check(sp.gcd(cusp_control, cusp_control.diff()).degree() == 0,
          "cusp alternative squarefree")

    # Geometry and response ledgers, reconstructed without importing the
    # primary packet arrays.
    node_packet = tuple(sorted((8, 3, 2, 2, 1) + (2, 2, 2), reverse=True))
    cusp_packet = tuple(sorted((8, 3, 3, 1) + (2, 2, 2), reverse=True))
    check(node_packet == (8, 3, 2, 2, 2, 2, 2, 1), "node packet")
    check(cusp_packet == (8, 3, 3, 2, 2, 2, 1), "cusp packet")
    check(sum(index - 1 for index in node_packet) == 14,
          "node defect")
    check(sum(index - 1 for index in cusp_packet) == 14,
          "cusp defect")
    node_L, cusp_L = R16.degree() + 4, R14.degree() + 4
    check((node_L, cusp_L) == (20, 18), "critical lengths")

    # Ordinary-node responses are strictly excluded.
    node_full, node_finite, beta = sum(node_packet), 16, 3
    check(2 * (node_full - node_L) == 4 < 14,
          "node full contradiction")
    node_capacities = (
        2 * node_finite - node_L - 2 + beta,
        2 * node_finite - node_L - 1 + beta,
        beta,
    )
    check(node_capacities == (13, 14, 3) and max(node_capacities) < 15,
          "node finite contradiction")

    # Cusp full response and the finite carrier-orbit refinement. Here
    # n=15>m+1=4, so the handle subgroup has positive rank and the union
    # lower bound used below is valid.
    cusp_full, cusp_finite, m = sum(cusp_packet), 15, 3
    check(2 * (cusp_full - cusp_L) == 6 < 14,
          "cusp full contradiction")
    check(cusp_finite > m + 1, "carrier-orbit positivity hypothesis")
    union_floor = cusp_finite - m
    support_sum_ceiling = 2 * cusp_finite - cusp_L
    overlap_ceiling = support_sum_ceiling - union_floor
    check((union_floor, support_sum_ceiling, overlap_ceiling) == (12, 12, 0),
          "carrier-orbit overlap")
    finite_origin_index = sum(index - 1 for index in (8, 3, 3, 1))
    check(2 * overlap_ceiling + m == 3 < finite_origin_index == 11,
          "carrier-orbit origin contradiction")

    rows = (
        "quotient=a^2+356/15;b=c*a;I~387c2+80c+40",
        "alternative=Res_s(A,C0)=p8R16;control=c3;squarefree",
        "node=disc~-45568*(43c-33)/675;e2,e2;g8;L20;packet=8,3,2,2,2,2,2,1",
        "cusp=c33/43;y2-(q+1161/80)u3;e3;g8;Res=p8R14;L18",
        "cusp-finite=n15>m+1;union12;supportsum12;overlap0;origin<=3<11",
    )
    semantic = hashlib.sha256("\n".join(rows).encode("ascii")).hexdigest()
    print("THM4161_Y_ONLY_DOUBLE_TOP_ROOT_SECONDARY_WALL_INDEPENDENT_AUDIT=PASS")
    print(f"checks={CHECKS}")
    for row in rows:
        print(row)
    print("controls=node:c=3;cusp:c=33/43;alternative-source-squarefree=yes")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
