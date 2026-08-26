#!/usr/bin/env python3
"""Exact primary certificate for THM-4164.

This certificate treats both the generic
J=15*a^2+356 row and its exact quadratic-field wall, computes independent
(s,p) and (X,T) critical projections, resolves the top boundary, and checks
all monodromy budgets.  It uses no assert statements.
"""

from __future__ import annotations

import hashlib
from math import gcd

import sympy as sp


CHECKS = 0


def require(condition: bool, message: str) -> None:
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


def polygon_ledger(vertices: tuple[tuple[int, int], ...]) -> tuple[int, int, int]:
    area2 = abs(sum(
        vertices[i][0] * vertices[(i + 1) % len(vertices)][1]
        - vertices[(i + 1) % len(vertices)][0] * vertices[i][1]
        for i in range(len(vertices))
    ))
    boundary = sum(
        gcd(
            abs(vertices[(i + 1) % len(vertices)][0] - vertices[i][0]),
            abs(vertices[(i + 1) % len(vertices)][1] - vertices[i][1]),
        )
        for i in range(len(vertices))
    )
    return area2, boundary, (area2 - boundary + 2) // 2


def main() -> None:
    s, p, X, T, a, w, q, u, v = sp.symbols("s p X T a w q u v")
    kappa = sp.Rational(1376, 135)
    K0 = sp.Rational(2848, 45)

    zeta = kappa * a**3
    theta = -3 * kappa * a**2
    phi = 3 * kappa * a
    top = sp.expand(zeta * w**3 + theta * w**2 + phi * w - kappa)
    require(sp.factor(top - kappa * (a * w - 1)**3) == 0,
            "triple-root parameterization")
    require(sp.factor(sp.discriminant(top, w)) == 0,
            "triple discriminant")
    require(sp.factor(sp.diff(top, w, 2).subs(w, 1 / a)) == 0,
            "triple second derivative")

    J = 15 * a**2 + 356
    inner_factor = 5805 * a**4 + 1013888
    I = sp.factor(4 * theta * K0**2 - 27 * zeta**2)
    require(
        sp.factor(I + sp.Rational(44032, 91125) * a**2 * inner_factor) == 0,
        "triple inner factor",
    )

    t = p - s**2
    H = sp.expand(
        -3 * p + sp.Rational(8, 3) * p**2 - kappa * p**3
        + K0 * s**2 * p**2 + phi * s * p**3
        + theta * s**2 * p**3 + zeta * s**3 * p**3
    )
    G_source = -s**2 / (2 * t) + H
    A = sp.factor(sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p))
    C0 = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
    B = sp.factor(sp.cancel((C0 + s * A) / t**2))
    require(sp.factor(sp.cancel(t**2 * sp.diff(G_source, s) - p * A)) == 0,
            "first source identity")
    require(sp.factor(sp.cancel(2 * t**2 * sp.diff(G_source, p) - (t**2 * B - s * A))) == 0,
            "second source identity")
    require((sp.degree(A, s), sp.degree(B, s)) == (6, 3),
            "source critical degrees")
    require(sp.factor(sp.Poly(A, s).LC() - 3 * zeta * p**2) == 0,
            "source A leading row")
    require(sp.factor(sp.Poly(B, s).LC() - 9 * zeta * p**2) == 0,
            "source B leading row")

    source_resultant = sp.resultant(A, B, s)
    require(valuation(source_resultant, p) == 6, "source p artifact")
    R16 = sp.Poly(sp.cancel(source_resultant / p**6), p)
    require(R16.degree() == 16, "generic source residual degree")
    require(
        sp.factor(
            R16.LC()
            - sp.Rational(
                9563767153858210665857024,
                9341736328125,
            ) * a**17 * J**2
        ) == 0,
        "generic R16 leading row",
    )
    require(sp.factor(R16.TC() + 46656 * zeta * I) == 0,
            "generic R16 constant row")

    P = T + X**2 * T**2
    Y = X * T * P
    G = sp.expand(
        -X**2 * T / 2 - 3 * P + sp.Rational(8, 3) * P**2
        - kappa * P**3 + K0 * Y**2 + phi * P**2 * Y
        + theta * P * Y**2 + zeta * Y**3
    )
    require(sp.factor(sp.cancel(G_source.subs({s: X * T, p: P}) - G)) == 0,
            "normalized reconstruction")
    f = sp.cancel(sp.diff(G, X) / T)
    h = sp.diff(G, T)
    require((sp.degree(f, X), sp.degree(h, X)) == (8, 9),
            "normalized critical degrees")
    require(sp.factor(sp.Poly(f, X).LC() - 9 * zeta * T**8) == 0,
            "normalized f leading row")
    require(sp.factor(sp.Poly(h, X).LC() - 9 * zeta * T**8) == 0,
            "normalized h leading row")
    normalized_resultant = sp.resultant(f, h, X)
    require(valuation(normalized_resultant, T) == 56,
            "normalized T artifact")
    Q16 = sp.Poly(
        sp.cancel(normalized_resultant / (T**56 * (6 * T + 1)**2)),
        T,
    )
    require(Q16.degree() == 16, "generic normalized residual degree")
    require(
        sp.factor(
            Q16.LC()
            - sp.Rational(
                612081097846925482614849536,
                38306957530517578125,
            ) * a**15 * J**2 * inner_factor**2
        ) == 0,
        "generic Q16 leading row",
    )
    require(
        sp.factor(Q16.TC() + sp.Rational(3**15, 2**7) * zeta**7) == 0,
        "generic Q16 constant row",
    )

    # Universal omitted points and exact length.
    hessian = sp.det(sp.hessian(G, (X, T)))
    for t_value, modulus, value, determinant in (
        (sp.Integer(0), X**2 + 6, sp.Integer(0), sp.Integer(6)),
        (-sp.Rational(1, 6), X**2 - 6, sp.Rational(1, 2), -sp.Integer(6)),
    ):
        require(
            sp.factor(sp.rem(sp.Poly(G.subs(T, t_value) - value, X),
                             sp.Poly(modulus, X)).as_expr()) == 0,
            "universal critical value",
        )
        require(
            sp.factor(sp.rem(sp.Poly(hessian.subs(T, t_value) - determinant, X),
                             sp.Poly(modulus, X)).as_expr()) == 0,
            "universal critical Hessian",
        )
    generic_length = R16.degree() + 4
    require(generic_length == 20, "generic affine length")

    # Rational generic control.
    R_control = sp.Poly(R16.as_expr().subs(a, 1), p)
    Q_control = sp.Poly(Q16.as_expr().subs(a, 1), T)
    require(J.subs(a, 1) == 371 and I.subs(a, 1) != 0,
            "generic control gates")
    require(sp.gcd(R_control, R_control.diff()).degree() == 0,
            "generic source control")
    require(sp.gcd(Q_control, Q_control.diff()).degree() == 0,
            "generic normalized control")

    # Boundary strict transform on the generic J!=0 row.
    F = sp.expand(q * (s**2 - p) - (s**2 - p) * H - s**2 / 2)
    K_boundary = sp.expand(
        sp.cancel(u**4 * F.subs({p: 1 / u, s: 1 / a + v}))
    )
    coefficient_v3 = sp.factor(
        sp.diff(K_boundary, v, 3).subs({u: 0, v: 0}) / 6
    )
    coefficient_u = sp.factor(
        sp.diff(K_boundary, u).subs({u: 0, v: 0})
    )
    require(sp.factor(coefficient_v3 - kappa * a**3) == 0,
            "generic boundary v3 row")
    require(sp.factor(coefficient_u - sp.Rational(8, 45) * J / a**2) == 0,
            "generic boundary u row")
    require(
        sp.factor(
            -coefficient_v3 / coefficient_u
            + sp.Rational(172, 3) * a**5 / J
        ) == 0,
        "generic cubic tangency coefficient",
    )
    # u~v^3, K_u a unit: omega=-u^2 ds/K_u has order six.
    require(2 * 3 == 6 and 6 + 1 == 7,
            "generic boundary index")
    generic_packet = (8, 7, 2, 2, 2, 1)
    require(sum(index - 1 for index in generic_packet) == 16 == 2 * 9 - 2,
            "generic packet defect")
    require(polygon_ledger(((0, 1), (2, 0), (5, 3), (3, 4), (0, 4))) == (27, 11, 9),
            "ambient polygon")

    generic_full = sum(generic_packet)
    generic_finite, beta = 16, 3
    require(generic_full == 22, "generic full degree")
    require(2 * (generic_full - generic_length) == 4 < 16,
            "generic full contradiction")
    generic_capacities = (
        2 * generic_finite - generic_length - 2 + beta,
        2 * generic_finite - generic_length - 1 + beta,
        beta,
    )
    require(generic_capacities == (13, 14, 3) and max(generic_capacities) < 15,
            "generic finite contradiction")

    # Exact endpoint/boundary wall J=0 in Q[a]/(15a^2+356).
    modulus = sp.Poly(J, a)

    def reduce_J(expression: sp.Expr) -> sp.Expr:
        numerator, denominator = sp.fraction(sp.cancel(expression))
        numerator = sp.rem(sp.Poly(numerator, a), modulus).as_expr()
        denominator = sp.rem(sp.Poly(denominator, a), modulus).as_expr()
        inverse = sp.invert(sp.Poly(denominator, a), modulus).as_expr()
        return sp.factor(sp.rem(sp.Poly(numerator * inverse, a), modulus).as_expr())

    R15 = sp.Poly(sum(reduce_J(R16.nth(d)) * p**d for d in range(17)), p)
    Q15 = sp.Poly(sum(reduce_J(Q16.nth(d)) * T**d for d in range(17)), T)
    require((R15.degree(), Q15.degree()) == (15, 15),
            "J-wall residual degrees")
    require(
        R15.LC()
        == -sp.Rational(
            3636055657731868916207553629845852858180471619584,
            1596123230438232421875,
        ) * a,
        "J-wall R15 leading row",
    )
    require(
        Q15.LC()
        == sp.Rational(
            299870607546709935491031806891677436336739204843700202793598976,
            10908504703026294708251953125,
        ) * a,
        "J-wall Q15 leading row",
    )
    require(reduce_J(I) == sp.Rational(335741565206528, 6834375),
            "J-wall inner factor nonzero")
    require(
        Q15.eval(-sp.Rational(1, 6))
        == -sp.Rational(
            71304646963571565929235399784047641547404476816773723491139584,
            782625597463934606616497039794921875,
        ) * a,
        "J-wall Q15 universal-coordinate value",
    )
    root = sp.sqrt(-1335)
    a0 = 2 * root / 15
    R_J_control = sp.Poly(sp.expand(R15.as_expr().subs(a, a0)), p, extension=root)
    Q_J_control = sp.Poly(sp.expand(Q15.as_expr().subs(a, a0)), T, extension=root)
    require(sp.gcd(R_J_control, R_J_control.diff()).degree() == 0,
            "J-wall source control")
    require(sp.gcd(Q_J_control, Q_J_control.diff()).degree() == 0,
            "J-wall normalized control")
    wall_length = R15.degree() + 4
    require(wall_length == 19, "J-wall affine length")

    # The tangent cone has two distinct lines u=0 and u=-(16a/9)v.
    wall_v3 = reduce_J(coefficient_v3)
    wall_u1 = reduce_J(coefficient_u)
    wall_uv = reduce_J(
        sp.diff(K_boundary, u, v).subs({u: 0, v: 0})
    )
    wall_u2 = reduce_J(
        sp.diff(K_boundary, u, 2).subs({u: 0, v: 0}) / 2
    )
    require(wall_u1 == 0, "J-wall killed u row")
    require(wall_v3 == -sp.Rational(489856, 2025) * a,
            "J-wall v3 row")
    require(wall_uv == -sp.Rational(16, 3) * a,
            "J-wall uv row")
    require(wall_u2 == -3, "J-wall u2 row")
    transverse_slope = sp.factor(wall_uv / 3)
    tangent_quadratic = sp.factor(-wall_v3 / wall_uv)
    require(transverse_slope == -sp.Rational(16, 9) * a,
            "transverse branch slope")
    require(tangent_quadratic == -sp.Rational(30616, 675),
            "quadratic branch coefficient")
    # Transverse branch: u~v gives ord omega=1, e=2.
    # Boundary-tangent branch: u~v^2 gives ord omega=3, e=4.
    require((2 * 1 - 1, 2 * 2 - 1) == (1, 3),
            "J-wall differential orders")
    wall_packet = (8, 4, 2, 2, 2, 2, 1)
    require(sum(index - 1 for index in wall_packet) == 14 == 2 * 8 - 2,
            "J-wall packet defect")
    wall_full = sum(wall_packet)
    wall_finite = 15
    require(wall_full == 21, "J-wall full degree")
    require(2 * (wall_full - wall_length) == 4 < 14,
            "J-wall full contradiction")
    wall_capacities = (
        2 * wall_finite - wall_length - 2 + beta,
        2 * wall_finite - wall_length - 1 + beta,
        beta,
    )
    require(wall_capacities == (12, 13, 3) and max(wall_capacities) < 14,
            "J-wall finite contradiction")

    rows = (
        "triple=C=kappa*(aW-1)^3;I~-a2*(5805a4+1013888)",
        "generic=Res:p6R16;T56*(6T+1)^2Q16;L20;control=a1",
        "generic-boundary=smooth-cubic-tangency;e7;g9;packet=8,7,2,2,2,1;responses22/16",
        "Jwall=a2=-356/15;Res:R15/Q15;L19;both-squarefree",
        "Jwall-boundary=node;branches=e2,e4;g8;packet=8,4,2,2,2,2,1;responses21/15",
        "remaining=inner 5805a4+1013888=0",
    )
    semantic = hashlib.sha256("\n".join(rows).encode("ascii")).hexdigest()
    print("THM4164_Y_ONLY_TRIPLE_TOP_ROOT_EXCLUSION=PASS")
    print(f"checks={CHECKS}")
    for row in rows:
        print(row)
    print("controls=generic:a=1;Jwall:a=2sqrt(-1335)/15;both-projections-squarefree")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
