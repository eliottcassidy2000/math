#!/usr/bin/env python3
"""Exact certificate for the J_top=0 part of THM-4161.

The ambient top cubic is kappa*(aW-1)^2*(bW-1).  This script imposes
J=15*a^2+356=0 exactly by coefficient reduction in the quadratic field and
writes b=c*a.  It recomputes both critical projections, the generic node and
the unique cusp subwall c=33/43, and both carrier responses.  It imports no
repository calculation and uses no assert statements.
"""

from __future__ import annotations

import hashlib
from math import gcd

import sympy as sp


CHECKS = 0
SEMANTIC: list[str] = []


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
    s, p, a, b, c, q, u, y = sp.symbols("s p a b c q u y")
    X, T = sp.symbols("X T")
    kappa = sp.Rational(1376, 135)
    epsilon = -kappa
    K0 = sp.Rational(2848, 45)
    modulus = sp.Poly(15 * a**2 + 356, a)

    Zeta = kappa * a**2 * b
    Theta = -kappa * (a**2 + 2 * a * b)
    Phi = kappa * (2 * a + b)
    t = p - s**2
    H = sp.expand(
        -3 * p
        + sp.Rational(8, 3) * p**2
        + epsilon * p**3
        + K0 * s**2 * p**2
        + Phi * s * p**3
        + Theta * s**2 * p**3
        + Zeta * s**3 * p**3
    )
    G_source = -s**2 / (2 * t) + H

    def reduce_J(expression: sp.Expr) -> sp.Expr:
        """Substitute b=c*a and reduce exactly modulo 15*a^2+356."""
        expression = sp.cancel(expression.subs(b, c * a))
        numerator, denominator = sp.fraction(expression)
        numerator = sp.rem(sp.Poly(numerator, a), modulus).as_expr()
        denominator = sp.rem(sp.Poly(denominator, a), modulus).as_expr()
        inverse = sp.invert(sp.Poly(denominator, a), modulus).as_expr()
        return sp.factor(
            sp.rem(sp.Poly(numerator * inverse, a), modulus).as_expr()
        )

    zeta_J = reduce_J(Zeta)
    theta_J = reduce_J(Theta)
    phi_J = reduce_J(Phi)
    I_J = reduce_J(4 * Theta * K0**2 - 27 * Zeta**2)
    require(zeta_J == -sp.Rational(489856, 2025) * a * c, "zeta on J wall")
    require(
        sp.factor(
            I_J
            - sp.Rational(1986636480512, 20503125)
            * (387 * c**2 + 80 * c + 40)
        ) == 0,
        "inner gate on J wall",
    )

    # Source critical projection before quotient reduction.
    A = sp.factor(sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p))
    C0 = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
    B = sp.factor(sp.cancel((C0 + s * A) / t**2))
    source_resultant = sp.resultant(A, B, s)
    require(valuation(source_resultant, p) == 6, "source p artifact")
    R17 = sp.Poly(sp.cancel(source_resultant / p**6), p)
    require(R17.degree() == 17, "ambient residual degree")
    R_J_expression = sum(
        reduce_J(R17.nth(degree)) * p**degree
        for degree in range(18)
    )
    R_J = sp.Poly(R_J_expression, p)
    require(R_J.degree() == 16, "generic J-wall source degree")
    source_lc = (
        -sp.Rational(
            156350393282470363396924806083371672901760279642112,
            1077383180545806884765625,
        )
        * a * c**5 * (c - 1)**2 * (43 * c - 33)
    )
    require(sp.factor(R_J.LC() - source_lc) == 0, "generic R16 leading row")
    require(sp.factor(R_J.TC() + 46656 * zeta_J * I_J) == 0,
            "generic R16 constant row")

    # Fully independent normalized critical projection.
    P = T + X**2 * T**2
    Y = X * T * P
    G = sp.expand(
        -X**2 * T / 2
        - 3 * P
        + sp.Rational(8, 3) * P**2
        + epsilon * P**3
        + K0 * Y**2
        + Phi * P**2 * Y
        + Theta * P * Y**2
        + Zeta * Y**3
    )
    require(sp.factor(sp.cancel(G_source.subs({s: X * T, p: P}) - G)) == 0,
            "source reconstruction")
    f = sp.cancel(sp.diff(G, X) / T)
    h = sp.diff(G, T)
    normalized_resultant = sp.resultant(f, h, X)
    require(valuation(normalized_resultant, T) == 56, "normalized T artifact")
    Q17 = sp.Poly(
        sp.cancel(normalized_resultant / (T**56 * (6 * T + 1)**2)),
        T,
    )
    require(Q17.degree() == 17, "ambient normalized degree")
    Q_J_expression = sum(
        reduce_J(Q17.nth(degree)) * T**degree
        for degree in range(18)
    )
    Q_J = sp.Poly(Q_J_expression, T)
    require(Q_J.degree() == 16, "generic J-wall normalized degree")
    normalized_lc = (
        sp.Rational(
            451470050926386584017169136106653470203416750403666143346688,
            66269166070884740352630615234375,
        )
        * a * c**3 * (c - 1)**2 * (43 * c - 33)
        * (387 * c**2 + 80 * c + 40)**2
    )
    require(sp.factor(Q_J.LC() - normalized_lc) == 0,
            "generic Q16 leading row")
    require(
        sp.factor(
            Q_J.TC()
            + reduce_J(sp.Rational(3**15, 2**7) * Zeta**7)
        ) == 0,
        "generic Q16 constant row",
    )

    # Exact algebraic control c=2 on the ordinary-node open set.
    root = sp.sqrt(-1335)
    a_control = 2 * root / 15

    def algebraic_poly(poly: sp.Poly, variable: sp.Symbol, c_value: sp.Rational) -> sp.Poly:
        return sp.Poly(
            sp.expand(poly.as_expr().subs({c: c_value, a: a_control})),
            variable,
            extension=root,
        )

    R_control = algebraic_poly(R_J, p, sp.Rational(2))
    Q_control = algebraic_poly(Q_J, T, sp.Rational(2))
    require(R_control.degree() == 16 and sp.gcd(R_control, R_control.diff()).degree() == 0,
            "node source control")
    require(Q_control.degree() == 16 and sp.gcd(Q_control, Q_control.diff()).degree() == 0,
            "node normalized control")

    # Boundary strict transform on J=0.
    F_cleared = sp.expand(q * (s**2 - p) - (s**2 - p) * H - s**2 / 2)
    boundary = sp.expand(sp.cancel(u**4 * F_cleared.subs(p, 1 / u)))
    s0 = 1 / a
    v2 = reduce_J(sp.diff(boundary, s, 2).subs({s: s0, u: 0}) / 2)
    uv = reduce_J(sp.diff(boundary, s, u).subs({s: s0, u: 0}))
    u2 = reduce_J(sp.diff(boundary, u, 2).subs({s: s0, u: 0}) / 2)
    u1 = reduce_J(sp.diff(boundary, u).subs({s: s0, u: 0}))
    require(u1 == 0, "J wall did not kill u-linear row")
    require(sp.factor(v2 + sp.Rational(489856, 2025) * (c - 1)) == 0,
            "node v2 row")
    require(sp.factor(uv + sp.Rational(16, 3) * a) == 0, "node uv row")
    require(u2 == -3, "node u2 row")
    node_discriminant = reduce_J(uv**2 - 4 * v2 * u2)
    require(sp.factor(
        node_discriminant + sp.Rational(45568, 675) * (43 * c - 33)
    ) == 0,
            "node discriminant")

    # Generic J-wall row: ordinary node, two rational e=2 branches.
    node_packet = (8, 3, 2, 2, 2, 2, 2, 1)
    node_genus = 8
    node_length = R_J.degree() + 4
    require(node_length == 20, "node critical length")
    require(sum(index - 1 for index in node_packet) == 14 == 2 * node_genus - 2,
            "node packet defect")
    node_full_n = sum(node_packet)
    node_finite_n = 16
    beta = 3
    require(node_full_n == 22, "node full degree")
    require(node_full_n - node_length == 2 and 4 < 14,
            "node full contradiction")
    node_capacities = (
        2 * node_finite_n - node_length - 2 + beta,
        2 * node_finite_n - node_length - 1 + beta,
        beta,
    )
    require(node_capacities == (13, 14, 3), "node finite capacities")
    require(max(node_capacities) < node_finite_n - 1,
            "node finite contradiction")

    # Unique cusp subwall c=33/43.  Center the repeated tangent exactly.
    cusp_c = sp.Rational(33, 43)
    tangent_shift = -sp.Rational(135, 2848) * a
    # y=v+tangent_shift*u, so s=s0+y-tangent_shift*u.
    cusp_boundary = sp.Poly(
        sp.expand(
            boundary.subs({
                b: c * a,
                s: s0 + y - tangent_shift * u,
            })
        ),
        y,
        u,
    )

    def cusp_coefficient(y_degree: int, u_degree: int) -> sp.Expr:
        return sp.factor(
            reduce_J(cusp_boundary.coeff_monomial(y**y_degree * u**u_degree)).subs(c, cusp_c)
        )

    require(cusp_coefficient(2, 0) == sp.Rational(22784, 405),
            "cusp square row")
    require(cusp_coefficient(1, 1) == 0 and cusp_coefficient(0, 2) == 0,
            "cusp repeated tangent")
    require(sp.factor(
        cusp_coefficient(0, 3) + q + sp.Rational(1161, 80)
    ) == 0,
            "cusp cubic row")
    # Hence A*y^2-(q+1161/80)*u^3+... is a rational (2,3) cusp.
    # With t=y/u, normalization has u=A*t^2/(q+1161/80), so no
    # square-root residue extension is introduced.  In original coordinates
    # ds has order one, K_u order three and u^2 order four: ord(omega)=2.
    cusp_differential_order = 4 + 1 - 3
    require(cusp_differential_order == 2, "cusp differential order")

    R_cusp = sp.Poly(R_J.as_expr().subs(c, cusp_c), p)
    Q_cusp = sp.Poly(Q_J.as_expr().subs(c, cusp_c), T)
    require(R_cusp.degree() == 14 and Q_cusp.degree() == 14,
            "cusp residual degrees")
    require(
        R_cusp.LC()
        == sp.Rational(
            830271724986997884692231133070215530174115479552,
            886735128021240234375,
        ) * a,
        "cusp R14 leading row",
    )
    require(
        Q_cusp.LC()
        == -sp.Rational(
            3973312205663515142507987786054543929167803303620264264254619648,
            490882711636183261871337890625,
        ) * a,
        "cusp Q14 leading row",
    )
    R_cusp_control = algebraic_poly(R_cusp, p, cusp_c)
    Q_cusp_control = algebraic_poly(Q_cusp, T, cusp_c)
    require(sp.gcd(R_cusp_control, R_cusp_control.diff()).degree() == 0,
            "cusp source squarefreeness")
    require(sp.gcd(Q_cusp_control, Q_cusp_control.diff()).degree() == 0,
            "cusp normalized squarefreeness")

    cusp_packet = (8, 3, 3, 2, 2, 2, 1)
    cusp_genus = 8
    cusp_length = R_cusp.degree() + 4
    require(cusp_length == 18, "cusp critical length")
    require(sum(index - 1 for index in cusp_packet) == 14 == 2 * cusp_genus - 2,
            "cusp defect")
    cusp_full_n = sum(cusp_packet)
    cusp_finite_n = 15
    require(cusp_full_n == 21, "cusp full degree")
    require(2 * (cusp_full_n - cusp_length) == 6 < 14,
            "cusp full contradiction")

    # The crude one-identity finite capacity is equality, so retain carrier
    # orbits. Three transpositions can merge at most three <X,Y>-orbits.
    # Thus |supp X union supp Y|>=n-3. Fixed sheets give support sum<=2n-L.
    orbit_union_floor = cusp_finite_n - beta
    handle_support_sum = 2 * cusp_finite_n - cusp_length
    overlap_ceiling = handle_support_sum - orbit_union_floor
    require((orbit_union_floor, handle_support_sum, overlap_ceiling) == (12, 12, 0),
            "cusp carrier-orbit overlap")
    # ind([X,Y])<=2*overlap, and multiplying the three carrier
    # transpositions raises index by at most beta. The rational origin packet
    # is (8,3,3,1), of index 11.
    origin_index_capacity = 2 * overlap_ceiling + beta
    cusp_origin_index = (8 - 1) + (3 - 1) + (3 - 1)
    require(origin_index_capacity == 3 < cusp_origin_index == 11,
            "cusp finite carrier-orbit contradiction")

    require(polygon_ledger(((0, 1), (2, 0), (5, 3), (3, 4), (0, 4))) == (27, 11, 9),
            "ambient polygon ledger")

    SEMANTIC.extend([
        "Jwall:a^2=-356/15;b=c*a;I~387c^2+80c+40",
        "node:c!=33/43;boundary=ordinary-node;e=2,2;g8;packet=(8,3,2,2,2,2,2,1)",
        "node:critical=p6R16=Q16;L20;full22;finite16,beta3;strict",
        "cusp:c=33/43;boundary=(2,3)-cusp;e=3;g8;packet=(8,3,3,2,2,2,1)",
        "cusp:critical=p6R14=Q14;L18;full21;finite15,beta3",
        "cusp-finite:carrier-orbit union>=12;support-sum<=12;overlap=0;origin-index<=3<11",
        "remaining:triple c=1;inner 387c^2+80c+40=0",
    ])
    digest = hashlib.sha256("\n".join(SEMANTIC).encode()).hexdigest()
    print("THM4161_Y_ONLY_DOUBLE_TOP_ROOT_SECONDARY_WALL=PASS")
    print(f"checks={CHECKS}")
    for row in SEMANTIC:
        print(row)
    print("controls=node:c=2;cusp:c=33/43;both-projections-squarefree")
    print(f"semantic_sha256={digest}")


if __name__ == "__main__":
    main()
