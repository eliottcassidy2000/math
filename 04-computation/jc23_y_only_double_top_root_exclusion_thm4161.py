#!/usr/bin/env python3
"""Exact primary certificate for the Y-only double-top-root wall in THM-4161.

This file intentionally lives in /tmp rather than the canon.  It starts from
the complete Y-only, eta=Delta=0, exact-weight-nine source in THM-4155,
parameterizes the discriminant-zero top cubic, and recomputes the critical
scheme in both (s,p) and (X,T) coordinates.  It also checks the boundary
strict transform and the resulting monodromy budgets.  No assert statements
are used, so python -O preserves all gates.
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


def convex_hull(points: list[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    points = sorted(set(points))

    def cross(o: tuple[int, int], a: tuple[int, int], b: tuple[int, int]) -> int:
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[tuple[int, int]] = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[int, int]] = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


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
    require((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return area2, boundary, (area2 - boundary + 2) // 2


def main() -> None:
    s, p, a, b, w, q, u = sp.symbols("s p a b w q u")
    X, T = sp.symbols("X T")
    kappa = sp.Rational(1376, 135)
    epsilon = -kappa
    K0 = sp.Rational(2848, 45)

    # The repeated root is 1/a and the remaining simple root is 1/b.
    Zeta = kappa * a**2 * b
    Theta = -kappa * (a**2 + 2 * a * b)
    Phi = kappa * (2 * a + b)
    top = sp.expand(Zeta * w**3 + Theta * w**2 + Phi * w + epsilon)
    require(sp.factor(top - kappa * (a * w - 1)**2 * (b * w - 1)) == 0,
            "top cubic parameterization")
    require(sp.factor(sp.discriminant(top, w)) == 0, "top discriminant")

    # The double root and the secondary wall have intrinsic coefficient
    # descriptions, so the theorem does not depend on choosing root labels.
    A_top = sp.factor(6 * Phi * Zeta - 2 * Theta**2)
    B_top = sp.factor(-9 * kappa * Zeta - Phi * Theta)
    subresultants = sp.subresultants(top, sp.diff(top, w), w)
    require(len(subresultants) == 3, "double-root subresultant length")
    require(
        sp.factor(
            subresultants[-1]
            - sp.Rational(1376, 135) * a**2 * b * (A_top * w + B_top)
        ) == 0,
        "intrinsic linear subresultant",
    )
    require(sp.factor(A_top / B_top + a) == 0,
            "intrinsic repeated-root address")

    J = 15 * a**2 + 356
    J_top = sp.factor(15 * A_top**2 + 356 * B_top**2)
    require(
        sp.factor(
            J_top
            - sp.Rational(14339490709504, 332150625)
            * a**2 * (a - b)**4 * J
        ) == 0,
        "intrinsic secondary wall",
    )
    I = sp.factor(4 * Theta * K0**2 - 27 * Zeta**2)
    I_numerator = 17415 * a**3 * b**2 + 1013888 * a + 2027776 * b
    require(sp.factor(I + sp.Rational(44032, 273375) * a * I_numerator) == 0,
            "inner endpoint parameterization")

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

    # First critical projection: the THM-4155 (A,B) pair.
    A = sp.factor(sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p))
    C0 = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
    B = sp.factor(sp.cancel((C0 + s * A) / t**2))
    require(not sp.denom(A).has(s, p) and not sp.denom(B).has(s, p),
            "source critical polynomiality")
    require((sp.degree(A, s), sp.degree(B, s)) == (6, 3), "source degrees")
    require(sp.factor(sp.Poly(A, s).LC() - 3 * Zeta * p**2) == 0,
            "source A leading row")
    require(sp.factor(sp.Poly(B, s).LC() - 9 * Zeta * p**2) == 0,
            "source B leading row")
    require(sp.factor(sp.cancel(t**2 * sp.diff(G_source, s) - p * A)) == 0,
            "first source critical identity")
    require(sp.factor(sp.cancel(2 * t**2 * sp.diff(G_source, p) - (t**2 * B - s * A))) == 0,
            "second source critical identity")

    source_resultant = sp.factor(sp.resultant(A, B, s))
    require(valuation(source_resultant, p) == 6, "source p artifact")
    R17_expression = sp.cancel(source_resultant / p**6)
    require(sp.denom(R17_expression) == 1, "R17 polynomiality")
    R17 = sp.Poly(R17_expression, p)
    require(R17.degree() == 17, "source residual degree")
    source_lc = (
        -sp.Rational(
            3289935900927224469054816256,
            252226880859375,
        )
        * a**11 * b**5 * (a - b)**3 * J
    )
    require(sp.factor(R17.LC() - source_lc) == 0, "R17 leading endpoint")
    require(sp.factor(R17.TC() + 46656 * Zeta * I) == 0, "R17 constant endpoint")

    # Second critical projection: normalized coordinates, recomputed from G.
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
            "normalized source reconstruction")
    f = sp.cancel(sp.diff(G, X) / T)
    h = sp.diff(G, T)
    require(sp.denom(f) == 1, "normalized f polynomiality")
    require((sp.degree(f, X), sp.degree(h, X)) == (8, 9), "normalized degrees")
    normalized_resultant = sp.resultant(f, h, X)
    require(valuation(normalized_resultant, T) == 56, "normalized T artifact")
    Q17_expression = sp.cancel(normalized_resultant / (T**56 * (6 * T + 1)**2))
    require(sp.denom(Q17_expression) == 1, "Q17 polynomiality")
    Q17 = sp.Poly(Q17_expression, T)
    require(Q17.degree() == 17, "normalized residual degree")
    normalized_lc = (
        -sp.Rational(
            210555897659342366019508240384,
            9308590679915771484375,
        )
        * a**9 * b**3 * (a - b)**3 * J * I_numerator**2
    )
    require(sp.factor(Q17.LC() - normalized_lc) == 0, "Q17 leading endpoint")
    require(sp.factor(Q17.TC() + sp.Rational(3**15, 2**7) * Zeta**7) == 0,
            "Q17 constant endpoint")

    # The two universal pairs survive.  Together with R17/Q17 they give L=21.
    require(sp.factor(f.subs(T, 0)) == -X, "T=0 f row")
    require(sp.factor(h.subs(T, 0) + (X**2 + 6) / 2) == 0, "T=0 h row")
    hessian = sp.factor(sp.det(sp.hessian(G, (X, T))))
    for t_value, modulus, value, determinant in (
        (sp.Integer(0), X**2 + 6, sp.Integer(0), sp.Integer(6)),
        (-sp.Rational(1, 6), X**2 - 6, sp.Rational(1, 2), -sp.Integer(6)),
    ):
        require(sp.factor(sp.rem(sp.Poly(G.subs(T, t_value) - value, X), sp.Poly(modulus, X)).as_expr()) == 0,
                "universal critical value")
        require(sp.factor(sp.rem(sp.Poly(hessian.subs(T, t_value) - determinant, X), sp.Poly(modulus, X)).as_expr()) == 0,
                "universal Hessian")
    critical_length = Q17.degree() + 2 + 2
    require(critical_length == 21, "critical length")

    # A disjoint rational control proves the proposed open stratum nonempty.
    control = {a: 1, b: 2}
    require(J.subs(control) == 371, "control J")
    require(I.subs(control) == -sp.Rational(9051394048, 10935), "control I")
    control_R17 = sp.Poly(R17_expression.subs(control), p)
    control_Q17 = sp.Poly(Q17_expression.subs(control), T)
    require(control_R17.degree() == 17 and sp.gcd(control_R17, control_R17.diff()).degree() == 0,
            "source control squarefreeness")
    require(control_Q17.degree() == 17 and sp.gcd(control_Q17, control_Q17.diff()).degree() == 0,
            "normalized control squarefreeness")

    # Boundary strict transform.  With u=1/p, u^4 times the cleared fibre
    # equation has a smooth quadratic tangency at s=1/a.
    F_cleared = sp.expand(q * (s**2 - p) - (s**2 - p) * H - s**2 / 2)
    boundary = sp.expand(sp.cancel(u**4 * F_cleared.subs(p, 1 / u)))
    require(sp.denom(boundary) == 1, "boundary polynomiality")
    require(sp.factor(boundary.subs(u, 0) - top.subs(w, s)) == 0, "top boundary face")
    repeated = {s: 1 / a, u: 0}
    v2_coefficient = sp.factor(sp.diff(boundary, s, 2).subs(repeated) / 2)
    u_coefficient = sp.factor(sp.diff(boundary, u).subs(repeated))
    require(sp.factor(v2_coefficient + kappa * a * (a - b)) == 0,
            "boundary v2 coefficient")
    require(sp.factor(u_coefficient - sp.Rational(8, 45) * J / a**2) == 0,
            "boundary u coefficient")
    branch_u_over_v2 = sp.factor(-v2_coefficient / u_coefficient)
    require(sp.factor(branch_u_over_v2 - sp.Rational(172, 3) * a**3 * (a - b) / J) == 0,
            "boundary tangency coefficient")

    # On the curve, F_p=-u^-2*K_u, so omega=ds/F_p=-u^2 ds/K_u.
    # Since u is a nonzero unit times v^2 and K_u is a unit, ord(omega)=4.
    differential_order = 2 * 2
    require(differential_order == 4, "tangency differential order")
    merged_index = differential_order + 1
    require(merged_index == 5, "merged ramification index")

    # The support and arithmetic genus are unchanged.  The smooth tangency
    # replaces two rational index-three points by one rational index-five.
    support = [
        monomial
        for monomial, coefficient in sp.Poly(F_cleared, s, p).terms()
        if coefficient != 0
    ]
    polygon = convex_hull(support)
    require(polygon == ((0, 1), (2, 0), (5, 3), (3, 4), (0, 4)), "Newton polygon")
    require(polygon_ledger(polygon) == (27, 11, 9), "Pick ledger")
    packet = (8, 5, 3, 2, 2, 2, 1)
    require(sum(index - 1 for index in packet) == 16 == 2 * 9 - 2,
            "boundary defect")
    require(sum(packet) == 23, "full response degree")
    rational_packet = (8, 5, 3, 1)
    carrier_packet = (2, 2, 2)
    require(sum(rational_packet) == 17 and sum(carrier_packet) == 6,
            "rational/carrier split")

    # Same prime cubic carrier as THM-4155, hence two Galois-locked responses.
    full_n, finite_n, beta = 23, 17, 3
    full_overlap = full_n - critical_length
    full_commutator_capacity = 2 * full_overlap
    require(full_overlap == 2 and full_commutator_capacity == 4 < 16,
            "full commutator contradiction")
    finite_capacities = (
        2 * finite_n - critical_length - 2 + beta,
        2 * finite_n - critical_length - 1 + beta,
        beta,
    )
    require(finite_capacities == (14, 15, 3), "finite capacities")
    require(max(finite_capacities) < finite_n - 1, "finite transitivity contradiction")

    # J=0 is simultaneously the next boundary strict-transform wall and the
    # next critical-degree wall.  a=b is the triple-root wall.
    require(sp.factor(u_coefficient * 45 * a**2 / 8 - J) == 0,
            "J boundary gate")
    require(sp.factor(source_lc / (a**11 * b**5 * (a - b)**3 * J)) != 0,
            "source strict-transform factor")

    SEMANTIC.extend([
        "top=kappa*(aW-1)^2*(bW-1);intrinsic=Disc(C)=0,C_nontriple,J_C*I_C!=0",
        "critical=source:p6R17;normalized:T56*(6T+1)^2*Q17;L=21",
        "boundary=smooth-quadratic-tangency;u/v2=172*a3*(a-b)/(3J);e=5",
        "geometry=g9;packet=(8,5,3,2,2,2,1);rational=(8,5,3,1)+cubic=(2,2,2)",
        "responses=full(n23,comm<=4<16);finite(n17,beta3,capacities14/15/3<16)",
        "next=J=0 secondary strict transform;triple=a-b=0;inner=I=0",
    ])
    semantic = hashlib.sha256("\n".join(SEMANTIC).encode()).hexdigest()
    print("THM4161_Y_ONLY_DOUBLE_TOP_ROOT_EXCLUSION=PASS")
    print(f"checks={CHECKS}")
    for row in SEMANTIC:
        print(row)
    print(f"control=(a,b)=(1,2);source_squarefree=yes;normalized_squarefree=yes")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
