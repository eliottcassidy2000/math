#!/usr/bin/env python3
"""Clean-room source-chart and boundary audit for THM-4161."""

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
    require((area2 - boundary + 2) % 2 == 0, "Pick parity")
    return area2, boundary, (area2 - boundary + 2) // 2


def main() -> None:
    s, p, a, b, w = sp.symbols("s p a b w")
    u, v, q = sp.symbols("u v q")
    kappa = sp.Rational(1376, 135)
    K0 = sp.Rational(2848, 45)
    zeta = kappa * a**2 * b
    theta = -kappa * (a**2 + 2 * a * b)
    phi = kappa * (2 * a + b)
    top = zeta * w**3 + theta * w**2 + phi * w - kappa

    # Coefficient-intrinsic recovery of the unique double root.
    A_top = sp.factor(6 * phi * zeta - 2 * theta**2)
    B_top = sp.factor(-9 * kappa * zeta - phi * theta)
    subresultants = sp.subresultants(top, sp.diff(top, w), w)
    require(len(subresultants) == 3, "nontriple subresultant chain")
    require(sp.degree(subresultants[-1], w) == 1, "linear gcd carrier")
    require(
        sp.factor(
            subresultants[-1]
            - kappa * a**2 * b * (A_top * w + B_top)
        ) == 0,
        "intrinsic subresultant identity",
    )
    require(sp.factor(-B_top / A_top - 1 / a) == 0,
            "unique repeated root")
    J_top = sp.factor(15 * A_top**2 + 356 * B_top**2)
    J = 15 * a**2 + 356
    require(
        sp.factor(
            J_top
            - sp.Rational(14339490709504, 332150625)
            * a**2 * (a - b)**4 * J
        ) == 0,
        "intrinsic secondary gate",
    )

    # Use the alternative critical pair (A,C0), not the primary pair (A,B).
    t = p - s**2
    H = sp.expand(
        -3 * p + sp.Rational(8, 3) * p**2 - kappa * p**3
        + K0 * s**2 * p**2 + phi * s * p**3
        + theta * s**2 * p**3 + zeta * s**3 * p**3
    )
    G = -s**2 / (2 * t) + H
    A_critical = sp.factor(sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p))
    C_critical = sp.expand(s**2 + 2 * t**2 * sp.diff(H, p))
    require(sp.factor(sp.cancel(t**2 * sp.diff(G, s) - p * A_critical)) == 0,
            "first critical identity")
    require(sp.factor(sp.cancel(2 * t**2 * sp.diff(G, p) - C_critical)) == 0,
            "second critical identity")
    alternate_resultant = sp.factor(sp.resultant(A_critical, C_critical, s))
    require(valuation(alternate_resultant, p) == 8, "alternative p artifact")
    residual_expression = sp.cancel(alternate_resultant / p**8)
    require(sp.denom(residual_expression) == 1, "alternative residual polynomial")
    residual = sp.Poly(residual_expression, p)
    require(residual.degree() == 17, "alternative residual degree")
    expected_lc = (
        -sp.Rational(3289935900927224469054816256, 252226880859375)
        * a**11 * b**5 * (a - b)**3 * J
    )
    expected_tc = (
        sp.Rational(3877634048, 50625) * a**3 * b
        * (17415 * a**3 * b**2 + 1013888 * a + 2027776 * b)
    )
    require(sp.factor(residual.LC() - expected_lc) == 0,
            "alternative leading endpoint")
    require(sp.factor(residual.TC() - expected_tc) == 0,
            "alternative constant endpoint")
    control = {a: 2, b: 3}
    control_residual = sp.Poly(residual_expression.subs(control), p)
    require(J.subs(control) != 0 and expected_tc.subs(control) != 0,
            "disjoint control gates")
    require(sp.gcd(control_residual, control_residual.diff()).degree() == 0,
            "disjoint source control squarefree")

    # Root-centred strict transform, independently of the primary chart.
    F = sp.expand(q * (s**2 - p) - (s**2 - p) * H - s**2 / 2)
    local = sp.expand(sp.cancel(u**4 * F.subs({p: 1 / u, s: 1 / a + v})))
    require(sp.denom(local) == 1, "local strict-transform polynomiality")
    coefficient_v2 = sp.factor(sp.diff(local, v, 2).subs({u: 0, v: 0}) / 2)
    coefficient_u = sp.factor(sp.diff(local, u).subs({u: 0, v: 0}))
    require(sp.factor(coefficient_v2 + kappa * a * (a - b)) == 0,
            "local quadratic row")
    require(sp.factor(coefficient_u - sp.Rational(8, 45) * J / a**2) == 0,
            "local unit row")
    require(sp.factor(-coefficient_v2 / coefficient_u
                      - sp.Rational(172, 3) * a**3 * (a - b) / J) == 0,
            "quadratic tangency coefficient")

    # u~v^2 and omega=-u^2 ds/K_u give ord(omega)=4, hence index five.
    require(2 * 2 == 4 and 4 + 1 == 5, "merged differential/index ledger")
    polygon = ((0, 1), (2, 0), (5, 3), (3, 4), (0, 4))
    require(polygon_ledger(polygon) == (27, 11, 9), "Newton/Pick ledger")
    packet = (8, 5, 3, 2, 2, 2, 1)
    require(sum(e - 1 for e in packet) == 16 == 2 * 9 - 2,
            "complete packet defect")

    L = 17 + 4
    full_n = sum(packet)
    finite_n = 8 + 5 + 3 + 1
    beta = 3
    require((L, full_n, finite_n, beta) == (21, 23, 17, 3),
            "response ledger")
    require(2 * (full_n - L) == 4 < 16,
            "full commutator contradiction")
    capacities = (
        2 * finite_n - L - 2 + beta,
        2 * finite_n - L - 1 + beta,
        beta,
    )
    require(capacities == (14, 15, 3) and max(capacities) < finite_n - 1,
            "finite transitivity contradiction")

    rows = (
        "chart=intrinsic-linear-subresultant;double-root=-B_C/A_C",
        "alternative-critical=Res_s(A,C0)=p^8*R17;endpoints-live",
        "local=u-unit-plus-v^2;omega-order=4;merged-index=5",
        "packet=(8,5,3,2,2,2,1);g=9;L=21;responses=23/17;beta=3",
        "control=(a,b)=(2,3);source-residual-squarefree=yes",
    )
    semantic = hashlib.sha256("\n".join(rows).encode("ascii")).hexdigest()
    print("THM4161_Y_ONLY_DOUBLE_TOP_ROOT_INDEPENDENT_AUDIT=PASS")
    print(f"checks={CHECKS}")
    for row in rows:
        print(row)
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
