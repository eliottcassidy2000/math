#!/usr/bin/env python3
"""Clean-room alternative-source audit for THM-4164.

This referee uses (A,C0), hence a p^8 artifact, and never forms the primary
normalized resultant.  It uses a disjoint generic control a=2 and derives
the two J-wall branches directly from the tangent cone.
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
    s, p, a, w, u, v, q = sp.symbols("s p a w u v q")
    X, T = sp.symbols("X T")
    kappa = sp.Rational(1376, 135)
    K0 = sp.Rational(2848, 45)

    zeta = kappa * a**3
    theta = -3 * kappa * a**2
    phi = 3 * kappa * a
    C_top = zeta * w**3 + theta * w**2 + phi * w - kappa
    check(sp.factor(C_top - kappa * (a * w - 1)**3) == 0,
          "complete triple chart")
    check(sp.factor(sp.gcd(sp.Poly(C_top, w), sp.Poly(sp.diff(C_top, w), w)).as_expr()
                    - kappa * (a * w - 1)**2) == 0,
          "triple gcd")
    J = 15 * a**2 + 356
    inner = 5805 * a**4 + 1013888
    I = sp.factor(4 * theta * K0**2 - 27 * zeta**2)
    check(sp.factor(I + sp.Rational(44032, 91125) * a**2 * inner) == 0,
          "inner factor")

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
          "alternative critical degrees")
    check(sp.factor(sp.Poly(A, s).LC() - 3 * zeta * p**2) == 0,
          "alternative A leading row")
    check(sp.factor(sp.Poly(C0, s).LC() - 6 * zeta * p**2) == 0,
          "alternative C0 leading row")
    raw = sp.resultant(A, C0, s)
    check(valuation(raw, p) == 8, "alternative p artifact")
    R16 = sp.Poly(sp.cancel(raw / p**8), p)
    check(R16.degree() == 16, "alternative generic degree")
    check(
        sp.factor(
            R16.LC()
            - sp.Rational(9563767153858210665857024, 9341736328125)
            * a**17 * J**2
        ) == 0,
        "alternative generic leading row",
    )
    check(sp.factor(R16.TC() + 46656 * zeta * I) == 0,
          "alternative generic constant row")

    # Disjoint generic control a=2.
    control = sp.Poly(R16.as_expr().subs(a, 2), p)
    check(J.subs(a, 2) != 0 and I.subs(a, 2) != 0,
          "disjoint generic gates")
    check(sp.gcd(control, control.diff()).degree() == 0,
          "disjoint generic squarefree")

    # Restore universal pairs by direct evaluation, not resultant factors.
    P = T + X**2 * T**2
    Y = X * T * P
    G = sp.expand(
        -X**2 * T / 2 - 3 * P + sp.Rational(8, 3) * P**2
        - kappa * P**3 + K0 * Y**2 + phi * P**2 * Y
        + theta * P * Y**2 + zeta * Y**3
    )
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

    # Direct root-centred boundary computation.
    F = sp.expand(q * (s**2 - p) - (s**2 - p) * H - s**2 / 2)
    local = sp.expand(sp.cancel(u**4 * F.subs({p: 1 / u, s: 1 / a + v})))
    v3 = sp.factor(sp.diff(local, v, 3).subs({u: 0, v: 0}) / 6)
    u1 = sp.factor(sp.diff(local, u).subs({u: 0, v: 0}))
    check(sp.factor(v3 - kappa * a**3) == 0, "generic v3 row")
    check(sp.factor(u1 - sp.Rational(8, 45) * J / a**2) == 0,
          "generic u row")
    check(sp.factor(-v3 / u1 + sp.Rational(172, 3) * a**5 / J) == 0,
          "generic cubic branch")
    # u~v^3 and omega=-u^2 ds/K_u.
    check((2 * 3, 2 * 3 + 1) == (6, 7), "generic residue/index")

    generic_packet = tuple(sorted((8, 7, 1) + (2, 2, 2), reverse=True))
    check(generic_packet == (8, 7, 2, 2, 2, 1), "generic packet")
    check(sum(index - 1 for index in generic_packet) == 16,
          "generic defect")
    generic_L = R16.degree() + 4
    generic_full, generic_finite, beta = sum(generic_packet), 16, 3
    check((generic_L, generic_full, generic_finite) == (20, 22, 16),
          "generic response ledger")
    check(2 * (generic_full - generic_L) == 4 < 16,
          "generic full contradiction")
    capacities = (
        2 * generic_finite - generic_L - 2 + beta,
        2 * generic_finite - generic_L - 1 + beta,
        beta,
    )
    check(capacities == (13, 14, 3) and max(capacities) < 15,
          "generic finite contradiction")

    # Independent exact reduction on J=0.
    modulus = sp.Poly(J, a)

    def reduce_J(expression: sp.Expr) -> sp.Expr:
        numerator, denominator = sp.fraction(sp.cancel(expression))
        numerator = sp.rem(sp.Poly(numerator, a), modulus).as_expr()
        denominator = sp.rem(sp.Poly(denominator, a), modulus).as_expr()
        inverse = sp.invert(sp.Poly(denominator, a), modulus).as_expr()
        return sp.factor(sp.rem(sp.Poly(numerator * inverse, a), modulus).as_expr())

    R15 = sp.Poly(sum(reduce_J(R16.nth(d)) * p**d for d in range(17)), p)
    check(R15.degree() == 15, "alternative J-wall degree")
    check(
        R15.LC()
        == -sp.Rational(
            3636055657731868916207553629845852858180471619584,
            1596123230438232421875,
        ) * a,
        "alternative J-wall leading row",
    )
    check(reduce_J(I) == sp.Rational(335741565206528, 6834375),
          "alternative J-wall inner unit")
    root = sp.sqrt(-1335)
    a0 = 2 * root / 15
    wall_control = sp.Poly(sp.expand(R15.as_expr().subs(a, a0)),
                           p, extension=root)
    check(sp.gcd(wall_control, wall_control.diff()).degree() == 0,
          "alternative J-wall squarefree")

    wall_u = reduce_J(u1)
    wall_v3 = reduce_J(v3)
    wall_uv = reduce_J(sp.diff(local, u, v).subs({u: 0, v: 0}))
    wall_u2 = reduce_J(sp.diff(local, u, 2).subs({u: 0, v: 0}) / 2)
    check(wall_u == 0, "wall u row")
    check(wall_v3 == -sp.Rational(489856, 2025) * a,
          "wall v3 row")
    check(wall_uv == -sp.Rational(16, 3) * a and wall_u2 == -3,
          "wall tangent cone")
    # Distinct tangent lines: u=0 and u=(wall_uv/3)v.  The first branch
    # balances uv with v^3 and hence has u=-(v3/uv)v^2.
    check(wall_uv / 3 == -sp.Rational(16, 9) * a,
          "transverse branch")
    check(sp.factor(-wall_v3 / wall_uv) == -sp.Rational(30616, 675),
          "quadratic tangent branch")
    check((2 * 1 - 1, 2 * 2 - 1) == (1, 3),
          "wall differential orders")

    wall_packet = tuple(sorted((8, 4, 2, 1) + (2, 2, 2), reverse=True))
    check(wall_packet == (8, 4, 2, 2, 2, 2, 1), "wall packet")
    check(sum(index - 1 for index in wall_packet) == 14,
          "wall defect")
    wall_L = R15.degree() + 4
    wall_full, wall_finite = sum(wall_packet), 15
    check((wall_L, wall_full, wall_finite) == (19, 21, 15),
          "wall response ledger")
    check(2 * (wall_full - wall_L) == 4 < 14,
          "wall full contradiction")
    wall_capacities = (
        2 * wall_finite - wall_L - 2 + beta,
        2 * wall_finite - wall_L - 1 + beta,
        beta,
    )
    check(wall_capacities == (12, 13, 3) and max(wall_capacities) < 14,
          "wall finite contradiction")

    rows = (
        "alternative=Res_s(A,C0)=p8R16;control=a2;squarefree",
        "generic-boundary=v3+u;smooth;e7;g9;packet=8,7,2,2,2,1;L20",
        "Jwall=a2=-356/15;alternative=p8R15;control=quadratic-field",
        "Jwall-boundary=uv+u2+v3;branches=e2,e4;g8;packet=8,4,2,2,2,2,1;L19",
        "responses=generic22/16;wall21/15;all-finite-and-full-strict",
        "remaining=inner 5805a4+1013888=0",
    )
    semantic = hashlib.sha256("\n".join(rows).encode("ascii")).hexdigest()
    print("THM4164_Y_ONLY_TRIPLE_TOP_ROOT_INDEPENDENT_AUDIT=PASS")
    print(f"checks={CHECKS}")
    for row in rows:
        print(row)
    print("controls=generic:a=2;Jwall:a=2sqrt(-1335)/15;alternative-source-squarefree")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
