#!/usr/bin/env python3
"""Exact weighted-face referee for THM-2741.

For each possible highest nonzero odd Faber seed in the full split degree-22
family, this script reconstructs the three quartic Faber observables, checks
the weighted Newton separation from every permitted lower column, and proves
that the restricted homogeneous response face is nonzero on all six G2
tangent lines.  It is an exact symbolic companion, not a JC(2) solver.
"""

from __future__ import annotations

from math import gcd

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def faber_coefficients(
    degree: int, p: sp.Expr, q: sp.Expr, r: sp.Expr, extra: int = 3
) -> list[sp.Expr]:
    exponent = sp.Rational(degree, 4)
    coefficients = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + extra + 1):
        value = sum(
            quartic[step]
            * ((exponent + 1) * step - index)
            * coefficients[index - step]
            for step in range(2, min(4, index) + 1)
        ) / index
        coefficients.append(sp.factor(value))
    return coefficients


def observables(
    degree: int, d: sp.Symbol, q: sp.Symbol, s: sp.Symbol
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    coefficients = faber_coefficients(degree, 2 * d, q, d**2 - s)
    phi = sp.factor(4 * coefficients[degree + 1])
    psi = sp.factor(4 * coefficients[degree + 2])
    response = sp.factor(
        4 * coefficients[degree + 3] + 2 * d * coefficients[degree + 1]
    )
    return phi, psi, response


def local_order(expression: sp.Expr, q: sp.Symbol, s: sp.Symbol) -> int:
    terms = sp.Poly(sp.expand(expression), q, s).terms()
    return min(sum(monomial) for monomial, _ in terms)


def local_face(expression: sp.Expr, q: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    order = local_order(expression, q, s)
    return sp.factor(
        sum(
            coefficient * q**monomial[0] * s**monomial[1]
            for monomial, coefficient in sp.Poly(sp.expand(expression), q, s).terms()
            if sum(monomial) == order
        )
    )


def main() -> None:
    d, q, s = sp.symbols("d q s")
    odd_degrees = list(range(21, 0, -2))
    even_degrees = [14, 10, 6, 2]
    allowed_degrees = odd_degrees + even_degrees
    rows = {degree: observables(degree, d, q, s) for degree in allowed_degrees + [22]}

    phi22, psi22, response22 = (entry.subs(d, 1) for entry in rows[22])
    phi_face = local_face(phi22, q, s)
    psi_face = local_face(psi22, q, s)
    response_face = local_face(response22, q, s)
    cross_face = q * s * (q**2 - 3 * s**2) * (3 * q**2 - s**2)
    g2_face = (q**2 - s**2) * (q**4 - 14 * q**2 * s**2 + s**4)
    require(local_order(phi22, q, s) == 6, "Phi22 order changed")
    require(local_order(psi22, q, s) == 6, "Psi22 order changed")
    require(local_order(response22, q, s) == 6, "R22 order changed")
    require(
        sp.expand(phi_face + sp.Rational(231, 128) * cross_face) == 0,
        "Phi22 face changed",
    )
    require(
        sp.expand(psi_face + sp.Rational(231, 256) * g2_face) == 0,
        "Psi22 face changed",
    )
    require(
        sp.expand(response_face - sp.Rational(231, 256) * cross_face) == 0,
        "R22 face changed",
    )

    print("DEGREE22 HIGHEST-ODD RESPONSE POLE ATLAS")
    print("initial_system=Psi22_6=0;Phi22_6+a_j*Phi_j(1,0,0)*h^r=0")
    print("j,r,g,coarse_branches,response_face_coefficient,response_pole_order,min_weight_margin")

    for selected in odd_degrees:
        r = 22 - selected
        common_gcd = gcd(r, 6)
        phi_selected, psi_selected, response_selected = rows[selected]
        phi_zero = sp.factor(phi_selected.subs({d: 1, q: 0, s: 0}))
        response_zero = sp.factor(response_selected.subs({d: 1, q: 0, s: 0}))
        require(phi_zero != 0, f"Phi_{selected}(1,0,0) vanished")
        require(response_zero != 0, f"R_{selected}(1,0,0) vanished")
        require(
            local_order(psi_selected.subs(d, 1), q, s) >= 1,
            f"Psi_{selected} acquired a constant term",
        )

        restricted_response = sp.factor(
            response_face - (response_zero / phi_zero) * phi_face
        )
        quotient = sp.factor(restricted_response / cross_face)
        require(sp.denom(quotient) != 0 and quotient != 0, "response coefficient vanished")
        require(
            sp.expand(restricted_response - quotient * cross_face) == 0,
            "response face is not the expected cross arrangement",
        )
        require(
            sp.gcd(sp.Poly(restricted_response, q, s), sp.Poly(psi_face, q, s))
            == 1,
            f"response face meets a Psi tangent for j={selected}",
        )

        # Use weights v(h)=6/g and v(q)=v(s)=r/g, with the common 1/g
        # suppressed.  Both top faces and the chosen odd Phi column have
        # weight 6r.  Every other permitted contribution must be strictly
        # larger.  Higher odd columns are zero by the definition of
        # 'highest odd seed'.
        initial_weight = 6 * r
        margins: list[int] = []
        for degree in allowed_degrees:
            if degree % 2 == 1 and degree > selected:
                continue
            gap = 22 - degree
            phi, psi, response = (entry.subs(d, 1) for entry in rows[degree])
            for name, expression in (
                ("Phi", phi),
                ("Psi", psi),
                ("R", response),
            ):
                weight = 6 * gap + r * local_order(expression, q, s)
                if degree == selected and name in {"Phi", "R"}:
                    require(weight == initial_weight, "chosen face weight changed")
                    continue
                require(
                    weight > initial_weight,
                    f"{name}_{degree} intrudes for highest odd j={selected}",
                )
                margins.append(weight - initial_weight)
        margins.extend([6 * 23 - initial_weight, 6 * 24 - initial_weight])
        require(min(margins) > 0, "flux constant intrudes into initial face")

        coarse_branches = 3 * common_gcd
        response_pole_order = (150 - 6 * r) // common_gcd
        require(coarse_branches >= 3, "coarse branch floor failed")
        require(response_pole_order > 0, "response ceased to have a pole")
        print(
            f"{selected},{r},{common_gcd},{coarse_branches},"
            f"{quotient},{response_pole_order},{min(margins)}"
        )

    print("all_eleven_odd_seeds=three_or_nine_coarse_response_poles")
    print("source_response_capacity=at_most_one_pole_by_rational_primitive")
    print("scope=INTEGRAL_MEMBERS_ONLY_REDUCIBLE_NONREDUCED_AND_JC2_REMAIN_OPEN")
    print("PASS")


if __name__ == "__main__":
    main()
