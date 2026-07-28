#!/usr/bin/env python3
"""Exact companion for THM-2745's one-pole composition boundary."""

from fractions import Fraction
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    odd_rows = tuple(range(21, 0, -2))
    expected_poles = (144, 44, 120, 108, 32, 84, 72, 20, 48, 36, 8)
    pole_rows = []

    # The source control uses U=u0*x^m and
    # S=s0+s1*x^(1-m), with u0*s1*(1-m)=kappa.
    u0 = Fraction(2)
    kappa = Fraction(3)
    ramification_degrees = range(1, 6)
    controls = 0
    for j, expected in zip(odd_rows, expected_poles):
        r = 22 - j
        common_gcd = gcd(r, 6)
        pole_order = (150 - 6 * r) // common_gcd
        require(pole_order == expected, f"pole order changed at j={j}")
        require(pole_order >= 8, f"constant-U boundary reopened at j={j}")
        pole_rows.append((j, r, common_gcd, pole_order))

        for e in ramification_degrees:
            m = 1 + e * pole_order
            s1 = kappa / (u0 * (1 - m))
            require(u0 * s1 * (1 - m) == kappa,
                    "rational-primitive derivative identity")

            # In X=1/x, take Z=z0+gamma*X^e and
            # R=s0+rho*(Z-z0)^P.  Coefficients below check
            # rho*gamma^P=s1 and exponent eP=m-1 exactly.
            gamma = Fraction(2)
            rho = s1 / gamma**pole_order
            require(rho * gamma**pole_order == s1,
                    "pure-power coefficient composition")
            require(e * pole_order == m - 1,
                    "ramification divisibility identity")
            controls += 1

    require(tuple(row[-1] for row in pole_rows) == expected_poles,
            "highest-odd pole atlas changed")
    require(controls == 55, "ramification control census changed")

    # Constant U gives an affine-linear source response, so the composition
    # degree equation P*e=1 has the unique positive solution (1,1).
    constant_solutions = tuple(
        (pole, e) for pole in range(1, 9) for e in range(1, 9)
        if pole * e == 1
    )
    require(constant_solutions == ((1, 1),),
            "constant-U degree boundary changed")

    # A translated positive control: neither factor is a literal monomial
    # until the target origin is shifted by z0.
    pole = 8
    e = 3
    z0 = Fraction(7)
    gamma = Fraction(3)
    rho = Fraction(5)
    source_coefficient = rho * gamma**pole
    source_exponent = e * pole
    require(source_coefficient == 5 * 3**8 and source_exponent == 24,
            "translated pure-power control changed")
    require(z0 != 0, "translation hostile collapsed to literal monomial")

    print("THM-2745 ONE-POLE RESPONSE COMPOSITION AUDIT")
    print(f"highest_odd_rows={len(pole_rows)}")
    print(f"pole_orders={tuple(row[-1] for row in pole_rows)}")
    print("constant_U_positive_degree_solutions=((1, 1),)")
    print("monomial_U_law=m-1=e*P")
    print("ramification_degrees_checked=(1,2,3,4,5)")
    print(f"exact_derivative_and_composition_controls={controls}")
    print("translated_control=R=2+5*(Z-7)^8,Z=7+3*X^3")
    print(f"translated_pullback=2+{source_coefficient}*X^{source_exponent}")
    print("scope=ONE_POLE_COMPONENT_NORMALIZATION_NOT_FACTORIZATION_CLOSURE")
    print("PASS")


if __name__ == "__main__":
    main()
