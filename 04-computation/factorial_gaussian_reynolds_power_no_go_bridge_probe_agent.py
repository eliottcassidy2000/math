#!/usr/bin/env python3
"""Exact hostile for the factorial/Gaussian torus bridge.

This is a scratch companion, not maintained canon.  It verifies two facts:

1. If U,V are independent standard real Gaussians, then
   ((U^2+V^2)/2)^k has moment k!.  Hence the diagonal/radial embedding of
   the factorial functional into a standard-real-Gaussian functional is
   exact, monomial by monomial.
2. For the polynomial P of THM-3290, the U(1) Reynolds projection does not
   commute with powers already at m=2.  The projected first power has a
   strictly positive Gaussian second moment although the projection of P^2
   still has Gaussian moment zero.

Only integer arithmetic and fractions are used.
"""

from fractions import Fraction
from math import comb, factorial


Poly = dict[tuple[int, ...], int]


def add(*polys: Poly) -> Poly:
    out: Poly = {}
    for poly in polys:
        for exponent, coefficient in poly.items():
            out[exponent] = out.get(exponent, 0) + coefficient
            if out[exponent] == 0:
                del out[exponent]
    return out


def scale(poly: Poly, scalar: int) -> Poly:
    return {e: scalar * c for e, c in poly.items() if scalar * c}


def multiply(left: Poly, right: Poly) -> Poly:
    out: Poly = {}
    for e_left, c_left in left.items():
        for e_right, c_right in right.items():
            exponent = tuple(a + b for a, b in zip(e_left, e_right))
            out[exponent] = out.get(exponent, 0) + c_left * c_right
            if out[exponent] == 0:
                del out[exponent]
    return out


def power(poly: Poly, exponent: int) -> Poly:
    if exponent < 0:
        raise ValueError("negative exponent")
    one = {(0,) * len(next(iter(poly))): 1}
    out = one
    base = poly
    e = exponent
    while e:
        if e & 1:
            out = multiply(out, base)
        base = multiply(base, base)
        e //= 2
    return out


def variable(number: int, index: int) -> Poly:
    exponent = [0] * number
    exponent[index] = 1
    return {tuple(exponent): 1}


def odd_double_factorial(exponent: int) -> int:
    """Return E[G^exponent] for a standard real Gaussian G."""
    if exponent % 2:
        return 0
    if exponent == 0:
        return 1
    out = 1
    for k in range(1, exponent, 2):
        out *= k
    return out


def radial_real_gaussian_moment(exponent: int) -> Fraction:
    """E[((U^2+V^2)/2)^exponent], expanded in independent U,V."""
    numerator = 0
    for j in range(exponent + 1):
        numerator += (
            comb(exponent, j)
            * odd_double_factorial(2 * j)
            * odd_double_factorial(2 * (exponent - j))
        )
    return Fraction(numerator, 2**exponent)


def reynolds_xy(poly: Poly) -> Poly:
    """Keep x^a y^a t^c and rename xy=q; output variables are (q,t)."""
    out: Poly = {}
    for (a, b, c), coefficient in poly.items():
        if a == b:
            exponent = (a, c)
            out[exponent] = out.get(exponent, 0) + coefficient
    return {e: c for e, c in out.items() if c}


def gaussian_qt_moment(poly: Poly) -> int:
    """Expectation for q=U^2+V^2 and independent t~N(0,1)."""
    out = 0
    for (a, b), coefficient in poly.items():
        out += (
            coefficient
            * 2**a
            * factorial(a)
            * odd_double_factorial(b)
        )
    return out


def canonical(poly: Poly) -> tuple[tuple[tuple[int, ...], int], ...]:
    return tuple(sorted(poly.items(), reverse=True))


def main() -> None:
    # Exact factorial/complex-radial dictionary.
    radial_checks = tuple(
        radial_real_gaussian_moment(k) == factorial(k) for k in range(13)
    )
    if not all(radial_checks):
        raise RuntimeError("factorial/Gaussian radial identity failed")

    # THM-3290 variables and polynomial.
    x = variable(3, 0)
    y = variable(3, 1)
    t = variable(3, 2)
    rho = add(power(t, 2), multiply(x, y))
    A = add(rho, power(x, 2))
    C = add(
        multiply(y, power(rho, 2)),
        scale(multiply(multiply(x, power(t, 2)), rho), -2),
        scale(multiply(power(x, 3), power(t, 2)), -1),
    )
    P = multiply(A, power(C, 2))
    weights = tuple(sorted({a - b for (a, b, _), coefficient in P.items() if coefficient}))
    if weights != (-2, 0, 2, 4, 6, 8):
        raise RuntimeError("THM-3290 weight support failed")

    R1 = reynolds_xy(P)
    R2 = reynolds_xy(power(P, 2))

    q = variable(2, 0)
    tau = variable(2, 1)
    q_plus_t2 = add(q, power(tau, 2))
    expected_R1 = multiply(
        multiply(q, add(q, scale(power(tau, 2), -4))),
        power(q_plus_t2, 4),
    )
    expected_R2 = multiply(
        multiply(power(q, 2), power(q_plus_t2, 8)),
        add(power(q, 2), scale(multiply(q, power(tau, 2)), -20), scale(power(tau, 4), 24)),
    )
    expected_defect = scale(
        multiply(
            multiply(power(q, 2), power(tau, 2)),
            multiply(power(q_plus_t2, 8), add(scale(q, 3), scale(power(tau, 2), -2))),
        ),
        -4,
    )

    defect = add(R2, scale(power(R1, 2), -1))
    if canonical(R1) != canonical(expected_R1):
        raise RuntimeError("R(P) formula failed")
    if canonical(R2) != canonical(expected_R2):
        raise RuntimeError("R(P^2) formula failed")
    if canonical(defect) != canonical(expected_defect):
        raise RuntimeError("Reynolds power-defect formula failed")

    m_R1 = gaussian_qt_moment(R1)
    m_R2 = gaussian_qt_moment(R2)
    m_R1_squared = gaussian_qt_moment(power(R1, 2))
    m_defect = gaussian_qt_moment(defect)
    if (m_R1, m_R2, m_R1_squared, m_defect) != (
        0,
        0,
        3_212_537_328_000,
        -3_212_537_328_000,
    ):
        raise RuntimeError("Gaussian hostile moments failed")

    # A quotient cannot carry the power operation: P and R(P) have the same
    # first projection, but their squared projections differ.
    same_first_projection = canonical(reynolds_xy(P)) == canonical(R1)
    different_second_projection = canonical(R2) != canonical(power(R1, 2))
    if not (same_first_projection and different_second_projection):
        raise RuntimeError("power-congruence hostile failed")

    print("status=FINITE-EXACT SCRATCH; not canon")
    print("factorial_gaussian_radial_identity=k! for k=0..12")
    print(f"thm3290_P_terms={len(P)}")
    print("thm3290_weight_set=-2,0,2,4,6,8")
    print("R(P)=q*(q-4*t^2)*(q+t^2)^4")
    print("R(P^2)=q^2*(q+t^2)^8*(q^2-20*q*t^2+24*t^4)")
    print("power_defect=-4*q^2*t^2*(q+t^2)^8*(3*q-2*t^2)")
    print(f"E_RP={m_R1}")
    print(f"E_RP2_projection={m_R2}")
    print(f"E_RP_squared={m_R1_squared}")
    print(f"E_power_defect={m_defect}")
    print("consequence=Reynolds projection preserves expectation but not powers")
    print("hostile=the THM-3290 null sequence does not descend via P -> R(P)")


if __name__ == "__main__":
    main()
