#!/usr/bin/env python3
"""Exact arithmetic companion for THM-3070.

The theorem itself is a valuation/associated-graded argument.  This companion
checks its exponent lattice, the balances used in the rational autonomous ODE
classification, polynomial-composition contradiction, two sharp Laurent
hostiles, and a leading-cancellation local survivor.
"""

from __future__ import annotations

import math

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def check_weight_lattice(p: int, cutoff: int = 90) -> None:
    """Check all nonnegative weight-zero/three solutions in a large box."""

    r = p + 3
    g = math.gcd(p, 3)
    r_reduced = r // g
    p_reduced = p // g

    actual_zero = {
        (i, j)
        for i in range(cutoff + 1)
        for j in range(cutoff + 1)
        if -p * i + r * j == 0
    }
    expected_zero = {
        (r_reduced * k, p_reduced * k)
        for k in range(cutoff + 1)
        if r_reduced * k <= cutoff and p_reduced * k <= cutoff
    }
    require(actual_zero == expected_zero, f"weight-zero lattice for p={p}")

    actual_three = {
        (i, j)
        for i in range(cutoff + 1)
        for j in range(cutoff + 1)
        if -p * i + r * j == 3
    }
    expected_three = {
        (1 + r_reduced * k, 1 + p_reduced * k)
        for k in range(cutoff + 1)
        if 1 + r_reduced * k <= cutoff and 1 + p_reduced * k <= cutoff
    }
    require(actual_three == expected_three, f"weight-three lattice for p={p}")


def main() -> None:
    u, s = sp.symbols("u s", nonzero=True)
    p, r = sp.symbols("p r", integer=True, positive=True)
    A = sp.Function("A")(u)
    B = sp.Function("B")(u)

    # The first possible coefficient of dx wedge dy for
    # x=A*s^-p and y=B*s^r.
    x_lead = A * s ** (-p)
    y_lead = B * s**r
    wedge_lead = sp.factor(
        sp.diff(x_lead, u) * sp.diff(y_lead, s)
        - sp.diff(x_lead, s) * sp.diff(y_lead, u)
    )
    expected_wedge = (
        r * sp.diff(A, u) * B + p * A * sp.diff(B, u)
    ) * s ** (r - p - 1)
    require(sp.simplify(wedge_lead - expected_wedge) == 0, "generic leading wedge")

    # Include both gcd regimes, especially p divisible by three.
    for exponent in range(1, 25):
        check_weight_lattice(exponent)

    # Universal integer balances in the rational ODE proof.  A polynomial C
    # would require n-1=(d+1)n.  A finite pole requires
    # m+1=m(d+1), whose sole positive solution is (m,d)=(1,1).
    polynomial_balances = [
        (n, d)
        for n in range(1, 80)
        for d in range(0, 80)
        if n - 1 == (d + 1) * n
    ]
    require(polynomial_balances == [], "no polynomial autonomous solution")
    pole_balances = [
        (m, d)
        for m in range(1, 80)
        for d in range(0, 80)
        if m + 1 == m * (d + 1)
    ]
    require(pole_balances == [(1, 1)], "unique pole/degree balance")
    zero_linear_balances = [n for n in range(1, 80) if n - 1 == n]
    zero_quadratic_balances = [n for n in range(1, 80) if n - 1 == 2 * n]
    require(zero_linear_balances == [], "no finite zero with nonzero linear term")
    require(zero_quadratic_balances == [], "no finite zero with pure quadratic term")
    infinity_linear_balances = [n for n in range(1, 80) if n + 1 == n]
    infinity_quadratic_balances = [n for n in range(1, 80) if n + 1 == 2 * n]
    require(infinity_linear_balances == [], "linear term impossible at infinity")
    require(infinity_quadratic_balances == [1], "one simple pole globally")

    # The unique rational normal form when H(z)=a*z.
    lam, a, u0 = sp.symbols("lam a u0", nonzero=True)
    C = -1 / (lam * a * (u - u0))
    H_of_C = a * C
    require(sp.factor(sp.diff(C, u) - lam * C * H_of_C) == 0, "ODE normal form")

    # No polynomial G(C(u)) can equal u.  After putting G(C)-u over its
    # reduced denominator of degree d, the numerator has degree d+1 and its
    # leading coefficient is the negative of the denominator's.
    z = sp.symbols("z")
    for degree in range(0, 11):
        coefficients = sp.symbols(f"g0:{degree + 1}")
        G = sum(coefficients[k] * z**k for k in range(degree + 1))
        numerator, denominator = sp.fraction(sp.cancel(G.subs(z, C) - u))
        numerator_poly = sp.Poly(numerator, u)
        require(numerator_poly.degree() == degree + 1, f"G(C)-u degree {degree}")
        denominator_poly = sp.Poly(denominator, u)
        require(
            sp.factor(numerator_poly.LC() + denominator_poly.LC()) == 0,
            f"G(C)-u leading coefficient {degree}",
        )
        require(denominator_poly.degree() == degree, f"G(C)-u denominator {degree}")

    # Hostile 1 (THM-3068): P is polynomial and Q is Laurent.  Its inverse
    # has p=1,r=4 and a nonzero leading wedge.
    x, y = sp.symbols("x y", nonzero=True)
    P_one = x**4 * y / 3
    Q_one = x ** -3
    jac_one = sp.factor(
        sp.diff(P_one, x) * sp.diff(Q_one, y)
        - sp.diff(P_one, y) * sp.diff(Q_one, x)
    )
    require(jac_one == 1, "P-polynomial/Q-Laurent hostile Jacobian")
    x_one = s ** -1
    y_one = 3 * u * s**4
    require(sp.simplify(P_one.subs({x: x_one, y: y_one}) - u) == 0, "hostile one P")
    require(sp.simplify(Q_one.subs(x, x_one) - s**3) == 0, "hostile one Q")

    # Hostile 2: Q is polynomial and P is Laurent.  This is the p=2,r=5
    # dual boundary, again with Jacobian one.
    P_two = 9 / (x**5 * y**2)
    Q_two = -(x**6 * y**3) / 27
    jac_two = sp.factor(
        sp.diff(P_two, x) * sp.diff(Q_two, y)
        - sp.diff(P_two, y) * sp.diff(Q_two, x)
    )
    require(jac_two == 1, "P-Laurent/Q-polynomial hostile Jacobian")
    x_two = u ** -1 * s ** -2
    y_two = -3 * u**2 * s**5
    require(sp.simplify(P_two.subs({x: x_two, y: y_two}) - u) == 0, "hostile two P")
    require(sp.simplify(Q_two.subs({x: x_two, y: y_two}) - s**3) == 0, "hostile two Q")

    # A formal leading-cancellation control.  It satisfies the local Keller
    # two-form exactly, but it is not asserted to globalize polynomially.
    x_cancel = u * s ** -1
    y_cancel = u ** -1 * s + sp.Rational(3, 4) * s**4
    wedge_cancel = sp.factor(
        sp.diff(x_cancel, u) * sp.diff(y_cancel, s)
        - sp.diff(x_cancel, s) * sp.diff(y_cancel, u)
    )
    require(wedge_cancel == 3 * s**2, "cancellation control wedge")
    leading_cancel = sp.factor(1 * sp.diff(u, u) * u ** -1 + 1 * u * sp.diff(u ** -1, u))
    require(leading_cancel == 0, "leading coefficient cancellation")

    print("theorem=THM-3070")
    print("status=PROVED_VERIFIED_EXACT_CANDIDATE")
    print("wedge_lead=(r*A_prime*B+p*A*B_prime)*s^(r-p-1)")
    print("nondegenerate_C3_order=2_forces=r-p=3")
    print("weight0=x^((p+3)/g) y^(p/g) powers;weight3=xy_times_weight0_powers")
    print("p_divisible_by_3_checked=p=3,6,9,12,15,18,21,24;g=3")
    print("autonomous_ODE=C_prime=(3*kappa^-1/g)*C*H(C)")
    print("rational_ODE_only=H(z)=a*z;C=-1/(lambda*a*(u-u0))")
    print("polynomial_P_contradiction=u=G(C);G_polynomial;C_tends_0_at_u_infinity")
    print("hostile1=P_polynomial,Q_Laurent,p=1,r=4,Jacobian=1")
    print("hostile2=P_Laurent,Q_polynomial,p=2,r=5,Jacobian=1")
    print("cancellation_control=p=1,r=1;A*B=1;dx_wedge_dy=du_wedge_d(s^3)")
    print("necessary_survivor=r<=p+2;A^r*B^p=constant;not_sufficient")
    print("scope=no_full_C3_A4_S4_or_JC2_exclusion")


if __name__ == "__main__":
    main()
