#!/usr/bin/env python3
"""Exact companion for THM-3074.

Checks the two-pole leading cancellation, unimodular toric charts, first
key-depth invoice, first-stage gcd lattice, and sharp local controls which
separate leading detection from later cancellation.
"""

from __future__ import annotations

import math

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def bezout_negative(a: int, b: int) -> tuple[int, int]:
    """Return c,d with a*c+b*d=-1 for coprime positive a,b."""

    # Use a tiny exact search to avoid any library return-convention issue.
    require(math.gcd(a, b) == 1, "Bezout inputs must be coprime")
    for c in range(-b, b + 1):
        numerator = -1 - a * c
        if numerator % b == 0:
            return c, numerator // b
    raise RuntimeError("Bezout search failed")


def bezout_positive(a: int, b: int) -> tuple[int, int]:
    """Return c,d with a*d-b*c=1 for coprime nonnegative a,b."""

    require(math.gcd(a, b) == 1, "Bezout inputs must be coprime")
    if a == 0:
        require(b == 1, "primitive zero case")
        return -1, 0
    for c in range(-a, a + 1):
        numerator = 1 + b * c
        if numerator % a == 0:
            return c, numerator // a
    raise RuntimeError("positive Bezout search failed")


def laurent_order(expr: sp.Expr, variable: sp.Symbol) -> int:
    numerator, denominator = sp.fraction(sp.cancel(expr))
    numerator_poly = sp.Poly(numerator, variable)
    denominator_poly = sp.Poly(denominator, variable)
    require(not numerator_poly.is_zero, "zero has no finite Laurent order")
    numerator_order = min(monomial[0] for monomial, coefficient in numerator_poly.terms() if coefficient)
    denominator_order = min(
        monomial[0] for monomial, coefficient in denominator_poly.terms() if coefficient
    )
    return numerator_order - denominator_order


def main() -> None:
    u, s = sp.symbols("u s", nonzero=True)
    p, q = sp.symbols("p q", integer=True, positive=True)
    A = sp.Function("A")(u)
    B = sp.Function("B")(u)

    # Both-coordinate escape: the first wedge term lies far below order two
    # and therefore has to cancel.
    x_lead = A * s ** (-p)
    y_lead = B * s ** (-q)
    wedge_lead = sp.factor(
        sp.diff(x_lead, u) * sp.diff(y_lead, s)
        - sp.diff(x_lead, s) * sp.diff(y_lead, u)
    )
    expected_wedge = (
        p * A * sp.diff(B, u) - q * sp.diff(A, u) * B
    ) * s ** (-p - q - 1)
    require(sp.simplify(wedge_lead - expected_wedge) == 0, "two-pole leading wedge")

    # Exact integer toric charts in both the two-pole and one-pole regimes.
    # The loops include non-coprime p,q and p,r, plus r=0.
    for p_value in range(1, 25):
        for q_value in range(1, 25):
            h = math.gcd(p_value, q_value)
            a = p_value // h
            b = q_value // h
            c, d = bezout_negative(a, b)
            require(a * c + b * d == -1, "two-pole Bezout identity")
            determinant = (-b) * d - a * c
            require(determinant == 1, "two-pole toric determinant")
            r_value = -p_value * c - q_value * d
            require(r_value == h, "two-pole complementary value")

        for r_source in range(0, 25):
            h = math.gcd(p_value, r_source)
            a = r_source // h
            b = p_value // h
            c, d = bezout_positive(a, b)
            require(a * d - b * c == 1, "one-pole Bezout identity")
            r_value = -p_value * c + r_source * d
            require(r_value == h, "one-pole complementary value")

    # First-stage lattice arithmetic.  A detected value three must lie in
    # h*Z+ell*Z; this is exactly gcd(h,ell)|3.  No converse is asserted.
    for h in range(1, 30):
        for ell in range(1, 40):
            divisor = math.gcd(h, ell)
            detected_three = any(
                j * h + n * ell == 3
                for j in range(-50, 51)
                for n in range(0, 51)
            )
            require(detected_three == (3 % divisor == 0), f"first lattice h={h},ell={ell}")

    # Equality-depth two-pole hostile: p=q=1, h=1, D=5, ell=5.
    # It is the punctured rational Keller map
    # P=x^4(y-x)/3, Q=x^-3, with polynomial P and Laurent Q.
    x, y = sp.symbols("x y", nonzero=True)
    x_equal = s ** -1
    y_equal = s ** -1 + 3 * u * s**4
    wedge_equal = sp.factor(
        sp.diff(x_equal, u) * sp.diff(y_equal, s)
        - sp.diff(x_equal, s) * sp.diff(y_equal, u)
    )
    require(wedge_equal == 3 * s**2, "equality-depth hostile wedge")
    M_equal = sp.cancel(y_equal / x_equal)
    R_equal = sp.cancel(1 / x_equal)
    require(sp.expand(M_equal - (1 + 3 * u * s**5)) == 0, "equality-depth M")
    require(R_equal == s, "equality-depth R")
    require(laurent_order(M_equal - 1, s) == 5, "equality depth ell=5")
    P_equal = x**4 * (y - x) / 3
    Q_equal = x ** -3
    jac_equal = sp.factor(
        sp.diff(P_equal, x) * sp.diff(Q_equal, y)
        - sp.diff(P_equal, y) * sp.diff(Q_equal, x)
    )
    require(jac_equal == 1, "equality hostile Jacobian")
    require(sp.simplify(P_equal.subs({x: x_equal, y: y_equal}) - u) == 0, "equality hostile P")
    require(sp.simplify(Q_equal.subs(x, x_equal) - s**3) == 0, "equality hostile Q")

    # Strict-depth hostile: p=q=2, h=2, D=7, ell=4.  Its first key is a
    # function of R and makes no wedge contribution; the later correction
    # supplies the exact order-two physical wedge.
    R_strict = u * s**2
    M_strict = 1 + R_strict**2 - 3 * u**3 * s**7
    x_strict = sp.cancel(1 / R_strict)
    y_strict = sp.cancel(M_strict / R_strict)
    wedge_strict = sp.factor(
        sp.diff(x_strict, u) * sp.diff(y_strict, s)
        - sp.diff(x_strict, s) * sp.diff(y_strict, u)
    )
    require(wedge_strict == 3 * s**2, "strict-depth hostile wedge")
    require(laurent_order(x_strict, s) == -2, "strict hostile x order")
    require(laurent_order(y_strict, s) == -2, "strict hostile y order")
    require(laurent_order(M_strict - 1, s) == 4, "strict hostile ell=4")
    first_m = u**2
    first_coefficient = sp.factor(
        2 * sp.diff(first_m, u) - 4 * first_m * sp.diff(u, u) / u
    )
    require(first_coefficient == 0, "strict first-key coefficient cancellation")
    require(math.gcd(2, 4) == 2 and 3 % 2 != 0, "off-lattice gcd control")

    # The polynomial T=x(y-x)-1 has valuation three despite the first
    # lattice 2Z.  Its value is exposed only after Z/R^2 loses its leading 1.
    T_later = sp.expand(x * (y - x) - 1)
    T_later_branch = sp.factor(T_later.subs({x: x_strict, y: y_strict}))
    require(T_later_branch == -3 * u * s**3, "later cancellation polynomial")
    require(laurent_order(T_later_branch, s) == 3, "later off-lattice value")
    ratio_later = sp.factor((M_strict - 1) / R_strict**2 - 1)
    require(ratio_later == -3 * u * s**3, "toric later-cancellation identity")

    # Cross-check the one-pole equality packet from THM-3070 in the same
    # first-key coordinates: p=r=h=1, D=ell=3, M=xy, R=y.
    x_one = u * s ** -1
    y_one = u ** -1 * s + sp.Rational(3, 4) * s**4
    M_one = sp.factor(x_one * y_one)
    R_one = y_one
    wedge_one = sp.factor(
        sp.diff(x_one, u) * sp.diff(y_one, s)
        - sp.diff(x_one, s) * sp.diff(y_one, u)
    )
    require(
        sp.expand(M_one - (1 + sp.Rational(3, 4) * u * s**3)) == 0,
        "one-pole M",
    )
    require(laurent_order(M_one - 1, s) == 3, "one-pole ell=3")
    require(laurent_order(R_one, s) == 1, "one-pole R order")
    require(wedge_one == 3 * s**2, "one-pole exact wedge")

    print("theorem=THM-3074")
    print("status=PROVED_VERIFIED_EXACT_CANDIDATE")
    print("two_pole_lead=(p*A*B_prime-q*A_prime*B)*s^(-p-q-1);must_cancel")
    print("two_pole_binomial=Y^(p/h)-c*X^(q/h);h=gcd(p,q)")
    print("one_pole_binomial=X^(r/h)*Y^(p/h)-c;inherited_from_THM3070")
    print("toric_chart=det(M,R)=1;val(M)=0;val(R)=h;Z=M/c-1;ell=val(Z)")
    print("depth_bound=1<=ell<=D;D=p+q+3_two_pole_or_p+3-r_one_pole")
    print("ell_eq_D:AB*(h*m_prime-ell*m*R0_prime/R0)=3*kappa^-1*tau")
    print("ell_lt_D:first_coefficient=0;m^h/R0^ell=constant;second_key_required")
    print("first_lattice=h*Z+ell*Z=gcd(h,ell)*Z")
    print("detected_target_value_3_implies=gcd(h,ell)_divides_3")
    print("off_lattice_hostile=h=2;ell=4;polynomial_x(y-x)-1_has_value_3_after_cancellation")
    print("equality_hostile=p=q=1;ell=D=5;P_polynomial;Q_Laurent;Jacobian=1")
    print("scope=coordinate_line_only;no_value_semigroup_or_full_C3_A4_S4_JC2_claim")


if __name__ == "__main__":
    main()
