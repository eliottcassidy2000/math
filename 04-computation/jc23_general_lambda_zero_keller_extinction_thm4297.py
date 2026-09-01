#!/usr/bin/env python3
"""Exact symbolic certificate for THM-4297.

Universe: the full exact-M=12 reduced (2,3) seam on
Lambda=U+W+Z=0 and U*Z*D!=0.  The script verifies the general A23 contact,
the central generic-point unit, transport of the complete repeated-face
critical ladder from THM-4292, and every valuation/form-order inequality used
to extinguish exceptional components.  The geometric extension and
proper-flat degree argument are proved in the theorem text.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise AssertionError(label)


def require_zero(value: sp.Expr, label: str) -> None:
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def coefficient(expr: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Expr:
    return sp.expand(expr).coeff(variable, degree)


def main() -> None:
    t, u, y, X, S, P, r = sp.symbols("t u y X S P r")
    U = sp.symbols("U", nonzero=True)
    W = sp.symbols("W")
    Z = -U - W
    delta = 2 * U + W
    c, alpha, upsilon, eta, Delta = sp.symbols(
        "c alpha_11 upsilon_5 eta Delta", nonzero=True
    )
    a0 = sp.Rational(7168, 135)
    k0 = sp.Rational(1376, 135)

    # The discriminant gate is exactly the local transverse coefficient.
    D = W**2 - 4 * U * Z
    require_zero(D - delta**2, "Lambda-wall discriminant square")
    top = U * r**6 + W * r**5 + Z * r**4
    require_zero(top - r**4 * (r - 1) * (U * r - Z), "top factorization")
    require_zero(top.subs(r, 1), "single contact support")
    require_zero(sp.diff(top, r).subs(r, 1) - delta, "contact derivative")

    # The raw central component has one length-twelve contact with R and is
    # transverse there because D!=0 means delta!=0.
    b = sp.symbols("b")
    central_toric = b**12 - U * r**6 - W * r**5 - Z * r**4
    require_zero(central_toric.subs(r, 1) - b**12, "contact algebra")
    require_zero(
        sp.diff(central_toric, r).subs({r: 1, b: 0}) + delta,
        "central branch derivative",
    )

    # Exact central order-nine unit test.  Modulo C, G_P=R*C_P.  On Lambda=0
    # the rational factor R is a unit on C, and C_P is generically a unit.
    central = 1 - U * P**6 - W * S**2 * P**5 - Z * S**4 * P**4
    central_p = sp.diff(central, P)
    require_zero(central.subs(P, S**2) - 1, "R restricted to C")
    domain = sp.QQ.frac_field(U, W, S)
    gcd = sp.gcd(
        sp.Poly(central, P, domain=domain),
        sp.Poly(central_p, P, domain=domain),
    )
    require(gcd.degree() == 0, "C_P is a central generic-point unit")
    uniform_resultant = sp.factor(
        sp.resultant(central.subs(S, 0), central_p.subs(S, 0), P)
    )
    require_zero(
        uniform_resultant - 46656 * U**6,
        "uniform S=0 central derivative resultant",
    )

    # Repetition forces c1=...=c5=0 but not the individual paired rows.
    # Substitute q=t^6*y.  Every W-dependent correction beyond delta first
    # occurs at t^6, after the complete critical ladder.
    r_scaled = 1 + t**6 * y
    hhat = (
        U * r_scaled**6
        + W * r_scaled**5
        + Z * r_scaled**4
        + alpha * t * (r_scaled**5 - r_scaled**4)
        + upsilon * t**2 * (r_scaled**5 - r_scaled**4)
        + eta * t**3 * (r_scaled**4 - r_scaled**3)
        + Delta * t**4 * (r_scaled**4 - r_scaled**3)
        + t**6 * (-k0 * r_scaled**3 + (c + k0) * r_scaled**2)
        + sp.Rational(8, 3) * t**8 * r_scaled**2
        - 3 * t**10 * r_scaled
    )
    quotient = sp.cancel(hhat / t**6)
    require(sp.denom(quotient) == 1, "Hhat/t^6 polynomial on Lambda=0")
    local = sp.expand(y * quotient - u**6 * y - sp.Rational(1, 2))
    local_mod_t6 = sp.rem(sp.Poly(local, t), sp.Poly(t**6, t)).as_expr()
    expected_general = (
        delta * y**2
        + c * y
        - sp.Rational(1, 2)
        + alpha * t * y**2
        + t**2 * (upsilon * y**2 + sp.Rational(8, 3) * y)
        + eta * t**3 * y**2
        + t**4 * (Delta * y**2 - 3 * y)
        - u**6 * y
    )
    require_zero(local_mod_t6 - expected_general, "general local expansion")

    # Balanced repetition is c^2+2*delta=0.  It turns the general equation
    # into exactly the THM-4292 critical polynomial through t^5.
    repeat_W = -2 * U - c**2 / 2
    repeated = sp.expand(local_mod_t6.subs(W, repeat_W))
    expected_repeated = (
        -sp.Rational(1, 2) * (c * y - 1) ** 2
        + alpha * t * y**2
        + t**2 * (upsilon * y**2 + sp.Rational(8, 3) * y)
        + eta * t**3 * y**2
        + t**4 * (Delta * y**2 - 3 * y)
        - u**6 * y
    )
    require_zero(repeated - expected_repeated, "transported repeated ladder")

    balanced_discriminant = sp.expand((c - X**6) ** 2 + 2 * delta)
    repeated_discriminant = sp.factor(balanced_discriminant.subs(W, repeat_W))
    require_zero(
        repeated_discriminant - X**6 * (X**6 - 2 * c),
        "balanced repeated discriminant",
    )
    nonzero_factor = X**6 - 2 * c
    require(
        sp.factor(sp.resultant(nonzero_factor, sp.diff(nonzero_factor, X), X)) != 0,
        "nonzero balanced roots are simple",
    )

    # Reconstruct all four moving-critical coefficients rather than merely
    # importing their values from THM-4292.
    quadratic_a = -c**2 / 2 + alpha * t + upsilon * t**2 + eta * t**3 + Delta * t**4
    quadratic_b = c + sp.Rational(8, 3) * t**2 - 3 * t**4
    critical = sp.series(
        -sp.Rational(1, 2) - quadratic_b**2 / (4 * quadratic_a),
        t,
        0,
        7,
    ).removeO().expand()
    C1 = sp.factor(coefficient(critical, t, 1))
    require_zero(C1 - alpha / c**2, "C1")
    critical_1 = sp.factor(critical.subs(alpha, 0))
    C2 = sp.factor(coefficient(critical_1, t, 2))
    require_zero(C2 - upsilon / c**2 - sp.Rational(8, 3) / c, "C2")
    critical_2 = sp.factor(critical_1.subs(upsilon, -sp.Rational(8, 3) * c))
    C3 = sp.factor(coefficient(critical_2, t, 3))
    require_zero(C3 - eta / c**2, "C3")
    critical_3 = sp.factor(critical_2.subs(eta, 0))
    C4 = sp.factor(coefficient(critical_3, t, 4))
    require_zero(
        C4 - (Delta + sp.Rational(32, 9) - 3 * c) / c**2,
        "C4 moving-critical correction",
    )

    # The maximal cancellation is still an allowed W!=0 point.
    c_star = sp.Rational(5152, 405)
    Delta_star = sp.Rational(4672, 135)
    U_star = sp.Integer(1)
    W_star = sp.factor(-2 * U_star - c_star**2 / 2)
    Z_star = sp.factor(-U_star - W_star)
    require_zero(c_star - (a0 - sp.Rational(7, 6) * Delta_star), "c relation")
    require_zero(Delta_star + sp.Rational(32, 9) - 3 * c_star, "C4 cancellation")
    require(W_star != 0 and Z_star != 0, "genuine W-nonzero hostile")
    require_zero((W_star**2 - 4 * U_star * Z_star) - c_star**4 / 4, "hostile D")

    # Exhaust every equality slope and positive Keller-form order.  These
    # are homogeneous in the ramification index, so primitive integer pairs
    # suffice for the equality rows.
    equality_rows: list[tuple[int, int, int, sp.Rational, sp.Rational]] = []
    for j in range(1, 5):
        s_value = 6 - j
        beta_value = 6 + j
        ratio = sp.Rational(beta_value, s_value)
        d_j = j * (s_value + beta_value)
        d_b = 6 * (beta_value - s_value)
        require(d_j == d_b, f"gap equality j={j}")
        form_order = sp.Rational(3 * s_value + 5 * beta_value) - sp.Rational(d_j, 2)
        require(form_order > 0, f"positive form order j={j}")
        equality_rows.append((j, s_value, beta_value, ratio, form_order))

    grid_checks = 0
    for s_value in range(1, 65):
        for beta_value in range(s_value + 1, 65):
            require(6 * s_value + 2 * beta_value > 0, "b-split order")
            require(6 * (beta_value - s_value) < 6 * (s_value + beta_value), "b before t6")
            for j in range(1, 5):
                twice_order = (6 - j) * s_value + (10 - j) * beta_value
                require(twice_order > 0, f"critical order j={j}")
                grid_checks += 1

    print("THM-4297 PRIMARY EXACT AUDIT: PASS")
    print(f"Lambda_wall: Z={Z}, D={sp.factor(D)}, delta={delta}")
    print("contact: algebra=b^12, derivative=-delta, analytic_type=A23")
    print(
        f"central: gcd_degree={gcd.degree()}, S0_resultant={uniform_resultant}, "
        "good_differential_order=9"
    )
    print("repeat_condition: c^2+2*delta=0")
    print(f"repeated_discriminant={repeated_discriminant}")
    print(f"critical_coefficients: C1={C1}; C2={C2}; C3={C3}; C4={C4}")
    print(f"W_nonzero_maximal_hostile: U={U_star}, W={W_star}, Z={Z_star}, c={c_star}")
    print("equality_rows: j,s,beta,beta/s,form_order")
    for row in equality_rows:
        print("  " + ",".join(str(value) for value in row))
    print(f"positive_order_grid_checks={grid_checks}")
    print("scope: Lambda=0, U*Z*D!=0, exact-M12 reduced seam, characteristic zero")


if __name__ == "__main__":
    main()
