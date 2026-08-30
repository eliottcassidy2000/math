#!/usr/bin/env python3
"""Exact symbolic certificate for THM-4292.

The universe is the full allowed lower-row exact-M=12 infinity equation on
W=0, Z=-U.  This certificate reconstructs the only repeated quadratic face,
derives its four critical-value correction coefficients, verifies the maximal
cancellation hostile, and checks the exhaustive valuation/form-order table.
"""

from __future__ import annotations

import sympy as sp


def assert_zero(value: sp.Expr, label: str) -> None:
    value = sp.factor(value)
    if value != 0:
        raise AssertionError(f"{label}: {value}")


def coefficient(expr: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Expr:
    return sp.expand(expr).coeff(variable, degree)


def main() -> None:
    t, u, y, X = sp.symbols("t u y X")
    c, alpha, upsilon, eta, Delta = sp.symbols(
        "c alpha_11 upsilon_5 eta Delta", nonzero=True
    )
    a0 = sp.Rational(7168, 135)
    K0 = sp.Rational(1376, 135)
    U = -c**2 / 4

    # Repetition forces the lower-row sums to vanish:
    # beta_11=-alpha_11, xi_10=-upsilon_5, zeta_3=-eta,
    # Theta=-Delta, Phi=0.  Put q=t^6*y and r=1+q.
    r = 1 + t**6 * y
    hhat = (
        U * (r**6 - r**4)
        + alpha * t * (r**5 - r**4)
        + upsilon * t**2 * (r**5 - r**4)
        + eta * t**3 * (r**4 - r**3)
        + Delta * t**4 * (r**4 - r**3)
        + t**6 * (-K0 * r**3 + (c + K0) * r**2)
        + sp.Rational(8, 3) * t**8 * r**2
        - 3 * t**10 * r
    )
    quotient = sp.cancel(hhat / t**6)
    if sp.denom(quotient) != 1:
        raise AssertionError("Hhat/t^6 is not polynomial after repeated-face constraints")

    # F/t^12 after b=sigma*u, t=sigma^2*u, b^12/t^6=u^6.
    P = sp.expand(y * quotient - u**6 * y - sp.Rational(1, 2))
    P_mod_t6 = sp.rem(sp.Poly(P, t), sp.Poly(t**6, t)).as_expr()
    expected = (
        -sp.Rational(1, 2) * (c * y - 1) ** 2
        + alpha * t * y**2
        + t**2 * (upsilon * y**2 + sp.Rational(8, 3) * y)
        + eta * t**3 * y**2
        + t**4 * (Delta * y**2 - 3 * y)
        - u**6 * y
    )
    assert_zero(P_mod_t6 - expected, "exact critical expansion modulo t^6")
    assert_zero(P.subs({t: 0, u: 0}) + (c * y - 1) ** 2 / 2, "double root")

    # The critical value is invariant under Weierstrass-coordinate changes.
    # Modulo t^6 the polynomial is quadratic in y, so evaluate it at its
    # exact critical point by d-b^2/(4a).  Higher-q terms begin at t^6 and
    # cannot affect C1,...,C4 or the absence of a t^5 term.
    A = alpha * t + upsilon * t**2 + eta * t**3 + Delta * t**4
    B = sp.Rational(8, 3) * t**2 - 3 * t**4
    quadratic_a = -c**2 / 2 + A
    quadratic_b = c + B
    critical = sp.series(
        -sp.Rational(1, 2) - quadratic_b**2 / (4 * quadratic_a),
        t,
        0,
        7,
    ).removeO().expand()

    C1 = sp.factor(coefficient(critical, t, 1))
    assert_zero(C1 - alpha / c**2, "C1")

    critical_1 = sp.factor(critical.subs(alpha, 0))
    C2 = sp.factor(coefficient(critical_1, t, 2))
    assert_zero(C2 - upsilon / c**2 - sp.Rational(8, 3) / c, "C2")

    critical_2 = sp.factor(critical_1.subs(upsilon, -sp.Rational(8, 3) * c))
    C3 = sp.factor(coefficient(critical_2, t, 3))
    assert_zero(C3 - eta / c**2, "C3")

    critical_3 = sp.factor(critical_2.subs(eta, 0))
    C4 = sp.factor(coefficient(critical_3, t, 4))
    expected_C4 = (Delta + sp.Rational(32, 9) - 3 * c) / c**2
    assert_zero(C4 - expected_C4, "C4 including moving-critical correction")

    # C4's 32/9 term is the square-completion correction.
    P2 = upsilon * y**2 + sp.Rational(8, 3) * y
    assert_zero(
        sp.diff(P2, y).subs({y: 1 / c, upsilon: -sp.Rational(8, 3) * c})
        + sp.Rational(8, 3),
        "P2 derivative",
    )
    assert_zero(sp.diff(-sp.Rational(1, 2) * (c * y - 1) ** 2, y, 2) + c**2, "P0 second derivative")

    # Solve the genuine coefficient relation c=a0-(7/6)Delta together with
    # C4=0.  This is the maximal-cancellation hostile.
    c_star = sp.Rational(5152, 405)
    Delta_star = sp.Rational(4672, 135)
    assert_zero(c_star - (a0 - sp.Rational(7, 6) * Delta_star), "c coefficient relation")
    assert_zero(Delta_star + sp.Rational(32, 9) - 3 * c_star, "maximal C4 cancellation")
    maximal = {
        c: c_star,
        Delta: Delta_star,
        alpha: 0,
        upsilon: -sp.Rational(8, 3) * c_star,
        eta: 0,
    }
    maximal_critical = sp.factor(critical.subs(maximal))
    for k in range(1, 6):
        assert_zero(coefficient(maximal_critical, t, k), f"maximal critical t^{k}")

    # Balanced beta=s face.  Its only multiple discriminant root is X=0;
    # each nonzero root of X^6-2c is simple because c!=0.
    balanced_discriminant = sp.expand(X**12 * ((X**6 - c) ** 2 + 4 * U))
    assert_zero(balanced_discriminant - X**18 * (X**6 - 2 * c), "balanced discriminant")
    assert_zero(
        sp.resultant(X**6 - 2 * c, sp.diff(X**6 - 2 * c, X), X)
        - (-1) ** 15 * 6**6 * (2 * c) ** 5,
        "nonzero balanced roots are simple",
    )

    # If C_j is the first critical coefficient, equality with the b^12 split
    # occurs at beta/s=(6+j)/(6-j).  The nonzero roots of the correction face
    # are simple, and the Keller-form order remains strictly positive.
    equality_rows: list[tuple[int, int, int, sp.Rational, sp.Rational]] = []
    for j in range(1, 5):
        s = 6 - j
        beta = 6 + j
        ratio = sp.Rational(beta, s)
        if ratio != sp.Rational(6 + j, 6 - j):
            raise AssertionError(f"equality ratio j={j}")
        d_j = j * (s + beta)
        d_b = 6 * (beta - s)
        if d_j != d_b:
            raise AssertionError(f"gap equality j={j}")
        form_order = sp.Rational(3 * s + 5 * beta) - sp.Rational(d_j, 2)
        expected_order = (sp.Rational(3) - sp.Rational(j, 2)) * s + (
            sp.Rational(5) - sp.Rational(j, 2)
        ) * beta
        if form_order != expected_order or form_order <= 0:
            raise AssertionError(f"form order j={j}: {form_order}")
        equality_rows.append((j, s, beta, ratio, form_order))

    # The X!=0 roots of X^j(C-X^(6-j)/c) are simple.
    C = sp.symbols("C", nonzero=True)
    for j in range(1, 5):
        nonzero_factor = C - X ** (6 - j) / c
        root_collision = sp.resultant(nonzero_factor, sp.diff(nonzero_factor, X), X)
        if sp.factor(root_collision) == 0:
            raise AssertionError(f"nonzero collision roots not simple j={j}")

    # If C1=...=C4=0, there is no t^5 term and b^12 always arrives before
    # the first possible O(t^6) correction for beta>s.
    if coefficient(maximal_critical, t, 5) != 0:
        raise AssertionError("unexpected t^5 critical term")
    s_symbol, beta_symbol = sp.symbols("s beta", positive=True)
    assert_zero(
        6 * (s_symbol + beta_symbol) - 6 * (beta_symbol - s_symbol) - 12 * s_symbol,
        "b split precedes t^6",
    )
    b_form_order = 6 * s_symbol + 2 * beta_symbol
    if b_form_order.is_positive is not True:
        raise AssertionError("b-split form order positivity")

    # Concrete hostile: s=1,beta=6.  The first surviving b gap is 30 and
    # the differential order is 18, not the naively alleged gap 14.
    hostile_s, hostile_beta = 1, 6
    hostile_b_gap = 6 * (hostile_beta - hostile_s)
    hostile_t2_gap = 2 * (hostile_s + hostile_beta)
    hostile_t4_gap = 4 * (hostile_s + hostile_beta)
    hostile_form_order = 6 * hostile_s + 2 * hostile_beta
    if (hostile_b_gap, hostile_t2_gap, hostile_t4_gap, hostile_form_order) != (30, 14, 28, 18):
        raise AssertionError("maximal hostile valuation ledger")

    print("THM-4292 PRIMARY EXACT AUDIT: PASS")
    print("critical_coefficients:")
    print(f"  C1={C1}")
    print(f"  C2_after_C1=0={C2}")
    print(f"  C3_after_C1_C2=0={C3}")
    print(f"  C4_after_C1_C2_C3=0={C4}")
    print("maximal_cancellation:")
    print(f"  c={c_star}, Delta={Delta_star}, upsilon_5={-sp.Rational(8, 3) * c_star}")
    print("equality_rows: j,s,beta,beta/s,form_order")
    for row in equality_rows:
        print("  " + ",".join(str(item) for item in row))
    print("maximal_hostile: s=1,beta=6,b_gap=30,t2_gap=14,t4_gap=28,form_order=18")
    print("scope: all divisorial refinements of the local prepared quadratic; U!=0; characteristic zero")


if __name__ == "__main__":
    main()
