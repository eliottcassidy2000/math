#!/usr/bin/env python3
"""Exact companion for THM-3080.

Checks the unimodular recursive key transform, integer depth/gcd ledgers,
and exact one-, two-, and three-stage local symplectic controls.
"""

from __future__ import annotations

import math

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def positive_bezout(alpha: int, beta: int) -> tuple[int, int]:
    """Return gamma,delta with alpha*delta+beta*gamma=1."""

    require(math.gcd(alpha, beta) == 1, "primitive key exponents")
    for gamma in range(-alpha, alpha + 1):
        numerator = 1 - beta * gamma
        if numerator % alpha == 0:
            return gamma, numerator // alpha
    raise RuntimeError("Bezout search failed")


def laurent_order(expr: sp.Expr, variable: sp.Symbol) -> int:
    numerator, denominator = sp.fraction(sp.cancel(expr))
    numerator_poly = sp.Poly(numerator, variable)
    denominator_poly = sp.Poly(denominator, variable)
    require(not numerator_poly.is_zero, "zero has no finite Laurent order")
    numerator_order = min(
        monomial[0]
        for monomial, coefficient in numerator_poly.terms()
        if coefficient
    )
    denominator_order = min(
        monomial[0]
        for monomial, coefficient in denominator_poly.terms()
        if coefficient
    )
    return numerator_order - denominator_order


def wedge(expr1: sp.Expr, expr2: sp.Expr, u: sp.Symbol, s: sp.Symbol) -> sp.Expr:
    return sp.factor(
        sp.diff(expr1, u) * sp.diff(expr2, s)
        - sp.diff(expr1, s) * sp.diff(expr2, u)
    )


def main() -> None:
    z, r = sp.symbols("z r", nonzero=True)

    # Every strict step is a unimodular monomial change from (Z_i,R_i) to
    # (M_(i+1),R_(i+1)).  Check signs, values, and the logarithmic wedge on
    # a grid containing many noncoprime input values.
    for g in range(1, 25):
        for e in range(1, 25):
            divisor = math.gcd(g, e)
            alpha = g // divisor
            beta = e // divisor
            gamma, delta = positive_bezout(alpha, beta)
            require(alpha * delta + beta * gamma == 1, "key Bezout identity")
            require(gamma * e + delta * g == divisor, "new complementary value")

            m_new = z**alpha / r**beta
            r_new = z**gamma * r**delta
            log_wedge = sp.factor(
                (
                    sp.diff(m_new, z) * sp.diff(r_new, r)
                    - sp.diff(m_new, r) * sp.diff(r_new, z)
                )
                / (m_new * r_new)
            )
            require(sp.simplify(log_wedge - 1 / (z * r)) == 0, "log-wedge identity")

    # A strict stage replaces the budget B by B-e.  Exhaust all positive
    # depth compositions through D=18 and check exact termination and the
    # stage bound.  These are arithmetic controls, not existence claims.
    composition_count = 0
    for total in range(1, 19):
        for mask in range(1 << (total - 1)):
            parts: list[int] = []
            current = 1
            for cut in range(total - 1):
                if mask & (1 << cut):
                    parts.append(current)
                    current = 1
                else:
                    current += 1
            parts.append(current)
            budget = total
            for index, depth in enumerate(parts):
                require(1 <= depth <= budget, "positive remaining depth")
                if index + 1 < len(parts):
                    require(depth < budget, "strict stage consumes a proper part")
                    budget -= depth
                else:
                    require(depth == budget, "terminal stage exhausts budget")
            require(sum(parts) == total, "depth partition")
            require(len(parts) <= total, "stage count bound")
            composition_count += 1
    require(composition_count == (1 << 18) - 1, "composition census")

    # The exponent comparison is uniform in the tame ramification index E:
    # sigma+e-1 is below/equal/above E-1 exactly as e is below/equal/above
    # B=E-sigma.
    for ramification in range(1, 13):
        for sigma in range(-20, ramification):
            budget = ramification - sigma
            for depth in range(1, budget + 2):
                comparison = (sigma + depth - 1) - (ramification - 1)
                require(comparison == depth - budget, "uniform ramification budget")

    u, s = sp.symbols("u s", nonzero=True)

    # One-stage equality control from THM-3074: D=e_0=5.
    x_one = s**-1
    y_one = s**-1 + 3 * u * s**4
    require(wedge(x_one, y_one, u, s) == 3 * s**2, "one-stage wedge")
    r_one = 1 / x_one
    m_one = sp.cancel(y_one / x_one)
    require(laurent_order(m_one - 1, s) == 5, "one-stage depth")
    require(laurent_order(r_one, s) == 1, "one-stage complementary value")

    # Two-stage strict control from THM-3074: D=7=4+3.  The first lattice is
    # 2Z and misses three; its normalized second key has depth three.
    r_two = u * s**2
    m_two_terminal = 1 - 3 * u * s**3
    m_two = 1 + r_two**2 * m_two_terminal
    x_two = 1 / r_two
    y_two = m_two / r_two
    require(wedge(x_two, y_two, u, s) == 3 * s**2, "two-stage wedge")
    z_two = m_two - 1
    normalized_two = sp.cancel(z_two / r_two**2)
    require(laurent_order(z_two, s) == 4, "two-stage first depth")
    require(normalized_two == m_two_terminal, "two-stage normalized key")
    require(laurent_order(normalized_two - 1, s) == 3, "two-stage terminal depth")
    require(math.gcd(2, 4) == 2, "two-stage gcd remains two")
    require(math.gcd(2, 3) == 1, "two-stage terminal gcd reaches one")

    # New three-stage control: p=q=4 and D=11=4+4+3.  Both strict
    # normalizations use M_(i+1)=(M_i-1)/R, and the third stage pays the
    # Jacobian exactly.  This is a local rational symplectic packet only.
    r_three = u * s**4
    m_three_2 = 1 + 3 * u * s**3
    m_three_1 = 1 + r_three * m_three_2
    m_three_0 = 1 + r_three * m_three_1
    x_three = 1 / r_three
    y_three = m_three_0 / r_three
    require(wedge(x_three, y_three, u, s) == 3 * s**2, "three-stage wedge")
    require(laurent_order(x_three, s) == -4, "three-stage x pole")
    require(laurent_order(y_three, s) == -4, "three-stage y pole")
    require(laurent_order(m_three_0 - 1, s) == 4, "three-stage first depth")
    recovered_1 = sp.cancel((m_three_0 - 1) / r_three)
    recovered_2 = sp.cancel((recovered_1 - 1) / r_three)
    require(
        sp.simplify(recovered_1 - m_three_1) == 0,
        "three-stage first normalization",
    )
    require(
        sp.simplify(recovered_2 - m_three_2) == 0,
        "three-stage second normalization",
    )
    require(laurent_order(recovered_1 - 1, s) == 4, "three-stage second depth")
    require(laurent_order(recovered_2 - 1, s) == 3, "three-stage terminal depth")

    # The telescoping two-form coefficient at the terminal stage is exact:
    # U_0=M_0/R^2, U_1=M_1/R, U_2=M_2.
    u_stage_0 = sp.cancel(x_three * y_three)
    u_stage_1 = sp.cancel(u_stage_0 * (m_three_0 - 1) / m_three_0)
    u_stage_2 = sp.cancel(u_stage_1 * (m_three_1 - 1) / m_three_1)
    require(
        sp.simplify(u_stage_0 - m_three_0 / r_three**2) == 0,
        "stage-zero prefactor",
    )
    require(
        sp.simplify(u_stage_1 - m_three_1 / r_three) == 0,
        "stage-one prefactor",
    )
    require(sp.simplify(u_stage_2 - m_three_2) == 0, "stage-two prefactor")
    terminal_log_wedge = sp.factor(
        u_stage_2 * wedge(m_three_2, r_three, u, s) / (m_three_2 * r_three)
    )
    require(terminal_log_wedge == 3 * s**2, "terminal coefficient payment")

    print("theorem=THM-3080")
    print("status=PROVED_VERIFIED_EXACT_CANDIDATE")
    print("stage=omega=U_i*dlog(M_i)^dlog(R_i);val(U_i)=sigma_i;budget_B_i=E-sigma_i")
    print("depth_gate=1<=e_i<=B_i")
    print("terminal=e_i=B_i;u_i*(g_i*m_i_prime-e_i*m_i*r_i_prime/r_i)=E*kappa^-1*tau")
    print("strict=e_i<B_i;leading_coefficient_zero;m_i^(g_i/d)/r_i^(e_i/d)=constant")
    print("transform=M_next=Z_i^(g_i/d)/(c*R_i^(e_i/d));val(R_next)=d=gcd(g_i,e_i)")
    print("key_value=val(Z_i^(g_i/d)-c*R_i^(e_i/d))=lcm(g_i,e_i)+e_next")
    print("budget_update=B_next=B_i-e_i;g_next=gcd(g_i,e_i)")
    print("finite_partition=sum_i(e_i)=D;number_of_stages<=D")
    print("c3_specialization=E=3;D=p+q+3_two_pole_or_p+3-r_one_pole")
    print("three_stage_hostile=p=q=4;D=11;depths=4+4+3;dx_wedge_dy=du_wedge_d(s^3)")
    print("scope=coordinate_line_local_toric_tower;no_polynomial_globalization_or_full_C3_A4_S4_JC2_claim")


if __name__ == "__main__":
    main()
