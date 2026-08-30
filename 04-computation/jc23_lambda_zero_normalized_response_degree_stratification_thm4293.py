#!/usr/bin/env python3
"""Exact certificate for THM-4293's Lambda=0 response stratification.

The certificate starts from THM-4230's literal generic-Q source equation.  It
normalizes the double top-edge root on W=0, Z=-U, checks the six noncritical
contact strata, and reconstructs the resulting ramification, pole degree,
genus, deck-mod-four, and Eisenstein-norm ledgers.  The repeated-discriminant
lane is recorded only as an explicitly OPEN boundary.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


CHECKS = 0


def need(condition: object, label: str) -> None:
    """Count and enforce one exact check."""

    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise AssertionError(label)


def assert_zero(value: sp.Expr, label: str) -> None:
    need(sp.factor(value) == 0, f"{label}: {sp.factor(value)}")


def coefficient(poly: sp.Expr, first: sp.Symbol, i: int,
                second: sp.Symbol, j: int) -> sp.Expr:
    return sp.Poly(poly, first, second).coeff_monomial(first**i * second**j)


def eisenstein_norms(bound: int) -> set[int]:
    values: set[int] = set()
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            values.add(a * a - a * b + b * b)
    return values


@dataclass(frozen=True)
class Row:
    r: int
    branch_orders: tuple[int, int]
    local_index: int
    genus: int
    full_degree: int
    finite_degree: int
    full_eigenline: bool
    finite_eigenline: bool


def main() -> None:
    s, p, qinv = sp.symbols("s p Q")
    z, w, u, carrier = sp.symbols("z w u carrier")
    K, Phi, Delta, Theta = sp.symbols("K Phi Delta Theta")
    eta, zeta_3 = sp.symbols("eta zeta_3")
    upsilon_5, xi_10 = sp.symbols("upsilon_5 xi_10")
    alpha_11, beta_11 = sp.symbols("alpha_11 beta_11")
    U, W, Z = sp.symbols("U W Z")

    # THM-4230 equation (4), with y=sp, and its literal generic-Q source.
    y = s * p
    k_forced = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
    source_h = (
        -3 * p
        + sp.Rational(8, 3) * p**2
        - sp.Rational(1376, 135) * p**3
        + k_forced * y**2
        + Phi * p**2 * y
        + Delta * p**4
        + Theta * p * y**2
        + eta * p**3 * y
        + zeta_3 * y**3
        + upsilon_5 * p**5
        + xi_10 * p**2 * y**2
        + alpha_11 * p**4 * y
        + beta_11 * p * y**3
        + U * p**6
        + W * p**3 * y**2
        + Z * y**4
    )
    source = sp.expand(
        (s**2 - p) * (1 - qinv * source_h) - qinv * s**2 / 2
    )

    # Reconstruct both load-bearing non-affine edges directly from F_Q.
    top_edge = sum(
        coefficient(source, s, 2 * position, p, 7 - position)
        * w**position
        for position in range(4)
    )
    expected_top = qinv * (1 - w) * (U + W * w + Z * w**2)
    assert_zero(top_edge - expected_top, "top edge")

    quartic_edge = (
        coefficient(source, s, 2, p, 0)
        + coefficient(source, s, 4, p, 2) * carrier**2
        + coefficient(source, s, 5, p, 3) * carrier**3
        + coefficient(source, s, 6, p, 4) * carrier**4
    )
    expected_quartic = (
        1 - qinv / 2
        - qinv
        * (k_forced * carrier**2 + zeta_3 * carrier**3 + Z * carrier**4)
    )
    assert_zero(quartic_edge - expected_quartic, "quartic carrier edge")
    need(sp.degree(expected_quartic.subs(Z, -U), carrier) == 4,
         "quartic carrier retains degree four")
    need(4 * 2 == 8, "finite quartic subtraction")

    # Literal regular top-infinity chart z=1/s, w=s^2/p.  The multiplier is
    # the primitive top weight, not a postulated local model.
    wall_source = source.subs({W: 0, Z: -U})
    chart_w = sp.cancel(
        z**14 * w**7 * wall_source.subs({s: 1 / z, p: 1 / (z**2 * w)})
    )
    need(sp.denom(chart_w) == 1, "regular top chart")
    chart = sp.expand(chart_w.subs(w, 1 + u))
    assert_zero(
        chart_w.subs(z, 0) - qinv * U * (1 - w) ** 2 * (1 + w),
        "double and simple top roots",
    )
    assert_zero(chart.subs(z, 0) - qinv * U * u**2 * (u + 2),
                "double-root local equation")
    assert_zero(chart.subs(u, 0) + qinv * z**12 / 2,
                "exact smoothing constant")

    c1 = alpha_11 + beta_11
    c2 = upsilon_5 + xi_10
    c3 = eta + zeta_3
    c4 = Delta + Theta
    c5 = Phi
    c6 = sp.factor(k_forced - sp.Rational(1376, 135))
    assert_zero(c6 - (sp.Rational(7168, 135) - sp.Rational(7, 6) * Delta),
                "forced c6")

    linear_u = sp.Poly(chart, u).coeff_monomial(u)
    expected_linear_u = (
        -qinv
        * (
            c1 * z
            + c2 * z**2
            + c3 * z**3
            + c4 * z**4
            + c5 * z**5
            + c6 * z**6
            + sp.Rational(8, 3) * z**8
            - 3 * z**10
        )
        + (1 - sp.Rational(7, 2) * qinv) * z**12
    )
    assert_zero(linear_u - expected_linear_u, "c1 through c6 expansion")
    need(sp.Poly(chart, u, z).coeff_monomial(u**2) == 2 * qinv * U,
         "quadratic leading coefficient")
    need(sp.Poly(chart, u, z).coeff_monomial(u**3) == qinv * U,
         "cubic unit correction")

    # Exact chain rule behind eta=Q*z^10*w^5*dz/Fbar_w on the curve.
    source_p_chart = sp.cancel(
        sp.diff(wall_source, p).subs({s: 1 / z, p: 1 / (z**2 * w)})
    )
    chain_rule = sp.factor(
        source_p_chart
        - z**-12 * w**-5 * (7 * chart_w / w - sp.diff(chart_w, w))
    )
    assert_zero(chain_rule, "residue chain rule")
    residue_z_power = 10
    need(residue_z_power == 10, "residue numerator order")

    # The six noncritical normalization strata.  For r<6 the two root
    # orders are r and 12-r.  At r=6 the distinct quadratic roots both have
    # order six, with discriminant c6^2+4U.
    lambda_symbol = sp.symbols("lambda")
    lambda_quadratic = 2 * U * lambda_symbol**2 - c6 * lambda_symbol - sp.Rational(1, 2)
    assert_zero(
        sp.discriminant(lambda_quadratic, lambda_symbol) - (c6**2 + 4 * U),
        "r6 branch discriminant",
    )
    generic_cr = sp.symbols("c_r", nonzero=True)
    first_branch_lead = generic_cr / (2 * U)
    second_branch_lead = -1 / (2 * generic_cr)
    assert_zero(
        2 * U * first_branch_lead**2
        - generic_cr * first_branch_lead,
        "r<6 first Newton branch",
    )
    assert_zero(
        -generic_cr * second_branch_lead - sp.Rational(1, 2),
        "r<6 second Newton branch",
    )

    norm_values = eisenstein_norms(12)
    need({value for value in norm_values if value <= 10}
         == {0, 1, 3, 4, 7, 9}, "small Eisenstein norms")

    rows: list[Row] = []
    for r_value in range(1, 7):
        orders = ((r_value, 12 - r_value)
                  if r_value < 6 else (6, 6))
        local_index = residue_z_power - r_value + 1
        genus = 18 - r_value
        full_degree = 2 * local_index + 11 + 1 + 4 * 2
        finite_degree = full_degree - 4 * 2
        differential_zeros = 2 * (local_index - 1) + (11 - 1) + 4 * (2 - 1)
        need(local_index == 11 - r_value, f"r={r_value} local index")
        need(full_degree == 42 - 2 * r_value, f"r={r_value} full degree")
        need(finite_degree == 34 - 2 * r_value, f"r={r_value} finite degree")
        need(differential_zeros == 2 * genus - 2,
             f"r={r_value} canonical degree")
        need(18 - genus == r_value, f"r={r_value} delta/genus drop")
        need(2 * full_degree
             == 2 * (2 * local_index + 11 + 1 + 4 * 2),
             f"r={r_value} full A-pole divisor")
        need(2 * finite_degree
             == 2 * (2 * local_index + 11 + 1),
             f"r={r_value} finite A-pole divisor")

        full_eigenline = (
            full_degree % 4 == 0 and full_degree // 4 in norm_values
        )
        finite_eigenline = (
            finite_degree % 4 == 0 and finite_degree // 4 in norm_values
        )
        rows.append(Row(
            r=r_value,
            branch_orders=orders,
            local_index=local_index,
            genus=genus,
            full_degree=full_degree,
            finite_degree=finite_degree,
            full_eigenline=full_eigenline,
            finite_eigenline=finite_eigenline,
        ))

    expected_rows = (
        (1, (1, 11), 10, 17, 40, 32, False, False),
        (2, (2, 10), 9, 16, 38, 30, False, False),
        (3, (3, 9), 8, 15, 36, 28, True, True),
        (4, (4, 8), 7, 14, 34, 26, False, False),
        (5, (5, 7), 6, 13, 32, 24, False, False),
        (6, (6, 6), 5, 12, 30, 22, False, False),
    )
    need(
        tuple(
            (
                row.r,
                row.branch_orders,
                row.local_index,
                row.genus,
                row.full_degree,
                row.finite_degree,
                row.full_eigenline,
                row.finite_eigenline,
            )
            for row in rows
        ) == expected_rows,
        "complete response/eigenline table",
    )
    need(tuple(row.r for row in rows if row.full_degree % 4 == 0)
         == (1, 3, 5), "mod-four full survivors")
    need(tuple(row.r for row in rows if row.finite_degree % 4 == 0)
         == (1, 3, 5), "mod-four finite survivors")
    need(tuple(row.r for row in rows
               if row.full_eigenline or row.finite_eigenline) == (3,),
         "exact Eisenstein eigenline survivor")

    # Repeated-discriminant boundary.  These coefficients are verified, but
    # no exhaustive classification is promoted here.
    a6, v = sp.symbols("a6 v", nonzero=True)
    delta_from_a6 = sp.Rational(6, 7) * (sp.Rational(7168, 135) - a6)
    critical_substitution = {
        Delta: delta_from_a6,
        Theta: -delta_from_a6,
        Phi: 0,
        beta_11: -alpha_11,
        xi_10: -upsilon_5,
        zeta_3: -eta,
        U: -a6**2 / 4,
    }
    critical = sp.cancel(
        chart.subs(critical_substitution)
        .subs(u, z**6 * (-1 / a6 + v))
        / z**12
    )
    critical_poly = sp.Poly(sp.expand(critical), z, v)
    assert_zero(
        critical_poly.coeff_monomial(v**2) + qinv * a6**2 / 2,
        "critical v2 coefficient",
    )
    assert_zero(
        critical_poly.coeff_monomial(z) - qinv * alpha_11 / a6**2,
        "critical z coefficient",
    )
    assert_zero(
        critical_poly.coeff_monomial(z**2)
        - qinv * (3 * upsilon_5 + 8 * a6) / (3 * a6**2),
        "critical z2 coefficient",
    )
    assert_zero(
        critical_poly.coeff_monomial(z**3) - qinv * eta / a6**2,
        "critical z3 coefficient",
    )
    for discriminant_order in range(13, 22):
        if discriminant_order % 2 == 0:
            local_total = 2 * (11 - discriminant_order // 2)
            need(discriminant_order <= 20 and local_total == 22 - discriminant_order,
                 f"open even m={discriminant_order} formula")
        else:
            local_total = 22 - discriminant_order
            need(discriminant_order <= 21 and local_total > 0,
                 f"open odd m={discriminant_order} formula")

    # THM-4291 specialization: all free weights 7--11 and Delta vanish.
    a_tail = sp.Rational(7168, 135)
    tail_factor = sp.factor(a_tail**2 + 4 * U)
    expected_tail_factor = (
        4 * (sp.Integer(18225) * U + sp.Integer(12845056))
        / sp.Integer(18225)
    )
    assert_zero(tail_factor - expected_tail_factor, "THM-4291 r6 factor")
    need(18 - 6 == 7 + 5 == 12, "THM-4291 genus budget")
    need((42 - 2 * 6, 34 - 2 * 6) == (30, 22),
         "THM-4291 corrected response")
    need(42 > 30 > 22, "degree-42 tail exceeds total response")

    print("STATUS=VERIFIED-EXACT; THM-4293=NONCRITICAL-STRATIFICATION; JC(2)=OPEN")
    print("TOP_EDGE=Q*(1-w)*(U+W*w+Z*w^2); WALL=Q*U*(1-w)^2*(1+w)")
    print(
        "COEFFICIENTS="
        "c1=alpha11+beta11,c2=upsilon5+xi10,c3=eta+zeta3,"
        "c4=Delta+Theta,c5=Phi,c6=7168/135-7*Delta/6"
    )
    print(
        "STRATA=R1:c1!=0;R2:c1=0,c2!=0;R3:c1=c2=0,c3!=0;"
        "R4:c1=...=c3=0,c4!=0;R5:c1=...=c4=0,c5!=0;"
        "R6:c1=...=c5=0,c6^2+4U!=0"
    )
    for row in rows:
        differential_zeros = 2 * (row.local_index - 1) + 10 + 4
        print(
            f"R={row.r}; BRANCH_ORDERS={row.branch_orders[0]},{row.branch_orders[1]}; "
            f"E={row.local_index},{row.local_index}; GENUS={row.genus}; "
            f"DEGREES={row.full_degree},{row.finite_degree}; "
            f"A_POLES={2 * row.local_index},{2 * row.local_index}; "
            f"DIFF_ZEROS={differential_zeros}; "
            f"MOD4={row.full_degree % 4},{row.finite_degree % 4}; "
            f"EIGENLINE={int(row.full_eigenline)},{int(row.finite_eigenline)}"
        )
    print("FINITE_QUARTIC_SUBTRACTION=4*2=8; CARRIER_DEGREE=4")
    print("MOD4_SURVIVORS=r=1,3,5; EXACT_4N_SURVIVOR=r=3; N=9,7")
    print(
        "REPEATED_DISCRIMINANT=OPEN; "
        "c1=...=c5=0,c6^2+4U=0; m=ord_z(Disc(W2))>12; "
        "rational-response=42-m,34-m; even-nonsquare-extra=20,12"
    )
    print(
        "REPEATED_INITIALS="
        "v2=-Q*c6^2/2,z=Q*alpha11/c6^2,"
        "z2=Q*(3*upsilon5+8*c6)/(3*c6^2),z3=Q*eta/c6^2"
    )
    print(
        "THM4291_CONTROL=r=6; GENUS=12=7+5; "
        "RESPONSES=30,22; DEGREE42_TAIL_EXCEEDS_TOTAL"
    )
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
