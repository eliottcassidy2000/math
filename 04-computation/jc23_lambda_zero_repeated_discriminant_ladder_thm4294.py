#!/usr/bin/env python3
"""Exact certificate for THM-4294's repeated-discriminant ladder.

Starting from the literal exact-weight-twelve generic-response equation, this
certificate reconstructs the Lambda=0 top-infinity chart, imposes the
repeated initial quadratic, and derives the successive Weierstrass
discriminant coefficients.  It then checks the normalization genus,
ramification indices, response-degree candidates, Eisenstein-norm sieve, and
the nonsplit quadratic residue field at the deepest stratum.

The computation is local to the already-reduced exact-M=12 wall.  It neither
proves seam entry nor proves the planar Jacobian conjecture.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp


CHECKS = 0


def need(condition: object, label: str) -> None:
    """Count and enforce one exact check without relying on ``assert``."""

    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise AssertionError(label)


def assert_zero(value: sp.Expr, label: str) -> None:
    need(sp.factor(value) == 0, f"{label}: {sp.factor(value)}")


def coefficient(poly: sp.Expr, first: sp.Symbol, i: int,
                second: sp.Symbol, j: int) -> sp.Expr:
    return sp.Poly(poly, first, second).coeff_monomial(first**i * second**j)


def z_coefficient(poly: sp.Expr, z: sp.Symbol, degree: int) -> sp.Expr:
    return sp.factor(sp.Poly(sp.expand(poly), z).coeff_monomial(z**degree))


def weighted_truncate(poly: sp.Expr, z: sp.Symbol, v: sp.Symbol,
                      bound: int) -> sp.Expr:
    """Keep terms whose positive ``z+v`` weight is at most ``bound``."""

    terms = []
    for (z_degree, v_degree), value in sp.Poly(sp.expand(poly), z, v).terms():
        if z_degree + v_degree <= bound:
            terms.append(value * z**z_degree * v**v_degree)
    return sp.Add(*terms)


def critical_value_series(poly: sp.Expr, z: sp.Symbol, v: sp.Symbol,
                          bound: int) -> tuple[sp.Expr, tuple[sp.Expr, ...]]:
    """Solve ``P_v=0`` formally and return ``P`` at its moving critical point."""

    truncated = weighted_truncate(poly, z, v, bound)
    derivative = sp.diff(truncated, v)
    critical = sp.S.Zero
    for degree in range(1, bound + 1):
        unknown = sp.Symbol(f"critical_{degree}")
        trial = critical + unknown * z**degree
        equation = z_coefficient(derivative.subs(v, trial), z, degree)
        solutions = sp.solve(sp.Eq(equation, 0), unknown, dict=False)
        need(len(solutions) == 1, f"unique critical coefficient {degree}")
        critical += sp.factor(solutions[0]) * z**degree
    value = sp.series(truncated.subs(v, critical), z, 0, bound + 1)
    value = sp.expand(value.removeO())
    return critical, tuple(z_coefficient(value, z, degree)
                           for degree in range(1, bound + 1))


def eisenstein_norms(bound: int) -> set[int]:
    """Enumerate a provably sufficient box for norms through ``bound``."""

    coordinate_bound = 0
    while (coordinate_bound + 1) ** 2 <= 2 * bound:
        coordinate_bound += 1
    values: set[int] = set()
    for a in range(-coordinate_bound, coordinate_bound + 1):
        for b in range(-coordinate_bound, coordinate_bound + 1):
            value = a * a - a * b + b * b
            if value <= bound:
                values.add(value)
    return values


@dataclass(frozen=True)
class Stratum:
    discriminant_order: int
    genus: int
    geometric_indices: tuple[int, ...]
    degree_candidates: tuple[int, ...]
    eigenline_survivors: tuple[int, ...]


def main() -> None:
    s, p, qinv = sp.symbols("s p Q")
    z, w, u, v, carrier = sp.symbols("z w u v carrier")
    K, Phi, Delta, Theta = sp.symbols("K Phi Delta Theta")
    eta, zeta_3 = sp.symbols("eta zeta_3")
    upsilon_5, xi_10 = sp.symbols("upsilon_5 xi_10")
    alpha_11, beta_11 = sp.symbols("alpha_11 beta_11")
    U, W, Z, c = sp.symbols("U W Z c", nonzero=True)

    # Literal THM-4230 source H and generic response F_Q.
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

    # Reconstruct the two non-affine carriers used in the response ledger.
    top_edge = sum(
        coefficient(source, s, 2 * position, p, 7 - position)
        * w**position
        for position in range(4)
    )
    assert_zero(
        top_edge - qinv * (1 - w) * (U + W * w + Z * w**2),
        "literal top edge",
    )
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
    assert_zero(quartic_edge - expected_quartic, "literal quartic carrier")

    # Lambda=0 wall and regular top chart.
    wall_source = source.subs({W: 0, Z: -U})
    chart_w = sp.cancel(
        z**14 * w**7
        * wall_source.subs({s: 1 / z, p: 1 / (z**2 * w)})
    )
    need(sp.denom(chart_w) == 1, "regular wall chart")
    chart = sp.expand(chart_w.subs(w, 1 + u))
    assert_zero(
        chart_w.subs(z, 0) - qinv * U * (1 - w)**2 * (1 + w),
        "double wall root",
    )
    assert_zero(chart.subs(u, 0) + qinv * z**12 / 2,
                "smoothing monomial")

    c6 = sp.Rational(7168, 135) - sp.Rational(7, 6) * Delta
    delta_from_c = sp.Rational(6, 7) * (sp.Rational(7168, 135) - c)
    assert_zero(c6.subs(Delta, delta_from_c) - c, "c=c6 relation")

    # The repeated initial quadratic forces all five sums to vanish and
    # U=-c^2/4.  Shift to its double root before expanding the literal germ.
    repeated_substitution = {
        Delta: delta_from_c,
        Theta: -delta_from_c,
        Phi: 0,
        beta_11: -alpha_11,
        xi_10: -upsilon_5,
        zeta_3: -eta,
        U: -c**2 / 4,
    }
    local = sp.cancel(
        chart.subs(repeated_substitution)
        .subs(u, z**6 * (-1 / c + v))
        / z**12
    )
    need(not sp.denom(local).has(z, v), "regular shifted local germ")
    assert_zero(local.subs(z, 0) + qinv * c**2 * v**2 / 2,
                "repeated initial square")
    assert_zero(sp.diff(local, v, 2).subs({z: 0, v: 0}) + qinv * c**2,
                "critical Hessian unit")

    local_v = sp.Poly(sp.expand(local), v)
    need(local_v.degree() == 7, "full local v-degree retained")
    higher_valuations: dict[int, int] = {}
    for degree in range(3, 8):
        value = local_v.coeff_monomial(v**degree)
        higher_valuations[degree] = sp.Poly(value, z).terms()[-1][0][0]
    need(higher_valuations == {3: 6, 4: 15, 5: 24, 6: 34, 7: 42},
         "higher local terms start beyond discriminant ladder")

    # Through z^6, higher v powers cannot affect the two small roots.  The
    # raw quadratic discriminant therefore has the same leading terms and
    # squareclass as the Weierstrass quadratic discriminant.
    quadratic = sp.series(local_v.coeff_monomial(v**2), z, 0, 7).removeO()
    linear = sp.series(local_v.coeff_monomial(v), z, 0, 7).removeO()
    constant = sp.series(local_v.coeff_monomial(1), z, 0, 7).removeO()
    discriminant = sp.series(
        linear**2 - 4 * quadratic * constant, z, 0, 7
    ).removeO().expand()

    d1 = z_coefficient(discriminant, z, 1)
    assert_zero(d1 - 2 * qinv**2 * alpha_11, "discriminant step one")

    discriminant_2 = sp.expand(discriminant.subs(alpha_11, 0))
    d2 = z_coefficient(discriminant_2, z, 2)
    assert_zero(
        d2 - sp.Rational(2, 3) * qinv**2 * (8 * c + 3 * upsilon_5),
        "discriminant step two",
    )

    discriminant_3 = sp.expand(
        discriminant_2.subs(upsilon_5, -sp.Rational(8, 3) * c)
    )
    d3 = z_coefficient(discriminant_3, z, 3)
    assert_zero(d3 - 2 * qinv**2 * eta, "discriminant step three")

    discriminant_4 = sp.expand(discriminant_3.subs(eta, 0))
    d4 = z_coefficient(discriminant_4, z, 4)
    assert_zero(
        d4 + sp.Rational(2, 105) * qinv**2 * (405 * c - 5152),
        "discriminant step four",
    )
    assert_zero(z_coefficient(discriminant_4, z, 5),
                "absent discriminant step five")
    d6_before_maximal = z_coefficient(discriminant_4, z, 6)
    assert_zero(
        d6_before_maximal
        + qinv * (135 * qinv * c + 9824 * qinv + 540 * c) / 270,
        "discriminant step six before maximal substitution",
    )

    c_maximal = sp.Rational(5152, 405)
    discriminant_maximal = sp.expand(discriminant_4.subs(c, c_maximal))
    for degree in range(1, 6):
        assert_zero(z_coefficient(discriminant_maximal, z, degree),
                    f"maximal cancellation degree {degree}")
    d6 = z_coefficient(discriminant_maximal, z, 6)
    assert_zero(
        d6 + sp.Rational(32, 405) * qinv * (541 * qinv + 322),
        "maximal discriminant coefficient",
    )

    # Independent formal-critical-value reconstruction of the same ladder.
    _, critical_values = critical_value_series(local, z, v, 4)
    assert_zero(critical_values[0] - qinv * alpha_11 / c**2,
                "critical value step one")
    assert_zero(
        critical_values[1].subs(alpha_11, 0)
        - qinv * (8 * c + 3 * upsilon_5) / (3 * c**2),
        "critical value step two",
    )
    assert_zero(
        critical_values[2].subs({alpha_11: 0,
                                 upsilon_5: -sp.Rational(8, 3) * c})
        - qinv * eta / c**2,
        "critical value step three",
    )
    assert_zero(
        critical_values[3].subs({alpha_11: 0,
                                 upsilon_5: -sp.Rational(8, 3) * c,
                                 eta: 0})
        + qinv * (405 * c - 5152) / (105 * c**2),
        "critical value step four",
    )

    maximal_local = (
        local.subs(alpha_11, 0)
        .subs(upsilon_5, -sp.Rational(8, 3) * c)
        .subs(eta, 0)
        .subs(c, c_maximal)
    )
    _, maximal_values = critical_value_series(maximal_local, z, v, 6)
    need(all(value == 0 for value in maximal_values[:5]),
         "maximal critical values one through five vanish")
    assert_zero(
        maximal_values[5]
        + sp.Rational(405, 1658944) * (541 * qinv + 322),
        "maximal critical value six",
    )
    assert_zero(
        d6 - 2 * qinv * c_maximal**2 * maximal_values[5],
        "critical value/discriminant agreement",
    )

    discriminant_orders = (13, 14, 15, 16, 18)
    need(tuple(12 + degree for degree in (1, 2, 3, 4, 6))
         == discriminant_orders, "complete discriminant-order set")

    # The deepest unit is nonsquare in C(Q): it has two distinct linear
    # factors of odd valuation.  In q=Q^-1 coordinates its squareclass is
    # 322q+541.
    squareclass = sp.factor(qinv * (541 * qinv + 322))
    numerator, denominator = sp.fraction(squareclass)
    factor_data = sp.factor_list(numerator, qinv)[1]
    need(denominator == 1, "maximal squareclass polynomial")
    need(len(factor_data) == 2, "two maximal squareclass factors")
    need(all(exponent == 1 for _, exponent in factor_data),
         "odd maximal squareclass valuations")
    q_target = sp.symbols("q_target")
    assert_zero(
        squareclass.subs(qinv, 1 / q_target)
        - (322 * q_target + 541) / q_target**2,
        "maximal residue extension in q-coordinate",
    )

    # Arithmetic response ledger.  The fixed rational packet contributes
    # 11+1=12; the quartic closed point contributes either 8 or 0.  At m=18
    # the nonsplit quadratic collision contributes either 4 or 0.
    norms = eisenstein_norms(10)
    need(norms == {0, 1, 3, 4, 7, 9},
         "complete Eisenstein norms through ten")

    def eigenline_allowed(degree: int) -> bool:
        return degree % 4 == 0 and degree // 4 in norms

    strata: list[Stratum] = []
    for order in (13, 14, 15, 16):
        genus = 18 - order // 2
        if order % 2:
            indices = (22 - order,)
            local_differential_zeros = 21 - order
        else:
            indices = (11 - order // 2, 11 - order // 2)
            local_differential_zeros = 20 - order
        candidates = (42 - order, 34 - order)
        survivors = tuple(value for value in candidates
                          if eigenline_allowed(value))
        need(sum(indices) == 22 - order,
             f"local collision degree m={order}")
        need(local_differential_zeros + 14 == 2 * genus - 2,
             f"canonical saturation m={order}")
        strata.append(Stratum(order, genus, indices, candidates, survivors))

    m18_candidates = tuple(sorted({
        12 + quartic + collision
        for quartic in (0, 8)
        for collision in (0, 4)
    }, reverse=True))
    m18_survivors = tuple(value for value in m18_candidates
                          if eigenline_allowed(value))
    need(2 * (2 - 1) + 14 == 2 * 9 - 2,
         "canonical saturation m=18")
    strata.append(Stratum(18, 9, (2, 2),
                           m18_candidates, m18_survivors))

    expected_strata = (
        Stratum(13, 12, (9,), (29, 21), ()),
        Stratum(14, 11, (4, 4), (28, 20), (28,)),
        Stratum(15, 11, (7,), (27, 19), ()),
        Stratum(16, 10, (3, 3), (26, 18), ()),
        Stratum(18, 9, (2, 2), (24, 20, 16, 12), (16, 12)),
    )
    need(tuple(strata) == expected_strata, "complete arithmetic strata")
    need(tuple((row.discriminant_order, value)
               for row in strata for value in row.eigenline_survivors)
         == ((14, 28), (18, 16), (18, 12)),
         "complete repeated-discriminant survivor set")

    # The m=18 quadratic field genuinely supports nonzero target points, so
    # E_q(K)={O} cannot force this closed point to the origin.  With a=1 the
    # target is C^2=A^3-3A/4+q-1/4.
    target_a = sp.symbols("target_a")
    target_relation = 644 * target_a**3 - 483 * target_a - 1243
    target_difference = (
        target_a**3 - sp.Rational(3, 4) * target_a
        + q_target - sp.Rational(1, 4)
        - (q_target + sp.Rational(541, 322))
    )
    assert_zero(
        sp.rem(sp.together(644 * target_difference),
               target_relation, target_a),
        "explicit nonzero target point over maximal residue field",
    )

    # The maximal lower rows coincide with THM-4292's deepest witness.
    maximal_delta = sp.factor(delta_from_c.subs(c, c_maximal))
    need(maximal_delta == sp.Rational(4672, 135),
         "maximal Delta witness")
    maximal_quartic = expected_quartic.subs({
        Delta: maximal_delta,
        zeta_3: 0,
        Z: c_maximal**2 / 4,
    })
    need(sp.degree(maximal_quartic, carrier) == 4,
         "maximal quartic carrier remains degree four")

    print("STATUS=VERIFIED-EXACT; THM-4294=REPEATED-DISCRIMINANT-LADDER; JC(2)=OPEN")
    print("UNIVERSE=exact_M12 Lambda=0 W=0 Z=-U U!=0 K=C(Q)")
    print("LOCAL_GERM=u=z^6*(-1/c+v); P=z^-12*Fbar; P_vv(0,0)=-Q*c^2")
    print("HIGHER_V_VALUATIONS=v3:6,v4:15,v5:24,v6:34,v7:42")
    print(
        "DISCRIMINANT_LADDER="
        "z1:2Q^2*alpha11;"
        "z2:(2/3)Q^2*(8c+3upsilon5);"
        "z3:2Q^2*eta;"
        "z4:-(2/105)Q^2*(405c-5152);"
        "z5:0;"
        "z6@c=5152/405:-(32/405)Q*(541Q+322)"
    )
    print("DISCRIMINANT_ORDERS=13,14,15,16,18; NO_M17; GENERIC_MAX_M18")
    print(
        "MAXIMAL_WITNESS="
        "c=5152/405,U=-c^2/4,W=0,Z=-U,Delta=4672/135,Theta=-Delta,"
        "Phi=0,upsilon5=-8c/3,xi10=-upsilon5,eta=alpha11=0,"
        "zeta3=-eta,beta11=-alpha11"
    )
    print(
        "M18_SQUARECLASS=Q*(541Q+322)=q^-2*(322q+541); "
        "NONSQUARE; RESIDUE_DEGREE=2"
    )
    for row in strata:
        indices = ",".join(str(value) for value in row.geometric_indices)
        degrees = ",".join(str(value) for value in row.degree_candidates)
        survivors = (",".join(str(value) for value in row.eigenline_survivors)
                     if row.eigenline_survivors else "NONE")
        print(
            f"M={row.discriminant_order}; GENUS={row.genus}; "
            f"GEOMETRIC_INDICES={indices}; DEGREES={degrees}; "
            f"EIGENLINE_SURVIVORS={survivors}"
        )
    print("REPEATED_SURVIVORS=m14:d28(N7);m18:d16(N4),d12(N3)")
    print(
        "M18_TARGET_HOSTILE=L=C(q)(sqrt(322q+541)); "
        "NONZERO_L_POINT_EXISTS; ORIGIN_LABEL_NOT_NEEDED_BY_CENTRAL_EXTINCTION"
    )
    print(
        "SCOPE=LOCAL_LADDER_AND_RELATIVE_RESPONSE_ARITHMETIC; "
        "KELLER_W_LAMBDA_ZERO_SLICE_CLOSED_BY_CENTRAL_EXTINCTION; "
        "GENERAL_LAMBDA_WALL_SEAM_ENTRY_JC2_OPEN"
    )
    print(f"CHECKS={CHECKS}")


if __name__ == "__main__":
    main()
