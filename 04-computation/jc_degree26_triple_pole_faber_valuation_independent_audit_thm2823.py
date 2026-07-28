from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


t, q, T, d, s = sp.symbols("t q T d s")
r, zeta = sp.symbols("r zeta")
c18, c14, c10, c6, c2 = sp.symbols("c18 c14 c10 c6 c2")


def faber_coefficients(degree: int, stop: int) -> list[sp.Expr]:
    """Coefficients of (1+2dt^2+qt^3+(d^2-s)t^4)^(degree/4).

    This uses the logarithmic-derivative recurrence A C'=alpha A' C,
    rather than a series engine or either THM-2823 companion.
    """

    alpha = sp.Rational(degree, 4)
    terms = {2: 2 * d, 3: q, 4: d**2 - s}
    coefficients = [sp.Integer(0)] * (stop + 1)
    coefficients[0] = sp.Integer(1)
    for index in range(1, stop + 1):
        coefficients[index] = sp.expand(
            sum(
                ((alpha + 1) * exponent - index)
                * coefficient
                * coefficients[index - exponent]
                for exponent, coefficient in terms.items()
                if index >= exponent
            )
            / index
        )
    return coefficients


def replace_q_square(expression: sp.Expr, divide_q: bool = False) -> sp.Expr:
    if divide_q:
        expression = sp.cancel(expression / q)
    polynomial = sp.Poly(sp.expand(expression), q)
    answer = sp.Integer(0)
    for (exponent,), coefficient in polynomial.terms():
        require(exponent % 2 == 0, "an odd q power survived parity reduction")
        answer += coefficient * T ** (exponent // 2)
    return sp.expand(answer)


def build_degree_26_bank():
    phi_raw = sp.Integer(0)
    psi_raw = sp.Integer(0)
    k_raw = sp.Integer(0)
    for degree, multiplier in (
        (26, sp.Integer(1)),
        (18, c18),
        (14, c14),
        (10, c10),
        (6, c6),
        (2, c2),
    ):
        coefficients = faber_coefficients(degree, degree + 3)
        phi_raw += multiplier * 4 * coefficients[degree + 1]
        psi_raw += multiplier * 4 * coefficients[degree + 2]
        k_raw += multiplier * (
            4 * coefficients[degree + 3] + 2 * d * coefficients[degree + 1]
        )
    phi_over_q = replace_q_square(phi_raw, divide_q=True)
    psi = replace_q_square(psi_raw)
    k = replace_q_square(k_raw, divide_q=True)
    h = sp.cancel((k + d * phi_over_q / 2) / T)
    require(sp.denom(h) == 1, "K+d Phi/(2q) is not divisible by T")
    return phi_over_q, psi, sp.expand(k), sp.expand(h)


def weighted_initial(
    polynomial: sp.Expr, weights: tuple[int, int, int]
) -> tuple[int, sp.Expr]:
    poly = sp.Poly(sp.expand(polynomial), T, d, s)
    records = []
    for exponents, coefficient in poly.terms():
        weight = sum(exponent * value for exponent, value in zip(exponents, weights))
        records.append((weight, exponents, coefficient))
    least = min(weight for weight, _, _ in records)
    initial = sum(
        coefficient * T ** exponents[0] * d ** exponents[1] * s ** exponents[2]
        for weight, exponents, coefficient in records
        if weight == least
    )
    return least, sp.expand(initial)


def primitive_polynomial(expression: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    polynomial = sp.Poly(expression, variable, domain=sp.QQ)
    denominator_lcm = sp.ilcm(
        *[coefficient.q for coefficient in polynomial.all_coeffs()]
    )
    integral = sp.Poly(polynomial.as_expr() * denominator_lcm, variable, domain=sp.ZZ)
    content, primitive = sp.polys.polytools.primitive(integral)
    require(content != 0, "zero primitive content")
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive.as_expr()


def generic_ratio_form(initial: sp.Expr, common_power: int) -> sp.Expr:
    polynomial = sp.Poly(initial, T, d, s)
    answer = sp.Integer(0)
    for (t_power, d_power, s_power), coefficient in polynomial.terms():
        require(t_power == d_power, "generic initial form is not a polynomial in Td")
        require(
            2 * t_power + s_power == common_power,
            "generic initial form has the wrong common s power",
        )
        answer += coefficient * r**t_power
    return sp.expand(answer)


def exceptional_ratio_form(initial: sp.Expr, common_power: int) -> sp.Expr:
    polynomial = sp.Poly(initial, T, d, s)
    answer = sp.Integer(0)
    for (t_power, d_power, s_power), coefficient in polynomial.terms():
        require((t_power - d_power) % 2 == 0, "half-integral zeta exponent")
        r_power = d_power
        zeta_power = (t_power - d_power) // 2
        recovered_s_power = s_power + 2 * r_power + 3 * zeta_power
        require(
            recovered_s_power == common_power,
            "exceptional initial form has inconsistent common scale",
        )
        answer += coefficient * r**r_power * zeta**zeta_power
    return sp.expand(answer)


def valuation_lane(a_order: int, b_order: int | None) -> str:
    if a_order >= 3:
        return "regular-a>=3"
    if a_order == 2 and (b_order is None or b_order >= 1):
        return "regular-a=2,b>=1"
    if a_order <= 1 and (b_order is None or b_order >= 2):
        return "regular-a<=1,b>=2"
    if (a_order, b_order) == (0, 1):
        return "exceptional-(0,1)"
    if (a_order, b_order) in ((0, 0), (1, 0), (2, 0), (1, 1)):
        return "generic-polar"
    raise RuntimeError(f"unclassified valuation lane {(a_order, b_order)}")


def audit_case_partition() -> None:
    for a_order in range(13):
        for b_order in tuple(range(13)) + (None,):
            valuation_lane(a_order, b_order)

    generic = ((0, 0), (1, 0), (2, 0), (1, 1))
    for a_order, b_order in generic:
        t_order = 2 * a_order - 3
        d_order = 2 * b_order - 3
        s_order = a_order + b_order - 3
        require(t_order + d_order == 2 * s_order < 0, "generic scale relation failed")

    require(valuation_lane(0, 1) == "exceptional-(0,1)", "exceptional lane lost")
    require((-3) + (-1) == 2 * (-2), "exceptional scale relation failed")


def audit_carrier_multiplicities() -> None:
    x = sp.symbols("x")
    power_v = x**5 * (x**3 - 1) ** 3
    power_m = x * (x**3 - 1) * (x**3 - sp.Rational(1, 2))
    cube_factor = x**3 - 1
    require(
        sp.rem(power_v, cube_factor**3, x) == 0
        and sp.rem(power_v, cube_factor**4, x) != 0,
        "power V does not have cube-root order three",
    )
    require(
        sp.rem(power_m, cube_factor, x) == 0
        and sp.rem(power_m, cube_factor**2, x) != 0,
        "power M does not have cube-root order one",
    )

    y = sp.symbols("y")
    cheb_v = (y**2 - sp.Rational(1, 4)) ** 4 * (y**2 - 1) ** 3
    cheb_m = (4 * y**3 - 3 * y) * (y**2 - sp.Rational(1, 4)) * (y**2 - 1)
    for point in (-1, 1):
        shifted = sp.symbols("shifted")
        v_shift = sp.Poly(sp.expand(cheb_v.subs(y, shifted + point)), shifted)
        m_shift = sp.Poly(sp.expand(cheb_m.subs(y, shifted + point)), shifted)
        v_min = min(exponent[0] for exponent, coefficient in v_shift.terms() if coefficient)
        m_min = min(exponent[0] for exponent, coefficient in m_shift.terms() if coefficient)
        require((v_min, m_min) == (3, 1), "Chebyshev local orders are not (3,1)")


def main() -> None:
    phi, psi, k, h = build_degree_26_bank()

    h_expected = (
        72 * c18 * T**2
        - 4032 * c18 * T * d * s
        + 6720 * c18 * s**3
        + 896 * c14 * T * d
        - 4480 * c14 * s**2
        + 2560 * c10 * s
        - 1024 * c6
        - 143 * T**3 * d
        + 5148 * T**2 * d**2 * s
        + 1287 * T**2 * s**2
        - 24024 * T * d * s**3
        + 12012 * s**5
    ) / 4096
    require(sp.expand(h - h_expected) == 0, "independent recurrence does not recover H_26")
    require(sp.expand(k + d * phi / 2 - T * h) == 0, "top-row decomposition failed")

    h_a3_weight, _ = weighted_initial(h, (3, -3, 0))
    h_a2_weight, _ = weighted_initial(h, (1, -1, 0))
    require(h_a3_weight >= 0, "a>=3 H bound failed")
    require(h_a2_weight >= 0, "a=2,b>=1 H bound failed")

    t4_coefficient = sp.Poly(phi, T, d, s).coeff_monomial(T**4)
    require(t4_coefficient == sp.Rational(143, 16384), "wrong T^4 coefficient")
    for weights, expected_weight in (
        ((-3, 0, -1), -12),
        ((-3, 0, 0), -12),
        ((-1, 0, 0), -4),
    ):
        least, initial = weighted_initial(phi, weights)
        require(least == expected_weight, "regular Phi weight is wrong")
        require(
            sp.expand(initial - t4_coefficient * T**4) == 0,
            "T^4 is not the unique regular-lane initial monomial",
        )

    f = r**3 - 21 * r**2 + 35 * r - 7
    g = 7 * r**3 - 35 * r**2 + 21 * r - 1
    generic_weights = (
        ((0, 0), (-3, -3, -3), -18, -21),
        ((1, 0), (-1, -3, -2), -12, -14),
        ((2, 0), (1, -3, -1), -6, -7),
        ((1, 1), (-1, -1, -1), -6, -7),
    )
    for lane, weights, phi_weight, psi_weight in generic_weights:
        least_phi, initial_phi = weighted_initial(phi, weights)
        least_psi, initial_psi = weighted_initial(psi, weights)
        require((least_phi, least_psi) == (phi_weight, psi_weight), f"{lane}: bad weights")
        reduced_phi = generic_ratio_form(initial_phi, 6)
        reduced_psi = generic_ratio_form(initial_psi / s, 6)
        require(
            sp.factor(reduced_phi / f) == -sp.Rational(429, 512),
            f"{lane}: Phi initial form differs",
        )
        require(
            sp.factor(reduced_psi / g) == sp.Rational(429, 512),
            f"{lane}: Psi initial form differs",
        )
    require(sp.resultant(f, g, r) == -(2**21), "generic cubic resultant differs")

    least_phi, exceptional_phi = weighted_initial(phi, (-3, -1, -2))
    least_psi, exceptional_psi = weighted_initial(psi, (-3, -1, -2))
    least_h, exceptional_h = weighted_initial(h, (-3, -1, -2))
    require((least_phi, least_psi, least_h) == (-12, -14, -10), "exceptional weights differ")

    P = (
        143 * zeta**2
        - 13728 * r**3
        - 20592 * zeta * r
        + 288288 * r**2
        + 48048 * zeta
        - 480480 * r
        + 96096
    )
    G = (
        1287 * zeta**2
        + 5148 * zeta * r**2
        - 96096 * r**3
        - 72072 * zeta * r
        + 480480 * r**2
        + 60060 * zeta
        - 288288 * r
        + 13728
    )
    J = 143 * ((9 - r) * zeta + 36 * r**2 - 168 * r + 84)
    require(
        sp.expand(exceptional_ratio_form(exceptional_phi, 6) - P / 16384) == 0,
        "exceptional Phi initial form differs",
    )
    require(
        sp.expand(exceptional_ratio_form(exceptional_psi, 7) + G / 16384) == 0,
        "exceptional Psi initial form differs",
    )
    require(
        sp.expand(exceptional_ratio_form(exceptional_h, 5) - J / 4096) == 0,
        "exceptional H initial form differs",
    )
    unit_basis = sp.groebner([P, G, J], zeta, r, order="lex", domain=sp.QQ)
    require(
        len(unit_basis.polys) == 1 and unit_basis.polys[0].as_expr() == 1,
        "exceptional ideal is not the unit ideal",
    )
    require(J.subs(r, 9) != 0, "r=9 was not excluded")

    zeta_solution = 12 * (3 * r**2 - 14 * r + 7) / (r - 9)
    p = 2 * r**5 + 3 * r**4 - 488 * r**3 + 2842 * r**2 - 6930 * r + 4011
    g5 = (
        13 * r**5
        - 182 * r**4
        + 1974 * r**3
        - 8776 * r**2
        + 13173 * r
        - 5130
    )
    p_substitution = sp.together(P.subs(zeta, zeta_solution)).as_numer_denom()[0]
    g_substitution = sp.together(G.subs(zeta, zeta_solution)).as_numer_denom()[0]
    require(sp.expand(p_substitution + 6864 * p) == 0, "P substitution quintic differs")
    require(sp.expand(g_substitution - 6864 * g5) == 0, "G substitution quintic differs")
    expected_resultant = -(2**36) * 3**2 * 31**9 * 37
    require(
        sp.resultant(p, g5, r) == expected_resultant,
        "exceptional quintic resultant differs",
    )

    a0, b0, v0 = sp.symbols("a0 b0 v0", nonzero=True)
    leading_T = a0**2 / v0
    leading_d = -(b0**2) / (4 * v0)
    leading_s = a0 * b0 / (2 * v0)
    forced_r = sp.cancel(leading_T * leading_d / leading_s**2)
    require(forced_r == -1, "source leading coefficients do not force r=-1")
    require((f.subs(r, -1), g.subs(r, -1)) == (-64, -64), "forced generic ratio survives")
    forced_zeta = sp.solve(J.subs(r, -1), zeta)[0]
    require(forced_zeta == -sp.Rational(144, 5), "forced exceptional zeta differs")
    require(
        P.subs({r: -1, zeta: forced_zeta}) == -sp.Rational(24490752, 25)
        and G.subs({r: -1, zeta: forced_zeta}) == -sp.Rational(50189568, 25),
        "forced-r exceptional shortcut does not obstruct",
    )

    audit_case_partition()
    audit_carrier_multiplicities()

    print("THM-2823 INDEPENDENT DEGREE-26 TRIPLE-POLE VALUATION AUDIT")
    print("method=log-derivative Faber recurrence + weighted Newton initials; no primary import")
    print("degree26_bank=Phi,Psi,K reconstructed; K=-(d/2)(Phi/q)+T*H_26")
    print("valuation_partition=3 regular regions + 4 generic polar lanes + 1 exceptional lane")
    print("regular_gates=H bounds 0,0; unique Phi/q monomial (143/16384)T^4")
    print("generic_initials=s^6*f(r),s^7*g(r); resultant=-2^21")
    print("exceptional_initials=P,G,J; Groebner ideal=<1>")
    print("exceptional_quintic_resultant=-2^36*3^2*31^9*37")
    print("stronger_source_relation=r=-1; f(-1)=g(-1)=-64; exceptional direct gate=nonzero")
    print("carrier_local_orders=power and Chebyshev both have (ord V,ord M)=(3,1)")
    print("ALL INDEPENDENT EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
