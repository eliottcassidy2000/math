#!/usr/bin/env python3
"""Exact companion for THM-2371's H2 common-root elimination.

The script works in THM-2357's gauge where the unique three-cycle is at
``y=1``.  It proves that the two hypotheses used by THM-2360,

    Res_y(p3,q5) != 0  and  q5(1) != 0,

are automatic on the ``H2*S4^2`` coefficient locus.  The proof uses three
small univariate subresultant certificates.  Each rational coprimality
claim is also certified by a nonzero Sylvester determinant over a prime
field, computed by a custom modular determinant routine.
"""

from __future__ import annotations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primitive_coefficient(
    polynomial: sp.Expr,
    y: sp.Symbol,
    power: int,
    parameter: sp.Symbol,
) -> sp.Poly:
    coefficient = sp.Poly(polynomial, y).coeff_monomial(y**power)
    numerator = sp.together(coefficient).as_numer_denom()[0]
    content, primitive = sp.Poly(numerator, parameter).primitive()
    require(content != 0, "subresultant coefficient unexpectedly vanished")
    return primitive


def subresultant_of_degree(
    polynomial: sp.Expr,
    y: sp.Symbol,
    degree: int,
) -> tuple[sp.Expr, list[int]]:
    sequence = sp.subresultants(polynomial, sp.diff(polynomial, y), y)
    profile = [sp.degree(item, y) for item in sequence]
    matches = [item for item in sequence if sp.degree(item, y) == degree]
    require(len(matches) == 1, f"subresultant degree {degree} is not unique")
    return matches[0], profile


def factor_order(polynomial: sp.Poly, factor: sp.Poly) -> tuple[int, sp.Poly]:
    order = 0
    quotient = polynomial
    while sp.rem(quotient, factor).is_zero:
        quotient = sp.exquo(quotient, factor)
        order += 1
    return order, quotient


def determinant_mod(matrix: list[list[int]], prime: int) -> int:
    values = [[entry % prime for entry in row] for row in matrix]
    determinant = 1
    size = len(values)
    for column in range(size):
        pivot = next(
            (
                row
                for row in range(column, size)
                if values[row][column] % prime != 0
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            values[column], values[pivot] = values[pivot], values[column]
            determinant = -determinant
        pivot_value = values[column][column] % prime
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, size):
            multiplier = values[row][column] * inverse % prime
            if multiplier:
                values[row] = [
                    (
                        values[row][index]
                        - multiplier * values[column][index]
                    )
                    % prime
                    for index in range(size)
                ]
    return determinant % prime


def sylvester_resultant_mod(
    first: sp.Poly,
    second: sp.Poly,
    variable: sp.Symbol,
    prime: int,
) -> int:
    first_mod = sp.Poly(first, variable, modulus=prime)
    second_mod = sp.Poly(second, variable, modulus=prime)
    first_degree = first_mod.degree()
    second_degree = second_mod.degree()
    require(
        first_degree == first.degree() and second_degree == second.degree(),
        "degree dropped in a modular certificate",
    )
    first_coefficients = [
        int(coefficient) % prime for coefficient in first_mod.all_coeffs()
    ]
    second_coefficients = [
        int(coefficient) % prime for coefficient in second_mod.all_coeffs()
    ]
    matrix: list[list[int]] = []
    for shift in range(second_degree):
        matrix.append(
            [0] * shift
            + first_coefficients
            + [0] * (second_degree - 1 - shift)
        )
    for shift in range(first_degree):
        matrix.append(
            [0] * shift
            + second_coefficients
            + [0] * (first_degree - 1 - shift)
        )
    custom = determinant_mod(matrix, prime)
    library = int(
        sp.resultant(first_mod, second_mod, variable)
    ) % prime
    require(custom == library, "custom and library modular resultants differ")
    return custom


def even_quotient(
    polynomial: sp.Poly,
    old_variable: sp.Symbol,
    new_variable: sp.Symbol,
) -> sp.Poly:
    result = sp.Integer(0)
    for (power,), coefficient in polynomial.terms():
        require(power % 2 == 0, "polynomial stopped being even")
        result += coefficient * new_variable ** (power // 2)
    return sp.Poly(result, new_variable).primitive()[1]


def reciprocal_quotient(
    polynomial: sp.Poly,
    old_variable: sp.Symbol,
    new_variable: sp.Symbol,
    half_degree: int,
) -> sp.Poly:
    expected_degree = 2 * half_degree
    require(polynomial.degree() == expected_degree, "wrong reciprocal degree")
    reciprocal = sp.expand(
        old_variable**expected_degree
        * polynomial.as_expr().subs(old_variable, 1 / old_variable)
    )
    require(
        reciprocal == polynomial.as_expr(),
        "obstruction stopped being reciprocal",
    )
    traces = [sp.Integer(2), new_variable]
    for _ in range(2, half_degree + 1):
        traces.append(sp.expand(new_variable * traces[-1] - traces[-2]))
    quotient = polynomial.coeff_monomial(old_variable**half_degree)
    for offset in range(1, half_degree + 1):
        quotient += (
            polynomial.coeff_monomial(old_variable ** (half_degree + offset))
            * traces[offset]
        )
    quotient_poly = sp.Poly(sp.expand(quotient), new_variable).primitive()[1]
    reconstructed = sp.cancel(
        old_variable**half_degree
        * quotient_poly.as_expr().subs(
            new_variable,
            old_variable + 1 / old_variable,
        )
        - polynomial.as_expr()
    )
    require(reconstructed == 0, "reciprocal quotient failed to reconstruct")
    return quotient_poly


def main() -> None:
    y, bvar, cvar, avar = sp.symbols("y B C a")

    dvar, wvar = sp.symbols("D W")
    quartic_p = (
        245 * y**4
        + 1890 * bvar * y**2
        - 24300 * bvar**2
        + 122472 * dvar
    )
    sextic_q = (
        539 * y**6
        + 11340 * bvar * y**4
        + 183708 * cvar * y**3
        + (72900 * bvar**2 - 367416 * dvar) * y**2
        + (2361960 * bvar * cvar + 2480058 * wvar) * y
    )
    d_incidence = (
        24300 * bvar**2 - 1890 * bvar - 245
    ) / sp.Integer(122472)
    w_incidence = (
        -sp.Rational(20, 21) * bvar * cvar
        - sp.Rational(2, 27) * cvar
        - (91 + 1215 * bvar) / sp.Integer(177147)
    )
    quartic_incidence = sp.expand(quartic_p.subs(dvar, d_incidence))
    sextic_incidence = sp.expand(
        sextic_q.subs({dvar: d_incidence, wvar: w_incidence})
    )
    p3 = 35 * (y + 1) * (54 * bvar + 7 * y**2 + 7)
    q5 = 7 * y * (
        1620 * bvar * y**2
        + 1620 * bvar * y
        + 2430 * bvar
        + 26244 * cvar * (y + 1)
        + 77 * (y**4 + y**3 + y**2 + y)
        + 182
    )
    require(
        sp.expand(quartic_incidence - (y - 1) * p3) == 0
        and sp.expand(sextic_incidence - (y - 1) * q5) == 0,
        "moving-root divisions stopped matching the degree-eighteen covariants",
    )
    residual = sp.expand(4 * (y - 1) * p3**3 + 49 * q5**2)
    full_mordell = sp.expand(
        4 * quartic_incidence**3 + 49 * sextic_incidence**2
    )
    require(
        sp.degree(p3, y) == 3
        and sp.Poly(p3, y).LC() == 245
        and sp.degree(q5, y) == 5
        and sp.Poly(q5, y).LC() == 539
        and sp.degree(residual, y) == 10
        and sp.Poly(residual, y).LC() == 73060029,
        "moving-root degrees or leading coefficients changed",
    )
    require(
        sp.expand(full_mordell - (y - 1) ** 2 * residual) == 0,
        "degree-ten residual stopped matching the full Mordell polynomial",
    )
    rvar, zvar_scale = sp.symbols("r z", nonzero=True)
    d_moving = (
        24300 * bvar**2
        - 1890 * bvar * rvar**2
        - 245 * rvar**4
    ) / sp.Integer(122472)
    w_moving = (
        -sp.Rational(20, 21) * bvar * cvar
        - sp.Rational(2, 27) * cvar * rvar**2
        - rvar**3 * (91 * rvar**2 + 1215 * bvar) / sp.Integer(177147)
    )
    p_moving = sp.cancel(
        quartic_p.subs(dvar, d_moving) / (y - rvar)
    )
    q_moving = sp.cancel(
        sextic_q.subs({dvar: d_moving, wvar: w_moving}) / (y - rvar)
    )
    normalized_p = p3.subs(
        {
            y: zvar_scale,
            bvar: bvar / rvar**2,
            cvar: cvar / rvar**3,
        }
    )
    normalized_q = q5.subs(
        {
            y: zvar_scale,
            bvar: bvar / rvar**2,
            cvar: cvar / rvar**3,
        }
    )
    require(
        sp.factor(
            p_moving
            - rvar**3 * normalized_p.subs(zvar_scale, y / rvar)
        )
        == 0
        and sp.factor(
            q_moving
            - rvar**5 * normalized_q.subs(zvar_scale, y / rvar)
        )
        == 0,
        "weighted moving-root covariance changed",
    )
    smooth = 245 + 2835 * bvar + 26244 * cvar
    require(
        p3.subs(y, 0) == 35 * (54 * bvar + 7)
        and q5.subs(y, 0) == 0
        and p3.subs(y, -1) == 0
        and q5.subs(y, -1) == -14 * (1215 * bvar + 91)
        and p3.subs(y, 1) == 140 * (27 * bvar + 7)
        and q5.subs(y, 1) == 14 * smooth,
        "distinguished root evaluations changed",
    )

    psi = (
        767400804 * bvar**4
        - 172777374 * bvar**3
        + 1750329 * bvar**2
        + 65086642152 * bvar * cvar**2
        + 1562436540 * bvar * cvar
        + 6260436 * bvar
        + 16874314632 * cvar**2
        + 117021996 * cvar
        + 405769
    )
    expected_resultant = (
        7061881225000
        * (54 * bvar + 7)
        * (1215 * bvar + 91)
        * psi
    )
    require(
        sp.factor(sp.resultant(p3, q5, y) - expected_resultant) == 0,
        "common-root resultant atlas changed",
    )
    factorization = sp.factor_list(psi, bvar, cvar)
    require(
        factorization == (1, [(psi, 1)]),
        "generic resultant component stopped being irreducible over Q",
    )
    require(
        sp.factor(sp.discriminant(psi, cvar))
        == -4821232752 * (54 * bvar + 7) ** 3 * (513 * bvar - 91) ** 2,
        "generic component discriminant changed",
    )

    # The q5(1)=0 singular-order wall.
    c_singular = -(245 + 2835 * bvar) / sp.Integer(26244)
    singular_residual = sp.factor(residual.subs(cvar, c_singular))
    t9 = sp.Poly(
        sp.cancel(singular_residual / (343 * (y - 1))),
        y,
    )
    expected_t9 = sp.Poly(
        213003 * y**9
        + 639009 * y**8
        + 7938 * (720 * bvar + 161) * y**7
        + 490 * (28755 * bvar + 3808) * y**6
        + 105 * (466560 * bvar**2 + 209790 * bvar + 19061) * y**5
        + 21 * (3936600 * bvar**2 + 1162350 * bvar + 83447) * y**4
        + (
            78732000 * bvar**3
            + 114434775 * bvar**2
            + 25401600 * bvar
            + 1539041
        )
        * y**3
        + 3
        * (
            78732000 * bvar**3
            + 40441275 * bvar**2
            + 6495930 * bvar
            + 333739
        )
        * y**2
        + 1500 * (54 * bvar + 7) ** 3 * y
        + 500 * (54 * bvar + 7) ** 3,
        y,
    )
    require(
        t9 == expected_t9
        and t9.degree() == 9
        and t9.LC() == 213003
        and sp.factor(t9.eval(1)) == 32000 * (27 * bvar + 7) ** 3,
        "singular-order forced-root quotient changed",
    )
    t9_sres3, t9_profile = subresultant_of_degree(
        t9.as_expr(),
        y,
        3,
    )
    t9_first = primitive_coefficient(t9_sres3, y, 3, bvar)
    t9_second = primitive_coefficient(t9_sres3, y, 2, bvar)
    require(
        t9_profile == list(range(9, -1, -1))
        and t9_first.degree() == 15
        and t9_second.degree() == 15
        and sp.gcd(t9_first, t9_second).degree() == 0,
        "singular-order subresultant obstruction changed",
    )
    singular_modulus = 11
    singular_mod_resultant = sylvester_resultant_mod(
        t9_first,
        t9_second,
        bvar,
        singular_modulus,
    )
    require(
        singular_mod_resultant == 6,
        "singular-order modular resultant changed",
    )

    exceptional_b = -sp.Rational(7, 27)
    exceptional_sextic = (
        621 * y**6
        + 3726 * y**5
        + 8721 * y**4
        + 10396 * y**3
        + 7536 * y**2
        + 3000 * y
        + 500
    )
    exceptional_residual = sp.factor(
        singular_residual.subs(bvar, exceptional_b)
    )
    require(
        exceptional_residual
        == 7**6 * (y - 1) ** 4 * exceptional_sextic,
        "singular-order exceptional factorization changed",
    )
    require(
        exceptional_sextic.subs(y, 1) == 34500,
        "exceptional sextic met the forced root",
    )
    exceptional_gcd = sp.gcd(
        sp.Poly(exceptional_residual, y),
        sp.Poly(sp.diff(exceptional_residual, y), y),
    ).monic()
    require(
        sp.expand(exceptional_gcd.as_expr() - (y - 1) ** 3) == 0,
        "singular-order exceptional gcd changed",
    )
    exceptional_mod_resultant = sylvester_resultant_mod(
        sp.Poly(exceptional_sextic, y),
        sp.Poly(sp.diff(exceptional_sextic, y), y),
        y,
        7,
    )
    require(
        exceptional_mod_resultant == 2,
        "exceptional sextic squarefreeness certificate changed",
    )

    # The fixed antipodal common-root component.
    antipodal_b = -sp.Rational(91, 1215)
    e8 = (
        22182741 * y**8
        - 214326 * y**6
        + 2946308904 * cvar * y**5
        - 11638431 * y**4
        - 1696359672 * cvar * y**3
        + (502096953744 * cvar**2 - 8207696) * y**2
        - 1344364
    )
    require(
        sp.factor(
            residual.subs(bvar, antipodal_b)
            - sp.Rational(2401, 729) * (y + 1) ** 2 * e8
        )
        == 0
        and sp.degree(e8, y) == 8
        and sp.Poly(e8, y).LC() == 22182741,
        "antipodal forced-root quotient changed",
    )
    e8_sres2, e8_profile = subresultant_of_degree(e8, y, 2)
    e8_first = primitive_coefficient(e8_sres2, y, 2, cvar)
    e8_second = primitive_coefficient(e8_sres2, y, 1, cvar)
    require(
        e8_profile == list(range(8, -1, -1))
        and sp.expand(e8_first.as_expr().subs(cvar, -cvar))
        == e8_first.as_expr()
        and sp.expand(e8_second.as_expr().subs(cvar, -cvar))
        == -e8_second.as_expr()
        and e8_second.eval(0) == 0,
        "antipodal parity reduction changed",
    )
    xvar = sp.symbols("X")
    antipodal_first = even_quotient(e8_first, cvar, xvar)
    e8_second_over_c = sp.exquo(e8_second, sp.Poly(cvar, cvar))
    antipodal_second = even_quotient(e8_second_over_c, cvar, xvar)
    require(
        antipodal_first.degree() == 5
        and antipodal_second.degree() == 4
        and antipodal_first.eval(0) != 0
        and int(
            sp.Poly(antipodal_first, xvar, modulus=17).eval(0)
        )
        % 17
        == 11
        and sp.gcd(antipodal_first, antipodal_second).degree() == 0,
        "antipodal reduced obstruction changed",
    )
    antipodal_modulus = 17
    antipodal_mod_resultant = sylvester_resultant_mod(
        antipodal_first,
        antipodal_second,
        xvar,
        antipodal_modulus,
    )
    require(
        antipodal_mod_resultant == 5,
        "antipodal modular resultant changed",
    )

    # The irreducible moving quadratic-root component.
    b_of_a = -sp.Rational(7, 54) * (avar**2 + 1)
    c_of_a = (
        7
        * (
            19 * avar**4
            + 19 * avar**3
            + 64 * avar**2
            + 19 * avar
            + 19
        )
        / (sp.Integer(26244) * (avar + 1))
    )
    q_at_a = sp.factor(q5.subs({y: avar, bvar: b_of_a}))
    require(
        q_at_a
        == -7
        * avar
        * (
            -26244 * cvar * avar
            - 26244 * cvar
            + 133 * avar**4
            + 133 * avar**3
            + 448 * avar**2
            + 133 * avar
            + 133
        )
        and sp.factor(psi.subs({bvar: b_of_a, cvar: c_of_a})) == 0,
        "moving common-root parametrization changed",
    )
    require(
        sp.factor(
            smooth.subs({bvar: b_of_a, cvar: c_of_a})
            - 7
            * (avar - 1) ** 2
            * (38 * avar**2 + 9 * avar + 3)
            / (2 * (avar + 1))
        )
        == 0,
        "moving component's singular-order intersections changed",
    )
    u8 = (
        621 * (avar + 1) ** 2 * y**8
        + 1242 * (avar + 1) ** 3 * y**7
        - 27
        * (avar + 1) ** 2
        * (11 * avar**2 - 92 * avar + 11)
        * y**6
        - 2
        * (avar + 1)
        * (
            709 * avar**4
            + 1006 * avar**3
            - 110 * avar**2
            + 1006 * avar
            + 709
        )
        * y**5
        - (avar + 1) ** 2
        * (
            139 * avar**4
            + 2836 * avar**3
            - 1415 * avar**2
            + 2836 * avar
            + 139
        )
        * y**4
        - 2
        * avar
        * (avar + 1)
        * (
            139 * avar**4
            + 1203 * avar**3
            - 1008 * avar**2
            + 1203 * avar
            + 139
        )
        * y**3
        + avar**2
        * (
            361 * avar**4
            + 1012 * avar**3
            + 5398 * avar**2
            + 1012 * avar
            + 361
        )
        * y**2
        + 1000 * avar**3 * (avar + 1) ** 3 * y
        + 500 * avar**4 * (avar + 1) ** 2
    )
    generic_residual = residual.subs({bvar: b_of_a, cvar: c_of_a})
    require(
        sp.factor(
            (avar + 1) ** 2 * generic_residual
            - 7**6 * (y - avar) ** 2 * u8
        )
        == 0
        and sp.degree(u8, y) == 8
        and sp.expand(
            sp.Poly(u8, y).LC() - 621 * (avar + 1) ** 2
        )
        == 0,
        "generic forced-root quotient changed",
    )
    u8_sres2, u8_profile = subresultant_of_degree(u8, y, 2)
    u8_first_raw = primitive_coefficient(u8_sres2, y, 2, avar)
    u8_second_raw = primitive_coefficient(u8_sres2, y, 1, avar)
    first_a_order, first_after_a = factor_order(
        u8_first_raw,
        sp.Poly(avar, avar),
    )
    first_minus_order, generic_first = factor_order(
        first_after_a,
        sp.Poly(avar + 1, avar),
    )
    second_a_order, second_after_a = factor_order(
        u8_second_raw,
        sp.Poly(avar, avar),
    )
    second_minus_order, generic_second = factor_order(
        second_after_a,
        sp.Poly(avar + 1, avar),
    )
    generic_first = generic_first.primitive()[1]
    generic_second = generic_second.primitive()[1]
    require(
        u8_profile == list(range(8, -1, -1))
        and (first_a_order, first_minus_order) == (2, 12)
        and (second_a_order, second_minus_order) == (3, 13)
        and generic_first.degree() == 36
        and generic_second.degree() == 34,
        "generic subresultant factor ledger changed",
    )
    zvar = sp.symbols("z")
    generic_first_quotient = reciprocal_quotient(
        generic_first,
        avar,
        zvar,
        18,
    )
    generic_second_quotient = reciprocal_quotient(
        generic_second,
        avar,
        zvar,
        17,
    )
    require(
        generic_first_quotient.degree() == 18
        and generic_second_quotient.degree() == 17
        and sp.gcd(
            generic_first_quotient,
            generic_second_quotient,
        ).degree()
        == 0,
        "generic reciprocal obstruction changed",
    )
    generic_modulus = 11
    generic_mod_resultant = sylvester_resultant_mod(
        generic_first_quotient,
        generic_second_quotient,
        zvar,
        generic_modulus,
    )
    require(
        generic_mod_resultant == 9,
        "generic modular resultant changed",
    )

    # Hostile controls: each raw wall is genuinely nonempty before imposing
    # the H2 square-class multiplicity.
    k_control = sp.Poly(
        residual.subs(
            {
                bvar: 0,
                cvar: -sp.Rational(245, 26244),
            }
        ),
        y,
    )
    antipodal_control = sp.Poly(
        residual.subs({bvar: antipodal_b, cvar: 0}),
        y,
    )
    a_control = sp.Rational(2)
    generic_control = sp.Poly(
        residual.subs(
            {
                bvar: b_of_a.subs(avar, a_control),
                cvar: c_of_a.subs(avar, a_control),
            }
        ),
        y,
    )
    control_gcd_degrees = (
        sp.gcd(k_control, k_control.diff()).degree(),
        sp.gcd(antipodal_control, antipodal_control.diff()).degree(),
        sp.gcd(generic_control, generic_control.diff()).degree(),
    )
    require(
        control_gcd_degrees == (0, 1, 1),
        "wall hostile controls changed",
    )
    positive_h2 = (y**2 + y + 1) * (y**4 + 2 * y + 3) ** 2
    positive_u8 = (y**2 + y + 1) * (y**3 + y + 2) ** 2
    require(
        sp.gcd(
            sp.Poly(positive_h2, y),
            sp.Poly(sp.diff(positive_h2, y), y),
        ).degree()
        >= 4
        and sp.gcd(
            sp.Poly(positive_u8, y),
            sp.Poly(sp.diff(positive_u8, y), y),
        ).degree()
        >= 3,
        "positive square-class controls changed",
    )

    print("THM-2371 degree-18 H2 common-root elimination exact companion")
    print("gauge: unique three-cycle at y=1")
    print("degrees/lc: p3=(3,245), q5=(5,539), R10=(10,73060029)")
    print(
        "resultant components: "
        "(54B+7)*(1215B+91)*Psi(B,C); Psi irreducible over Q"
    )
    print(
        "Psi discriminant in C: "
        "-4821232752*(54B+7)^3*(513B-91)^2"
    )
    print(
        "K=0 quotient: R10=343*(y-1)*T9, "
        "T9(1)=32000*(27B+7)^3"
    )
    print(
        "K=0 generic obstruction: degrees (15,15), "
        f"resultant mod {singular_modulus}={singular_mod_resultant}"
    )
    print(
        "K=0 exceptional B=-7/27: gcd degree 3; "
        f"sextic resultant mod 7={exceptional_mod_resultant}"
    )
    print(
        "antipodal quotient: R10=(2401/729)*(y+1)^2*E8"
    )
    print(
        "antipodal X=C^2 obstruction: degrees (5,4), "
        f"resultant mod {antipodal_modulus}={antipodal_mod_resultant}"
    )
    print(
        "generic component: "
        "B=-7*(a^2+1)/54, "
        "C=7*(19a^4+19a^3+64a^2+19a+19)/(26244*(a+1))"
    )
    print(
        "generic quotient: (a+1)^2*R10=7^6*(y-a)^2*U8, "
        "lc(U8)=621*(a+1)^2"
    )
    print(
        "generic raw factor orders: "
        "y2 -> a^2*(a+1)^12; y1 -> a^3*(a+1)^13"
    )
    print(
        "generic z=a+a^-1 obstruction: degrees (18,17), "
        f"resultant mod {generic_modulus}={generic_mod_resultant}"
    )
    print(f"hostile wall gcd degrees: {control_gcd_degrees}")
    print("VERDICT: Res(p3,q5)!=0 and q5(1)!=0 on every H2 point")
    print("THM-2360 cube descent is unconditional on the H2 stratum")
    print("OPEN: the coprime linear-times-cube coefficient locus and H4")


if __name__ == "__main__":
    main()
