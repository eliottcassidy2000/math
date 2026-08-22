#!/usr/bin/env python3
"""Exact companion for THM-3686.

This script audits the ``y=0`` restriction of the weighted quartic Keller
map from THM-3438.  It verifies the image equation, the source immersion and
an explicit Bezout identity, the cubic normalization presentation and its
smoothness, the retained collision on distinct normalization sheets, and
explicit three- and two-bracket decompositions of one.

Every truth gate uses ``require`` and therefore remains active under
``python -O``.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def bracket(first: sp.Expr, second: sp.Expr, x: sp.Symbol, z: sp.Symbol) -> sp.Expr:
    return sp.expand(sp.diff(first, x) * sp.diff(second, z) - sp.diff(first, z) * sp.diff(second, x))


def main() -> None:
    x, z = sp.symbols("x z")
    t = x**2 * z
    A = sp.expand(3 * z * (2 - t))
    B = sp.expand(2 * x * z * (2 - t))
    C = sp.expand(x * (1 - t))
    g = 1 - t

    image_relation = sp.expand(16 * A**3 * C**2 + 27 * B**4 - 36 * A * B**2)
    require(image_relation == 0, ("image relation", sp.factor(image_relation)))

    J_AB = sp.factor(bracket(A, B, x, z))
    J_AC = sp.factor(bracket(A, C, x, z))
    J_BC = sp.factor(bracket(B, C, x, z))
    expected_minors = (
        -12 * z * (t - 2) * (t - 1),
        -6 * (2 * t**2 - 4 * t + 1),
        -2 * x * (3 * t**2 - 6 * t + 2),
    )
    require(tuple(sp.expand(left - right) for left, right in zip((J_AB, J_AC, J_BC), expected_minors)) == (0, 0, 0), "immersion minors")

    bezout = sp.expand((-sp.Rational(1, 6) + 2 * t - t**2) * J_AC + 2 * x * z * (t - 2) * J_BC)
    require(bezout == 1, ("source Bezout identity", sp.factor(bezout)))

    # Abstract normalization presentation.
    aa, bb, cc, gg = sp.symbols("A B C g")
    normalization_relations = (
        bb * cc - 2 * gg * (1 - gg**2),
        3 * bb * gg - 2 * aa * cc,
        4 * aa * (1 - gg**2) - 3 * bb**2,
    )
    substitution = {aa: A, bb: B, cc: C, gg: g}
    require(tuple(sp.expand(relation.subs(substitution)) for relation in normalization_relations) == (0, 0, 0), "normalization relations")
    require(sp.expand(g**3 - g + B * C / 2) == 0, "integral cubic")
    require(sp.expand(3 * B * g - 2 * A * C) == 0, "birational sheet coordinate")

    variables = (aa, bb, cc, gg)
    relation_jacobian = sp.Matrix(
        [[sp.diff(relation, variable) for variable in variables] for relation in normalization_relations]
    )
    two_minors = [
        sp.expand(relation_jacobian.extract(rows, columns).det())
        for rows in combinations(range(3), 2)
        for columns in combinations(range(4), 2)
    ]
    three_minors = [
        sp.expand(relation_jacobian[:, columns].det())
        for columns in combinations(range(4), 3)
    ]
    relation_groebner = sp.groebner(normalization_relations, gg, cc, bb, aa, order="lex")
    require(all(relation_groebner.reduce(minor)[1] == 0 for minor in three_minors), "normalization rank exceeded two")
    smooth_groebner = sp.groebner((*normalization_relations, *two_minors), gg, cc, bb, aa, order="grevlex")
    require(list(smooth_groebner) == [1], ("normalization singular ideal", list(smooth_groebner)))
    require(2 * cc**2 in two_minors and 6 * bb**2 in two_minors, "rank forcing minors missing")

    collision_points = ((1, 0), (-1, 2))
    collision_rows = []
    for point in collision_points:
        values = {x: point[0], z: point[1]}
        collision_rows.append(tuple(sp.expand(F).subs(values) for F in (A, B, C, g)))
    require(tuple(collision_rows) == ((0, 0, 1, 1), (0, 0, 1, -1)), collision_rows)

    # The image hypersurface is singular on the two coordinate axes.  The
    # collision lies on the C-axis, and its cubic normalization fibre has
    # the three distinct sheet values 0,+1,-1.
    H = 16 * aa**3 * cc**2 + 27 * bb**4 - 36 * aa * bb**2
    H_gradient = tuple(sp.diff(H, variable) for variable in (aa, bb, cc))
    require(all(sp.expand(entry.subs({aa: 0, bb: 0})) == 0 for entry in (H, *H_gradient)), "C-axis singularity")
    require(all(sp.expand(entry.subs({bb: 0, cc: 0})) == 0 for entry in (H, *H_gradient)), "A-axis singularity")
    fibre_polynomial = sp.factor((gg**3 - gg + bb * cc / 2).subs({aa: 0, bb: 0, cc: 1}))
    require(fibre_polynomial == gg * (gg - 1) * (gg + 1), fibre_polynomial)

    # Canonical-form divisor on the normalization boundary.  Near D0 use
    # (C,g); near Dminus use (A,B), with dg obtained from f3.
    x_cg = cc / gg
    z_cg = gg**2 * (1 - gg) / cc**2
    omega_cg = sp.factor(
        sp.diff(x_cg, cc) * sp.diff(z_cg, gg) - sp.diff(x_cg, gg) * sp.diff(z_cg, cc)
    )
    require(omega_cg == -gg / cc**2, ("omega on C,g chart", omega_cg))

    g_A = (1 - gg**2) / (2 * aa * gg)
    g_B = -3 * bb / (4 * aa * gg)
    x_ab = 3 * bb / (2 * aa)
    z_ag = aa / (3 * (1 + gg))
    x_A = sp.diff(x_ab, aa)
    x_B = sp.diff(x_ab, bb)
    z_A = sp.diff(z_ag, aa) + sp.diff(z_ag, gg) * g_A
    z_B = sp.diff(z_ag, bb) + sp.diff(z_ag, gg) * g_B
    omega_ab = sp.factor(x_A * z_B - x_B * z_A)
    omega_ab_reduced = sp.factor(omega_ab.subs(bb**2, 4 * aa * (1 - gg**2) / 3))
    require(omega_ab_reduced == -1 / (4 * aa * gg), ("omega on A,B chart", omega_ab_reduced))

    # A short additive Poisson certificate in the target subalgebra.
    Q0 = -B**2 / 4 + A / 6 + A**2 * C**2 / 6
    three_brackets = (
        bracket(C, Q0, x, z),
        sp.Rational(7, 4) * bracket(B * C**2, A * B * C, x, z),
        -sp.Rational(13, 8) * bracket(A * C**3, B**2, x, z),
    )
    require(sp.expand(sum(three_brackets)) == 1, ("three-bracket identity", sp.factor(sum(three_brackets))))

    degree_three_four_brackets = (
        bracket(C, A / 6 - B**2 / 4, x, z),
        bracket(C**2, A**2 * C, x, z) / 12,
        -sp.Rational(13, 6) * bracket(C**3, A * B**2, x, z),
        5 * bracket(B * C**2, A * B * C, x, z),
    )
    require(sp.expand(sum(degree_three_four_brackets)) == 1, ("degree-three four-bracket identity", sp.factor(sum(degree_three_four_brackets))))

    # Rank-two compression of the coefficient tensor.  The two scalar
    # equations define nonempty coefficient data over C and force b != 0.
    a, b = sp.symbols("a b")
    coefficient_relations = (108 * a**2 + 2040 * a - 403, 23324 * a + 920 * b**2 - 4459)
    forced_a_at_b_zero = sp.Rational(4459, 23324)
    require(sp.expand(coefficient_relations[0].subs(a, forced_a_at_b_zero)) != 0, "coefficient system might allow b=0")

    alpha = (6 * a + 104) / 17
    beta = (84 * a - 91) / (102 * b)
    P1 = C + alpha * A * C**3
    Q1 = A / 6 - B**2 / 4 + (sp.Rational(1, 6) - 3 * a / 2) * A**2 * C**2 - 12 * b * A * B * C / 7
    P2 = B * C**2 + beta * A * C**3
    Q2 = 8 * b * A / 7 + 3 * b * B**2 / 28 - 15 * b * A**2 * C**2 / 14 + sp.Rational(7, 4) * A * B * C
    two_bracket_error = sp.together(bracket(P1, Q1, x, z) + bracket(P2, Q2, x, z) - 1)
    error_numerator, error_denominator = sp.fraction(two_bracket_error)
    require(error_denominator != 0, "zero denominator")
    coefficient_groebner = sp.groebner(coefficient_relations, b, a, order="lex", domain=sp.QQ)
    error_poly = sp.Poly(sp.expand(error_numerator), x, z)
    coefficient_remainders = [
        coefficient_groebner.reduce(coefficient)[1] for coefficient in error_poly.coeffs()
    ]
    require(all(remainder == 0 for remainder in coefficient_remainders), ("two-bracket remainder", coefficient_remainders))

    # Verify the rank-two factorization used to obtain P_i,Q_i.
    coefficient_matrix = sp.Matrix(
        [
            [sp.Rational(1, 6), -sp.Rational(1, 4), sp.Rational(1, 6) - 3 * a / 2, -12 * b / 7],
            [8 * b / 7, 3 * b / 28, -15 * b / 14, sp.Rational(7, 4)],
            [a, -sp.Rational(13, 8), 0, b],
        ]
    )
    rank_error = coefficient_matrix.row(2) - alpha * coefficient_matrix.row(0) - beta * coefficient_matrix.row(1)
    rank_remainders = []
    for entry in rank_error:
        numerator = sp.fraction(sp.together(entry))[0]
        rank_remainders.append(coefficient_groebner.reduce(numerator)[1])
    require(all(remainder == 0 for remainder in rank_remainders), ("rank-two factorization", rank_remainders))

    # Constant-linear compression inside the four displayed functions is
    # impossible.  Represent the coefficient field in the Q-basis
    # (1,a,b,ab), solve for a general alternating two-form, and check that
    # its unique solution has nonzero Pluecker scalar.
    four_functions = (P1, Q1, sp.expand(102 * b * P2), Q2)
    pairs = tuple(combinations(range(4), 2))
    wedge_parameters = sp.symbols("w0:24")
    field_basis = (sp.Integer(1), a, b, a * b)
    reduced_field_wedges: list[sp.Expr] = []
    for left, right in pairs:
        source_wedge = bracket(four_functions[left], four_functions[right], x, z)
        for basis_element in field_basis:
            reduced_source_wedge = 0
            for (power_x, power_z), coefficient in sp.Poly(sp.expand(basis_element * source_wedge), x, z).terms():
                remainder = coefficient_groebner.reduce(sp.expand(coefficient))[1]
                reduced_source_wedge += remainder * x**power_x * z**power_z
            reduced_field_wedges.append(sp.expand(reduced_source_wedge))
    general_wedge_error = sp.expand(
        sum(parameter * source_wedge for parameter, source_wedge in zip(wedge_parameters, reduced_field_wedges)) - 1
    )
    wedge_equations = sp.Poly(general_wedge_error, x, z, a, b).coeffs()
    wedge_matrix, wedge_vector = sp.linear_eq_to_matrix(wedge_equations, wedge_parameters)
    require(wedge_matrix.rank() == 24, ("four-function wedge solution not unique", wedge_matrix.rank()))
    require(wedge_matrix.row_join(wedge_vector).rank() == 24, "four-function wedge system inconsistent")

    inverse_102b = -sp.Rational(144, 31213) * a * b - sp.Rational(46708, 530621) * b
    inverse_check_numerator = sp.fraction(sp.together(1 / (102 * b) - inverse_102b))[0]
    require(coefficient_groebner.reduce(sp.expand(inverse_check_numerator))[1] == 0, "inverse of 102b")
    expected_wedges = (sp.Integer(1), 0, 0, 0, 0, inverse_102b)
    expected_parameter_values: list[sp.Expr] = []
    for coefficient in expected_wedges:
        polynomial = sp.Poly(coefficient, a, b)
        expected_parameter_values.extend(
            (
                polynomial.coeff_monomial(1),
                polynomial.coeff_monomial(a),
                polynomial.coeff_monomial(b),
                polynomial.coeff_monomial(a * b),
            )
        )
    expected_vector = sp.Matrix(expected_parameter_values)
    require(wedge_matrix * expected_vector == wedge_vector, "four-function wedge candidate failed")
    pluecker = sp.expand(expected_wedges[0] * expected_wedges[5] - expected_wedges[1] * expected_wedges[4] + expected_wedges[2] * expected_wedges[3])
    require(pluecker == inverse_102b and pluecker != 0, ("Pluecker scalar", pluecker))

    print("theorem=THM-3686-y0-collision-normalization-and-bracket-anatomy")
    print("restriction=A=3z(2-x^2z);B=2xz(2-x^2z);C=x(1-x^2z)")
    print("kernel=16A^3C^2+27B^4-36AB^2")
    print(f"immersion_minors={(J_AB, J_AC, J_BC)}")
    print("source_bezout=1=(-1/6+2t-t^2){A,C}+2xz(t-2){B,C};t=x^2z")
    print("normalization_generator=g=1-x^2z=2AC/(3B);cubic=g^3-g+BC/2")
    print("normalization_relations=BC-2g(1-g^2);3Bg-2AC;4A(1-g^2)-3B^2")
    print("normalization_jacobian_rank=2_everywhere;smooth_normal_surface")
    print(f"collision_normalization_rows={tuple(collision_rows)};fibre_g={fibre_polynomial}")
    print("canonical_form=-(g/C^2)dC^dg_on_D0_chart=-(1/(4Ag))dA^dB_on_Dminus_chart;divisor=D0")
    print("three_bracket=1={C,-B^2/4+A/6+A^2C^2/6}+(7/4){BC^2,ABC}-(13/8){AC^3,B^2}")
    print("degree3_four_bracket=1={C,A/6-B^2/4}+{C^2,A^2C}/12-(13/6){C^3,AB^2}+5{BC^2,ABC}")
    print("two_bracket_coefficients=108a^2+2040a-403=0;23324a+920b^2-4459=0;b!=0")
    print("two_bracket_support=P1,Q1,P2,Q2_target_degree_at_most_4;identity=verified_mod_coefficient_field")
    print("four_function_linear_compression=unique_w01=1,w23=1/(102b),cross=0;Pluecker=1/(102b)!=0")
    print("scope=bracket_width_at_most_2;one_bracket_Darboux_pair_open;homogeneous_all_degree_blocked_in_theorem")
    print("commands=python3 -B 04-computation/jacobian_y0_collision_normalization_bracket_anatomy_thm3686.py;python3 -B -O 04-computation/jacobian_y0_collision_normalization_bracket_anatomy_thm3686.py")


if __name__ == "__main__":
    main()
