#!/usr/bin/env python3
"""Primary exact certificate for THM-4399, inspired by arXiv:2608.23777.

First continue the complete fixed residual-weight-at-most-14 source-normal
family of THM-4395 through row eleven.  Quotienting by every row-ten tangent
response leaves one primitive source hypersurface F=0; projected P_2/P_3
depth is then automatic with an affine-eight terminal fibre.

Second imitate the paper's triangular Hamiltonian correction: adjoin one
source term lambda18*p**3*y**4.  It has t-valuation eleven and residual weight
eighteen, so it changes no earlier row.  Its even leading coefficient is a
nonzero constant direction in the row-eleven compatibility quotient and
globally absorbs F.  The complete valuation-eleven monomial table records the
odd/even boundary of this mechanism.  This is a finite-row session scout, not
a theorem companion or a complete weight-at-most-18 calculation.
"""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_weight14_row10_global_affine_absorption_thm4395 as W10  # noqa: E402


R9 = W10.R9
R8 = R9.R8
x = R9.x
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def exact(value: sp.Expr) -> sp.Expr:
    return R9.exact(value)


def primitive(value: sp.Expr) -> sp.Expr:
    return R9.primitive(value)


def proportional(left: sp.Expr, right: sp.Expr) -> bool:
    return R9.proportional(left, right)


def row_ten_state():
    a9, c9, theta9, bracket_subs, grows, gate, c14_graph = W10.build_row_nine_state()
    g10 = exact(grows[10].subs(bracket_subs).subs(gate).subs(R9.c14, c14_graph))
    select10, bracket10, matrix10, pivots10 = R9.solve_bracket_fibre(
        a9, c9, 10, g10, theta9[:7]
    )
    assert (matrix10.shape, matrix10.rank(), pivots10) == (
        (11, 7), 7, tuple(range(7))
    )

    xi_equation = (
        13365000 * R8.Phi**2
        + 15035625 * R8.Phi * R8.eta
        + 6014250 * R9.c42
        + 50787000 * R9.c70
        + 57672000 * R8.xi10
        - 964604821504
    )
    xi_condition = [value for value in bracket10 if proportional(value, xi_equation)]
    assert len(xi_condition) == 1
    remaining = [value for value in bracket10 if value is not xi_condition[0]]
    assert len(remaining) == 1
    xi_graph = exact(sp.solve(xi_equation, R8.xi10)[0])

    selected_a9 = R9.impose(a9, select10, {R8.xi10: xi_graph})
    selected_c9 = R9.impose(c9, select10, {R8.xi10: xi_graph})
    a10_trial, c10_trial, theta10 = R9.append_tangent(
        selected_a9, selected_c9, 10, "w14r11_theta10"
    )
    depth10, _, _ = R9.projected_depth_residuals(a10_trial, c10_trial, 10)
    terminal10, depth10_conditions, matrix_depth10, pivots_depth10 = (
        R9.eliminate_linear_fibre(depth10, theta10)
    )
    assert (matrix_depth10.shape, matrix_depth10.rank(), pivots_depth10) == (
        (36, 11), 3, (8, 9, 10)
    )
    beta_equation = -91 * R8.Phi + 15 * R8.beta11 + 18 * R8.eta
    assert len(depth10_conditions) == 1
    assert proportional(depth10_conditions[0], beta_equation)
    beta_graph = exact(sp.solve(beta_equation, R8.beta11)[0])

    bracket_remaining = (
        -104916222000 * R8.Phi**2
        + 122625090000 * R8.Phi * R8.alpha11
        + 19707603750 * R8.Phi * R8.beta11
        + 20802470625 * R8.Phi * R9.c51
        - 246422138625 * R8.Phi * R8.eta
        + 20802470625 * R8.alpha11 * R8.eta
        - 89131914000 * R9.c42
        - 194981256000 * R9.c70
        + 61312545000 * R8.eta**2
        + 2707389207937024
    )
    assert proportional(
        remaining[0].subs(R8.xi10, xi_graph), bracket_remaining
    )
    final_bracket = exact(bracket_remaining.subs(R8.beta11, beta_graph))
    assert sp.diff(final_bracket, R9.c42) == -89131914000
    c42_graph = exact(sp.solve(final_bracket, R9.c42)[0])
    xi_final = exact(xi_graph.subs(R9.c42, c42_graph))
    c14_final = exact(
        c14_graph.subs(R8.xi10, xi_final).subs(R9.c42, c42_graph)
    )
    graph10 = {
        R8.beta11: beta_graph,
        R9.c42: c42_graph,
        R8.xi10: xi_final,
        R9.c14: c14_final,
    }
    a10 = R9.impose(a10_trial, terminal10, graph10)
    c10 = R9.impose(c10_trial, terminal10, graph10)
    g11 = exact(grows[11].subs(bracket_subs).subs(gate).subs(graph10))
    free10 = tuple(
        symbol for index, symbol in enumerate(theta10) if index not in pivots_depth10
    )
    return a10, c10, free10, g11, grows, bracket_subs, gate, graph10


def main() -> None:
    a10, c10, free10, g11, grows, bracket_subs, gate, graph10 = row_ten_state()
    check(len(free10) == 8, "row-ten depth leaves eight tangent responses")

    # Complete weight-at-most-14 row-eleven bracket quotient.
    difference11 = exact(g11 - R8.predicted_G(11, a10, c10))
    equations11 = [R9.xcoeff(difference11, degree) for degree in range(12)]
    check(
        exact(difference11 - sum(equations11[j] * x**j for j in range(12))) == 0,
        "row-eleven coefficient universe is exhaustive",
    )
    matrix_direct, rhs_direct = sp.linear_eq_to_matrix(equations11, free10)
    left_null = matrix_direct.T.nullspace()
    check(
        (matrix_direct.shape, matrix_direct.rank(), len(left_null)) == ((12, 8), 8, 4),
        "row-eleven direct quotient dimensions",
    )
    select11, conditions11, matrix11, pivots11 = R9.solve_bracket_fibre(
        a10, c10, 11, g11, free10
    )
    check(
        (matrix11.shape, matrix11.rank(), pivots11, len(conditions11))
        == ((12, 8), 8, tuple(range(8)), 1),
        "row-eleven bracket has one source condition",
    )
    F = primitive(conditions11[0])
    source_variables = [
        R8.Phi,
        R8.eta,
        R8.alpha11,
        R9.c51,
        R9.c23,
        R9.c70,
    ]
    fpoly = sp.Poly(F, *source_variables, domain=sp.QQ)
    check(fpoly.content() == 1, "F is primitive over the integers")
    check(
        len(fpoly.terms()) == 36 and fpoly.total_degree() == 4,
        "F term count and total degree",
    )
    check(
        tuple(fpoly.degree(variable) for variable in source_variables)
        == (4, 4, 2, 2, 1, 2),
        "F multidegree",
    )
    check(sp.factor(F) == F, "F has no nontrivial rational factor")
    c23_response = sp.diff(F, R9.c23)
    check(
        c23_response == 21736146783278091456000000 * R8.Phi,
        "old c23 response vanishes exactly on Phi=0",
    )
    check(
        all(
            sp.diff(F, variable) == 0
            or bool(sp.diff(F, variable).free_symbols)
            for variable in source_variables
        ),
        "old source graph has no constant pivot",
    )

    # The selected residual ideal is literally (F), not just set-theoretically
    # supported on F=0.
    reduced11 = [exact(value.subs(select11)) for value in equations11]
    nonzero_reduced = [value for value in reduced11 if value != 0]
    reduced_ratios = [exact(value / F) for value in nonzero_reduced]
    check(
        nonzero_reduced
        and all(value.is_Rational and value != 0 for value in reduced_ratios),
        "every residual is a rational multiple of F",
    )
    check(
        any(value != 0 for value in reduced_ratios),
        "selected residual ideal contains a unit multiple of F over Q",
    )

    quotient0 = sp.Matrix([(left.T * rhs_direct)[0] for left in left_null])
    quotient_ratios = [exact(value / F) for value in quotient0]
    check(
        all(value.is_Rational for value in quotient_ratios),
        "four quotient coordinates are constant multiples of F",
    )
    active = [index for index, value in enumerate(quotient_ratios) if value != 0]
    check(active == [0], "only one primitive quotient coordinate is active")

    # Depth is automatic after row-eleven bracket selection.  This is over the
    # source polynomial ring: no localization or source graph is imposed.
    a10_selected = R9.impose(a10, select11)
    c10_selected = R9.impose(c10, select11)
    a11_trial, c11_trial, theta11 = R9.append_tangent(
        a10_selected, c10_selected, 11, "w14r11_theta11"
    )
    depth11, amatrix11, cmatrix11 = R9.projected_depth_residuals(
        a11_trial, c11_trial, 11
    )
    terminal11, depth11_conditions, depth11_matrix, depth11_pivots = (
        R9.eliminate_linear_fibre(depth11, theta11)
    )
    check(
        (amatrix11.shape, amatrix11.rank()) == ((102, 228), 77),
        "row-eleven P2 universe",
    )
    check(
        (cmatrix11.shape, cmatrix11.rank()) == ((114, 361), 94),
        "row-eleven P3 universe",
    )
    check(
        (depth11_matrix.shape, depth11_matrix.rank(), depth11_pivots)
        == ((45, 12), 4, (8, 9, 10, 11)),
        "row-eleven joint depth fibre",
    )
    check(
        len(depth11) == 45 and depth11_conditions == [],
        "row-eleven depth has no source condition",
    )
    check(
        all(exact(value.subs(terminal11)) == 0 for value in depth11),
        "all base row-eleven depth equations vanish",
    )
    check(
        len(theta11) - depth11_matrix.rank() == 8,
        "base row-eleven terminal fibre is affine eight",
    )

    # Full t-valuation-eleven source table.  These are all p^a*y^b with
    # a+2b=11.  Their first-row powers x^b expose the Student parity diagonal.
    valuation_monomials = [R8.p ** (11 - 2 * b) * R8.y**b for b in range(6)]
    weights = tuple(2 * (11 - 2 * b) + 3 * b for b in range(6))
    check(
        weights == (22, 21, 20, 19, 18, 17),
        "complete valuation-eleven residual weights",
    )
    for b, monomial in enumerate(valuation_monomials):
        check(
            all(exact(R9.tcoeff(monomial, row)) == 0 for row in range(11)),
            f"valuation-eleven monomial {b} is absent through row ten",
        )
        check(
            exact(R9.tcoeff(monomial, 11) - x**b) == 0,
            f"valuation-eleven monomial {b} leading row",
        )

    moment11 = R9.R9.primitive_student_row(11)
    student_ratios = tuple(sp.Rational(moment11[b], moment11[0]) for b in range(6))
    expected_ratios = (
        sp.Integer(1),
        sp.Integer(0),
        sp.Rational(2, 7),
        sp.Integer(0),
        sp.Rational(36, 133),
        sp.Integer(0),
    )
    check(
        student_ratios == expected_ratios,
        "row-eleven Student parity/transversality diagonal",
    )

    active_index = active[0]
    active_scale = quotient_ratios[active_index]
    response_constants = []
    for b in range(6):
        response_vector = sp.Matrix([-left[b] for left in left_null])
        check(
            all(
                response_vector[index] == 0
                for index in range(4)
                if index != active_index
            ),
            f"valuation-eleven monomial {b} creates no hidden quotient coordinate",
        )
        response_constants.append(exact(response_vector[active_index] / active_scale))
    response_constants = tuple(response_constants)
    expected_constants = (
        sp.Integer(25358837913824440032000000),
        sp.Integer(0),
        sp.Integer(7245382261092697152000000),
        sp.Integer(0),
        sp.Integer(6864046352614134144000000),
        sp.Integer(0),
    )
    check(
        response_constants == expected_constants,
        "full valuation-eleven quotient response table",
    )
    check(
        tuple(exact(value / response_constants[0]) for value in response_constants)
        == student_ratios,
        "quotient response equals the Student parity diagonal",
    )

    # The minimum-weight transverse late channel is p^3*y^4 (weight 18).
    # The lower-weight p*y^5 channel is odd and exactly invisible.
    lambda18, mu17 = sp.symbols("lambda18 mu17")
    correction_constant = response_constants[4]
    corrected_g11 = exact(g11 + lambda18 * x**4)
    select11_corrected, corrected_conditions, corrected_matrix, corrected_pivots = (
        R9.solve_bracket_fibre(a10, c10, 11, corrected_g11, free10)
    )
    check(
        (
            corrected_matrix.shape,
            corrected_matrix.rank(),
            corrected_pivots,
            len(corrected_conditions),
        )
        == ((12, 8), 8, tuple(range(8)), 1),
        "late correction row-eleven bracket dimensions",
    )
    corrected_F = primitive(corrected_conditions[0])
    check(
        corrected_F == F + correction_constant * lambda18,
        "late p3y4 correction has the exact constant pivot",
    )
    lambda_graph = exact(sp.solve(corrected_F, lambda18)[0])
    check(lambda_graph == -F / correction_constant, "global lambda18 graph")
    check(
        sp.fraction(sp.together(lambda_graph))[1] == correction_constant,
        "constant lambda18 denominator",
    )

    corrected_difference = exact(corrected_g11 - R8.predicted_G(11, a10, c10))
    corrected_equations = [
        R9.xcoeff(corrected_difference, degree) for degree in range(12)
    ]
    check(
        all(
            exact(value.subs(select11_corrected).subs(lambda18, lambda_graph)) == 0
            for value in corrected_equations
        ),
        "all twelve corrected bracket coefficients vanish",
    )

    corrected_a10 = R9.impose(a10, select11_corrected, {lambda18: lambda_graph})
    corrected_c10 = R9.impose(c10, select11_corrected, {lambda18: lambda_graph})
    corrected_a11, corrected_c11, corrected_theta11 = R9.append_tangent(
        corrected_a10, corrected_c10, 11, "w18late_theta11"
    )
    corrected_depth, corrected_amatrix, corrected_cmatrix = (
        R9.projected_depth_residuals(corrected_a11, corrected_c11, 11)
    )
    (
        corrected_terminal,
        corrected_depth_conditions,
        corrected_depth_matrix,
        corrected_depth_pivots,
    ) = R9.eliminate_linear_fibre(corrected_depth, corrected_theta11)
    check(
        (
            corrected_amatrix.shape,
            corrected_amatrix.rank(),
            corrected_cmatrix.shape,
            corrected_cmatrix.rank(),
        )
        == ((102, 228), 77, (114, 361), 94),
        "corrected row-eleven depth universes",
    )
    check(
        (
            corrected_depth_matrix.shape,
            corrected_depth_matrix.rank(),
            corrected_depth_pivots,
        )
        == ((45, 12), 4, (8, 9, 10, 11)),
        "corrected row-eleven joint depth fibre",
    )
    check(
        corrected_depth_conditions == [] and len(corrected_depth) == 45,
        "corrected row-eleven depth is automatic",
    )
    check(
        all(
            exact(value.subs(corrected_terminal).subs(lambda18, lambda_graph)) == 0
            for value in corrected_depth
        ),
        "all forty-five corrected depth equations vanish",
    )
    check(
        len(corrected_theta11) - corrected_depth_matrix.rank() == 8,
        "corrected terminal fibre is affine eight",
    )

    odd_g11 = exact(g11 + mu17 * x**5)
    _, odd_conditions, odd_matrix, _ = R9.solve_bracket_fibre(
        a10, c10, 11, odd_g11, free10
    )
    check(
        odd_matrix.rank() == 8 and len(odd_conditions) == 1,
        "odd weight-seventeen hostile quotient dimensions",
    )
    check(primitive(odd_conditions[0]) == F, "p*y^5 is exactly bracket-blind")

    # Hostiles: the old family misses row eleven at a dense point and at the
    # Phi=eta=0 corner; the constant correction graph works at both.
    dense = {
        R8.Phi: 1,
        R8.eta: 2,
        R8.alpha11: 3,
        R9.c51: 5,
        R9.c23: 7,
        R9.c70: 11,
    }
    corner = {
        R8.Phi: 0,
        R8.eta: 0,
        R8.alpha11: 1,
        R9.c51: 2,
        R9.c23: 3,
        R9.c70: 4,
    }
    for label, control in (("dense", dense), ("Phi_eta_zero", corner)):
        fvalue = exact(F.subs(control))
        lambda_value = exact(lambda_graph.subs(control))
        check(
            fvalue != 0 and lambda_value != 0,
            f"{label} hostile misses the uncorrected hypersurface",
        )
        corrected_control = {**control, lambda18: lambda_value}
        check(
            all(
                exact(value.subs(select11_corrected).subs(corrected_control)) == 0
                for value in corrected_equations
            ),
            f"{label} corrected bracket control",
        )
        check(
            all(
                exact(value.subs(corrected_terminal).subs(corrected_control)) == 0
                for value in corrected_depth
            ),
            f"{label} corrected depth control",
        )

    print("THM4399_source_normal_row11_late_weight18_response_primary")
    print("imports=THM4395_primary_state;field=characteristic_zero")
    print("universe_base=complete_fixed_residual_weight_at_most_14")
    print("row11_bracket=shape(12,8),rank8,pivots(0..7),primitive_conditions1")
    print("F_terms=36;F_total_degree=4;F_multidegree=(4,4,2,2,1,2);F_QQ_irreducible=yes")
    print("F=" + sp.sstr(F))
    print("old_c23_response=21736146783278091456000000*Phi;constant_source_pivot=none")
    print("selected_bracket_residual_ideal=(F)")
    print("base_row11_depth=P2(102x228,rank77),P3(114x361,rank94),joint(45x12,rank4,pivots8,9,10,11),conditions0,fibre=A8")
    print("valuation11_table_columns=b,monomial,weight,leading_x_power,Student_ratio,quotient_response")
    for b in range(6):
        print(
            "valuation11_row="
            + sp.sstr(
                (
                    b,
                    f"p^{11 - 2 * b}*y^{b}",
                    weights[b],
                    b,
                    student_ratios[b],
                    response_constants[b],
                )
            )
        )
    print("minimum_transverse_late_channel=lambda18*p^3*y^4;residual_weight=18;t_valuation=11")
    print("corrected_row11_equation=F+6864046352614134144000000*lambda18")
    print("lambda18_graph=-F/6864046352614134144000000;parameter_localization=none")
    print("corrected_literal_bracket_coefficients=12/12_zero")
    print("corrected_row11_depth=45/45_zero;terminal_fibre=A8")
    print("weight17_hostile=mu17*p*y^5;Student_ratio=0;row11_equation_remains_F")
    print(
        "dense_control_F="
        + sp.sstr(F.subs(dense))
        + ";lambda18="
        + sp.sstr(lambda_graph.subs(dense))
        + ";corrected_bracket_depth=PASS"
    )
    print(
        "Phi_eta_zero_control_F="
        + sp.sstr(F.subs(corner))
        + ";lambda18="
        + sp.sstr(lambda_graph.subs(corner))
        + ";corrected_bracket_depth=PASS"
    )
    print("scope=finite_row11_fixed_chart;late_single_channel_not_complete_weight18;entry_and_all_row_lift_open")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
