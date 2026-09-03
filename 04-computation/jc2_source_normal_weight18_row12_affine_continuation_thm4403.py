#!/usr/bin/env python3
"""Primary exact certificate for THM-4403.

This certificate imports the frozen row-eleven construction, rebuilds its
selected state, and then tests the last channel of the residual-weight-18 face. Under
the source-normal substitution, p^3*y^4 first enters row eleven and y^6 first
enters row twelve. It determines the complete row-twelve bracket quotient
after the row-eleven correction, solves the constant-pivot source graph, and
tests the projected P_2/P_3 depth equations. It does not declare the complete
weight-at-most-18 family.
"""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_row11_late_weight18_response_thm4399 as Q11  # noqa: E402


R9 = Q11.R9
R8 = Q11.R8
x = R9.x


def exact(value: sp.Expr) -> sp.Expr:
    return R9.exact(value)


def primitive(value: sp.Expr) -> sp.Expr:
    return R9.primitive(value)


def build_row_eleven_state():
    (
        a10,
        c10,
        free10,
        g11,
        grows,
        bracket_subs,
        gate,
        graph10,
    ) = Q11.row_ten_state()

    _, base_conditions, _, _ = R9.solve_bracket_fibre(a10, c10, 11, g11, free10)
    if len(base_conditions) != 1:
        raise AssertionError("row-eleven base obstruction count changed")
    F = primitive(base_conditions[0])

    lambda18 = sp.symbols("lambda18")
    g11_late = exact(g11 + lambda18 * x**4)
    select11, conditions11, matrix11, pivots11 = R9.solve_bracket_fibre(
        a10, c10, 11, g11_late, free10
    )
    if (matrix11.shape, matrix11.rank(), pivots11, len(conditions11)) != (
        (12, 8),
        8,
        tuple(range(8)),
        1,
    ):
        raise AssertionError("row-eleven late bracket geometry changed")
    corrected_F = primitive(conditions11[0])
    correction_constant = sp.diff(corrected_F, lambda18)
    if correction_constant == 0 or correction_constant.free_symbols:
        raise AssertionError("row-eleven correction is not a constant pivot")
    lambda_graph = exact(sp.solve(corrected_F, lambda18)[0])

    a10_selected = R9.impose(a10, select11, {lambda18: lambda_graph})
    c10_selected = R9.impose(c10, select11, {lambda18: lambda_graph})
    a11_trial, c11_trial, theta11 = R9.append_tangent(
        a10_selected, c10_selected, 11, "w18r12_theta11"
    )
    depth11, _, _ = R9.projected_depth_residuals(a11_trial, c11_trial, 11)
    terminal11, conditions_depth11, matrix_depth11, pivots_depth11 = (
        R9.eliminate_linear_fibre(depth11, theta11)
    )
    if conditions_depth11:
        raise AssertionError("row-eleven late depth acquired a source condition")
    if (matrix_depth11.shape, matrix_depth11.rank(), pivots_depth11) != (
        (45, 12),
        4,
        (8, 9, 10, 11),
    ):
        raise AssertionError("row-eleven late depth geometry changed")
    a11 = R9.impose(a11_trial, terminal11, {lambda18: lambda_graph})
    c11 = R9.impose(c11_trial, terminal11, {lambda18: lambda_graph})
    free11 = tuple(
        symbol for index, symbol in enumerate(theta11) if index not in pivots_depth11
    )
    return (
        a11,
        c11,
        free11,
        F,
        lambda18,
        lambda_graph,
        grows,
        bracket_subs,
        gate,
        graph10,
    )


def main() -> None:
    (
        a11,
        c11,
        free11,
        F,
        lambda18,
        lambda_graph,
        grows,
        bracket_subs,
        gate,
        graph10,
    ) = build_row_eleven_state()
    if len(free11) != 8:
        raise AssertionError("row-eleven terminal fibre is not affine eight")

    nu18 = sp.symbols("nu18")
    p3y4_row12 = exact(R9.tcoeff(R8.p**3 * R8.y**4, 12))
    y6_row12 = exact(R9.tcoeff(R8.y**6, 12))
    if p3y4_row12 != 7 * x**6 or y6_row12 != x**6:
        raise AssertionError("weight-18 row-twelve source flag changed")
    if any(exact(R9.tcoeff(R8.y**6, row)) != 0 for row in range(12)):
        raise AssertionError("y^6 enters before row twelve")

    g12_base = exact(
        grows[12].subs(bracket_subs).subs(gate).subs(graph10)
    )
    g12_late = exact(
        g12_base
        + lambda_graph * p3y4_row12
        + nu18 * y6_row12
    )
    select12, conditions12, matrix12, pivots12 = R9.solve_bracket_fibre(
        a11, c11, 12, g12_late, free11
    )

    print("THM4403_weight18_row12_affine_continuation_primary")
    print(f"row12_bracket_matrix={matrix12.shape};rank={matrix12.rank()};pivots={pivots12}")
    print(f"row12_condition_count={len(conditions12)}")
    for index, condition in enumerate(conditions12):
        value = primitive(condition)
        print(f"condition_{index}_terms={len(sp.Poly(value, *sorted(value.free_symbols, key=str), domain=sp.QQ).terms())}")
        print(f"condition_{index}_dnu18={sp.sstr(exact(sp.diff(value, nu18)))}")
        print(f"condition_{index}={sp.sstr(value)}")

    # Audit the entire exact-valuation-12 diagonal.  These seven monomials
    # all leave rows <=11 unchanged, but their x^b leading rows need not be
    # equivalent after the four depth-selected tangent directions have been
    # spent.  This catches responses invisible to the scalar Student moment.
    rho = sp.symbols("rho12_0:7")
    valuation12 = tuple(R8.p ** (12 - 2 * b) * R8.y**b for b in range(7))
    weights12 = tuple(24 - b for b in range(7))
    leading12 = tuple(exact(R9.tcoeff(monomial, 12)) for monomial in valuation12)
    if leading12 != tuple(x**b for b in range(7)):
        raise AssertionError("valuation-twelve leading diagonal changed")
    if any(
        exact(R9.tcoeff(monomial, row)) != 0
        for monomial in valuation12
        for row in range(12)
    ):
        raise AssertionError("valuation-twelve monomial enters an earlier row")
    g12_all = exact(
        g12_base
        + lambda_graph * p3y4_row12
        + sum(rho[b] * leading12[b] for b in range(7))
    )
    _, conditions_all, matrix_all, pivots_all = R9.solve_bracket_fibre(
        a11, c11, 12, g12_all, free11
    )
    if (matrix_all, pivots_all) != (matrix12, pivots12):
        raise AssertionError("valuation-twelve response changed tangent matrix")
    response = sp.Matrix(
        [
            [exact(sp.diff(condition, rho[b])) for b in range(7)]
            for condition in conditions_all
        ]
    )
    print(f"valuation12_weights={weights12}")
    print(f"valuation12_condition_count={len(conditions_all)}")
    print(f"valuation12_response_rank={response.rank()}")
    print(f"valuation12_response_matrix={sp.sstr(response)}")
    for b in range(7):
        print(
            f"valuation12_channel=b{b},monomial=p^{12 - 2*b}*y^{b},"
            f"weight={weights12[b]},response={sp.sstr(tuple(response[:, b]))}"
        )
    transverse_pairs = []
    for left in range(7):
        for right in range(left + 1, 7):
            minor = response[:, (left, right)]
            if minor.rows == 2 and minor.det() != 0:
                transverse_pairs.append((max(weights12[left], weights12[right]), left, right))
    if transverse_pairs:
        best = min(transverse_pairs)
        print(f"minimum_two_channel_pair=max_weight{best[0]},b{best[1]},b{best[2]}")

    if len(conditions12) == 1 and sp.diff(conditions12[0], nu18) != 0:
        nu_graph = exact(sp.solve(conditions12[0], nu18)[0])
        print(f"nu18_graph={sp.sstr(nu_graph)}")
        a11_selected = R9.impose(a11, select12, {nu18: nu_graph})
        c11_selected = R9.impose(c11, select12, {nu18: nu_graph})
        a12_trial, c12_trial, theta12 = R9.append_tangent(
            a11_selected, c11_selected, 12, "w18r12_theta12"
        )
        depth12, amatrix12, cmatrix12 = R9.projected_depth_residuals(
            a12_trial, c12_trial, 12
        )
        terminal12, depth_conditions12, depth_matrix12, depth_pivots12 = (
            R9.eliminate_linear_fibre(depth12, theta12)
        )
        print(
            "row12_depth_universe="
            f"P2{amatrix12.shape}/rank{amatrix12.rank()},"
            f"P3{cmatrix12.shape}/rank{cmatrix12.rank()}"
        )
        print(
            "row12_depth_fibre="
            f"shape{depth_matrix12.shape},rank{depth_matrix12.rank()},"
            f"pivots{depth_pivots12},conditions{len(depth_conditions12)},"
            f"dimension{len(theta12)-depth_matrix12.rank()}"
        )
        for index, condition in enumerate(depth_conditions12):
            print(f"depth_condition_{index}={sp.sstr(primitive(condition))}")

        if len(depth_conditions12) != 1:
            raise AssertionError("row-twelve depth condition count changed")
        depth_condition = primitive(depth_conditions12[0])
        c70_derivative = exact(sp.diff(depth_condition, R9.c70))
        if c70_derivative == 0 or c70_derivative.free_symbols:
            raise AssertionError("row-twelve depth has no constant c70 pivot")
        c70_graph = exact(sp.solve(depth_condition, R9.c70)[0])
        print(f"c70_graph={sp.sstr(c70_graph)}")
        bracket_residual = exact(g12_late - R8.predicted_G(12, a11, c11))
        bracket_equations = [R9.xcoeff(bracket_residual, degree) for degree in range(13)]
        print(
            "row12_bracket_literal_zero="
            + str(
                all(
                    exact(value.subs(select12).subs(nu18, nu_graph)) == 0
                    for value in bracket_equations
                )
            )
        )
        print(
            "row12_depth_literal_zero_after_terminal="
            + str(
                all(
                    exact(value.subs(terminal12).subs(nu18, nu_graph)) == 0
                    for value in depth12
                )
            )
        )

    if len(conditions12) == 2:
        primitive_conditions = [primitive(value) for value in conditions12]
        pivot_candidates = [
            value for value in primitive_conditions if sp.diff(value, nu18) != 0
        ]
        side_candidates = [
            value for value in primitive_conditions if sp.diff(value, nu18) == 0
        ]
        if len(pivot_candidates) != 1 or len(side_candidates) != 1:
            raise AssertionError("row-twelve condition triangularity changed")
        pivot_condition = pivot_candidates[0]
        side_condition = side_candidates[0]
        side_derivative = exact(sp.diff(side_condition, R9.c23))
        if side_derivative == 0 or side_derivative.free_symbols:
            raise AssertionError("persistent side condition has no constant c23 pivot")
        c23_graph = exact(sp.solve(side_condition, R9.c23)[0])
        pivot_on_side = exact(pivot_condition.subs(R9.c23, c23_graph))
        nu_graph = exact(sp.solve(pivot_on_side, nu18)[0])
        print(f"persistent_side_condition={sp.sstr(side_condition)}")
        print(f"c23_graph={sp.sstr(c23_graph)}")
        print(f"nu18_graph={sp.sstr(nu_graph)}")

        a11_selected = R9.impose(
            a11, select12, {R9.c23: c23_graph}, {nu18: nu_graph}
        )
        c11_selected = R9.impose(
            c11, select12, {R9.c23: c23_graph}, {nu18: nu_graph}
        )
        a12_trial, c12_trial, theta12 = R9.append_tangent(
            a11_selected, c11_selected, 12, "w18r12_theta12"
        )
        depth12, amatrix12, cmatrix12 = R9.projected_depth_residuals(
            a12_trial, c12_trial, 12
        )
        terminal12, depth_conditions12, depth_matrix12, depth_pivots12 = (
            R9.eliminate_linear_fibre(depth12, theta12)
        )
        print(
            "row12_depth_universe="
            f"P2{amatrix12.shape}/rank{amatrix12.rank()},"
            f"P3{cmatrix12.shape}/rank{cmatrix12.rank()}"
        )
        print(
            "row12_depth_fibre="
            f"shape{depth_matrix12.shape},rank{depth_matrix12.rank()},"
            f"pivots{depth_pivots12},conditions{len(depth_conditions12)},"
            f"dimension{len(theta12)-depth_matrix12.rank()}"
        )
        for index, condition in enumerate(depth_conditions12):
            print(f"depth_condition_{index}={sp.sstr(primitive(condition))}")

        if len(depth_conditions12) != 1:
            raise AssertionError("row-twelve depth condition count changed")
        depth_condition = primitive(depth_conditions12[0])
        c70_derivative = exact(sp.diff(depth_condition, R9.c70))
        if c70_derivative == 0 or c70_derivative.free_symbols:
            raise AssertionError("row-twelve depth has no constant c70 pivot")
        c70_graph = exact(sp.solve(depth_condition, R9.c70)[0])
        print(f"c70_graph={sp.sstr(c70_graph)}")

        bracket_residual = exact(g12_late - R8.predicted_G(12, a11, c11))
        bracket_equations = [R9.xcoeff(bracket_residual, degree) for degree in range(13)]
        bracket_graphs = ({R9.c23: c23_graph}, {nu18: nu_graph})
        print(
            "row12_bracket_literal_zero="
            + str(
                all(
                    exact(R9.impose([value.subs(select12)], *bracket_graphs)[0]) == 0
                    for value in bracket_equations
                )
            )
        )
        depth_zero = all(
            exact(value.subs(terminal12).subs(R9.c70, c70_graph)) == 0
            for value in depth12
        )
        print("row12_depth_literal_zero_after_terminal_and_c70=" + str(depth_zero))
        if not depth_zero:
            raise AssertionError("row-twelve depth graph does not solve every equation")

        final_lambda_graph = exact(
            lambda_graph.subs(R9.c23, c23_graph).subs(R9.c70, c70_graph)
        )
        final_nu_graph = exact(nu_graph.subs(R9.c70, c70_graph))
        graph_expressions = (
            c23_graph,
            c70_graph,
            final_lambda_graph,
            final_nu_graph,
        )
        if any(sp.fraction(sp.together(value))[1].free_symbols for value in graph_expressions):
            raise AssertionError("one row-twelve source graph requires localization")
        print("source_projection=global_A4_over(Phi,eta,alpha11,c51)")
        print("solved_source_coordinates=(c23,c70,lambda18,nu18)")
        print("all_source_graph_denominators=constant")

        dense = {R8.Phi: 1, R8.eta: 2, R8.alpha11: 3, R9.c51: 5}
        corner = {R8.Phi: 0, R8.eta: 0, R8.alpha11: 1, R9.c51: 2}
        for label, control in (("dense", dense), ("Phi_eta_zero", corner)):
            control_graph = {
                R9.c23: exact(c23_graph.subs(control)),
                R9.c70: exact(c70_graph.subs(control)),
                lambda18: exact(final_lambda_graph.subs(control)),
                nu18: exact(final_nu_graph.subs(control)),
            }
            if any(value == 0 for value in control_graph.values()):
                raise AssertionError(f"{label} control lost a solved coordinate")
            substitutions = {**control, **control_graph}
            if exact(
                F.subs(substitutions)
                + sp.Integer(6864046352614134144000000)
                * substitutions[lambda18]
            ) != 0:
                raise AssertionError(f"{label} row-eleven correction failed")
            if not all(
                exact(
                    R9.impose(
                        [value.subs(select12)],
                        {R9.c23: c23_graph},
                        {R9.c70: c70_graph},
                        {nu18: final_nu_graph},
                    )[0].subs(control)
                )
                == 0
                for value in bracket_equations
            ):
                raise AssertionError(f"{label} row-twelve bracket failed")
            if not all(
                exact(value.subs(terminal12).subs(R9.c70, c70_graph).subs(control))
                == 0
                for value in depth12
            ):
                raise AssertionError(f"{label} row-twelve depth failed")
            print(
                f"{label}_control="
                f"c23:{sp.sstr(control_graph[R9.c23])},"
                f"c70:{sp.sstr(control_graph[R9.c70])},"
                f"lambda18:{sp.sstr(control_graph[lambda18])},"
                f"nu18:{sp.sstr(control_graph[nu18])};PASS"
            )
    print("scope=two_channel_weight18_tail;not_complete_weight18;not_JC2")
    print("result=PASS")


if __name__ == "__main__":
    main()
