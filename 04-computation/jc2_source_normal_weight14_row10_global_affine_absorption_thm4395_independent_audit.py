#!/usr/bin/env python3
"""Independent exact audit for THM-4395.

The script imports only the audited THM-4308/4315 row operators.  It does
not import the THM-4395 primary candidate.  Starting from the literal complete
residual-weight-at-most-14 source, it independently reconstructs the bracket
and projected P_2/P_3 depth fibres through row ten.
"""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_bracket_hasse_rows8_thm4308 as R8  # noqa: E402
import jc2_source_normal_student_stein_row9_thm4315 as R9  # noqa: E402


x, t = R8.x, R8.t
c51, c23, c70, c42, c14 = sp.symbols("c51 c23 c70 c42 c14")
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def exact(value: sp.Expr) -> sp.Expr:
    return sp.cancel(sp.together(sp.expand(value)))


def xcoeff(value: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(value).coeff(x, degree)


def tcoeff(value: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(value).coeff(t, degree)


def primitive(value: sp.Expr) -> sp.Expr:
    numerator, _ = sp.fraction(sp.together(exact(value)))
    variables = tuple(sorted(numerator.free_symbols, key=str))
    if not variables:
        return numerator
    return sp.Poly(numerator, *variables, domain=sp.QQ).primitive()[1].as_expr()


def proportional(left: sp.Expr, right: sp.Expr) -> bool:
    variables = tuple(sorted(left.free_symbols | right.free_symbols, key=str))
    lp = sp.Poly(sp.fraction(sp.together(exact(left)))[0], *variables, domain=sp.QQ)
    rp = sp.Poly(sp.fraction(sp.together(exact(right)))[0], *variables, domain=sp.QQ)
    return lp.monic() == rp.monic()


def apply(rows: list[sp.Expr], *maps: dict[sp.Symbol, sp.Expr]) -> list[sp.Expr]:
    answer = rows
    for mapping in maps:
        answer = [exact(row.subs(mapping)) for row in answer]
    return answer


def resolve(value: sp.Expr, mapping: dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    answer = exact(value)
    for _ in range(5):
        updated = exact(answer.subs(mapping, simultaneous=False))
        if updated == answer:
            return answer
        answer = updated
    return answer


def eliminate_linear(
    equations: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> tuple[
    dict[sp.Symbol, sp.Expr],
    list[sp.Expr],
    list[sp.Expr],
    sp.Matrix,
    tuple[int, ...],
]:
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    pivots = tuple(matrix.rref()[1])
    row_pivots = tuple(matrix.T.rref()[1])
    if pivots:
        square = matrix.extract(row_pivots, pivots)
        check(square.rows == square.cols and square.det() != 0, "linear pivot minor")
        values = square.inv() * rhs.extract(row_pivots, (0,))
        substitutions = {
            variables[column]: exact(values[index])
            for index, column in enumerate(pivots)
        }
    else:
        substitutions = {}
    raw = [exact(value.subs(substitutions)) for value in equations]
    conditions: list[sp.Expr] = []
    for value in raw:
        if value == 0:
            continue
        candidate = primitive(value)
        if all(
            exact(candidate - old) != 0 and exact(candidate + old) != 0
            for old in conditions
        ):
            conditions.append(candidate)
    return substitutions, conditions, raw, matrix, pivots


def depth_equations(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
) -> tuple[list[sp.Expr], sp.Matrix, sp.Matrix, int]:
    acoords, amatrix = R9.depth_matrix(2, row)
    ccoords, cmatrix = R9.depth_matrix(3, row)
    avec = sp.Matrix([xcoeff(arows[n], degree) for n, degree in acoords])
    cvec = sp.Matrix([xcoeff(crows[n], degree) for n, degree in ccoords])
    aeq = [sp.expand((left.T * avec)[0]) for left in amatrix.T.nullspace()]
    ceq = [sp.expand((left.T * cvec)[0]) for left in cmatrix.T.nullspace()]
    return aeq + ceq, amatrix, cmatrix, len(aeq)


def bracket_equations(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    source: sp.Expr,
) -> list[sp.Expr]:
    difference = exact(source - R8.predicted_G(row, arows, crows))
    check(
        difference == 0 or sp.Poly(difference, x, domain=sp.EX).degree() <= row,
        f"row-{row} coefficient universe",
    )
    equations = [xcoeff(difference, degree) for degree in range(row + 1)]
    check(exact(difference - sum(equations[j] * x**j for j in range(row + 1))) == 0,
          f"row-{row} coefficient exhaustion")
    return equations


def append_tangent(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    prefix: str,
) -> tuple[list[sp.Expr], list[sp.Expr], tuple[sp.Symbol, ...]]:
    abase, cbase = R8.particular_row(row, R8.B_row(row, arows, crows))
    symbols = tuple(sp.symbols(f"{prefix}_0:{row + 1}"))
    theta = sum(symbols[j] * x**j for j in range(row + 1))
    return (
        arows + [sp.expand(abase + theta * sp.diff(R8.A0, x))],
        crows + [sp.expand(cbase + theta * sp.diff(R8.C0, x))],
        symbols,
    )


def constant_denominator(value: sp.Expr) -> int:
    denominator = sp.fraction(sp.together(value))[1]
    check(not denominator.free_symbols and denominator != 0, "constant graph denominator")
    return int(denominator)


def main() -> None:
    face13 = c51 * R8.p**5 * R8.y + c23 * R8.p**2 * R8.y**3
    face14 = c70 * R8.p**7 + c42 * R8.p**4 * R8.y**2 + c14 * R8.p * R8.y**4
    source = sp.expand(R8.H + face13 + face14)
    grows = {row: tcoeff(source, row) for row in range(4, 15)}

    indices13 = [(i, j) for j in range(15) for i in range(15) if 2 * i + 3 * j == 13]
    indices14 = [(i, j) for j in range(15) for i in range(15) if 2 * i + 3 * j == 14]
    check(indices13 == [(5, 1), (2, 3)], "complete weight-13 face")
    check(indices14 == [(7, 0), (4, 2), (1, 4)], "complete weight-14 face")
    check(
        sp.Matrix([[1, 0], [6, 1]]).det() == 1
        and sp.Matrix([[1, 0, 0], [7, 1, 0], [21, 6, 1]]).det() == 1,
        "integral Pascal sidecars",
    )

    expected13 = {
        7: c51 * x,
        8: (6 * c51 + c23) * x**3,
        9: (15 * c51 + 5 * c23) * x**5,
        10: (20 * c51 + 10 * c23) * x**7,
    }
    expected14 = {
        7: c70,
        8: (7 * c70 + c42) * x**2,
        9: (21 * c70 + 6 * c42 + c14) * x**4,
        10: (35 * c70 + 15 * c42 + 5 * c14) * x**6,
    }
    check(
        all(exact(tcoeff(face13, row) - expected13[row]) == 0 for row in range(7, 11)),
        "literal weight-13 rows through ten",
    )
    check(
        all(exact(tcoeff(face14, row) - expected14[row]) == 0 for row in range(7, 11)),
        "literal weight-14 rows through ten",
    )

    # Independent bracket recursion through row eight.
    arows = [R8.A0, R8.A1, R8.A2, R8.A3]
    crows = [R8.C0, R8.C1, R8.C2, R8.C3]
    bracket_subs: dict[sp.Symbol, sp.Expr] = {}
    response_symbols = {5: R8.upsilon5, 6: R8.U, 7: R8.W, 8: R8.Z}
    for n in range(4, 8):
        abase, cbase = R8.particular_row(n, R8.B_row(n, arows, crows))
        m = n + 1
        difference = sp.expand(
            grows[m].subs(bracket_subs)
            - R8.predicted_G(m, arows + [abase], crows + [cbase])
        )
        moment = R9.primitive_student_row(m)
        obstruction = exact(sum(moment[j] * xcoeff(difference, j) for j in range(m + 1)))
        answers = sp.solve(obstruction, response_symbols[m])
        check(len(answers) == 1, f"row-{m} scalar response")
        bracket_subs[response_symbols[m]] = exact(answers[0])
        target = sp.expand(difference.subs(response_symbols[m], answers[0]))
        theta = R8.tangent_solve(m, target)
        arows.append(sp.expand(abase + theta * sp.diff(R8.A0, x)))
        crows.append(sp.expand(cbase + theta * sp.diff(R8.C0, x)))
        check(
            exact(R8.predicted_G(m, arows, crows) - grows[m].subs(bracket_subs)) == 0,
            f"row-{m} bracket reconstruction",
        )

    theta8 = tuple(sp.symbols("w4395_theta8_0:9"))
    tangent8 = sum(theta8[j] * x**j for j in range(9))
    abase8, cbase8 = R8.particular_row(8, R8.B_row(8, arows, crows))
    arows.append(sp.expand(abase8 + tangent8 * sp.diff(R8.A0, x)))
    crows.append(sp.expand(cbase8 + tangent8 * sp.diff(R8.C0, x)))

    depth8, ap8, cp8, split8 = depth_equations(arows, crows, 8)
    terminal8, conditions8, raw_depth8, matrix_depth8, pivots8 = eliminate_linear(depth8, theta8)
    check((ap8.shape, ap8.rank()) == ((63, 131), 51), "row-8 P2 universe")
    check((cp8.shape, cp8.rank()) == ((72, 204), 63), "row-8 P3 universe")
    check((matrix_depth8.shape, matrix_depth8.rank(), pivots8) == ((21, 9), 2, (7, 8)),
          "row-8 joint terminal")
    old_gate = (
        896 - 15 * R8.Delta,
        3 * R8.Phi + 2 * R8.zeta3,
        3030 * R8.Delta + 225 * R8.Theta - 182528,
    )
    check(len(conditions8) == 3 and all(any(proportional(a, b) for a in conditions8) for b in old_gate),
          "row-8 old Hasse ideal")
    gate = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * R8.Phi,
    }
    check(all(exact(value.subs(gate)) == 0 for value in raw_depth8), "row-8 depth coefficients vanish")

    # Row nine: select the affine-seven bracket fibre, then globally solve c14.
    a8 = apply(arows, gate, terminal8, gate)
    c8 = apply(crows, gate, terminal8, gate)
    g9 = exact(grows[9].subs(bracket_subs).subs(gate))
    select9, conditions9, raw9, matrix9, pivots9 = eliminate_linear(
        bracket_equations(a8, c8, 9, g9), theta8[:7]
    )
    check((matrix9.shape, matrix9.rank(), pivots9) == ((10, 7), 7, tuple(range(7))),
          "row-9 bracket fibre")
    check(len(conditions9) == 1, "row-9 unique condition")
    e9 = (
        613527750 * R8.Phi**2
        - 511211250 * R8.Phi * R8.alpha11
        - 3154140000 * R8.Phi * R8.eta
        - 255605625 * R8.eta**2
        + 6736896000 * R8.xi10
        - 46483785515008
    )
    face_linear = 3339765000 * c70 + 898128000 * c42 + 216513000 * c14
    row9_equation = sp.expand(e9 - face_linear)
    check(proportional(conditions9[0], row9_equation), "row-9 source equation")
    check(sp.diff(conditions9[0], c51) == sp.diff(conditions9[0], c23) == 0,
          "odd channels absent at row nine")
    c14_answers = sp.solve(row9_equation, c14)
    check(len(c14_answers) == 1, "global c14 graph")
    c14_graph = exact(c14_answers[0])
    check(constant_denominator(c14_graph) == 216513000, "c14 denominator")
    c14_subs = {c14: c14_graph}

    selected_a8 = apply(a8, select9, c14_subs)
    selected_c8 = apply(c8, select9, c14_subs)
    a9_trial, c9_trial, theta9 = append_tangent(
        selected_a8, selected_c8, 9, "w4395_theta9"
    )
    depth9, ap9, cp9, split9 = depth_equations(a9_trial, c9_trial, 9)
    terminal9, conditions_depth9, raw_depth9, matrix_depth9, pivots_depth9 = eliminate_linear(
        depth9, theta9
    )
    check((ap9.shape, ap9.rank()) == ((75, 160), 59), "row-9 P2 universe")
    check((cp9.shape, cp9.rank()) == ((85, 251), 73), "row-9 P3 universe")
    amat9, _ = sp.linear_eq_to_matrix(depth9[:split9], theta9)
    cmat9, _ = sp.linear_eq_to_matrix(depth9[split9:], theta9)
    check((amat9.rank(), amat9.rref()[1]) == (3, (7, 8, 9)), "row-9 P2 terminal")
    check((cmat9.rank(), cmat9.rref()[1]) == (2, (8, 9)), "row-9 P3 terminal")
    check((matrix_depth9.shape, matrix_depth9.rank(), pivots_depth9) == ((28, 10), 3, (7, 8, 9)),
          "row-9 joint terminal")
    check(not conditions_depth9 and all(value == 0 for value in raw_depth9), "row-9 depth automatic")

    # Row ten bracket: after c14, one condition globally solves xi10.
    a9 = apply(a9_trial, terminal9)
    c9 = apply(c9_trial, terminal9)
    g10 = exact(grows[10].subs(bracket_subs).subs(gate).subs(c14_subs))
    free9 = tuple(theta9[j] for j in range(len(theta9)) if j not in pivots_depth9)
    select10, conditions10, raw10, matrix10, pivots10 = eliminate_linear(
        bracket_equations(a9, c9, 10, g10), free9
    )
    check((matrix10.shape, matrix10.rank(), pivots10) == ((11, 7), 7, tuple(range(7))),
          "row-10 bracket fibre")
    check(len(conditions10) == 2, "row-10 has two source conditions")
    expected_xi_equation = (
        13365000 * R8.Phi**2
        + 15035625 * R8.Phi * R8.eta
        + 6014250 * c42
        + 50787000 * c70
        + 57672000 * R8.xi10
        - 964604821504
    )
    xi_conditions = [
        value
        for value in conditions10
        if R8.beta11 not in value.free_symbols and c51 not in value.free_symbols
    ]
    check(len(xi_conditions) == 1, "unique triangular xi10 equation")
    xi_condition = xi_conditions[0]
    check(proportional(xi_condition, expected_xi_equation), "exact xi10 equation")
    check(not sp.diff(xi_condition, R8.xi10).free_symbols,
          "xi10 has constant pivot coefficient")
    xi_answers = sp.solve(xi_condition, R8.xi10)
    check(len(xi_answers) == 1, "global xi10 graph")
    xi_graph = exact(xi_answers[0])
    xi_denominator = constant_denominator(xi_graph)
    check(xi_denominator == 57672000, "xi10 denominator")
    xi_subs = {R8.xi10: xi_graph}
    remaining_bracket = [exact(value.subs(xi_subs)) for value in conditions10 if value != xi_condition]
    remaining_bracket = [value for value in remaining_bracket if value != 0]
    check(len(remaining_bracket) == 1, "one bracket equation remains after xi10")
    expected_remaining_bracket = (
        -104916222000 * R8.Phi**2
        + 122625090000 * R8.Phi * R8.alpha11
        + 19707603750 * R8.Phi * R8.beta11
        + 20802470625 * R8.Phi * c51
        - 246422138625 * R8.Phi * R8.eta
        + 20802470625 * R8.alpha11 * R8.eta
        - 89131914000 * c42
        - 194981256000 * c70
        + 61312545000 * R8.eta**2
        + 2707389207937024
    )
    check(
        proportional(remaining_bracket[0], expected_remaining_bracket),
        "exact remaining bracket equation",
    )

    # Row-ten depth is independent of the remaining bracket equation.
    selected_a9 = apply(a9, select10, xi_subs)
    selected_c9 = apply(c9, select10, xi_subs)
    a10_trial, c10_trial, theta10 = append_tangent(
        selected_a9, selected_c9, 10, "w4395_theta10"
    )
    depth10, ap10, cp10, split10 = depth_equations(a10_trial, c10_trial, 10)
    terminal10, conditions_depth10, raw_depth10, matrix_depth10, pivots_depth10 = eliminate_linear(
        depth10, theta10
    )
    check((ap10.shape, ap10.rank()) == ((88, 193), 68), "row-10 P2 universe")
    check((cp10.shape, cp10.rank()) == ((99, 304), 83), "row-10 P3 universe")
    amat10, _ = sp.linear_eq_to_matrix(depth10[:split10], theta10)
    cmat10, _ = sp.linear_eq_to_matrix(depth10[split10:], theta10)
    check((amat10.rank(), amat10.rref()[1]) == (3, (8, 9, 10)), "row-10 P2 terminal")
    check((cmat10.rank(), cmat10.rref()[1]) == (3, (8, 9, 10)), "row-10 P3 terminal")
    check((matrix_depth10.shape, matrix_depth10.rank(), pivots_depth10) == ((36, 11), 3, (8, 9, 10)),
          "row-10 joint terminal")
    check(len(theta10) - len(pivots_depth10) == 8, "row-10 terminal affine-eight fibre")
    check(len(conditions_depth10) == 1, "unique row-10 depth condition")
    beta_equation = -91 * R8.Phi + 15 * R8.beta11 + 18 * R8.eta
    check(proportional(conditions_depth10[0], beta_equation), "row-10 beta depth equation")
    beta_answers = sp.solve(beta_equation, R8.beta11)
    check(len(beta_answers) == 1, "global beta11 graph")
    beta_graph = exact(beta_answers[0])
    check(constant_denominator(beta_graph) == 15, "beta11 denominator")
    beta_subs = {R8.beta11: beta_graph}

    # The second bracket equation becomes a constant-coefficient c42 graph.
    c42_equation = primitive(remaining_bracket[0].subs(beta_subs))
    c42_coefficient = sp.diff(c42_equation, c42)
    check(c42_coefficient != 0 and not c42_coefficient.free_symbols,
          "c42 has constant nonzero coefficient")
    c42_answers = sp.solve(c42_equation, c42)
    check(len(c42_answers) == 1, "global c42 graph")
    c42_graph = exact(c42_answers[0])
    c42_denominator = constant_denominator(c42_graph)
    check(c42_coefficient == -89131914000 and c42_denominator == 89131914000,
          "c42 coefficient and denominator")
    c42_subs = {c42: c42_graph}

    # Resolve the triangular graph because xi10 and c14 were solved before c42.
    xi_final = exact(xi_graph.subs(c42_subs))
    c14_final = exact(c14_graph.subs({R8.xi10: xi_final}).subs(c42_subs))
    final_graphs = {
        R8.xi10: xi_final,
        R8.beta11: beta_graph,
        c42: c42_graph,
        c14: c14_final,
    }
    xi_resolved_denominator = constant_denominator(xi_final)
    c14_resolved_denominator = constant_denominator(c14_final)
    check(
        (xi_resolved_denominator, c14_resolved_denominator)
        == (316913472000, 10829527551000),
        "exact resolved graph denominators",
    )
    source_parameters = {R8.Phi, R8.eta, R8.alpha11, c51, c23, c70}
    check(
        set().union(*(value.free_symbols for value in final_graphs.values()))
        <= source_parameters,
        "global affine-six source graph",
    )
    solved_symbols = set(final_graphs)
    check(all(not (value.free_symbols & solved_symbols) for value in final_graphs.values()),
          "triangular graph is resolved")
    check(
        all(not sp.fraction(sp.together(value))[1].free_symbols for value in final_graphs.values()),
        "resolved graph denominators remain constant",
    )

    # Necessity and sufficiency: every literal bracket/depth coefficient dies.
    check(all(resolve(value, final_graphs) == 0 for value in raw9),
          "all row-9 bracket coefficients vanish")
    check(all(resolve(value, final_graphs) == 0 for value in raw_depth9),
          "all row-9 depth coefficients vanish")
    check(all(resolve(value, final_graphs) == 0 for value in raw10),
          "all row-10 bracket coefficients vanish")
    check(all(resolve(value, final_graphs) == 0 for value in raw_depth10),
          "all row-10 depth coefficients vanish")

    # Branch-free Phi=0 control and a dense rational hostile to D_0=0.
    free_symbols = (R8.Phi, R8.eta, R8.alpha11, c51, c23, c70)

    def specialize_graphs(values: tuple[int, ...]) -> dict[sp.Symbol, sp.Expr]:
        base = dict(zip(free_symbols, map(sp.Integer, values)))
        specialized = {symbol: exact(value.subs(base)) for symbol, value in final_graphs.items()}
        return {**base, **specialized}

    boundary = specialize_graphs((0, 2, 3, 5, 7, 11))
    check(all(value.is_finite for value in boundary.values()), "Phi=0 graphs are finite")
    check(all(exact(value.subs(boundary)) == 0 for value in raw9 + raw_depth9 + raw10 + raw_depth10),
          "Phi=0 literal control")

    dense = specialize_graphs((1, 2, 3, 5, 7, 11))
    check(all(dense[symbol] != 0 for symbol in (c51, c23, c70, c42, c14)),
          "dense weight-13/14 coefficients nonzero")
    check(all(exact(value.subs(dense)) == 0 for value in raw9 + raw_depth9 + raw10 + raw_depth10),
          "dense literal row-ten control")
    old_elliptic = (
        7231154026500 * R8.Phi**3
        + 50541940696500 * R8.Phi**2 * R8.eta
        + 6793915500000 * R8.Phi * R8.eta**2
        - 631918028977864704 * R8.Phi
        + 353642000625 * R8.eta**3
        - 91584545734393856 * R8.eta
    )
    old_elliptic_value = exact(old_elliptic.subs(dense))
    check(old_elliptic_value != 0, "dense control is hostile to old elliptic carrier")

    print("THM-4395 weight-14 row-ten global affine absorption independent audit")
    print("imports=audited_THM4308_R8_and_THM4315_R9; THM4395_primary_import=no")
    print("universe=complete_fixed_residual_weight_at_most_14")
    print("row9=c14_global_graph; row9_depth=automatic")
    print("row10_bracket=two_conditions; xi10_global_graph=yes")
    print("row10_depth=-91*Phi+15*beta11+18*eta; beta11_global_graph=yes")
    print(f"row10_remaining_bracket=c42_global_graph coefficient={sp.sstr(c42_coefficient)}")
    print(
        "triangular_graph_denominators="
        f"c14:216513000,xi10:{xi_denominator},beta11:15,c42:{c42_denominator};"
        f"resolved_xi10:{xi_resolved_denominator},resolved_c14:{c14_resolved_denominator};"
        "parameter_localization=none"
    )
    print("row10_universe=P2(88x193,rank68),P3(99x304,rank83),joint(36x11,rank3,pivots8,9,10)")
    print("source_projection=A6(Phi,eta,alpha11,c51,c23,c70);terminal_fibre=A8")
    print("Phi_zero_control=PASS")
    print(
        "dense_control="
        + sp.sstr(tuple(dense[symbol] for symbol in (*free_symbols, R8.xi10, R8.beta11, c42, c14)))
        + f";old_elliptic_D={sp.sstr(old_elliptic_value)}_nonzero"
    )
    print("field=characteristic_zero finite_field_used=no")
    print("scope=fixed_source_normal_chart_through_row10_only;row11_and_entry_open")
    print(f"checks={CHECKS} result=PASS")


if __name__ == "__main__":
    main()
