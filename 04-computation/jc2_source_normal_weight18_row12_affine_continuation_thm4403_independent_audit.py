#!/usr/bin/env python3
"""Clean-room exact audit for THM-4403.

Only the audited THM-4308/4315 operator implementations are imported.  The
script reconstructs the complete residual-weight-at-most-14 source graph
through row ten, independently reconstructs the late p^3*y^4 row-eleven
response, and then adjoins y^6 at row twelve.  It deliberately does not
import the companion row-12 diagnostic or either row-eleven candidate.
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
lambda18, nu18 = sp.symbols("lambda18 nu18")
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
    for _ in range(12):
        updated = exact(answer.subs(mapping, simultaneous=False))
        if updated == answer:
            return answer
        answer = updated
    raise AssertionError("substitution graph did not resolve")


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
        if all(not proportional(candidate, old) for old in conditions):
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
    check(
        exact(difference - sum(equations[j] * x**j for j in range(row + 1))) == 0,
        f"row-{row} coefficient exhaustion",
    )
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
    # The complete weight <=14 source, plus precisely the two late weight-18
    # channels under audit.  Their exact t-valuations are 11 and 12.
    face13 = c51 * R8.p**5 * R8.y + c23 * R8.p**2 * R8.y**3
    face14 = c70 * R8.p**7 + c42 * R8.p**4 * R8.y**2 + c14 * R8.p * R8.y**4
    late11 = lambda18 * R8.p**3 * R8.y**4
    late12 = nu18 * R8.y**6
    source = sp.expand(R8.H + face13 + face14 + late11 + late12)
    grows = {row: tcoeff(source, row) for row in range(4, 13)}

    check(tcoeff(late11, 10) == 0 and tcoeff(late11, 11) == lambda18 * x**4,
          "p3y4 exact valuation eleven")
    check(tcoeff(late12, 11) == 0 and tcoeff(late12, 12) == nu18 * x**6,
          "y6 exact valuation twelve")
    check(all(lambda18 not in grows[row].free_symbols for row in range(4, 11)),
          "late row-eleven firewall")
    check(all(nu18 not in grows[row].free_symbols for row in range(4, 12)),
          "late row-twelve firewall")

    # Reconstruct rows four through eight directly from the audited bracket
    # operator, without importing any weight-14 or late-response certificate.
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
        check(exact(R8.predicted_G(m, arows, crows) - grows[m].subs(bracket_subs)) == 0,
              f"row-{m} bracket reconstruction")

    theta8 = tuple(sp.symbols("audit12_theta8_0:9"))
    tangent8 = sum(theta8[j] * x**j for j in range(9))
    abase8, cbase8 = R8.particular_row(8, R8.B_row(8, arows, crows))
    arows.append(sp.expand(abase8 + tangent8 * sp.diff(R8.A0, x)))
    crows.append(sp.expand(cbase8 + tangent8 * sp.diff(R8.C0, x)))
    depth8, ap8, cp8, _ = depth_equations(arows, crows, 8)
    terminal8, conditions8, raw_depth8, matrix_depth8, pivots8 = eliminate_linear(depth8, theta8)
    check((ap8.shape, ap8.rank(), cp8.shape, cp8.rank())
          == ((63, 131), 51, (72, 204), 63), "row-eight depth universes")
    check((matrix_depth8.shape, matrix_depth8.rank(), pivots8)
          == ((21, 9), 2, (7, 8)), "row-eight terminal selection")
    gate = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * R8.Phi,
    }
    old_gate = (896 - 15 * R8.Delta, 3 * R8.Phi + 2 * R8.zeta3,
                3030 * R8.Delta + 225 * R8.Theta - 182528)
    check(len(conditions8) == 3
          and all(any(proportional(a, b) for a in conditions8) for b in old_gate),
          "row-eight Hasse ideal")
    check(all(exact(value.subs(gate)) == 0 for value in raw_depth8),
          "row-eight depth sufficiency")

    # Row nine: select the old affine-seven tangent and solve the complete
    # even weight-14 response c14.  Then its depth equations are automatic.
    a8 = apply(arows, gate, terminal8, gate)
    c8 = apply(crows, gate, terminal8, gate)
    g9 = resolve(grows[9], {**bracket_subs, **gate})
    select9, conditions9, raw9, matrix9, pivots9 = eliminate_linear(
        bracket_equations(a8, c8, 9, g9), theta8[:7]
    )
    check((matrix9.shape, matrix9.rank(), pivots9)
          == ((10, 7), 7, tuple(range(7))), "row-nine bracket selector")
    check(len(conditions9) == 1, "row-nine unique source condition")
    e9 = (
        613527750 * R8.Phi**2 - 511211250 * R8.Phi * R8.alpha11
        - 3154140000 * R8.Phi * R8.eta - 255605625 * R8.eta**2
        + 6736896000 * R8.xi10 - 46483785515008
    )
    row9_expected = e9 - (3339765000 * c70 + 898128000 * c42 + 216513000 * c14)
    check(proportional(conditions9[0], row9_expected), "row-nine source equation")
    c14_graph = exact(sp.solve(row9_expected, c14)[0])
    check(constant_denominator(c14_graph) == 216513000, "c14 graph denominator")
    c14_subs = {c14: c14_graph}
    selected_a8 = apply(a8, select9, c14_subs)
    selected_c8 = apply(c8, select9, c14_subs)
    a9_trial, c9_trial, theta9 = append_tangent(selected_a8, selected_c8, 9,
                                                "audit12_theta9")
    depth9, ap9, cp9, _ = depth_equations(a9_trial, c9_trial, 9)
    terminal9, depth_conditions9, raw_depth9, matrix_depth9, pivots_depth9 = eliminate_linear(
        depth9, theta9
    )
    check((ap9.shape, ap9.rank(), cp9.shape, cp9.rank())
          == ((75, 160), 59, (85, 251), 73), "row-nine depth universes")
    check((matrix_depth9.shape, matrix_depth9.rank(), pivots_depth9)
          == ((28, 10), 3, (7, 8, 9)), "row-nine depth selector")
    check(not depth_conditions9 and all(value == 0 for value in raw_depth9),
          "row-nine depth automatic")

    # Row ten independently reconstructs THM-4395's four global graphs.
    a9 = apply(a9_trial, terminal9)
    c9 = apply(c9_trial, terminal9)
    g10 = resolve(grows[10], {**bracket_subs, **gate, **c14_subs})
    free9 = tuple(theta9[j] for j in range(10) if j not in pivots_depth9)
    select10, conditions10, raw10, matrix10, pivots10 = eliminate_linear(
        bracket_equations(a9, c9, 10, g10), free9
    )
    check((matrix10.shape, matrix10.rank(), pivots10)
          == ((11, 7), 7, tuple(range(7))), "row-ten bracket selector")
    xi_expected = (
        13365000 * R8.Phi**2 + 15035625 * R8.Phi * R8.eta
        + 6014250 * c42 + 50787000 * c70 + 57672000 * R8.xi10
        - 964604821504
    )
    xi_condition = next(value for value in conditions10 if proportional(value, xi_expected))
    xi_graph = exact(sp.solve(xi_condition, R8.xi10)[0])
    check(constant_denominator(xi_graph) == 57672000, "xi10 graph denominator")
    xi_subs = {R8.xi10: xi_graph}
    remaining10 = [exact(value.subs(xi_subs)) for value in conditions10
                   if value != xi_condition]
    remaining10 = [value for value in remaining10 if value != 0]
    check(len(remaining10) == 1, "row-ten second bracket condition")
    selected_a9 = apply(a9, select10, xi_subs)
    selected_c9 = apply(c9, select10, xi_subs)
    a10_trial, c10_trial, theta10 = append_tangent(selected_a9, selected_c9, 10,
                                                   "audit12_theta10")
    depth10, ap10, cp10, _ = depth_equations(a10_trial, c10_trial, 10)
    terminal10, depth_conditions10, raw_depth10, matrix_depth10, pivots_depth10 = eliminate_linear(
        depth10, theta10
    )
    check((ap10.shape, ap10.rank(), cp10.shape, cp10.rank())
          == ((88, 193), 68, (99, 304), 83), "row-ten depth universes")
    check((matrix_depth10.shape, matrix_depth10.rank(), pivots_depth10)
          == ((36, 11), 3, (8, 9, 10)), "row-ten depth selector")
    beta_expected = -91 * R8.Phi + 15 * R8.beta11 + 18 * R8.eta
    check(len(depth_conditions10) == 1
          and proportional(depth_conditions10[0], beta_expected),
          "row-ten beta equation")
    beta_graph = exact(sp.solve(beta_expected, R8.beta11)[0])
    beta_subs = {R8.beta11: beta_graph}
    c42_equation = primitive(remaining10[0].subs(beta_subs))
    check(sp.diff(c42_equation, c42) == -89131914000,
          "row-ten c42 constant pivot")
    c42_graph = exact(sp.solve(c42_equation, c42)[0])
    c42_subs = {c42: c42_graph}
    xi_final = exact(xi_graph.subs(c42_subs))
    c14_final = exact(c14_graph.subs({R8.xi10: xi_final}).subs(c42_subs))
    row10_graphs = {
        R8.beta11: beta_graph,
        c42: c42_graph,
        R8.xi10: xi_final,
        c14: c14_final,
    }
    check(set().union(*(value.free_symbols for value in row10_graphs.values()))
          <= {R8.Phi, R8.eta, R8.alpha11, c51, c23, c70},
          "reconstructed affine-six row-ten graph")
    check(all(resolve(value, row10_graphs) == 0
              for value in raw9 + raw_depth9 + raw10 + raw_depth10),
          "literal rows nine and ten vanish")

    # Row eleven.  The late p^3*y^4 channel should be transverse to the one
    # primitive source obstruction left after the eight tangent pivots.
    inherited_graphs = {**bracket_subs, **gate, **row10_graphs}
    a10 = apply(a10_trial, terminal10, row10_graphs)
    c10 = apply(c10_trial, terminal10, row10_graphs)
    g11 = resolve(grows[11], inherited_graphs)
    free10 = tuple(theta10[j] for j in range(11) if j not in pivots_depth10)
    equations11 = bracket_equations(a10, c10, 11, g11)
    select11, conditions11, raw11, matrix11, pivots11 = eliminate_linear(
        equations11, free10
    )
    check((matrix11.shape, matrix11.rank(), pivots11)
          == ((12, 8), 8, tuple(range(8))), "row-eleven bracket selector")
    check(len(conditions11) == 1, "row-eleven source obstruction is principal")
    full11 = next(value for value in conditions11 if lambda18 in value.free_symbols)
    base11 = exact(full11.subs(lambda18, 0))
    derivative_c23 = exact(sp.diff(base11, c23))
    expected_dF = sp.Integer(21736146783278091456000000) * R8.Phi
    scale11 = exact(expected_dF / derivative_c23)
    check(not scale11.free_symbols and scale11 != 0,
          "row-eleven normalization is constant")
    equation11 = sp.expand(scale11 * full11)
    F = sp.expand(equation11.subs(lambda18, 0))
    response11 = sp.Integer(6864046352614134144000000)
    check(sp.diff(F, c23) == expected_dF, "row-eleven F c23 derivative")
    check(sp.diff(equation11, lambda18) == response11,
          "p3y4 row-eleven response coefficient")
    lambda_graph = exact(-F / response11)
    lambda_subs = {lambda18: lambda_graph}
    check(all(resolve(value, lambda_subs) == 0 for value in raw11),
          "late row-eleven response kills full bracket residual")
    check(sp.Poly(F, R8.Phi, R8.eta, R8.alpha11, c51, c23, c70).total_degree() == 4,
          "row-eleven primitive quartic degree")
    check(sp.Poly(F, c23).degree() == 1 and exact(F.subs(R8.Phi, 0)) != 0,
          "row-eleven quartic is primitive linear in c23")

    selected_a10 = apply(a10, select11, lambda_subs)
    selected_c10 = apply(c10, select11, lambda_subs)
    a11_trial, c11_trial, theta11 = append_tangent(selected_a10, selected_c10, 11,
                                                   "audit12_theta11")
    depth11, ap11, cp11, _ = depth_equations(a11_trial, c11_trial, 11)
    terminal11, depth_conditions11, raw_depth11, matrix_depth11, pivots_depth11 = eliminate_linear(
        depth11, theta11
    )
    check((ap11.shape, ap11.rank(), cp11.shape, cp11.rank())
          == ((102, 228), 77, (114, 361), 94), "row-eleven depth universes")
    check((matrix_depth11.shape, matrix_depth11.rank(), pivots_depth11)
          == ((45, 12), 4, (8, 9, 10, 11)), "row-eleven depth selector")
    check(not depth_conditions11 and all(value == 0 for value in raw_depth11),
          "row-eleven depth automatic")

    # Row twelve.  The eight surviving row-eleven tangents leave two source
    # conditions: a nu18 graph and the independent linear equation L.
    a11 = apply(a11_trial, terminal11)
    c11 = apply(c11_trial, terminal11)
    g12 = resolve(grows[12], {**inherited_graphs, **lambda_subs})
    free11 = tuple(theta11[j] for j in range(12) if j not in pivots_depth11)
    equations12 = bracket_equations(a11, c11, 12, g12)
    select12, conditions12, raw12, matrix12, pivots12 = eliminate_linear(
        equations12, free11
    )
    check((matrix12.shape, matrix12.rank(), pivots12)
          == ((13, 8), 8, tuple(range(8))), "row-twelve bracket selector")
    check(len(conditions12) == 2, "row-twelve has exactly two source conditions")
    nu_candidates = [value for value in conditions12 if nu18 in value.free_symbols]
    check(nu_candidates, "row-twelve nu18 condition present")
    nu_raw = nu_candidates[0]
    nu_coefficient_raw = sp.diff(nu_raw, nu18)
    nu_response = sp.Integer(29652680243293059502080000000)
    scale12 = exact(nu_response / nu_coefficient_raw)
    check(not scale12.free_symbols and scale12 != 0,
          "row-twelve nu normalization is constant")
    nu_equation = sp.expand(scale12 * nu_raw)
    check(sp.diff(nu_equation, nu18) == nu_response,
          "y6 row-twelve response coefficient")
    nu_graph = exact(-nu_equation.subs(nu18, 0) / nu_response)
    nu_subs = {nu18: nu_graph}
    after_nu = [resolve(value, nu_subs) for value in raw12]
    after_nu_nonzero = [value for value in after_nu if value != 0]
    L = -23839 * R8.Phi - 675 * R8.alpha11 - 675 * c23 + 4266 * R8.eta
    check(after_nu_nonzero and all(proportional(value, L) for value in after_nu_nonzero),
          "row-twelve residual ideal after nu is L")
    c23_graph = exact(sp.solve(L, c23)[0])
    check(sp.diff(L, c23) == -675 and constant_denominator(c23_graph) == 675,
          "row-twelve c23 constant graph")
    c23_subs = {c23: c23_graph}
    check(all(resolve(value, {**nu_subs, **c23_subs}) == 0 for value in raw12),
          "row-twelve bracket sufficiency")

    # Audit all exact-valuation-12 monomials p^(12-2b)y^b.  Their leading
    # rows are x^b.  Compute their literal cokernel response vectors after
    # the same terminal elimination, without choosing a quotient basis.
    response_columns: list[sp.Matrix] = []
    response_symbols12 = sp.symbols("rho12_0:7")
    for b, rho in enumerate(response_symbols12):
        exponent_p = 12 - 2 * b
        monomial = R8.p**exponent_p * R8.y**b
        check(tcoeff(monomial, 11) == 0 and tcoeff(monomial, 12) == x**b,
              f"valuation-twelve monomial b={b}")
        perturbed = bracket_equations(a11, c11, 12, g12 + rho * x**b)
        _, _, perturbed_raw, perturbed_matrix, perturbed_pivots = eliminate_linear(
            perturbed, free11
        )
        check(perturbed_matrix == matrix12 and perturbed_pivots == pivots12,
              f"response b={b} uses same tangent quotient")
        response_columns.append(sp.Matrix([exact(sp.diff(value, rho))
                                           for value in perturbed_raw]))
    response_matrix = sp.Matrix.hstack(*response_columns)
    check(response_matrix.rank() == 1, "exact-valuation-twelve response rank one")
    check(all(response_columns[b].is_zero_matrix for b in (1, 3, 5)),
          "all odd exact-valuation-twelve responses invisible")
    check(not response_columns[6].is_zero_matrix, "y6 response is visible")
    L_vector = sp.Matrix([
        exact(value / L) if value != 0 else sp.Integer(0)
        for value in after_nu
    ])
    check(sp.Matrix.hstack(response_columns[6], L_vector).rank() == 2,
          "same-row response line is transverse to L residual")

    selected_a11 = apply(a11, select12, nu_subs, c23_subs)
    selected_c11 = apply(c11, select12, nu_subs, c23_subs)
    a12_trial, c12_trial, theta12 = append_tangent(selected_a11, selected_c11, 12,
                                                   "audit12_theta12")
    depth12, ap12, cp12, _ = depth_equations(a12_trial, c12_trial, 12)
    terminal12, depth_conditions12, raw_depth12, matrix_depth12, pivots_depth12 = eliminate_linear(
        depth12, theta12
    )
    check((ap12.shape, ap12.rank(), cp12.shape, cp12.rank())
          == ((117, 267), 87, (130, 424), 105), "row-twelve depth universes")
    check((matrix_depth12.shape, matrix_depth12.rank(), pivots_depth12)
          == ((55, 13), 4, (9, 10, 11, 12)), "row-twelve depth selector")
    G = (
        -95023209310875 * R8.Phi**2
        - 157540468455000 * R8.Phi * R8.alpha11
        - 21622415197500 * R8.Phi * c51
        - 4375612912500 * R8.Phi * R8.eta
        - 21622415197500 * R8.alpha11 * R8.eta
        - 368130186720000 * c70
        - 78770234227500 * R8.eta**2
        + 178582220348850176
    )
    check(len(depth_conditions12) == 1
          and all(proportional(value, G) for value in depth_conditions12),
          "row-twelve unique depth condition G")
    c70_graph = exact(sp.solve(G, c70)[0])
    check(sp.diff(G, c70) == -368130186720000
          and constant_denominator(c70_graph) == 368130186720000,
          "row-twelve c70 constant graph")
    c70_subs = {c70: c70_graph}
    check(all(resolve(value, c70_subs) == 0 for value in raw_depth12),
          "row-twelve depth sufficiency")
    check(len(theta12) - len(pivots_depth12) == 9,
          "row-twelve terminal affine-nine fibre")

    # Resolve the entire triangular source graph to four free coordinates.
    full_graphs = {
        **row10_graphs,
        lambda18: lambda_graph,
        nu18: nu_graph,
        c23: c23_graph,
        c70: c70_graph,
    }
    resolved_graphs = {symbol: resolve(value, full_graphs)
                       for symbol, value in full_graphs.items()}
    base = {R8.Phi, R8.eta, R8.alpha11, c51}
    check(set().union(*(value.free_symbols for value in resolved_graphs.values())) <= base,
          "global affine-four source graph")
    check(all(not sp.fraction(sp.together(value))[1].free_symbols
              for value in resolved_graphs.values()),
          "all source graph denominators constant")

    # Dense and boundary controls evaluate every retained equation, while a
    # hostile off-L source confirms no same-row monomial response can repair L.
    all_residuals = (raw9 + raw_depth9 + raw10 + raw_depth10 + raw11
                     + raw_depth11 + raw12 + raw_depth12)

    def control(values: tuple[int, int, int, int]) -> dict[sp.Symbol, sp.Expr]:
        specialization = dict(zip((R8.Phi, R8.eta, R8.alpha11, c51),
                                  map(sp.Integer, values)))
        for symbol, value in resolved_graphs.items():
            specialization[symbol] = exact(value.subs(specialization))
        return specialization

    dense = control((1, 2, 3, 5))
    boundary = control((0, 0, 3, 5))
    check(all(exact(value.subs(dense)) == 0 for value in all_residuals),
          "dense literal rows-nine-through-twelve control")
    check(all(exact(value.subs(boundary)) == 0 for value in all_residuals),
          "Phi-eta-zero boundary control")
    hostile = {R8.Phi: 0, R8.eta: 0, R8.alpha11: 1, c23: 0}
    check(exact(L.subs(hostile)) == -675,
          "off-L hostile is genuinely nonzero")
    check(sp.Matrix.hstack(response_columns[6], L_vector).rank() == 2,
          "off-L hostile survives arbitrary same-row response")

    response_ratios: list[sp.Expr | str] = []
    pivot_entry = next(entry for entry in response_columns[6] if entry != 0)
    for column in response_columns:
        if column.is_zero_matrix:
            response_ratios.append(sp.Integer(0))
        else:
            entry = next(column[j] for j in range(column.rows)
                         if response_columns[6][j] != 0)
            response_ratios.append(exact(entry / pivot_entry))

    print("THM-4403 two-channel weight-18 row-twelve independent audit")
    print("imports=audited_THM4308_R8_and_THM4315_R9;late_candidates_import=no")
    print("universe=complete_weight_le_14_plus_lambda18_p3y4_plus_nu18_y6")
    print("row10_source=A6(Phi,eta,alpha11,c51,c23,c70);terminal=A8")
    print("row11_bracket=12x8_rank8;primitive_quartic_F=yes")
    print(f"row11_dF_dc23={sp.sstr(expected_dF)}")
    print(f"row11_lambda_response={response11};parameter_localization=none")
    print("row11_depth=P2(102x228,rank77),P3(114x361,rank94),joint(45x12,rank4,pivots8,9,10,11)")
    print("row12_bracket=13x8_rank8;source_conditions=nu_graph_and_L")
    print(f"row12_nu_response={nu_response}")
    print(f"row12_L={sp.sstr(L)};c23_pivot=-675")
    print("valuation12_monomials=(p12,p10y,p8y2,p6y3,p4y4,p2y5,y6)")
    print(f"valuation12_response_rank={response_matrix.rank()};ratios_to_y6={sp.sstr(tuple(response_ratios))}")
    print("valuation12_odd_responses=zero;L_payable_by_same_row_source=no")
    print(f"row12_G={sp.sstr(G)};c70_pivot=-368130186720000")
    print("row12_depth=P2(117x267,rank87),P3(130x424,rank105),joint(55x13,rank4,pivots9,10,11,12)")
    print("source_projection=A4(Phi,eta,alpha11,c51);terminal_fibre=A9")
    print("dense_control=PASS;Phi_eta_zero_control=PASS;off_L_hostile=PASS")
    print("field=characteristic_zero;finite_field_used=no")
    print("scope=restricted_two_channel_weight18_tail_through_row12_only;complete_weights15_to18_and_entry_open")
    print(f"checks={CHECKS} result=PASS")


if __name__ == "__main__":
    main()
