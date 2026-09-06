#!/usr/bin/env python3
"""Clean-room exact audit of the restricted source-normal row-13 tail.

The only imported mathematics is the audited THM-4308/4315 row operator.
This script reconstructs the complete residual-weight-at-most-14 source and
the two selected weight-18 corrections through row twelve, then independently
tests the first-visible valuation-13 source diagonal.  It deliberately does
not import or inspect the companion row-13 scout.
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
lambda18, nu18, kappa20 = sp.symbols("lambda18 nu18 kappa20")
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


def primitive(value: sp.Expr, positive_in: sp.Symbol | None = None) -> sp.Expr:
    numerator, _ = sp.fraction(sp.together(exact(value)))
    variables = tuple(sorted(numerator.free_symbols, key=str))
    if variables:
        numerator = sp.Poly(numerator, *variables, domain=sp.QQ).primitive()[1].as_expr()
    if positive_in is not None and sp.diff(numerator, positive_in) < 0:
        numerator = -numerator
    return sp.expand(numerator)


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
    for _ in range(20):
        updated = exact(answer.subs(mapping, simultaneous=False))
        if updated == answer:
            return answer
        answer = updated
    raise AssertionError("substitution graph did not resolve")


def eliminate_linear(
    equations: list[sp.Expr], variables: tuple[sp.Symbol, ...]
) -> tuple[dict[sp.Symbol, sp.Expr], list[sp.Expr], list[sp.Expr], sp.Matrix, tuple[int, ...]]:
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    pivots = tuple(matrix.rref()[1])
    row_pivots = tuple(matrix.T.rref()[1])
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    if pivots:
        square = matrix.extract(row_pivots, pivots)
        check(square.rows == square.cols and square.det() != 0, "linear pivot minor")
        values = square.inv() * rhs.extract(row_pivots, (0,))
        substitutions = {
            variables[column]: exact(values[index])
            for index, column in enumerate(pivots)
        }
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
    arows: list[sp.Expr], crows: list[sp.Expr], row: int
) -> tuple[list[sp.Expr], sp.Matrix, sp.Matrix]:
    acoords, amatrix = R9.depth_matrix(2, row)
    ccoords, cmatrix = R9.depth_matrix(3, row)
    avec = sp.Matrix([xcoeff(arows[n], degree) for n, degree in acoords])
    cvec = sp.Matrix([xcoeff(crows[n], degree) for n, degree in ccoords])
    aeq = [sp.expand((left.T * avec)[0]) for left in amatrix.T.nullspace()]
    ceq = [sp.expand((left.T * cvec)[0]) for left in cmatrix.T.nullspace()]
    return aeq + ceq, amatrix, cmatrix


def bracket_equations(
    arows: list[sp.Expr], crows: list[sp.Expr], row: int, source: sp.Expr
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
    arows: list[sp.Expr], crows: list[sp.Expr], row: int, prefix: str
) -> tuple[list[sp.Expr], list[sp.Expr], tuple[sp.Symbol, ...]]:
    abase, cbase = R8.particular_row(row, R8.B_row(row, arows, crows))
    symbols = tuple(sp.symbols(f"{prefix}_0:{row + 1}"))
    theta = sum(symbols[j] * x**j for j in range(row + 1))
    return (
        arows + [sp.expand(abase + theta * sp.diff(R8.A0, x))],
        crows + [sp.expand(cbase + theta * sp.diff(R8.C0, x))],
        symbols,
    )


def constant_denominator(value: sp.Expr) -> bool:
    denominator = sp.fraction(sp.together(value))[1]
    return not denominator.free_symbols and denominator != 0


def main() -> None:
    # Complete inherited weight <=14 source, the audited selected weight-18
    # tail, and precisely one new valuation-13 channel of minimum candidate
    # weight 20.
    face13 = c51 * R8.p**5 * R8.y + c23 * R8.p**2 * R8.y**3
    face14 = c70 * R8.p**7 + c42 * R8.p**4 * R8.y**2 + c14 * R8.p * R8.y**4
    late11 = lambda18 * R8.p**3 * R8.y**4
    late12 = nu18 * R8.y**6
    late13 = kappa20 * R8.p * R8.y**6
    source = sp.expand(R8.H + face13 + face14 + late11 + late12 + late13)
    grows = {row: tcoeff(source, row) for row in range(4, 14)}
    check(tcoeff(late13, 12) == 0 and tcoeff(late13, 13) == kappa20 * x**6,
          "p*y6 exact valuation thirteen")
    check(all(kappa20 not in grows[row].free_symbols for row in range(4, 13)),
          "row-thirteen source firewall")

    # Rows four through eight from the literal THM-4308 source/operator.
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

    theta8 = tuple(sp.symbols("ind13_theta8_0:9"))
    tangent8 = sum(theta8[j] * x**j for j in range(9))
    abase8, cbase8 = R8.particular_row(8, R8.B_row(8, arows, crows))
    arows.append(sp.expand(abase8 + tangent8 * sp.diff(R8.A0, x)))
    crows.append(sp.expand(cbase8 + tangent8 * sp.diff(R8.C0, x)))
    depth8, ap8, cp8 = depth_equations(arows, crows, 8)
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

    # Row nine, including the complete weight-14 c14 response.
    a8 = apply(arows, gate, terminal8, gate)
    c8 = apply(crows, gate, terminal8, gate)
    g9 = resolve(grows[9], {**bracket_subs, **gate})
    select9, conditions9, raw9, matrix9, pivots9 = eliminate_linear(
        bracket_equations(a8, c8, 9, g9), theta8[:7]
    )
    check((matrix9.shape, matrix9.rank(), pivots9)
          == ((10, 7), 7, tuple(range(7))), "row-nine bracket selector")
    e9 = (
        613527750 * R8.Phi**2 - 511211250 * R8.Phi * R8.alpha11
        - 3154140000 * R8.Phi * R8.eta - 255605625 * R8.eta**2
        + 6736896000 * R8.xi10 - 46483785515008
    )
    row9_expected = e9 - (3339765000 * c70 + 898128000 * c42 + 216513000 * c14)
    check(len(conditions9) == 1 and proportional(conditions9[0], row9_expected),
          "row-nine source equation")
    c14_graph = exact(sp.solve(row9_expected, c14)[0])
    selected_a8 = apply(a8, select9, {c14: c14_graph})
    selected_c8 = apply(c8, select9, {c14: c14_graph})
    a9_trial, c9_trial, theta9 = append_tangent(selected_a8, selected_c8, 9, "ind13_theta9")
    depth9, ap9, cp9 = depth_equations(a9_trial, c9_trial, 9)
    terminal9, depth_conditions9, raw_depth9, matrix_depth9, pivots_depth9 = eliminate_linear(depth9, theta9)
    check((ap9.shape, ap9.rank(), cp9.shape, cp9.rank())
          == ((75, 160), 59, (85, 251), 73), "row-nine depth universes")
    check((matrix_depth9.shape, matrix_depth9.rank(), pivots_depth9)
          == ((28, 10), 3, (7, 8, 9)), "row-nine depth selector")
    check(not depth_conditions9 and all(value == 0 for value in raw_depth9),
          "row-nine depth automatic")

    # Row ten: independently recover the four constant-pivot graphs.
    a9 = apply(a9_trial, terminal9)
    c9 = apply(c9_trial, terminal9)
    g10 = resolve(grows[10], {**bracket_subs, **gate, c14: c14_graph})
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
    remaining10 = [exact(value.subs(R8.xi10, xi_graph)) for value in conditions10
                   if value != xi_condition]
    remaining10 = [value for value in remaining10 if value != 0]
    check(len(remaining10) == 1, "row-ten second bracket condition")
    selected_a9 = apply(a9, select10, {R8.xi10: xi_graph})
    selected_c9 = apply(c9, select10, {R8.xi10: xi_graph})
    a10_trial, c10_trial, theta10 = append_tangent(selected_a9, selected_c9, 10, "ind13_theta10")
    depth10, ap10, cp10 = depth_equations(a10_trial, c10_trial, 10)
    terminal10, depth_conditions10, raw_depth10, matrix_depth10, pivots_depth10 = eliminate_linear(depth10, theta10)
    check((ap10.shape, ap10.rank(), cp10.shape, cp10.rank())
          == ((88, 193), 68, (99, 304), 83), "row-ten depth universes")
    check((matrix_depth10.shape, matrix_depth10.rank(), pivots_depth10)
          == ((36, 11), 3, (8, 9, 10)), "row-ten depth selector")
    beta_expected = -91 * R8.Phi + 15 * R8.beta11 + 18 * R8.eta
    check(len(depth_conditions10) == 1 and proportional(depth_conditions10[0], beta_expected),
          "row-ten beta equation")
    beta_graph = exact(sp.solve(beta_expected, R8.beta11)[0])
    c42_equation = primitive(remaining10[0].subs(R8.beta11, beta_graph))
    check(sp.diff(c42_equation, c42) == -89131914000, "row-ten c42 pivot")
    c42_graph = exact(sp.solve(c42_equation, c42)[0])
    xi_final = exact(xi_graph.subs(c42, c42_graph))
    c14_final = exact(c14_graph.subs(R8.xi10, xi_final).subs(c42, c42_graph))
    row10_graphs = {
        R8.beta11: beta_graph,
        c42: c42_graph,
        R8.xi10: xi_final,
        c14: c14_final,
    }
    check(set().union(*(value.free_symbols for value in row10_graphs.values()))
          <= {R8.Phi, R8.eta, R8.alpha11, c51, c23, c70},
          "row-ten affine-six graph")

    # Row eleven: reconstruct the first selected weight-18 correction.
    inherited10 = {**bracket_subs, **gate, **row10_graphs}
    a10 = apply(a10_trial, terminal10, row10_graphs)
    c10 = apply(c10_trial, terminal10, row10_graphs)
    g11 = resolve(grows[11], inherited10)
    free10 = tuple(theta10[j] for j in range(11) if j not in pivots_depth10)
    select11, conditions11, raw11, matrix11, pivots11 = eliminate_linear(
        bracket_equations(a10, c10, 11, g11), free10
    )
    check((matrix11.shape, matrix11.rank(), pivots11)
          == ((12, 8), 8, tuple(range(8))), "row-eleven bracket selector")
    check(len(conditions11) == 1, "row-eleven principal residual")
    full11 = conditions11[0]
    base11 = exact(full11.subs(lambda18, 0))
    normalization11 = exact(
        (sp.Integer(21736146783278091456000000) * R8.Phi) / sp.diff(base11, c23)
    )
    check(not normalization11.free_symbols and normalization11 != 0,
          "row-eleven constant normalization")
    equation11 = sp.expand(normalization11 * full11)
    F = sp.expand(equation11.subs(lambda18, 0))
    response11 = sp.diff(equation11, lambda18)
    check(response11 == 6864046352614134144000000, "row-eleven response")
    lambda_graph = exact(-F / response11)
    selected_a10 = apply(a10, select11, {lambda18: lambda_graph})
    selected_c10 = apply(c10, select11, {lambda18: lambda_graph})
    a11_trial, c11_trial, theta11 = append_tangent(selected_a10, selected_c10, 11, "ind13_theta11")
    depth11, ap11, cp11 = depth_equations(a11_trial, c11_trial, 11)
    terminal11, depth_conditions11, raw_depth11, matrix_depth11, pivots_depth11 = eliminate_linear(depth11, theta11)
    check((ap11.shape, ap11.rank(), cp11.shape, cp11.rank())
          == ((102, 228), 77, (114, 361), 94), "row-eleven depth universes")
    check((matrix_depth11.shape, matrix_depth11.rank(), pivots_depth11)
          == ((45, 12), 4, (8, 9, 10, 11)), "row-eleven depth selector")
    check(not depth_conditions11 and all(value == 0 for value in raw_depth11),
          "row-eleven depth automatic")

    # Row twelve: recover the nu graph, the independent L condition, and the
    # later c70 depth graph exactly as in the proved chain.
    a11 = apply(a11_trial, terminal11)
    c11 = apply(c11_trial, terminal11)
    inherited11 = {**inherited10, lambda18: lambda_graph}
    g12 = resolve(grows[12], inherited11)
    free11 = tuple(theta11[j] for j in range(12) if j not in pivots_depth11)
    select12, conditions12, raw12, matrix12, pivots12 = eliminate_linear(
        bracket_equations(a11, c11, 12, g12), free11
    )
    check((matrix12.shape, matrix12.rank(), pivots12)
          == ((13, 8), 8, tuple(range(8))), "row-twelve bracket selector")
    check(len(conditions12) == 2, "row-twelve two source conditions")
    nu_raw = next(value for value in conditions12 if nu18 in value.free_symbols)
    normalization12 = exact(sp.Integer(29652680243293059502080000000) / sp.diff(nu_raw, nu18))
    check(not normalization12.free_symbols and normalization12 != 0,
          "row-twelve constant normalization")
    nu_equation = sp.expand(normalization12 * nu_raw)
    nu_response = sp.diff(nu_equation, nu18)
    nu_graph = exact(-nu_equation.subs(nu18, 0) / nu_response)
    after_nu = [resolve(value, {nu18: nu_graph}) for value in raw12]
    L = -23839 * R8.Phi - 675 * R8.alpha11 - 675 * c23 + 4266 * R8.eta
    check(all(value == 0 or proportional(value, L) for value in after_nu),
          "row-twelve residual L")
    c23_graph = exact(sp.solve(L, c23)[0])
    selected_a11 = apply(a11, select12, {nu18: nu_graph}, {c23: c23_graph})
    selected_c11 = apply(c11, select12, {nu18: nu_graph}, {c23: c23_graph})
    a12_trial, c12_trial, theta12 = append_tangent(selected_a11, selected_c11, 12, "ind13_theta12")
    depth12, ap12, cp12 = depth_equations(a12_trial, c12_trial, 12)
    terminal12, depth_conditions12, raw_depth12, matrix_depth12, pivots_depth12 = eliminate_linear(depth12, theta12)
    check((ap12.shape, ap12.rank(), cp12.shape, cp12.rank())
          == ((117, 267), 87, (130, 424), 105), "row-twelve depth universes")
    check((matrix_depth12.shape, matrix_depth12.rank(), pivots_depth12)
          == ((55, 13), 4, (9, 10, 11, 12)), "row-twelve depth selector")
    G = (
        -95023209310875 * R8.Phi**2
        -157540468455000 * R8.Phi * R8.alpha11
        -21622415197500 * R8.Phi * c51
        -4375612912500 * R8.Phi * R8.eta
        -21622415197500 * R8.alpha11 * R8.eta
        -368130186720000 * c70
        -78770234227500 * R8.eta**2
        +178582220348850176
    )
    check(len(depth_conditions12) == 1 and proportional(depth_conditions12[0], G),
          "row-twelve unique depth condition")
    c70_graph = exact(sp.solve(G, c70)[0])
    check(all(resolve(value, {c70: c70_graph}) == 0 for value in raw_depth12),
          "row-twelve depth sufficiency")

    # Resolve the inherited A4 graph before any row-thirteen calculation.
    graphs12 = {
        **row10_graphs,
        lambda18: lambda_graph,
        nu18: nu_graph,
        c23: c23_graph,
        c70: c70_graph,
    }
    resolved12 = {symbol: resolve(value, graphs12) for symbol, value in graphs12.items()}
    base = {R8.Phi, R8.eta, R8.alpha11, c51}
    check(set().union(*(value.free_symbols for value in resolved12.values())) <= base,
          "row-twelve global A4 graph")
    check(all(constant_denominator(value) for value in resolved12.values()),
          "row-twelve graph has no localization")

    # Row thirteen.  The projected row-twelve depth fibre supplies nine
    # tangents.  Five formal cokernel coordinates remain; the inherited source
    # occupies exactly one line, which must be identified intrinsically.
    a12 = apply(a12_trial, terminal12, {c70: c70_graph})
    c12 = apply(c12_trial, terminal12, {c70: c70_graph})
    all_inherited = {**bracket_subs, **gate, **resolved12}
    g13 = resolve(grows[13], all_inherited)
    free12 = tuple(theta12[j] for j in range(13) if j not in pivots_depth12)
    equations13 = bracket_equations(a12, c12, 13, g13)
    select13, conditions13, raw13, matrix13, pivots13 = eliminate_linear(equations13, free12)
    check((matrix13.shape, matrix13.rank(), pivots13)
          == ((14, 9), 9, tuple(range(9))), "row-thirteen bracket selector")
    check(matrix13.rows - matrix13.rank() == 5, "row-thirteen raw cokernel five")
    check(len(conditions13) == 1, "row-thirteen selected residual line")

    moment13 = tuple(R9.primitive_student_row(13))
    difference13 = exact(g13 - R8.predicted_G(13, a12, c12))
    student13 = exact(sum(moment13[j] * xcoeff(difference13, j) for j in range(14)))
    R13 = primitive(student13, positive_in=kappa20)
    check(proportional(conditions13[0], R13), "row-thirteen residual is Student coordinate")
    response13 = sp.diff(R13, kappa20)
    check(response13 != 0 and not response13.free_symbols,
          "row-thirteen correction is a constant pivot")
    base13 = sp.expand(R13.subs(kappa20, 0))
    variables13 = (R8.Phi, R8.eta, R8.alpha11, c51)
    poly13 = sp.Poly(base13, *variables13, domain=sp.QQ)
    check(poly13.total_degree() == 4 and len(poly13.terms()) == 29,
          "row-thirteen 29-term quartic")
    kappa_graph = exact(-base13 / response13)
    check(constant_denominator(kappa_graph), "row-thirteen global constant graph")
    check(all(resolve(value, {kappa20: kappa_graph}) == 0 for value in raw13),
          "row-thirteen bracket sufficiency")

    # Audit the full exact-valuation-13 diagonal, without using a chosen
    # quotient basis.  Raw post-pivot residual vectors determine its rank.
    rho13 = sp.symbols("rho13_0:7")
    monomials13 = tuple(R8.p ** (13 - 2*b) * R8.y**b for b in range(7))
    weights13 = tuple(26 - b for b in range(7))
    response_columns: list[sp.Matrix] = []
    for b, rho in enumerate(rho13):
        monomial = monomials13[b]
        check(all(tcoeff(monomial, row) == 0 for row in range(13))
              and tcoeff(monomial, 13) == x**b,
              f"valuation-thirteen monomial b={b}")
        perturbed = bracket_equations(a12, c12, 13, g13 + rho * x**b)
        _, _, perturbed_raw, perturbed_matrix, perturbed_pivots = eliminate_linear(
            perturbed, free12
        )
        check(perturbed_matrix == matrix13 and perturbed_pivots == pivots13,
              f"valuation-thirteen tangent quotient b={b}")
        response_columns.append(sp.Matrix([
            exact(sp.diff(value, rho)) for value in perturbed_raw
        ]))
    response_matrix = sp.Matrix.hstack(*response_columns)
    check(response_matrix.rank() == 1, "valuation-thirteen response rank one")
    check(all(response_columns[b].is_zero_matrix for b in (1, 3, 5)),
          "valuation-thirteen odd channels vanish")
    check(all(not response_columns[b].is_zero_matrix for b in (0, 2, 4, 6)),
          "valuation-thirteen even channels visible")
    pivot_entry_index = next(j for j in range(response_columns[6].rows)
                             if response_columns[6][j] != 0)
    pivot_entry = response_columns[6][pivot_entry_index]
    ratios13 = tuple(
        sp.Integer(0) if column.is_zero_matrix
        else exact(column[pivot_entry_index] / pivot_entry)
        for column in response_columns
    )
    check(ratios13 == tuple(exact(sp.Integer(value) / moment13[6]) for value in moment13[:7]),
          "valuation-thirteen ratios are Student moments")
    check(weights13[6] == 20 and all(
        response_columns[b].is_zero_matrix for b in range(7) if weights13[b] < 20
    ), "minimum visible valuation-thirteen weight is twenty")

    # Projected depth at row thirteen must be automatic after the source graph.
    selected_a12 = apply(a12, select13, {kappa20: kappa_graph})
    selected_c12 = apply(c12, select13, {kappa20: kappa_graph})
    a13_trial, c13_trial, theta13 = append_tangent(selected_a12, selected_c12, 13, "ind13_theta13")
    depth13, ap13, cp13 = depth_equations(a13_trial, c13_trial, 13)
    terminal13, depth_conditions13, raw_depth13, matrix_depth13, pivots_depth13 = eliminate_linear(
        depth13, theta13
    )
    check((ap13.shape, ap13.rank(), cp13.shape, cp13.rank())
          == ((133, 308), 97, (147, 491), 117), "row-thirteen depth universes")
    check((matrix_depth13.shape, matrix_depth13.rank(), pivots_depth13)
          == ((66, 14), 5, (9, 10, 11, 12, 13)), "row-thirteen depth selector")
    check(not depth_conditions13 and all(value == 0 for value in raw_depth13),
          "row-thirteen depth automatic")
    check(len(theta13) - len(pivots_depth13) == 9,
          "row-thirteen terminal affine-nine fibre")

    # Entire triangular source graph and two literal controls.
    all_graphs = {**graphs12, kappa20: kappa_graph}
    resolved13 = {symbol: resolve(value, all_graphs) for symbol, value in all_graphs.items()}
    check(set().union(*(value.free_symbols for value in resolved13.values())) <= base,
          "row-thirteen source remains A4")
    check(all(constant_denominator(value) for value in resolved13.values()),
          "row-thirteen no source localization")
    residuals = raw9 + raw_depth9 + raw10 + raw_depth10 + raw11 + raw_depth11
    residuals += raw12 + raw_depth12 + raw13 + raw_depth13

    def audit_control(values: tuple[int, int, int, int]) -> None:
        point = dict(zip((R8.Phi, R8.eta, R8.alpha11, c51), map(sp.Integer, values)))
        check(all(resolve(value, {**bracket_subs, **gate, **resolved13}).subs(point) == 0
                  for value in residuals), f"literal control {values}")

    audit_control((1, 2, 3, 5))
    audit_control((0, 0, 3, 5))
    hostile = exact(R13.subs(kappa20, kappa_graph + 1))
    check(hostile == response13, "graph-plus-one hostile leaves exact response")

    print("restricted source-normal row-thirteen clean-room independent audit")
    print("imports=audited_THM4308_and_THM4315_operators;row13_scout_import=no")
    print("universe=complete_weight_le_14_plus_lambda18_p3y4_plus_nu18_y6_plus_kappa20_py6")
    print("inheritance=row12_source_A4(Phi,eta,alpha11,c51);terminal_A9")
    print("row13_bracket=14x9_rank9;raw_cokernel=5;selected_residual_dimension=1")
    print(f"row13_student_primitive={sp.sstr(moment13)}")
    print(f"row13_residual_terms={len(poly13.terms())};total_degree={poly13.total_degree()}")
    print(f"row13_residual={sp.sstr(base13)}")
    print(f"row13_kappa20_response={response13}")
    print(f"row13_kappa20_graph={sp.sstr(kappa_graph)}")
    print(f"valuation13_weights={weights13}")
    print(f"valuation13_response_rank={response_matrix.rank()};ratios_to_py6={sp.sstr(ratios13)}")
    print("valuation13_odd_responses=zero;minimum_visible_weight=20;channel=p*y^6")
    print("row13_depth=P2(133x308,rank97),P3(147x491,rank117),joint(66x14,rank5,pivots9,10,11,12,13)")
    print("source_projection=A4(Phi,eta,alpha11,c51);terminal_fibre=A9")
    print("dense_control=PASS;Phi_eta_zero_control=PASS;graph_plus_one_hostile=PASS")
    print("field=characteristic_zero;finite_field_used=no")
    print("scope=restricted_selected_tail_through_row13;not_complete_weight20;not_entry;not_JC2")
    print(f"checks={CHECKS} result=PASS")


if __name__ == "__main__":
    main()
