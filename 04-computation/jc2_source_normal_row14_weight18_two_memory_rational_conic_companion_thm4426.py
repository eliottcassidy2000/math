#!/usr/bin/env python3
"""Clean-room two-memory rational-conic companion to THM-4426.

Only the audited THM-4308/4315 source-normal row operators are imported; no
row-fourteen scratch scout is read or imported.  The certificate reconstructs
rows 4--14 after restoring both early weight-eighteen coordinates

    z=[p^9]H,  h18=[p^6*y^2]H.

It derives the global two-coordinate unpaid bracket conic and a rational
section over every characteristic-zero field.  On the exact boundary
Phi=eta=0, alpha11=1 it computes the projected-depth conic, proves that its
affine curve is a rational G_m, and verifies every bracket/depth coefficient
on the rational parameterization.  This is finite-row partial-source work,
not a full B_2, termination, chart-entry, JC(2), or DC(2) claim.
"""

from __future__ import annotations

import hashlib
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
z, h18, lambda18, nu18, kappa20, rho22 = sp.symbols(
    "z h18 lambda18 nu18 kappa20 rho22"
)
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
        numerator = (
            sp.Poly(numerator, *variables, domain=sp.QQ)
            .primitive()[1]
            .as_expr()
        )
    numerator = sp.expand(numerator)
    if positive_in is not None:
        coefficient = sp.diff(numerator, positive_in)
        if not coefficient.free_symbols and coefficient < 0:
            numerator = -numerator
    return sp.expand(numerator)


def proportional(left: sp.Expr, right: sp.Expr) -> bool:
    variables = tuple(sorted(left.free_symbols | right.free_symbols, key=str))
    lp = sp.Poly(
        sp.fraction(sp.together(exact(left)))[0], *variables, domain=sp.QQ
    )
    rp = sp.Poly(
        sp.fraction(sp.together(exact(right)))[0], *variables, domain=sp.QQ
    )
    return lp.monic() == rp.monic()


def apply(rows: list[sp.Expr], *maps: dict[sp.Symbol, sp.Expr]) -> list[sp.Expr]:
    answer = rows
    for mapping in maps:
        answer = [exact(row.subs(mapping)) for row in answer]
    return answer


def resolve(value: sp.Expr, mapping: dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    answer = exact(value)
    for _ in range(24):
        updated = exact(answer.subs(mapping, simultaneous=False))
        if updated == answer:
            return answer
        answer = updated
    raise AssertionError("substitution graph did not resolve")


def eliminate_linear(
    equations: list[sp.Expr], variables: tuple[sp.Symbol, ...]
) -> tuple[
    dict[sp.Symbol, sp.Expr],
    list[sp.Expr],
    list[sp.Expr],
    sp.Matrix,
    tuple[int, ...],
    tuple[int, ...],
]:
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
    return substitutions, conditions, raw, matrix, pivots, row_pivots


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
    return denominator != 0 and not denominator.free_symbols


def denominator(value: sp.Expr) -> int:
    result = sp.fraction(sp.together(value))[1]
    check(result.is_Integer and result > 0, "positive integral denominator")
    return int(result)


def expression_hash(value: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(value)).encode("ascii")).hexdigest()


def remainder_mod_q(value: sp.Expr, q: sp.Expr) -> sp.Expr:
    numerator, denominator_q = sp.fraction(sp.together(exact(value)))
    check(exact(denominator_q) != 0, "quotient denominator nonzero")
    numerator_remainder = sp.rem(sp.Poly(numerator, z, domain=sp.QQ), sp.Poly(q, z, domain=sp.QQ)).as_expr()
    denominator_remainder = sp.rem(
        sp.Poly(denominator_q, z, domain=sp.QQ), sp.Poly(q, z, domain=sp.QQ)
    ).as_expr()
    # All denominators in this certificate are constant.  Retaining the
    # general coprimality check makes the exact field argument explicit.
    check(sp.gcd(sp.Poly(denominator_remainder, z), sp.Poly(q, z)).degree() == 0,
          "denominator is a unit in the quadratic quotient")
    return exact(numerator_remainder / denominator_remainder)


def main() -> None:
    face13 = c51 * R8.p**5 * R8.y + c23 * R8.p**2 * R8.y**3
    face14 = c70 * R8.p**7 + c42 * R8.p**4 * R8.y**2 + c14 * R8.p * R8.y**4
    memory9 = z * R8.p**9
    memory10 = h18 * R8.p**6 * R8.y**2
    late11 = lambda18 * R8.p**3 * R8.y**4
    late12 = nu18 * R8.y**6
    late13 = kappa20 * R8.p * R8.y**6
    paid14 = rho22 * R8.p**2 * R8.y**6
    source = sp.expand(
        R8.H + face13 + face14 + memory9 + memory10
        + late11 + late12 + late13 + paid14
    )
    grows = {row: tcoeff(source, row) for row in range(4, 15)}
    expected_p9 = tuple(
        sp.binomial(9, row - 9) * x ** (2 * (row - 9))
        for row in range(9, 15)
    )
    check(
        tuple(tcoeff(R8.p**9, row) for row in range(9, 15)) == expected_p9,
        "p9 binomial transport rows nine through fourteen",
    )
    check(
        all(z not in grows[row].free_symbols for row in range(4, 9)),
        "p9 leaves rows four through eight unchanged",
    )

    # Rows four through eight, reconstructed literally from the two audited
    # frozen operators.
    arows = [R8.A0, R8.A1, R8.A2, R8.A3]
    crows = [R8.C0, R8.C1, R8.C2, R8.C3]
    bracket_subs: dict[sp.Symbol, sp.Expr] = {}
    response_symbols = {5: R8.upsilon5, 6: R8.U, 7: R8.W, 8: R8.Z}
    for n in range(4, 8):
        abase, cbase = R8.particular_row(n, R8.B_row(n, arows, crows))
        row = n + 1
        difference = sp.expand(
            grows[row].subs(bracket_subs)
            - R8.predicted_G(row, arows + [abase], crows + [cbase])
        )
        moment = R9.primitive_student_row(row)
        obstruction = exact(
            sum(moment[j] * xcoeff(difference, j) for j in range(row + 1))
        )
        answers = sp.solve(obstruction, response_symbols[row])
        check(len(answers) == 1, f"row-{row} scalar response")
        bracket_subs[response_symbols[row]] = exact(answers[0])
        target = sp.expand(difference.subs(response_symbols[row], answers[0]))
        theta = R8.tangent_solve(row, target)
        arows.append(sp.expand(abase + theta * sp.diff(R8.A0, x)))
        crows.append(sp.expand(cbase + theta * sp.diff(R8.C0, x)))
        check(
            exact(
                R8.predicted_G(row, arows, crows)
                - grows[row].subs(bracket_subs)
            )
            == 0,
            f"row-{row} bracket reconstruction",
        )

    theta8 = tuple(sp.symbols("memref_theta8_0:9"))
    tangent8 = sum(theta8[j] * x**j for j in range(9))
    abase8, cbase8 = R8.particular_row(8, R8.B_row(8, arows, crows))
    arows.append(sp.expand(abase8 + tangent8 * sp.diff(R8.A0, x)))
    crows.append(sp.expand(cbase8 + tangent8 * sp.diff(R8.C0, x)))
    depth8, ap8, cp8 = depth_equations(arows, crows, 8)
    terminal8, conditions8, raw_depth8, matrix_depth8, pivots8, _ = eliminate_linear(
        depth8, theta8
    )
    check(
        (ap8.shape, ap8.rank(), cp8.shape, cp8.rank())
        == ((63, 131), 51, (72, 204), 63),
        "row-eight depth universes",
    )
    check(
        (matrix_depth8.shape, matrix_depth8.rank(), pivots8)
        == ((21, 9), 2, (7, 8)),
        "row-eight terminal selection",
    )
    gate = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * R8.Phi,
    }
    old_gate = (
        896 - 15 * R8.Delta,
        3 * R8.Phi + 2 * R8.zeta3,
        3030 * R8.Delta + 225 * R8.Theta - 182528,
    )
    check(
        len(conditions8) == 3
        and all(any(proportional(a, b) for a in conditions8) for b in old_gate),
        "row-eight Hasse ideal",
    )
    check(
        all(exact(value.subs(gate)) == 0 for value in raw_depth8),
        "row-eight depth sufficiency",
    )

    # Row nine: z first appears and is retained; c14 absorbs the condition.
    a8 = apply(arows, gate, terminal8, gate)
    c8 = apply(crows, gate, terminal8, gate)
    g9 = resolve(grows[9], {**bracket_subs, **gate})
    eq9 = bracket_equations(a8, c8, 9, g9)
    select9, conditions9, raw9, matrix9, pivots9, _ = eliminate_linear(
        eq9, theta8[:7]
    )
    check(
        (matrix9.shape, matrix9.rank(), pivots9, len(conditions9))
        == ((10, 7), 7, tuple(range(7)), 1),
        "row-nine bracket selector with p9 memory",
    )
    equation9 = primitive(conditions9[0], positive_in=c14)
    check(
        sp.diff(equation9, c14) == 216513000,
        "row-nine c14 constant pivot",
    )
    c14_graph = exact(sp.solve(equation9, c14)[0])
    check(constant_denominator(c14_graph), "row-nine c14 constant denominator")
    selected_a8 = apply(a8, select9, {c14: c14_graph})
    selected_c8 = apply(c8, select9, {c14: c14_graph})
    a9_trial, c9_trial, theta9 = append_tangent(
        selected_a8, selected_c8, 9, "memref_theta9"
    )
    depth9, ap9, cp9 = depth_equations(a9_trial, c9_trial, 9)
    terminal9, depth_conditions9, raw_depth9, matrix_depth9, pivots_depth9, _ = eliminate_linear(
        depth9, theta9
    )
    check(
        (ap9.shape, ap9.rank(), cp9.shape, cp9.rank())
        == ((75, 160), 59, (85, 251), 73),
        "row-nine depth universes",
    )
    check(
        (matrix_depth9.shape, matrix_depth9.rank(), pivots_depth9)
        == ((28, 10), 3, (7, 8, 9)),
        "row-nine depth selector",
    )
    check(
        not depth_conditions9 and all(value == 0 for value in raw_depth9),
        "row-nine depth automatic",
    )

    # Row ten: the same pivot pattern survives, while the solved graphs now
    # retain z.  Conditions are identified by their z=0 canonical limits.
    a9 = apply(a9_trial, terminal9)
    c9 = apply(c9_trial, terminal9)
    g10 = resolve(grows[10], {**bracket_subs, **gate, c14: c14_graph})
    free9 = tuple(theta9[j] for j in range(10) if j not in pivots_depth9)
    eq10 = bracket_equations(a9, c9, 10, g10)
    select10, conditions10, raw10, matrix10, pivots10, _ = eliminate_linear(
        eq10, free9
    )
    check(
        (matrix10.shape, matrix10.rank(), pivots10)
        == ((11, 7), 7, tuple(range(7))),
        "row-ten bracket selector with p9 memory",
    )
    xi_expected = (
        13365000 * R8.Phi**2
        + 15035625 * R8.Phi * R8.eta
        + 6014250 * c42
        + 50787000 * c70
        + 57672000 * R8.xi10
        - 964604821504
    )
    xi_candidates = [
        value for value in conditions10
        if proportional(value.subs({z: 0, h18: 0}), xi_expected)
    ]
    check(len(xi_candidates) == 1, "row-ten canonical xi limit")
    xi_equation = primitive(xi_candidates[0], positive_in=R8.xi10)
    check(
        sp.diff(xi_equation, R8.xi10) == 57672000,
        "row-ten xi constant pivot",
    )
    xi_graph = exact(sp.solve(xi_equation, R8.xi10)[0])
    remaining10 = [
        exact(value.subs(R8.xi10, xi_graph))
        for value in conditions10 if value != xi_candidates[0]
    ]
    remaining10 = [value for value in remaining10 if value != 0]
    check(len(remaining10) == 1, "row-ten second bracket condition")
    selected_a9 = apply(a9, select10, {R8.xi10: xi_graph})
    selected_c9 = apply(c9, select10, {R8.xi10: xi_graph})
    a10_trial, c10_trial, theta10 = append_tangent(
        selected_a9, selected_c9, 10, "memref_theta10"
    )
    depth10, ap10, cp10 = depth_equations(a10_trial, c10_trial, 10)
    terminal10, depth_conditions10, raw_depth10, matrix_depth10, pivots_depth10, _ = eliminate_linear(
        depth10, theta10
    )
    check(
        (ap10.shape, ap10.rank(), cp10.shape, cp10.rank())
        == ((88, 193), 68, (99, 304), 83),
        "row-ten depth universes",
    )
    check(
        (matrix_depth10.shape, matrix_depth10.rank(), pivots_depth10)
        == ((36, 11), 3, (8, 9, 10)),
        "row-ten depth selector",
    )
    beta_expected = -91 * R8.Phi + 15 * R8.beta11 + 18 * R8.eta
    beta_candidates = [
        value for value in depth_conditions10
        if proportional(value.subs({z: 0, h18: 0}), beta_expected)
    ]
    check(len(beta_candidates) == 1, "row-ten canonical beta limit")
    beta_equation = primitive(beta_candidates[0], positive_in=R8.beta11)
    check(
        sp.diff(beta_equation, R8.beta11) == 15,
        "row-ten beta constant pivot",
    )
    beta_graph = exact(sp.solve(beta_equation, R8.beta11)[0])
    c42_equation = primitive(
        remaining10[0].subs(R8.beta11, beta_graph), positive_in=c42
    )
    check(sp.diff(c42_equation, c42) == 89131914000, "row-ten c42 pivot")
    c42_graph = exact(sp.solve(c42_equation, c42)[0])
    xi_final = exact(xi_graph.subs(c42, c42_graph))
    c14_final = resolve(c14_graph, {R8.xi10: xi_final, c42: c42_graph})
    row10_graphs = {
        R8.beta11: beta_graph,
        c42: c42_graph,
        R8.xi10: xi_final,
        c14: c14_final,
    }
    check(
        all(constant_denominator(value) for value in row10_graphs.values()),
        "row-ten graphs have constant denominators",
    )

    # Row eleven: lambda retains its proved constant response.
    inherited10 = {**bracket_subs, **gate, **row10_graphs}
    a10 = apply(a10_trial, terminal10, row10_graphs)
    c10 = apply(c10_trial, terminal10, row10_graphs)
    g11 = resolve(grows[11], inherited10)
    free10 = tuple(theta10[j] for j in range(11) if j not in pivots_depth10)
    eq11 = bracket_equations(a10, c10, 11, g11)
    select11, conditions11, raw11, matrix11, pivots11, _ = eliminate_linear(
        eq11, free10
    )
    check(
        (matrix11.shape, matrix11.rank(), pivots11, len(conditions11))
        == ((12, 8), 8, tuple(range(8)), 1),
        "row-eleven bracket selector with p9 memory",
    )
    raw_equation11 = conditions11[0]
    scale11 = exact(
        sp.Integer(6864046352614134144000000)
        / sp.diff(raw_equation11, lambda18)
    )
    check(not scale11.free_symbols and scale11 != 0, "row-eleven normalization")
    equation11 = sp.expand(scale11 * raw_equation11)
    check(
        sp.diff(equation11, lambda18) == 6864046352614134144000000,
        "row-eleven lambda constant pivot",
    )
    lambda_graph = exact(sp.solve(equation11, lambda18)[0])
    selected_a10 = apply(a10, select11, {lambda18: lambda_graph})
    selected_c10 = apply(c10, select11, {lambda18: lambda_graph})
    a11_trial, c11_trial, theta11 = append_tangent(
        selected_a10, selected_c10, 11, "memref_theta11"
    )
    depth11, ap11, cp11 = depth_equations(a11_trial, c11_trial, 11)
    terminal11, depth_conditions11, raw_depth11, matrix_depth11, pivots_depth11, _ = eliminate_linear(
        depth11, theta11
    )
    check(
        (ap11.shape, ap11.rank(), cp11.shape, cp11.rank())
        == ((102, 228), 77, (114, 361), 94),
        "row-eleven depth universes",
    )
    check(
        (matrix_depth11.shape, matrix_depth11.rank(), pivots_depth11)
        == ((45, 12), 4, (8, 9, 10, 11)),
        "row-eleven depth selector",
    )
    check(
        not depth_conditions11 and all(value == 0 for value in raw_depth11),
        "row-eleven depth automatic",
    )

    # Row twelve: nu pays one condition, c23 the transported side condition,
    # and c70 the projected-depth condition.  Every pivot remains constant.
    a11 = apply(a11_trial, terminal11)
    c11 = apply(c11_trial, terminal11)
    inherited11 = {**inherited10, lambda18: lambda_graph}
    g12 = resolve(grows[12], inherited11)
    free11 = tuple(theta11[j] for j in range(12) if j not in pivots_depth11)
    eq12 = bracket_equations(a11, c11, 12, g12)
    select12, conditions12, raw12, matrix12, pivots12, _ = eliminate_linear(
        eq12, free11
    )
    check(
        (matrix12.shape, matrix12.rank(), pivots12, len(conditions12))
        == ((13, 8), 8, tuple(range(8)), 2),
        "row-twelve bracket selector with p9 memory",
    )
    nu_candidates = [value for value in conditions12 if nu18 in value.free_symbols]
    check(len(nu_candidates) == 1, "row-twelve unique nu equation")
    scale12 = exact(
        sp.Integer(29652680243293059502080000000)
        / sp.diff(nu_candidates[0], nu18)
    )
    check(not scale12.free_symbols and scale12 != 0, "row-twelve normalization")
    nu_equation = sp.expand(scale12 * nu_candidates[0])
    check(
        sp.diff(nu_equation, nu18) == 29652680243293059502080000000,
        "row-twelve nu constant pivot",
    )
    nu_graph = exact(sp.solve(nu_equation, nu18)[0])
    after_nu = [resolve(value, {nu18: nu_graph}) for value in raw12]
    side12 = [value for value in after_nu if value != 0]
    side12_primitives: list[sp.Expr] = []
    for value in side12:
        candidate = primitive(value, positive_in=c23)
        if all(not proportional(candidate, old) for old in side12_primitives):
            side12_primitives.append(candidate)
    check(len(side12_primitives) == 1, "row-twelve transported side class")
    equation_c23 = side12_primitives[0]
    check(sp.diff(equation_c23, c23) == 675, "row-twelve c23 constant pivot")
    c23_graph = exact(sp.solve(equation_c23, c23)[0])
    selected_a11 = apply(a11, select12, {nu18: nu_graph}, {c23: c23_graph})
    selected_c11 = apply(c11, select12, {nu18: nu_graph}, {c23: c23_graph})
    a12_trial, c12_trial, theta12 = append_tangent(
        selected_a11, selected_c11, 12, "memref_theta12"
    )
    depth12, ap12, cp12 = depth_equations(a12_trial, c12_trial, 12)
    terminal12, depth_conditions12, raw_depth12, matrix_depth12, pivots_depth12, _ = eliminate_linear(
        depth12, theta12
    )
    check(
        (ap12.shape, ap12.rank(), cp12.shape, cp12.rank())
        == ((117, 267), 87, (130, 424), 105),
        "row-twelve depth universes",
    )
    check(
        (matrix_depth12.shape, matrix_depth12.rank(), pivots_depth12, len(depth_conditions12))
        == ((55, 13), 4, (9, 10, 11, 12), 1),
        "row-twelve depth selector",
    )
    equation_c70 = primitive(depth_conditions12[0], positive_in=c70)
    check(
        sp.diff(equation_c70, c70) == 368130186720000,
        "row-twelve c70 constant pivot",
    )
    c70_graph = exact(sp.solve(equation_c70, c70)[0])
    check(
        all(resolve(value, {c70: c70_graph}) == 0 for value in raw_depth12),
        "row-twelve depth sufficiency",
    )

    graphs12 = {
        **row10_graphs,
        lambda18: lambda_graph,
        nu18: nu_graph,
        c23: c23_graph,
        c70: c70_graph,
    }
    resolved12 = {symbol: resolve(value, graphs12) for symbol, value in graphs12.items()}
    base5 = {R8.Phi, R8.eta, R8.alpha11, c51, z, h18}
    check(
        set().union(*(value.free_symbols for value in resolved12.values())) <= base5,
        "row-twelve source graph lies over affine five",
    )
    check(
        all(constant_denominator(value) for value in resolved12.values()),
        "row-twelve source graph has no localization",
    )

    # Row thirteen: kappa again gives the proved constant pivot and depth is
    # automatic, now over A^5 with z retained.
    a12 = apply(a12_trial, terminal12, {c70: c70_graph})
    c12 = apply(c12_trial, terminal12, {c70: c70_graph})
    inherited12 = {**bracket_subs, **gate, **resolved12}
    g13 = resolve(grows[13], inherited12)
    free12 = tuple(theta12[j] for j in range(13) if j not in pivots_depth12)
    eq13 = bracket_equations(a12, c12, 13, g13)
    select13, conditions13, raw13, matrix13, pivots13, _ = eliminate_linear(
        eq13, free12
    )
    check(
        (matrix13.shape, matrix13.rank(), pivots13, len(conditions13))
        == ((14, 9), 9, tuple(range(9)), 1),
        "row-thirteen bracket selector with p9 memory",
    )
    scale13 = exact(
        sp.Integer(146361421124462229507072000000000)
        / sp.diff(conditions13[0], kappa20)
    )
    check(not scale13.free_symbols and scale13 != 0, "row-thirteen normalization")
    equation13 = sp.expand(scale13 * conditions13[0])
    check(
        sp.diff(equation13, kappa20)
        == 146361421124462229507072000000000,
        "row-thirteen kappa constant pivot",
    )
    kappa_graph = exact(sp.solve(equation13, kappa20)[0])
    selected_a12 = apply(a12, select13, {kappa20: kappa_graph})
    selected_c12 = apply(c12, select13, {kappa20: kappa_graph})
    a13_trial, c13_trial, theta13 = append_tangent(
        selected_a12, selected_c12, 13, "memref_theta13"
    )
    depth13, ap13, cp13 = depth_equations(a13_trial, c13_trial, 13)
    terminal13, depth_conditions13, raw_depth13, matrix_depth13, pivots_depth13, _ = eliminate_linear(
        depth13, theta13
    )
    check(
        (ap13.shape, ap13.rank(), cp13.shape, cp13.rank())
        == ((133, 308), 97, (147, 491), 117),
        "row-thirteen depth universes",
    )
    check(
        (matrix_depth13.shape, matrix_depth13.rank(), pivots_depth13)
        == ((66, 14), 5, (9, 10, 11, 12, 13)),
        "row-thirteen depth selector",
    )
    check(
        not depth_conditions13 and all(value == 0 for value in raw_depth13),
        "row-thirteen depth automatic",
    )

    graphs13 = {**graphs12, kappa20: kappa_graph}
    resolved13 = {symbol: resolve(value, graphs13) for symbol, value in graphs13.items()}
    check(
        set().union(*(value.free_symbols for value in resolved13.values())) <= base5,
        "row-thirteen source graph remains affine-five base",
    )
    check(
        all(constant_denominator(value) for value in resolved13.values()),
        "row-thirteen source graph has no localization",
    )
    a13 = apply(a13_trial, terminal13)
    c13 = apply(c13_trial, terminal13)
    a13 = [resolve(value, resolved13) for value in a13]
    c13 = [resolve(value, resolved13) for value in c13]
    free13 = tuple(theta13[j] for j in range(14) if j not in pivots_depth13)
    check(len(free13) == 9, "row-thirteen terminal A9")

    # Row fourteen: form the direct 15x9 coefficient quotient.  Normalize
    # its unpaid class so that z=0 is literally THM-4415's frozen J14.
    g14 = resolve(grows[14], {**bracket_subs, **gate, **resolved13})
    eq14 = bracket_equations(a13, c13, 14, g14)
    select14, conditions14, raw14, matrix14, pivots14, _ = eliminate_linear(
        eq14, free13
    )
    check(
        (matrix14.shape, matrix14.rank(), pivots14, len(conditions14))
        == ((15, 9), 9, tuple(range(9)), 2),
        "row-fourteen bracket quotient with memory",
    )
    unpaid_candidates = [value for value in conditions14 if rho22 not in value.free_symbols]
    paid_candidates = [value for value in conditions14 if rho22 in value.free_symbols]
    check(
        len(unpaid_candidates) == len(paid_candidates) == 1,
        "one unpaid and one rho-paid class",
    )
    unpaid_raw = unpaid_candidates[0]
    origin4 = {
        R8.Phi: 0, R8.eta: 0, R8.alpha11: 0,
        c51: 0, z: 0, h18: 0,
    }
    frozen_j14_origin = sp.Integer(
        689373657131467757187366756679487062016
    )
    scale_unpaid = exact(frozen_j14_origin / unpaid_raw.subs(origin4))
    check(not scale_unpaid.free_symbols and scale_unpaid != 0, "J18 normalization")
    j18 = sp.expand(scale_unpaid * unpaid_raw)
    j14 = sp.expand(j18.subs({z: 0, h18: 0}))
    check(
        expression_hash(j14)
        == "d15ce156dea6b7883378afebf23262d7950e351b845b0894986cbff39938884f",
        "literal THM4415 J14 hash",
    )
    j18_poly = sp.Poly(j18, z, domain=sp.EX)
    check(j18_poly.degree() == 2, "unpaid class quadratic in z")
    quadratic_coefficient = exact(j18_poly.coeff_monomial(z**2))
    linear_coefficient = exact(j18_poly.coeff_monomial(z))
    check(
        quadratic_coefficient == 643145476450643616480975000000,
        "unpaid quadratic coefficient",
    )
    expected_linear = (
        -211041194570570137414961250000 * R8.Phi**2
        + 1125125025346225769232150000000 * R8.Phi * R8.alpha11
        + 127589692647756992119425000000 * R8.Phi * c51
        + 373391507810486191201875000000 * R8.Phi * R8.eta
        + 127589692647756992119425000000 * R8.alpha11 * R8.eta
        + 562562512673112884616075000000 * R8.eta**2
        + 63406745144628813691929829048320000
    )
    check(
        linear_coefficient.subs(h18, 0) == expected_linear,
        "unpaid linear coefficient restricts to the p9 slice",
    )

    paid_raw = paid_candidates[0]
    frozen_rho_pivot = sp.Integer(
        1536794921806853409824256000000000
    )
    scale_paid = exact(frozen_rho_pivot / sp.diff(paid_raw, rho22))
    check(not scale_paid.free_symbols and scale_paid != 0, "paid normalization")
    paid_equation = sp.expand(scale_paid * paid_raw)
    check(
        sp.diff(paid_equation, rho22) == frozen_rho_pivot,
        "constant rho pivot",
    )
    rho_graph_global = exact(-paid_equation.subs(rho22, 0) / frozen_rho_pivot)
    check(constant_denominator(rho_graph_global), "global rho constant denominator")

    two_parameter = sp.Poly(j18, z, h18, domain=sp.EX)
    j_h2 = exact(two_parameter.coeff_monomial(h18**2))
    j_hz = exact(two_parameter.coeff_monomial(h18*z))
    j_z2 = exact(two_parameter.coeff_monomial(z**2))
    j_h = exact(two_parameter.coeff_monomial(h18))
    j_z = exact(two_parameter.coeff_monomial(z))
    j_zero = exact(two_parameter.coeff_monomial(1))
    global_u = 801*h18 + 27826*z
    global_v = 46174707*h18 + 1363633922*z
    global_k = sp.Rational(16949646393750000, 1)
    check(
        j_h2*h18**2+j_hz*h18*z+j_z2*z**2
        == sp.expand(global_k*global_u*global_v),
        "global homogeneous conic rationally splits",
    )
    inverse_uv = sp.solve(
        (sp.Eq(sp.Symbol("u"), global_u), sp.Eq(sp.Symbol("v"), global_v)),
        (h18, z), dict=True,
    )
    check(len(inverse_uv) == 1, "global packet coordinate inverse")
    u_symbol, v_symbol = sp.symbols("u v")
    h_uv = exact(inverse_uv[0][h18].subs({sp.Symbol("u"): u_symbol, sp.Symbol("v"): v_symbol}))
    z_uv = exact(inverse_uv[0][z].subs({sp.Symbol("u"): u_symbol, sp.Symbol("v"): v_symbol}))
    alpha_global = exact(j_h*sp.diff(h_uv, u_symbol)+j_z*sp.diff(z_uv, u_symbol))
    beta_global = exact(j_h*sp.diff(h_uv, v_symbol)+j_z*sp.diff(z_uv, v_symbol))
    check(
        j18 == sp.expand(
            global_k*global_u*global_v
            + alpha_global*global_u + beta_global*global_v + j_zero
        ),
        "global bilinear normal form",
    )
    gamma_global = exact(j_zero-alpha_global*beta_global/global_k)
    tau = sp.symbols("tau", nonzero=True)
    u_global_parameter = exact(tau-beta_global/global_k)
    v_global_parameter = exact(-gamma_global/(global_k*tau)-alpha_global/global_k)
    h_global_parameter = exact(h_uv.subs({u_symbol: u_global_parameter, v_symbol: v_global_parameter}))
    z_global_parameter = exact(z_uv.subs({u_symbol: u_global_parameter, v_symbol: v_global_parameter}))
    check(
        exact(j18.subs({h18: h_global_parameter, z: z_global_parameter})) == 0,
        "global bracket conic rational parameterization",
    )
    rho_global_parameter = exact(
        rho_graph_global.subs({h18: h_global_parameter, z: z_global_parameter})
    )
    global_section = {
        h18: exact(h_global_parameter.subs(tau, 1)),
        z: exact(z_global_parameter.subs(tau, 1)),
    }
    global_section[rho22] = exact(rho_global_parameter.subs(tau, 1))
    check(
        all(constant_denominator(value) for value in global_section.values()),
        "global bracket section has only constant denominators",
    )
    check(
        exact(j18.subs(global_section)) == 0
        and exact(paid_equation.subs(global_section)) == 0,
        "global bracket section pays both conditions",
    )
    global_bracket_after = [
        resolve(value, {**select14, **global_section}) for value in eq14
    ]
    for index, value in enumerate(global_bracket_after):
        check(value == 0, f"global bracket section coefficient {index}")
    print("two_parameter_weight18_bracket_development")
    print(f"rank_profile=row14_{matrix14.shape}_rank{matrix14.rank()}_conditions{len(conditions14)}")
    print(f"degrees_z_h_total={two_parameter.degree(z)},{two_parameter.degree(h18)},{two_parameter.total_degree()}")
    print(f"term_count_z_h={len(two_parameter.terms())}")
    print("J18zh=" + sp.sstr(j18))
    print("J18zh_homogeneous_factor=" + sp.sstr(global_k*global_u*global_v))
    print(f"global_uv_inverse=h:{sp.sstr(h_uv)};z:{sp.sstr(z_uv)}")
    print("global_alpha_u=" + sp.sstr(alpha_global))
    print("global_beta_v=" + sp.sstr(beta_global))
    print("global_gamma=" + sp.sstr(gamma_global))
    print("global_parameter=U:tau;V:-gamma/(k*tau);u:U-beta/k;v:V-alpha/k")
    print(
        "global_tau1_section_hashes="
        + sp.sstr(tuple((str(key), expression_hash(value)) for key, value in global_section.items()))
    )
    print(f"rho_pivot={frozen_rho_pivot};rho_graph_denominator={denominator(rho_graph_global)}")
    print(f"checks={CHECKS};result=PASS")

    # Exact boundary Phi=eta=0, alpha11=1.  J18 becomes linear in c51 and
    # therefore gives a constant-denominator graph over the z-line.
    boundary = {R8.Phi: 0, R8.eta: 0, R8.alpha11: 1}
    j_boundary_raw = sp.expand(j18.subs(boundary))
    boundary_scale = exact(
        sp.Integer(45455611113114234086880000000)
        / sp.diff(j_boundary_raw, c51)
    )
    check(
        not boundary_scale.free_symbols and boundary_scale != 0,
        "two-parameter boundary normalization",
    )
    j_boundary = sp.expand(boundary_scale * j_boundary_raw)
    d_c51 = exact(sp.diff(j_boundary, c51))
    check(
        d_c51 == 45455611113114234086880000000,
        "boundary c51 pivot",
    )
    boundary_poly = sp.Poly(j_boundary, z, h18, c51, domain=sp.QQ)
    check(
        boundary_poly.degree(c51) == 1
        and boundary_poly.degree(z) == 2
        and boundary_poly.degree(h18) == 2,
        "boundary unpaid bidegree",
    )
    j_boundary_h0 = sp.expand(j_boundary.subs(h18, 0))
    a0 = exact(j_boundary_h0.coeff(z, 2))
    b0 = exact(j_boundary_h0.coeff(z, 1))
    c0 = exact(j_boundary_h0.subs({z: 0, c51: 0}))
    check(
        (a0, b0, c0)
        == (
            sp.Integer(10049148069541306507515234375),
            sp.Integer(990730392884825213936403578880000),
            sp.Integer(10771463517407623861844320362236985344),
        ),
        "boundary unpaid coefficients",
    )
    c51_graph = exact(sp.solve(j_boundary, c51)[0])
    check(constant_denominator(c51_graph), "boundary c51 constant denominator")
    rho_graph = resolve(rho_graph_global, {**boundary, c51: c51_graph})
    check(constant_denominator(rho_graph), "boundary rho constant denominator")

    # All fifteen bracket coefficients vanish modulo the bracket equation.
    bracket_after = [
        resolve(value, {**select14, **boundary, c51: c51_graph, rho22: rho_graph})
        for value in eq14
    ]
    for index, value in enumerate(bracket_after):
        check(value == 0, f"boundary bracket coefficient {index}")

    selected_a13 = apply(
        a13, select14, boundary, {c51: c51_graph}, {rho22: rho_graph}
    )
    selected_c13 = apply(
        c13, select14, boundary, {c51: c51_graph}, {rho22: rho_graph}
    )
    a14_trial, c14_trial, theta14 = append_tangent(
        selected_a13, selected_c13, 14, "memref_theta14"
    )
    depth14, ap14, cp14 = depth_equations(a14_trial, c14_trial, 14)
    terminal14, depth_conditions14, raw_depth14, matrix_depth14, pivots_depth14, row_pivots14 = eliminate_linear(
        depth14, theta14
    )
    check(
        (ap14.shape, ap14.rank(), cp14.shape, cp14.rank())
        == ((150, 353), 108, (165, 564), 129),
        "row-fourteen depth universes",
    )
    check(
        (matrix_depth14.shape, matrix_depth14.rank(), pivots_depth14, len(depth_conditions14))
        == ((78, 15), 5, (10, 11, 12, 13, 14), 1),
        "row-fourteen depth selector on boundary",
    )
    q_raw = depth_conditions14[0]
    scale_q = exact(
        a0 / sp.Poly(q_raw, z, h18, domain=sp.QQ).coeff_monomial(z**2)
    )
    check(not scale_q.free_symbols and scale_q != 0, "depth quadratic normalization")
    q = sp.expand(scale_q * q_raw)
    q_poly = sp.Poly(q, z, h18, domain=sp.QQ)
    check(q_poly.total_degree() == 2, "depth source condition quadratic")
    q_h0 = sp.expand(q.subs(h18, 0))
    q_h0_poly = sp.Poly(q_h0, z, domain=sp.QQ)
    c1 = exact(q_h0.subs(z, 0))
    check(
        (q_h0_poly.coeff_monomial(z**2), q_h0_poly.coeff_monomial(z), c1)
        == (
            a0,
            b0,
            sp.Integer(10771463883409470380030782972892985344),
        ),
        "boundary depth quadratic coefficients",
    )
    bracket_without_c51 = sp.expand(j_boundary.subs(c51, 0))
    depth_minus_bracket = sp.expand(q - bracket_without_c51)
    print("two_parameter_boundary_depth")
    print(f"boundary_depth_rank={matrix_depth14.shape}_rank{matrix_depth14.rank()}_pivots{pivots_depth14}_conditions{len(depth_conditions14)}")
    print("Qzh=" + sp.sstr(q))
    print("Qzh_factor=" + sp.sstr(sp.factor(q)))
    print("Qzh_minus_bracket_noc51=" + sp.sstr(depth_minus_bracket))
    print("Qzh_h_discriminant_factor=" + sp.sstr(sp.factor(sp.discriminant(q, h18))))
    print("Qzh_z_discriminant_factor=" + sp.sstr(sp.factor(sp.discriminant(q, z))))
    q_h2 = q_poly.coeff_monomial(h18**2)
    q_hz = q_poly.coeff_monomial(h18*z)
    q_z2 = q_poly.coeff_monomial(z**2)
    homogeneous = sp.expand(q_h2*h18**2 + q_hz*h18*z + q_z2*z**2)
    expected_homogeneous = sp.Rational(1059352899609375, 4) * (
        801*h18 + 27826*z
    ) * (
        46174707*h18 + 1363633922*z
    )
    check(homogeneous == sp.expand(expected_homogeneous),
          "rational split at infinity")

    u_expr = 801*h18 + 27826*z
    v_expr = 46174707*h18 + 1363633922*z
    k_uv = sp.Rational(1059352899609375, 4)
    alpha_u = sp.Rational(668569947360751846653880320000, 113)
    beta_v = sp.Rational(68455987735926217512960000, 113)
    q_constant = sp.Integer(10771463883409470380030782972892985344)
    check(
        q == sp.expand(k_uv*u_expr*v_expr + alpha_u*u_expr + beta_v*v_expr + q_constant),
        "bilinear normal form",
    )
    shift_u = exact(beta_v/k_uv)
    shift_v = exact(alpha_u/k_uv)
    rectangle_rhs = exact(-(q_constant-alpha_u*beta_v/k_uv)/k_uv)
    check(
        q == sp.expand(k_uv*((u_expr+shift_u)*(v_expr+shift_v)-rectangle_rhs)),
        "translated rectangular hyperbola",
    )

    parameter = sp.symbols("s", nonzero=True)
    u_parameter = exact(parameter-shift_u)
    v_parameter = exact(rectangle_rhs/parameter-shift_v)
    inverse = sp.solve(
        (sp.Eq(u_expr, u_parameter), sp.Eq(v_expr, v_parameter)),
        (h18, z), dict=True,
    )
    check(len(inverse) == 1, "invertible rational packet coordinates")
    h_parameter = exact(inverse[0][h18])
    z_parameter = exact(inverse[0][z])
    check(exact(q.subs({h18: h_parameter, z: z_parameter})) == 0,
          "rational parameterization lies on depth conic")

    constant_gap = exact(depth_minus_bracket)
    fixed_c51 = exact(constant_gap/d_c51)
    check(fixed_c51 == sp.Rational(1087, 135),
          "combined ideal fixes c51")
    check(
        exact(c51_graph.subs({h18: h_parameter, z: z_parameter})-fixed_c51) == 0,
        "parameterized c51 graph",
    )
    rho_parameter = exact(rho_graph.subs({h18: h_parameter, z: z_parameter}))

    pivot_minor = matrix_depth14.extract(row_pivots14, pivots_depth14).det()
    pivot_parameter = exact(pivot_minor.subs({h18: h_parameter, z: z_parameter}))
    check(pivot_parameter == sp.Rational(1, 32),
          "constant terminal pivot on rational conic")
    check(len(theta14)-len(pivots_depth14) == 10,
          "terminal fibre A10")
    depth_after = [exact(value.subs(terminal14)) for value in depth14]
    depth_parameter = [
        exact(value.subs({h18: h_parameter, z: z_parameter}))
        for value in depth_after
    ]
    for index, value in enumerate(depth_parameter):
        check(value == 0, f"rational conic depth coefficient {index}")

    point_parameter = shift_u
    h_point = exact(h_parameter.subs(parameter, point_parameter))
    z_point = exact(z_parameter.subs(parameter, point_parameter))
    rho_point = exact(rho_parameter.subs(parameter, point_parameter))
    check(
        exact(q.subs({h18: h_point, z: z_point})) == 0,
        "explicit rational point",
    )
    old_slice_hostile = exact(q.subs({h18: 0, z: 1}))
    check(old_slice_hostile != 0, "old h0 z1 hostile survives")

    gamma_primitive = primitive(gamma_global)
    check(sp.factor(gamma_primitive) == gamma_primitive,
          "global conic degeneration divisor is QQ-irreducible")
    graph_order = (
        ("c14", resolved13[c14]),
        ("beta11", resolved13[R8.beta11]),
        ("c42", resolved13[c42]),
        ("xi10", resolved13[R8.xi10]),
        ("lambda18", resolved13[lambda18]),
        ("nu18", resolved13[nu18]),
        ("c23", resolved13[c23]),
        ("c70", resolved13[c70]),
        ("kappa20", resolved13[kappa20]),
    )
    graph_data = []
    for name, value in graph_order:
        numerator = sp.fraction(sp.together(value))[0]
        polynomial = sp.Poly(numerator, z, h18, domain=sp.EX)
        graph_data.append(
            (name, polynomial.degree(z), polynomial.degree(h18),
             polynomial.total_degree(), denominator(value))
        )

    print("Qzh_homogeneous_factor=" + sp.sstr(expected_homogeneous))
    print(f"bilinear_coordinates=u:{sp.sstr(u_expr)};v:{sp.sstr(v_expr)}")
    print(f"bilinear_form=k*u*v+alpha*u+beta*v+c;k={k_uv};alpha={alpha_u};beta={beta_v};c={q_constant}")
    print(f"rectangle_shifts=u:{shift_u};v:{shift_v};rhs:{rectangle_rhs}")
    print(f"rational_parameter=s_nonzero;u=s-{shift_u};v={rectangle_rhs}/s-{shift_v}")
    print(f"rational_parameter_h={sp.sstr(h_parameter)}")
    print(f"rational_parameter_z={sp.sstr(z_parameter)}")
    print(f"combined_depth_bracket_c51={fixed_c51};rho_parameter={sp.sstr(rho_parameter)}")
    print(f"terminal_pivot_on_curve={pivot_parameter};terminal_fibre=A10")
    print(f"explicit_rational_point_s={point_parameter};h={h_point};z={z_point};c51={fixed_c51};rho={rho_point}")
    print(f"old_h0_z1_depth_hostile={old_slice_hostile}")
    print("prefix_rank_profile=row9_B10x7_r7_D28x10_r3;row10_B11x7_r7_D36x11_r3;row11_B12x8_r8_D45x12_r4;row12_B13x8_r8_D55x13_r4;row13_B14x9_r9_D66x14_r5")
    print("prefix_source_projection=A6(Phi,eta,alpha11,c51,z,h18);row13_terminal=A9")
    print("prefix_graph_z_h_total_degrees_and_denominators=" + sp.sstr(tuple(graph_data)))
    print(f"J18zh_hash={expression_hash(j18)};gamma_hash={expression_hash(gamma_global)};Qzh_hash={expression_hash(q)}")
    print(f"global_degenerate_divisor=gamma_zero;gamma_primitive={sp.sstr(gamma_primitive)}")
    print("global_parameter_domain=tau_nonzero;when_gamma_nonzero_it_covers_the_full_affine_conic;when_gamma_zero_it_covers_only_V0_with_U_nonzero_and_misses_U0")
    print("boundary_parameter_domain=s_nonzero_is_the_entire_affine_conic;excluded_s0_and_s_infinity_are_the_two_projective_points_at_infinity")
    print("coefficientwise_verification=15_of_15_bracket_and_78_of_78_depth_zero")
    print("field=characteristic_zero;finite_field_used=no")
    print("scope=finite_row14_partial_weight18_source;global_bracket_and_boundary_depth_only;not_full_B2_or_termination_or_chart_entry_or_JC2_or_DC2")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
