#!/usr/bin/env python3
"""Clean-room exact audit for the proposed THM-4399 row-eleven repair.

Only the proved THM-4308 bracket operators and THM-4315 Student/depth
operators are imported.  In particular this file neither imports the
THM-4395 row-ten implementation nor the THM-4399 candidate scout.  It
reconstructs the complete residual-weight-at-most-14 source, redoes the
row-nine and row-ten eliminations, computes the row-eleven compatibility
polynomial, and audits the late ``lambda18*p**3*y**4`` response.

The scope is finite row eleven in the fixed source-normal chart.  The late
single channel is not the complete weight-at-most-18 source family.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_bracket_hasse_rows8_thm4308 as R8  # noqa: E402
import jc2_source_normal_student_stein_row9_thm4315 as R9  # noqa: E402


x, t = R8.x, R8.t
c51, c23 = sp.symbols("c51 c23")
c70, c42, c14 = sp.symbols("c70 c42 c14")
lambda18 = sp.symbols("lambda18")
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
    numerator, _ = sp.fraction(exact(value))
    variables = tuple(sorted(numerator.free_symbols, key=str))
    if not variables:
        return numerator
    result = sp.Poly(numerator, *variables, domain=sp.QQ).primitive()[1].as_expr()
    return sp.expand(result)


def apply(rows: list[sp.Expr], *maps: dict[sp.Symbol, sp.Expr]) -> list[sp.Expr]:
    answer = rows
    for mapping in maps:
        answer = [exact(row.subs(mapping)) for row in answer]
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
    """Solve pivot coordinates and return the residual compatibility ideal."""

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
        if candidate == 0:
            continue
        if all(
            exact(candidate - old) != 0 and exact(candidate + old) != 0
            for old in conditions
        ):
            conditions.append(candidate)
    return substitutions, conditions, raw, matrix, pivots


@lru_cache(maxsize=None)
def depth_data(depth: int, row: int) -> tuple[tuple[tuple[int, int], ...], sp.Matrix, tuple[sp.Matrix, ...]]:
    coordinates, matrix = R9.depth_matrix(depth, row)
    left_null = tuple(matrix.T.nullspace())
    return tuple(coordinates), matrix, left_null


def depth_equations(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
) -> tuple[list[sp.Expr], sp.Matrix, sp.Matrix, int]:
    acoords, amatrix, anull = depth_data(2, row)
    ccoords, cmatrix, cnull = depth_data(3, row)
    avec = sp.Matrix([xcoeff(arows[n], degree) for n, degree in acoords])
    cvec = sp.Matrix([xcoeff(crows[n], degree) for n, degree in ccoords])
    aequations = [exact((left.T * avec)[0]) for left in anull]
    cequations = [exact((left.T * cvec)[0]) for left in cnull]
    return aequations + cequations, amatrix, cmatrix, len(aequations)


def bracket_equations(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    source_row: sp.Expr,
) -> list[sp.Expr]:
    difference = exact(source_row - R8.predicted_G(row, arows, crows))
    check(
        difference == 0 or sp.Poly(difference, x, domain=sp.EX).degree() <= row,
        f"row-{row} bracket degree cap",
    )
    return [xcoeff(difference, degree) for degree in range(row + 1)]


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


def normalize_by_coefficient(
    value: sp.Expr,
    monomial: sp.Expr,
    target: sp.Integer,
    variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    coefficient = sp.Poly(value, *variables, domain=sp.QQ).coeff_monomial(monomial)
    check(coefficient != 0, f"normalizing coefficient for {monomial}")
    return exact(value * target / coefficient)


F_EXPECTED = sp.expand(
    712092531095885035687500 * R8.Phi**4
    + 821011731481169865000000 * R8.Phi**3 * R8.alpha11
    + 139278775876269887812500 * R8.Phi**3 * c51
    - 316295855547474501562500 * R8.Phi**3 * R8.eta
    + 236647708617429900000000 * R8.Phi**2 * R8.alpha11**2
    + 80291186852342287500000 * R8.Phi**2 * R8.alpha11 * c51
    - 43058909665599049687500 * R8.Phi**2 * R8.alpha11 * R8.eta
    + 6810413170511176171875 * R8.Phi**2 * c51**2
    - 30932285940138480468750 * R8.Phi**2 * c51 * R8.eta
    + 3733889920481087640000000 * R8.Phi**2 * c70
    + 445628783011794805546875 * R8.Phi**2 * R8.eta**2
    - 76479730269188998356067968000 * R8.Phi**2
    + 80291186852342287500000 * R8.Phi * R8.alpha11**2 * R8.eta
    + 13620826341022352343750 * R8.Phi * R8.alpha11 * c51 * R8.eta
    + 2152506377265652800000000 * R8.Phi * R8.alpha11 * c70
    + 205715422677291419531250 * R8.Phi * R8.alpha11 * R8.eta**2
    - 44128098153451944185391360000 * R8.Phi * R8.alpha11
    + 21736146783278091456000000 * R8.Phi * c23
    + 365157331857566100000000 * R8.Phi * c51 * c70
    + 40145593426171143750000 * R8.Phi * c51 * R8.eta**2
    - 7301988813385081335237120000 * R8.Phi * c51
    - 829255929072250500000000 * R8.Phi * c70 * R8.eta
    - 91168842770934468750000 * R8.Phi * R8.eta**3
    + 16987951185108221412403200000 * R8.Phi * R8.eta
    + 6810413170511176171875 * R8.alpha11**2 * R8.eta**2
    + 12679418956912220016000000 * R8.alpha11**2
    + 365157331857566100000000 * R8.alpha11 * c70 * R8.eta
    + 40145593426171143750000 * R8.alpha11 * R8.eta**3
    - 7301988813385081335237120000 * R8.alpha11 * R8.eta
    + 25358837913824440032000000 * c51 * R8.eta
    + 4894705859649350400000000 * c70**2
    + 1076253188632826400000000 * c70 * R8.eta**2
    - 200365448846185933658357760000 * c70
    + 59161927154357475000000 * R8.eta**4
    - 22073830342778447233850880000 * R8.eta**2
    + 2077158990491529775184308229636096
)


def main() -> None:
    # Reconstruct the complete residual-weight <=14 face literally.
    indices13 = [(i, j) for j in range(14) for i in range(14) if 2 * i + 3 * j == 13]
    indices14 = [(i, j) for j in range(15) for i in range(15) if 2 * i + 3 * j == 14]
    check(indices13 == [(5, 1), (2, 3)], "complete weight-13 face")
    check(indices14 == [(7, 0), (4, 2), (1, 4)], "complete weight-14 face")
    face13 = sp.expand(c51 * R8.p**5 * R8.y + c23 * R8.p**2 * R8.y**3)
    face14 = sp.expand(c70 * R8.p**7 + c42 * R8.p**4 * R8.y**2 + c14 * R8.p * R8.y**4)
    source = sp.expand(R8.H + face13 + face14)
    grows = {row: tcoeff(source, row) for row in range(4, 12)}
    expected13 = {
        7: c51 * x,
        8: (6 * c51 + c23) * x**3,
        9: (15 * c51 + 5 * c23) * x**5,
        10: (20 * c51 + 10 * c23) * x**7,
        11: (15 * c51 + 10 * c23) * x**9,
    }
    expected14 = {
        7: c70,
        8: (7 * c70 + c42) * x**2,
        9: (21 * c70 + 6 * c42 + c14) * x**4,
        10: (35 * c70 + 15 * c42 + 5 * c14) * x**6,
        11: (35 * c70 + 20 * c42 + 10 * c14) * x**8,
    }
    for row in range(7, 12):
        check(exact(tcoeff(face13, row) - expected13[row]) == 0, f"weight-13 row {row}")
        check(exact(tcoeff(face14, row) - expected14[row]) == 0, f"weight-14 row {row}")

    # Rebuild rows 4--7 and the full row-eight tangent from THM-4308's
    # audited determinant/Student operators.
    arows = [R8.A0, R8.A1, R8.A2, R8.A3]
    crows = [R8.C0, R8.C1, R8.C2, R8.C3]
    bracket_subs: dict[sp.Symbol, sp.Expr] = {}
    response_symbols = {5: R8.upsilon5, 6: R8.U, 7: R8.W, 8: R8.Z}
    for n in range(4, 8):
        abase, cbase = R8.particular_row(n, R8.B_row(n, arows, crows))
        trial_a, trial_c = arows + [abase], crows + [cbase]
        row = n + 1
        difference = exact(grows[row].subs(bracket_subs) - R8.predicted_G(row, trial_a, trial_c))
        moment = R9.primitive_student_row(row)
        obstruction = exact(sum(moment[j] * xcoeff(difference, j) for j in range(row + 1)))
        answers = sp.solve(obstruction, response_symbols[row])
        check(len(answers) == 1, f"row-{row} scalar response")
        bracket_subs[response_symbols[row]] = exact(answers[0])
        target = exact(difference.subs(response_symbols[row], answers[0]))
        theta = R8.tangent_solve(row, target)
        arows.append(sp.expand(abase + theta * sp.diff(R8.A0, x)))
        crows.append(sp.expand(cbase + theta * sp.diff(R8.C0, x)))
        check(
            exact(R8.predicted_G(row, arows, crows) - grows[row].subs(bracket_subs)) == 0,
            f"row-{row} reconstructed",
        )

    arows, crows, theta8 = append_tangent(arows, crows, 8, "audit4399_theta8")
    depth8, ap8, cp8, _ = depth_equations(arows, crows, 8)
    terminal8, conditions8, _, dmatrix8, dpivots8 = eliminate_linear(depth8, theta8)
    check((ap8.shape, ap8.rank()) == ((63, 131), 51), "row-8 P2 universe")
    check((cp8.shape, cp8.rank()) == ((72, 204), 63), "row-8 P3 universe")
    check(dmatrix8.rank() == 2 and dpivots8 == (7, 8), "row-8 terminal pair")
    gate = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * R8.Phi,
    }
    check(len(conditions8) == 3 and all(exact(value.subs(gate)) == 0 for value in conditions8), "old Hasse gate")
    a8 = apply(arows, gate, terminal8, gate)
    c8 = apply(crows, gate, terminal8, gate)

    # Row nine: independently recover the global c14 graph.
    g9 = exact(grows[9].subs(bracket_subs).subs(gate))
    select9, conditions9, _, bmatrix9, bpivots9 = eliminate_linear(
        bracket_equations(a8, c8, 9, g9), theta8[:7]
    )
    check(bmatrix9.shape == (10, 7) and bmatrix9.rank() == 7, "row-9 bracket selector")
    check(bpivots9 == tuple(range(7)) and len(conditions9) == 1, "row-9 one condition")
    old_e9 = (
        613527750 * R8.Phi**2
        - 511211250 * R8.Phi * R8.alpha11
        - 3154140000 * R8.Phi * R8.eta
        - 255605625 * R8.eta**2
        + 6736896000 * R8.xi10
        - 46483785515008
    )
    e9_expected = old_e9 - 3339765000 * c70 - 898128000 * c42 - 216513000 * c14
    e9 = normalize_by_coefficient(
        conditions9[0], R8.Phi**2, sp.Integer(613527750),
        (R8.Phi, R8.alpha11, R8.eta, R8.xi10, c70, c42, c14),
    )
    check(exact(e9 - e9_expected) == 0, "row-9 compatibility equation")
    c14_graph = exact(sp.solve(e9, c14)[0])
    check(sp.denom(sp.together(c14_graph)) == 216513000, "row-9 c14 constant denominator")

    a9base, c9base, theta9 = append_tangent(apply(a8, select9), apply(c8, select9), 9, "audit4399_theta9")
    depth9, ap9, cp9, _ = depth_equations(a9base, c9base, 9)
    terminal9, conditions_d9, _, dmatrix9, dpivots9 = eliminate_linear(depth9, theta9)
    check((ap9.shape, ap9.rank()) == ((75, 160), 59), "row-9 P2 universe")
    check((cp9.shape, cp9.rank()) == ((85, 251), 73), "row-9 P3 universe")
    check(dmatrix9.rank() == 3 and dpivots9 == (7, 8, 9), "row-9 terminal triple")
    check(not conditions_d9, "row-9 depth automatic")
    a9 = apply(a9base, terminal9, {c14: c14_graph})
    c9 = apply(c9base, terminal9, {c14: c14_graph})

    # Row ten: recover the xi bracket graph, the beta depth graph, and the
    # remaining constant-pivot c42 graph without using THM-4395 code.
    g10 = exact(grows[10].subs(bracket_subs).subs(gate).subs(c14, c14_graph))
    select10, conditions10, _, bmatrix10, bpivots10 = eliminate_linear(
        bracket_equations(a9, c9, 10, g10), theta9[:7]
    )
    check(bmatrix10.shape == (11, 7) and bmatrix10.rank() == 7, "row-10 bracket selector")
    check(bpivots10 == tuple(range(7)) and len(conditions10) == 2, "row-10 two conditions")
    xi_conditions = [
        value
        for value in conditions10
        if sp.diff(value, R8.xi10) != 0
        and all(sp.diff(value, variable) == 0 for variable in (R8.alpha11, R8.beta11, c51))
    ]
    check(len(xi_conditions) == 1, "unique pure row-10 xi condition")
    xi_raw = xi_conditions[0]
    e10_xi_expected = (
        13365000 * R8.Phi**2
        + 15035625 * R8.Phi * R8.eta
        + 6014250 * c42
        + 50787000 * c70
        + 57672000 * R8.xi10
        - 964604821504
    )
    e10_xi = normalize_by_coefficient(
        xi_raw, R8.xi10, sp.Integer(57672000),
        (R8.Phi, R8.eta, R8.xi10, c42, c70),
    )
    check(exact(e10_xi - e10_xi_expected) == 0, "row-10 xi equation")
    xi_graph = exact(sp.solve(e10_xi, R8.xi10)[0])
    remaining10 = [exact(value.subs(R8.xi10, xi_graph)) for value in conditions10 if value != xi_raw]
    remaining10 = [value for value in remaining10 if value != 0]
    check(len(remaining10) == 1, "one row-10 equation after xi")
    b10_expected = (
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
    b10 = normalize_by_coefficient(
        remaining10[0], c42, sp.Integer(-89131914000),
        (R8.Phi, R8.eta, R8.alpha11, R8.beta11, c51, c42, c70),
    )
    check(exact(b10 - b10_expected) == 0, "row-10 second bracket equation")

    a10base, c10base, theta10 = append_tangent(
        apply(a9, select10, {R8.xi10: xi_graph}),
        apply(c9, select10, {R8.xi10: xi_graph}),
        10,
        "audit4399_theta10",
    )
    depth10, ap10, cp10, _ = depth_equations(a10base, c10base, 10)
    terminal10, conditions_d10, _, dmatrix10, dpivots10 = eliminate_linear(depth10, theta10)
    check((ap10.shape, ap10.rank()) == ((88, 193), 68), "row-10 P2 universe")
    check((cp10.shape, cp10.rank()) == ((99, 304), 83), "row-10 P3 universe")
    check(dmatrix10.rank() == 3 and dpivots10 == (8, 9, 10), "row-10 terminal triple")
    check(len(conditions_d10) == 1, "row-10 unique depth condition")
    beta_condition = normalize_by_coefficient(
        conditions_d10[0], R8.beta11, sp.Integer(15),
        (R8.Phi, R8.eta, R8.beta11),
    )
    check(exact(beta_condition - (-91 * R8.Phi + 15 * R8.beta11 + 18 * R8.eta)) == 0, "row-10 beta depth graph")
    beta_graph = exact(sp.solve(beta_condition, R8.beta11)[0])
    b10_depth = exact(b10.subs(R8.beta11, beta_graph))
    b10_depth_expected = (
        14643240750 * R8.Phi**2
        + 122625090000 * R8.Phi * R8.alpha11
        + 20802470625 * R8.Phi * c51
        - 270071263125 * R8.Phi * R8.eta
        + 20802470625 * R8.alpha11 * R8.eta
        - 89131914000 * c42
        - 194981256000 * c70
        + 61312545000 * R8.eta**2
        + 2707389207937024
    )
    check(exact(b10_depth - b10_depth_expected) == 0, "row-10 beta-substituted equation")
    c42_graph = exact(sp.solve(b10_depth, c42)[0])
    check(sp.denom(sp.together(c42_graph)) == 89131914000, "row-10 c42 constant denominator")

    # Finalize all old graphs.  The sequential substitutions matter because
    # c14 depends on xi and xi depends on c42.
    a10 = apply(a10base, terminal10, {R8.beta11: beta_graph}, {c42: c42_graph})
    c10 = apply(c10base, terminal10, {R8.beta11: beta_graph}, {c42: c42_graph})
    g11 = exact(
        grows[11]
        .subs(bracket_subs)
        .subs(gate)
        .subs(c14, c14_graph)
        .subs(R8.xi10, xi_graph)
        .subs(R8.beta11, beta_graph)
        .subs(c42, c42_graph)
    )

    # Row eleven: compute the primitive compatibility polynomial F.
    equations11 = bracket_equations(a10, c10, 11, g11)
    select11, conditions11, raw11, bmatrix11, bpivots11 = eliminate_linear(equations11, theta10[:8])
    check(bmatrix11.shape == (12, 8) and bmatrix11.rank() == 8, "row-11 bracket selector")
    check(bpivots11 == tuple(range(8)) and len(conditions11) == 1, "row-11 principal compatibility")
    fvars = (R8.Phi, R8.eta, R8.alpha11, c51, c23, c70)
    F = normalize_by_coefficient(
        conditions11[0], R8.Phi**4, sp.Integer(712092531095885035687500), fvars
    )
    check(exact(F - F_EXPECTED) == 0, "row-11 F exact coefficient audit")
    fpoly = sp.Poly(F, *fvars, domain=sp.QQ)
    check(len(fpoly.terms()) == 36 and fpoly.total_degree() == 4, "F term count and total degree")
    check(tuple(fpoly.degree(variable) for variable in fvars) == (4, 4, 2, 2, 1, 2), "F multidegree")
    check(fpoly.is_irreducible, "F irreducible over QQ")

    # Base projected depth at row eleven is automatic and leaves A^8.
    a11base, c11base, theta11 = append_tangent(apply(a10, select11), apply(c10, select11), 11, "audit4399_theta11")
    depth11, ap11, cp11, _ = depth_equations(a11base, c11base, 11)
    terminal11, conditions_d11, raw_d11, dmatrix11, dpivots11 = eliminate_linear(depth11, theta11)
    check((ap11.shape, ap11.rank()) == ((102, 228), 77), "row-11 P2 universe")
    check((cp11.shape, cp11.rank()) == ((114, 361), 94), "row-11 P3 universe")
    check(len(depth11) == 45, "row-11 depth equation count")
    check(dmatrix11.rank() == 4 and dpivots11 == (8, 9, 10, 11), "row-11 terminal quartet")
    check(not conditions_d11 and 12 - dmatrix11.rank() == 8, "row-11 depth automatic A8")
    check(all(exact(value.subs(terminal11)) == 0 for value in raw_d11), "row-11 all base depth equations")

    # Every t-valuation-eleven monomial has leading row x^b.  Recompute the
    # quotient response simultaneously instead of assuming the Student row.
    mus = tuple(sp.symbols("late11_0:6"))
    late_row = sum(mus[b] * x**b for b in range(6))
    late_equations = bracket_equations(a10, c10, 11, g11 + late_row)
    late_select, late_conditions, _, late_matrix, late_pivots = eliminate_linear(late_equations, theta10[:8])
    check(late_matrix == bmatrix11 and late_pivots == bpivots11, "late row keeps selector matrix")
    check(len(late_conditions) == 1, "late row has one compatibility equation")
    late_normalized = normalize_by_coefficient(
        late_conditions[0], R8.Phi**4, sp.Integer(712092531095885035687500), fvars + mus
    )
    base_response = sp.Integer(25358837913824440032000000)
    moment11 = R9.primitive_student_row(11)
    response_table: list[tuple[int, int, int, int, sp.Rational, sp.Expr]] = []
    response_linear = sp.Integer(0)
    for b in range(6):
        i = 11 - 2 * b
        weight = 2 * i + 3 * b
        ratio = sp.Rational(moment11[b], moment11[0])
        response = exact(sp.diff(late_normalized, mus[b]))
        check(response == base_response * ratio, f"valuation-11 response b={b}")
        response_linear += response * mus[b]
        response_table.append((b, i, b, weight, ratio, response))
    check(exact(late_normalized - F - response_linear) == 0, "full valuation-11 response equation")

    correction = sp.expand(lambda18 * R8.p**3 * R8.y**4)
    check(all(tcoeff(correction, row) == 0 for row in range(11)), "late correction leaves rows through ten")
    check(exact(tcoeff(correction, 11) - lambda18 * x**4) == 0, "late correction row-eleven leading term")
    lambda_response = response_table[4][-1]
    check(lambda_response == 6864046352614134144000000, "lambda18 constant response")
    lambda_graph = exact(-F / lambda_response)

    corrected_equations = bracket_equations(a10, c10, 11, g11 + lambda18 * x**4)
    corrected_select, corrected_conditions, _, corrected_matrix, corrected_pivots = eliminate_linear(
        corrected_equations, theta10[:8]
    )
    check(corrected_matrix == bmatrix11 and corrected_pivots == bpivots11, "corrected selector matrix")
    check(len(corrected_conditions) == 1, "corrected one equation")
    corrected_condition = normalize_by_coefficient(
        corrected_conditions[0], R8.Phi**4, sp.Integer(712092531095885035687500), fvars + (lambda18,)
    )
    check(exact(corrected_condition - F - lambda_response * lambda18) == 0, "corrected row-11 equation")
    check(
        all(exact(value.subs(corrected_select).subs(lambda18, lambda_graph)) == 0 for value in corrected_equations),
        "all twelve corrected bracket coefficients",
    )

    corrected_a11base, corrected_c11base, corrected_theta11 = append_tangent(
        apply(a10, corrected_select, {lambda18: lambda_graph}),
        apply(c10, corrected_select, {lambda18: lambda_graph}),
        11,
        "audit4399_corrected_theta11",
    )
    corrected_depth, _, _, _ = depth_equations(corrected_a11base, corrected_c11base, 11)
    corrected_terminal, corrected_depth_conditions, corrected_raw_depth, corrected_depth_matrix, corrected_depth_pivots = eliminate_linear(
        corrected_depth, corrected_theta11
    )
    check(corrected_depth_matrix.rank() == 4 and corrected_depth_pivots == (8, 9, 10, 11), "corrected depth quartet")
    check(not corrected_depth_conditions and 12 - corrected_depth_matrix.rank() == 8, "corrected depth A8")
    check(
        all(exact(value.subs(corrected_terminal)) == 0 for value in corrected_raw_depth),
        "all forty-five corrected depth equations",
    )

    # Boundary and hostile controls.  All graphs use constant denominators;
    # in particular Phi=0 and eta=0 are retained.
    graph_denominators = [
        sp.denom(sp.together(c14_graph)),
        sp.denom(sp.together(xi_graph)),
        sp.denom(sp.together(beta_graph)),
        sp.denom(sp.together(c42_graph)),
        sp.denom(sp.together(lambda_graph)),
    ]
    check(all(not denominator.free_symbols for denominator in graph_denominators), "all graph denominators constant")
    dense = {R8.Phi: 1, R8.eta: 2, R8.alpha11: 3, c51: 5, c23: 7, c70: 11}
    dense_f = F.subs(dense)
    check(dense_f == 2074612777259675413868913342459971, "dense hostile F")
    boundary = {R8.Phi: 0, R8.eta: 0, R8.alpha11: 1, c51: 2, c23: 3, c70: 4}
    boundary_f = F.subs(boundary)
    check(boundary_f == 2076357619690857742751501214596096, "Phi=eta=0 boundary F")
    check(exact((F + lambda_response * lambda18).subs(dense).subs(lambda18, lambda_graph.subs(dense))) == 0, "dense corrected equation")
    check(exact((F + lambda_response * lambda18).subs(boundary).subs(lambda18, lambda_graph.subs(boundary))) == 0, "boundary corrected equation")
    check(response_table[5][-1] == 0, "weight-17 odd hostile has zero response")

    print("THM-4399 row-eleven late weight-eighteen response independent audit")
    print("imports=audited_THM4308_R8_and_THM4315_R9;THM4395_primary=no;THM4399_scout=no")
    print("field=characteristic_zero;universe_base=complete_fixed_residual_weight_at_most_14")
    print("row9_graph=c14;denominator=216513000;localization=none")
    print("row10_graphs=xi10,beta11,c42;denominators=57672000,15,89131914000;localization=none")
    print("row11_bracket=shape(12,8),rank8,pivots(0..7),primitive_conditions1")
    print(f"F_terms={len(fpoly.terms())};F_total_degree={fpoly.total_degree()};F_multidegree={tuple(fpoly.degree(v) for v in fvars)};F_QQ_irreducible=yes")
    print(f"F={sp.sstr(F)}")
    print("base_row11_depth=P2(102x228,rank77),P3(114x361,rank94),joint(45x12,rank4,pivots8,9,10,11),conditions0,fibre=A8")
    print("valuation11_table_columns=b,monomial,weight,leading_x_power,Student_ratio,quotient_response")
    for b, i, j, weight, ratio, response in response_table:
        print(f"valuation11_row=({b}, p^{i}*y^{j}, {weight}, {b}, {ratio}, {response})")
    print("minimum_transverse_late_channel=lambda18*p^3*y^4;residual_weight=18;t_valuation=11")
    print(f"corrected_row11_equation=F+{lambda_response}*lambda18")
    print(f"lambda18_graph=-F/{lambda_response};parameter_localization=none")
    print("corrected_literal_bracket_coefficients=12/12_zero")
    print("corrected_row11_depth=45/45_zero;terminal_fibre=A8")
    print("weight17_hostile=mu17*p*y^5;Student_ratio=0;row11_equation_remains_F")
    print(f"dense_control_F={dense_f};lambda18={lambda_graph.subs(dense)};corrected_bracket_depth=PASS")
    print(f"Phi_eta_zero_control_F={boundary_f};lambda18={lambda_graph.subs(boundary)};corrected_bracket_depth=PASS")
    print("scope=finite_row11_fixed_chart;late_single_channel_not_complete_weight18;entry_and_all_row_lift_open")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
