#!/usr/bin/env python3
"""Independent exact audit for THM-4390's complete weight-14 row-nine face.

This audit imports the proved THM-4308/4315 row operators, but no THM-4390
primary certificate.  It reconstructs the literal source face, its integral
Pascal sidecar, the bracket rows through row nine, and the projected P_2/P_3
depth test at row nine.
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
c51, c23 = sp.symbols("c51 c23")
c70, c42, c14 = sp.symbols("c70 c42 c14")
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
    return sp.Poly(numerator, *variables, domain=sp.QQ).primitive()[1].as_expr()


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
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    pivots = tuple(matrix.rref()[1])
    row_pivots = tuple(matrix.T.rref()[1])
    if pivots:
        square = matrix.extract(row_pivots, pivots)
        check(square.rows == square.cols and square.det() != 0, "pivot minor")
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
) -> tuple[list[sp.Expr], sp.Matrix, sp.Matrix]:
    acoords, amatrix = R9.depth_matrix(2, row)
    ccoords, cmatrix = R9.depth_matrix(3, row)
    avec = sp.Matrix([xcoeff(arows[n], degree) for n, degree in acoords])
    cvec = sp.Matrix([xcoeff(crows[n], degree) for n, degree in ccoords])
    equations = [sp.expand((left.T * avec)[0]) for left in amatrix.T.nullspace()]
    equations += [sp.expand((left.T * cvec)[0]) for left in cmatrix.T.nullspace()]
    return equations, amatrix, cmatrix


def bracket_equations(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    source: sp.Expr,
) -> list[sp.Expr]:
    difference = exact(source - R8.predicted_G(row, arows, crows))
    check(
        difference == 0 or sp.Poly(difference, x, domain=sp.EX).degree() <= row,
        f"row-{row} degree cap",
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


def main() -> None:
    # Complete Diophantine faces at weights 13 and 14, with their unimodular
    # THM-4298 coordinates.  Together with R8.H this is the complete declared
    # residual-weight-at-most-14 source family.
    indices13 = [
        (i, j)
        for j in range(14)
        for i in range(14)
        if 2 * i + 3 * j == 13
    ]
    check(indices13 == [(5, 1), (2, 3)], "complete weight-13 face")
    indices = [
        (i, j)
        for j in range(15)
        for i in range(15)
        if 2 * i + 3 * j == 14
    ]
    check(indices == [(7, 0), (4, 2), (1, 4)], "complete weight-14 face")
    pascal = sp.Matrix([[1, 0, 0], [7, 1, 0], [21, 6, 1]])
    inverse = sp.Matrix([[1, 0, 0], [-7, 1, 0], [21, -6, 1]])
    check(pascal.det() == 1 and pascal * inverse == sp.eye(3), "integral Pascal inverse")
    pascal13 = sp.Matrix([[1, 0], [6, 1]])
    inverse13 = sp.Matrix([[1, 0], [-6, 1]])
    check(pascal13.det() == 1 and pascal13 * inverse13 == sp.eye(2), "weight-13 Pascal inverse")
    observer = (sp.Integer(1), sp.Rational(2, 5), sp.Rational(36, 85))
    for (row, degree), expected in zip(((7, 0), (8, 2), (9, 4)), observer):
        moment = R9.primitive_student_row(row)
        check(sp.Rational(moment[degree], moment[0]) == expected, f"row-{row} Student diagonal")
    check(sp.prod(observer) != 0, "even face observer is lossless")
    for row, degree in ((7, 1), (8, 3), (9, 5)):
        moment = R9.primitive_student_row(row)
        check(moment[degree] == 0, f"row-{row} odd weight-13 observer vanishes")
    face13 = sp.expand(c51 * R8.p**5 * R8.y + c23 * R8.p**2 * R8.y**3)
    face14 = sp.expand(c70 * R8.p**7 + c42 * R8.p**4 * R8.y**2 + c14 * R8.p * R8.y**4)
    source = sp.expand(R8.H + face13 + face14)
    grows = {row: tcoeff(source, row) for row in range(4, 15)}
    expected13 = {
        7: c51 * x,
        8: (6 * c51 + c23) * x**3,
        9: (15 * c51 + 5 * c23) * x**5,
        10: (20 * c51 + 10 * c23) * x**7,
        11: (15 * c51 + 10 * c23) * x**9,
        12: (6 * c51 + 5 * c23) * x**11,
        13: (c51 + c23) * x**13,
        14: sp.Integer(0),
    }
    expected14 = {
        7: c70,
        8: (7 * c70 + c42) * x**2,
        9: (21 * c70 + 6 * c42 + c14) * x**4,
        10: (35 * c70 + 15 * c42 + 5 * c14) * x**6,
        11: (35 * c70 + 20 * c42 + 10 * c14) * x**8,
        12: (21 * c70 + 15 * c42 + 10 * c14) * x**10,
        13: (7 * c70 + 6 * c42 + 5 * c14) * x**12,
        14: (c70 + c42 + c14) * x**14,
    }
    check(
        all(exact(tcoeff(face13, row) - delta) == 0 for row, delta in expected13.items()),
        "literal weight-13 row transform",
    )
    check(
        all(exact(tcoeff(face14, row) - delta) == 0 for row, delta in expected14.items()),
        "literal weight-14 row transform",
    )
    check(
        all(
            exact(grows[row] - tcoeff(R8.H, row) - expected13[row] - expected14[row]) == 0
            for row in range(7, 15)
        ),
        "complete weight-at-most-14 source rows",
    )

    # Rebuild bracket rows using the literal new source.  The c70/c42 face
    # directions are allowed to alter selected row lifts; only the source
    # compatibility ideal is compared with the old gate.
    arows = [R8.A0, R8.A1, R8.A2, R8.A3]
    crows = [R8.C0, R8.C1, R8.C2, R8.C3]
    bracket_subs: dict[sp.Symbol, sp.Expr] = {}
    response_symbols = {5: R8.upsilon5, 6: R8.U, 7: R8.W, 8: R8.Z}
    for n in range(4, 8):
        abase, cbase = R8.particular_row(n, R8.B_row(n, arows, crows))
        trial_a = arows + [abase]
        trial_c = crows + [cbase]
        m = n + 1
        difference = sp.expand(
            grows[m].subs(bracket_subs) - R8.predicted_G(m, trial_a, trial_c)
        )
        moment = R9.primitive_student_row(m)
        obstruction = exact(
            sum(moment[j] * xcoeff(difference, j) for j in range(m + 1))
        )
        answers = sp.solve(obstruction, response_symbols[m])
        check(len(answers) == 1, f"row-{m} response")
        bracket_subs[response_symbols[m]] = exact(answers[0])
        target = sp.expand(difference.subs(response_symbols[m], answers[0]))
        theta = R8.tangent_solve(m, target)
        arows.append(sp.expand(abase + theta * sp.diff(R8.A0, x)))
        crows.append(sp.expand(cbase + theta * sp.diff(R8.C0, x)))
        check(
            exact(R8.predicted_G(m, arows, crows) - grows[m].subs(bracket_subs)) == 0,
            f"row-{m} bracket reconstruction",
        )
    check(sp.diff(bracket_subs[R8.W], c70) == -sp.Rational(13, 6), "row-7 c70 response")
    check(sp.diff(bracket_subs[R8.Z], c70) == -sp.Rational(494, 81), "row-8 c70 response")
    check(sp.diff(bracket_subs[R8.Z], c42) == -sp.Rational(13, 18), "row-8 c42 response")
    check(sp.diff(bracket_subs[R8.W], c51) == 0, "row-7 odd channel has zero scalar response")
    check(
        sp.diff(bracket_subs[R8.Z], c51) == sp.diff(bracket_subs[R8.Z], c23) == 0,
        "row-8 odd channels have zero scalar response",
    )

    theta8 = tuple(sp.symbols("w14_audit_theta8_0:9"))
    tangent8 = sum(theta8[j] * x**j for j in range(9))
    abase8, cbase8 = R8.particular_row(8, R8.B_row(8, arows, crows))
    arows.append(sp.expand(abase8 + tangent8 * sp.diff(R8.A0, x)))
    crows.append(sp.expand(cbase8 + tangent8 * sp.diff(R8.C0, x)))

    depth8, ap8, cp8 = depth_equations(arows, crows, 8)
    terminal8, conditions8, _, matrix8, pivots8 = eliminate_linear(depth8, theta8)
    check((ap8.shape, ap8.rank()) == ((63, 131), 51), "row-8 P2 universe")
    check((cp8.shape, cp8.rank()) == ((72, 204), 63), "row-8 P3 universe")
    check(matrix8.rank() == 2 and pivots8 == (7, 8), "row-8 terminal pivots")
    old_gate = [
        896 - 15 * R8.Delta,
        3 * R8.Phi + 2 * R8.zeta3,
        3030 * R8.Delta + 225 * R8.Theta - 182528,
    ]
    check(
        len(conditions8) == 3
        and all(
            any(
                sp.Poly(value, *sorted(value.free_symbols, key=str), domain=sp.QQ).monic()
                == sp.Poly(expected, *sorted(expected.free_symbols, key=str), domain=sp.QQ).monic()
                for expected in old_gate
            )
            for value in conditions8
        ),
        "row-8 Hasse ideal unchanged",
    )
    gate = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * R8.Phi,
    }
    check(all(exact(value.subs(gate)) == 0 for value in conditions8), "old Hasse gate solves row 8")

    a8 = apply(arows, gate, terminal8, gate)
    c8 = apply(crows, gate, terminal8, gate)
    g9 = exact(grows[9].subs(bracket_subs).subs(gate))
    select9, conditions9, raw9, matrix9, pivots9 = eliminate_linear(
        bracket_equations(a8, c8, 9, g9),
        theta8[:7],
    )
    check(matrix9.rank() == 7 and pivots9 == tuple(range(7)), "row-9 affine-seven selection")
    check(len(conditions9) == 1, "unique row-9 source condition")
    old_e9 = (
        613527750 * R8.Phi**2
        - 511211250 * R8.Phi * R8.alpha11
        - 3154140000 * R8.Phi * R8.eta
        - 255605625 * R8.eta**2
        + 6736896000 * R8.xi10
        - 46483785515008
    )
    face_linear = 3339765000 * c70 + 898128000 * c42 + 216513000 * c14
    condition = conditions9[0]
    coefficient_old = sp.Poly(condition, R8.Phi).coeff_monomial(R8.Phi**2)
    normalized_condition = exact(condition * sp.Rational(613527750, coefficient_old))
    check(exact(normalized_condition - old_e9 + face_linear) == 0, "row-9 sign and coefficients")
    check(
        sp.diff(normalized_condition, c51) == sp.diff(normalized_condition, c23) == 0,
        "weight-13 channels vanish from row-9 condition",
    )
    check(sp.diff(normalized_condition, c14) == -216513000, "constant c14 pivot coefficient")

    # Row-nine projected depth is tested after the selected bracket lift but
    # before solving the source condition.  No residual condition remains.
    a8_selected = apply(a8, select9)
    c8_selected = apply(c8, select9)
    a9_trial, c9_trial, theta9 = append_tangent(
        a8_selected, c8_selected, 9, "w14_audit_theta9"
    )
    depth9, ap9, cp9 = depth_equations(a9_trial, c9_trial, 9)
    _, depth9_conditions, _, matrix_depth9, pivots_depth9 = eliminate_linear(depth9, theta9)
    check((ap9.shape, ap9.rank()) == ((75, 160), 59), "row-9 P2 universe")
    check((cp9.shape, cp9.rank()) == ((85, 251), 73), "row-9 P3 universe")
    a_depth_count = len(ap9.T.nullspace())
    amat9, _ = sp.linear_eq_to_matrix(depth9[:a_depth_count], theta9)
    cmat9, _ = sp.linear_eq_to_matrix(depth9[a_depth_count:], theta9)
    check(amat9.rank() == 3 and amat9.rref()[1] == (7, 8, 9), "row-9 P2 terminal")
    check(cmat9.rank() == 2 and cmat9.rref()[1] == (8, 9), "row-9 P3 terminal")
    check(matrix_depth9.rank() == 3 and pivots_depth9 == (7, 8, 9), "row-9 depth terminal")
    check(not depth9_conditions, "row-9 projected depth automatic")

    c14_answers = sp.solve(normalized_condition, c14)
    check(len(c14_answers) == 1, "c14 global solve")
    c14_graph = exact(c14_answers[0])
    denominator = sp.denom(sp.together(c14_graph))
    check(denominator == 216513000, "c14 solve has constant denominator 216513000")
    check(
        sp.diff(c14_graph, c51) == sp.diff(c14_graph, c23) == 0,
        "c14 graph independent of weight-13 channels",
    )
    check(exact(normalized_condition.subs(c14, c14_graph)) == 0, "c14 solve sufficiency")

    print("THM-4390 weight-14 row-nine face absorption independent audit")
    print("imports=audited_THM4308_R8_and_THM4315_R9; THM4390_primary_import=no")
    print("universe=complete_fixed_residual_weight_at_most_14")
    print("weight13_face=(p^5*y,p^2*y^3) channels=(c51,6*c51+c23) complete=yes")
    print("weight13_odd_Student_diagonal=(0,0,0) row9_derivatives=(0,0)")
    print("weight14_face=(p^7,p^4*y^2,p*y^4) complete=yes")
    print("pascal_channels=(c70,7*c70+c42,21*c70+6*c42+c14) determinant=1")
    print("pascal_inverse=(h0,h1-7*h0,h2-6*h1+21*h0)")
    print("Student_diagonal=(1,2/5,36/85) lossless=yes")
    print("rows_through_8=old_Hasse_gate_exactly")
    print("row9_condition=E9-(3339765000*c70+898128000*c42+216513000*c14)")
    print("row9_projected_P2_P3_depth=automatic terminal_rank=3 pivots=(7,8,9)")
    print("c14_absorption=global_affine_graph pivot_coefficient=-216513000 boundary_loss=none")
    print("field=characteristic_zero finite_field_used=no")
    print(f"checks={CHECKS} result=PASS")


if __name__ == "__main__":
    main()
