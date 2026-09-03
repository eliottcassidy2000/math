#!/usr/bin/env python3
"""Primary exact certificate for THM-4390.

Adjoin the complete residual-weight-thirteen and -fourteen faces

    c51*p**5*y + c23*p**2*y**3
    + c70*p**7 + c42*p**4*y**2 + c14*p*y**4

to the fixed THM-4308 source-normal chart.  This certificate reconstructs
the literal source rows, bracket recursion, and projected P_2/P_3 depth
through row nine.  It imports only the audited THM-4308/4315 row operators.
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
H14 = sp.expand(
    R8.H
    + c51 * R8.p**5 * R8.y
    + c23 * R8.p**2 * R8.y**3
    + c70 * R8.p**7
    + c42 * R8.p**4 * R8.y**2
    + c14 * R8.p * R8.y**4
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


def primitive(value: sp.Expr) -> sp.Expr:
    numerator, _ = sp.fraction(exact(value))
    variables = tuple(sorted(numerator.free_symbols, key=str))
    if not variables:
        return sp.signsimp(numerator)
    return sp.Poly(numerator, *variables, domain=sp.QQ).primitive()[1].as_expr()


def proportional(left: sp.Expr, right: sp.Expr) -> bool:
    variables = tuple(sorted(left.free_symbols | right.free_symbols, key=str))
    left_poly = sp.Poly(sp.fraction(exact(left))[0], *variables, domain=sp.QQ)
    right_poly = sp.Poly(sp.fraction(exact(right))[0], *variables, domain=sp.QQ)
    return left_poly.monic() == right_poly.monic()


def impose(rows: list[sp.Expr], *maps: dict[sp.Symbol, sp.Expr]) -> list[sp.Expr]:
    answer = rows
    for mapping in maps:
        answer = [exact(row.subs(mapping)) for row in answer]
    return answer


def projected_depth_residuals(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
) -> tuple[list[sp.Expr], sp.Matrix, sp.Matrix]:
    acoords, amatrix = R9.depth_matrix(2, row)
    ccoords, cmatrix = R9.depth_matrix(3, row)
    avec = sp.Matrix([xcoeff(arows[n], degree) for n, degree in acoords])
    cvec = sp.Matrix([xcoeff(crows[n], degree) for n, degree in ccoords])
    residuals = [sp.expand((left.T * avec)[0]) for left in amatrix.T.nullspace()]
    residuals += [sp.expand((left.T * cvec)[0]) for left in cmatrix.T.nullspace()]
    return residuals, amatrix, cmatrix


def eliminate_linear_fibre(
    residuals: list[sp.Expr],
    fibre_symbols: tuple[sp.Symbol, ...],
) -> tuple[dict[sp.Symbol, sp.Expr], list[sp.Expr], sp.Matrix, tuple[int, ...]]:
    matrix, rhs = sp.linear_eq_to_matrix(residuals, fibre_symbols)
    pivots = tuple(matrix.rref()[1])
    row_pivots = tuple(matrix.T.rref()[1])
    pivot_matrix = matrix.extract(row_pivots, pivots)
    check(pivot_matrix.rows == pivot_matrix.cols, "linear pivot block is square")
    check(pivot_matrix.det() != 0, "linear pivot block is invertible")
    solution_vector = pivot_matrix.inv() * rhs.extract(row_pivots, (0,))
    substitutions = {
        fibre_symbols[column]: exact(solution_vector[index])
        for index, column in enumerate(pivots)
    }
    reduced = [exact(value.subs(substitutions)) for value in residuals]
    conditions: list[sp.Expr] = []
    for value in reduced:
        if value == 0:
            continue
        candidate = primitive(value)
        if all(not proportional(candidate, old) for old in conditions):
            conditions.append(candidate)
    return substitutions, conditions, matrix, pivots


def append_tangent(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    prefix: str,
) -> tuple[list[sp.Expr], list[sp.Expr], tuple[sp.Symbol, ...]]:
    bvalue = R8.B_row(row, arows, crows)
    abase, cbase = R8.particular_row(row, bvalue)
    symbols = tuple(sp.symbols(f"{prefix}_0:{row + 1}"))
    theta = sum(symbols[j] * x**j for j in range(row + 1))
    anew = sp.expand(abase + theta * sp.diff(R8.A0, x))
    cnew = sp.expand(cbase + theta * sp.diff(R8.C0, x))
    return arows + [anew], crows + [cnew], symbols


def solve_bracket_fibre(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    source: sp.Expr,
    fibre_symbols: tuple[sp.Symbol, ...],
) -> tuple[dict[sp.Symbol, sp.Expr], list[sp.Expr], sp.Matrix, tuple[int, ...]]:
    difference = exact(source - R8.predicted_G(row, arows, crows))
    equations = [xcoeff(difference, degree) for degree in range(row + 1)]
    check(exact(difference - sum(equations[d] * x**d for d in range(row + 1))) == 0,
          f"row {row} coefficient universe is exhaustive")
    return eliminate_linear_fibre(equations, fibre_symbols)


def build_through_row_eight() -> tuple[
    list[sp.Expr],
    list[sp.Expr],
    dict[sp.Symbol, sp.Expr],
    dict[int, sp.Expr],
]:
    grows = {row: tcoeff(H14, row) for row in range(4, 15)}
    arows = [R8.A0, R8.A1, R8.A2, R8.A3]
    crows = [R8.C0, R8.C1, R8.C2, R8.C3]
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    response_symbols = {5: R8.upsilon5, 6: R8.U, 7: R8.W, 8: R8.Z}
    check(exact(R8.predicted_G(4, arows, crows) - grows[4]) == 0,
          "inherited row four")

    for n in range(4, 8):
        abase, cbase = R8.particular_row(n, R8.B_row(n, arows, crows))
        target_row = n + 1
        difference = sp.expand(
            grows[target_row].subs(substitutions)
            - R8.predicted_G(target_row, arows + [abase], crows + [cbase])
        )
        moment = R9.primitive_student_row(target_row)
        obstruction = exact(
            sum(moment[j] * xcoeff(difference, j) for j in range(target_row + 1))
        )
        response = response_symbols[target_row]
        answers = sp.solve(obstruction, response)
        check(len(answers) == 1, f"row {target_row} response is unique")
        substitutions[response] = exact(answers[0])
        target = sp.expand(difference.subs(response, substitutions[response]))
        theta = R8.tangent_solve(target_row, target)
        arows.append(sp.expand(abase + theta * sp.diff(R8.A0, x)))
        crows.append(sp.expand(cbase + theta * sp.diff(R8.C0, x)))
        check(
            exact(R8.predicted_G(target_row, arows, crows)
                  - grows[target_row].subs(substitutions)) == 0,
            f"row {target_row} bracket reconstruction",
        )

    _, _, theta8 = append_tangent(arows, crows, 8, "w14_theta8")
    b8 = R8.B_row(8, arows, crows)
    a8base, c8base = R8.particular_row(8, b8)
    tangent8 = sum(theta8[j] * x**j for j in range(9))
    arows.append(sp.expand(a8base + tangent8 * sp.diff(R8.A0, x)))
    crows.append(sp.expand(c8base + tangent8 * sp.diff(R8.C0, x)))
    return arows, crows, substitutions, grows


def main() -> None:
    pairs13 = [(a, b) for a in range(7) for b in range(5) if 2 * a + 3 * b == 13]
    pairs = [(a, b) for a in range(8) for b in range(5) if 2 * a + 3 * b == 14]
    check(pairs13 == [(2, 3), (5, 1)], "complete weight-thirteen face")
    check(pairs == [(1, 4), (4, 2), (7, 0)], "complete weight-fourteen face")

    expected13 = {
        row: sp.expand(
            (
                sp.binomial(6, row - 7) * c51
                + sp.binomial(5, row - 8) * c23
            )
            * x ** (2 * row - 13)
        )
        for row in range(7, 14)
    }
    expected14 = {
        row: sp.expand(
            (
                sp.binomial(7, row - 7) * c70
                + sp.binomial(6, row - 8) * c42
                + sp.binomial(5, row - 9) * c14
            )
            * x ** (2 * row - 14)
        )
        for row in range(7, 15)
    }
    for row in range(4, 15):
        delta = exact(tcoeff(H14 - R8.H, row))
        expected = (expected13.get(row, 0) + expected14.get(row, 0))
        check(delta == expected, f"source delta row {row}")

    odd_h7 = c51
    odd_h8 = 6 * c51 + c23
    check(exact(c51 - odd_h7) == 0, "odd channel inverse c51")
    check(exact(c23 - (odd_h8 - 6 * odd_h7)) == 0, "odd channel inverse c23")

    h7 = c70
    h8 = 7 * c70 + c42
    h9 = 21 * c70 + 6 * c42 + c14
    check(exact(c70 - h7) == 0, "channel inverse c70")
    check(exact(c42 - (h8 - 7 * h7)) == 0, "channel inverse c42")
    check(exact(c14 - (h9 - 6 * h8 + 21 * h7)) == 0, "channel inverse c14")

    # THM-4328's row-specific Student observer is diagonal and nonzero on
    # the three even channels (m,j)=(7,0),(8,2),(9,4).
    observer = (
        sp.Integer(1),
        sp.Rational(2, 5),
        sp.Rational(36, 85),
    )
    for (row, degree), value in zip(((7, 0), (8, 2), (9, 4)), observer):
        moment = R9.primitive_student_row(row)
        normalized = sp.Rational(moment[degree], moment[0])
        check(normalized == value, f"Student observer row {row}")
    check(sp.prod(observer) != 0, "even-face observer is lossless")
    for row, degree in ((7, 1), (8, 3)):
        moment = R9.primitive_student_row(row)
        check(moment[degree] == 0, f"odd face invisible at row {row}")

    arows, crows, bracket_subs, grows = build_through_row_eight()
    check(sp.diff(bracket_subs[R8.W], c70) == -sp.Rational(13, 6),
          "row-seven c70 response")
    check(sp.diff(bracket_subs[R8.Z], c70) == -sp.Rational(494, 81),
          "row-eight c70 response")
    check(sp.diff(bracket_subs[R8.Z], c42) == -sp.Rational(13, 18),
          "row-eight c42 response")
    check(all(sp.diff(bracket_subs[symbol], variable) == 0
              for symbol in (R8.W, R8.Z) for variable in (c51, c23)),
          "odd face does not alter selected row-seven/eight responses")

    theta8 = tuple(sp.symbols("w14_theta8_0:9"))
    depth8, amatrix8, cmatrix8 = projected_depth_residuals(arows, crows, 8)
    terminal8, depth8_conditions, matrix8, pivots8 = eliminate_linear_fibre(depth8, theta8)
    check((amatrix8.shape, amatrix8.rank()) == ((63, 131), 51), "row-eight P2 universe")
    check((cmatrix8.shape, cmatrix8.rank()) == ((72, 204), 63), "row-eight P3 universe")
    check((matrix8.shape, matrix8.rank(), pivots8) == ((21, 9), 2, (7, 8)),
          "row-eight terminal system")
    old_depth8 = (
        896 - 15 * R8.Delta,
        3 * R8.Phi + 2 * R8.zeta3,
        3030 * R8.Delta + 225 * R8.Theta - 182528,
    )
    check(len(depth8_conditions) == 3, "three row-eight conditions")
    check(all(any(proportional(new, old) for new in depth8_conditions) for old in old_depth8),
          "row-eight Hasse gate unchanged")
    gate = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * R8.Phi,
    }
    check(all(exact(condition.subs(gate)) == 0 for condition in depth8_conditions),
          "row-eight gate solves depth")
    check(all(exact(value.subs(terminal8).subs(gate)) == 0 for value in depth8),
          "all row-eight depth equations vanish")

    a8 = impose(arows, gate, terminal8, gate)
    c8rows = impose(crows, gate, terminal8, gate)
    g9 = exact(grows[9].subs(bracket_subs).subs(gate))
    row9_select, row9_conditions, row9_matrix, row9_pivots = solve_bracket_fibre(
        a8, c8rows, 9, g9, theta8[:7]
    )
    check((row9_matrix.shape, row9_matrix.rank(), row9_pivots) == ((10, 7), 7, tuple(range(7))),
          "row-nine bracket selector")
    check(len(row9_conditions) == 1, "unique row-nine condition")

    face_linear = 3339765000 * c70 + 898128000 * c42 + 216513000 * c14
    row9_equation = sp.expand(R9.E9_GATE - face_linear)
    check(proportional(row9_conditions[0], row9_equation), "row-nine face equation")
    check(all(sp.diff(row9_conditions[0], variable) == 0 for variable in (c51, c23)),
          "weight-thirteen face remains absent from row-nine equation")
    check(sp.diff(row9_equation, c14) == -216513000, "c14 transverse coefficient")
    c14_answers = sp.solve(row9_equation, c14)
    check(len(c14_answers) == 1, "global c14 graph")
    c14_graph = exact(c14_answers[0])
    check(sp.fraction(sp.together(c14_graph))[1] == 216513000,
          "constant c14 denominator")

    row9_difference = exact(g9 - R8.predicted_G(9, a8, c8rows))
    check(
        all(
            exact(xcoeff(row9_difference, degree).subs(row9_select).subs(c14, c14_graph))
            == 0
            for degree in range(10)
        ),
        "all row-nine bracket coefficients vanish on the c14 graph",
    )

    selected_a8 = impose(a8, row9_select, {c14: c14_graph})
    selected_c8 = impose(c8rows, row9_select, {c14: c14_graph})
    a9, c9, theta9 = append_tangent(selected_a8, selected_c8, 9, "w14_theta9")
    depth9, amatrix9, cmatrix9 = projected_depth_residuals(a9, c9, 9)
    terminal9, depth9_conditions, matrix9, pivots9 = eliminate_linear_fibre(depth9, theta9)
    check((amatrix9.shape, amatrix9.rank()) == ((75, 160), 59), "row-nine P2 universe")
    check((cmatrix9.shape, cmatrix9.rank()) == ((85, 251), 73), "row-nine P3 universe")
    check((matrix9.shape, matrix9.rank(), pivots9) == ((28, 10), 3, (7, 8, 9)),
          "row-nine terminal depth")
    check(depth9_conditions == [], "row-nine depth automatic")
    check(all(exact(value.subs(terminal9)) == 0 for value in depth9),
          "all row-nine depth equations vanish")

    control_base = {
        R8.Phi: 1,
        R8.eta: 2,
        R8.alpha11: 3,
        R8.xi10: 5,
        c70: 7,
        c42: 11,
    }
    control_c14 = exact(c14_graph.subs(control_base))
    control = {**control_base, c14: control_c14}
    check(exact(row9_equation.subs(control)) == 0, "rational carrier control")
    check(exact(R9.E9_GATE.subs(control)) == exact(face_linear.subs(control)) != 0,
          "hostile to inherited row-nine gate")

    print("THM-4390 weight-fourteen row-nine face absorption")
    print("imports=audited_THM4308_R8_and_THM4315_R9")
    print("field=characteristic_zero")
    print("universe=complete_fixed_residual_weight_at_most_14")
    print("weight13_face=c51*p^5*y+c23*p^2*y^3;odd_row9_invisible=yes")
    print("weight14_face=c70*p^7+c42*p^4*y^2+c14*p*y^4;complete=yes")
    print("channels=(c70,7*c70+c42,21*c70+6*c42+c14);unimodular=yes")
    print("Student_diagonal=(1,2/5,36/85);lossless=yes")
    print("row8_gate=(Delta=896/15,Theta=512/75,zeta3=-3*Phi/2);unchanged=yes")
    print("row8_depth=P2(63x131,rank51),P3(72x204,rank63),terminal(21x9,rank2)")
    print("row9_equation=E9-(3339765000*c70+898128000*c42+216513000*c14)")
    print("row9_c14_graph=" + sp.sstr(c14_graph))
    print("row9_depth=P2(75x160,rank59),P3(85x251,rank73),terminal(28x10,rank3);automatic=yes")
    print(
        "rational_control=(Phi,eta,alpha11,xi10,c70,c42,c14)="
        + sp.sstr((1, 2, 3, 5, 7, 11, control_c14))
        + ";old_E9_nonzero=yes"
    )
    print("scope=fixed_source_normal_chart_through_row9_only;row10_and_entry_open")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
