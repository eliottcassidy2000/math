#!/usr/bin/env python3
"""Primary exact certificate for proposed THM-4366.

This certificate restricts the audited source-normal row-eight family to
U=0.  It treats Phi=0 separately, continues the Phi!=0 branch through the
row-ten bracket and both projected depth modules, identifies the new
L_(10,11,3) P_2 diagonal consumer, and kills the resulting six source points
by the row-eleven bracket.  Everything is finite-row and relative to the
source-normal residual-weight-at-most-twelve universe; there is no JC(2)
claim.
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


CHECKS = 0
x = R8.x
P = R8.Phi
eta = R8.eta
alpha = R8.alpha11
beta = R8.beta11
Y = sp.symbols("Y")

XI_U = sp.Rational(237757952, 54675)

F0 = sp.Poly(
    5646560625 * eta**2 + 379697122115584,
    eta,
    domain=sp.QQ,
)

ETA10 = -8 * (4100625 * P**2 - 219011178496) / (36905625 * P)
ALPHA10 = (
    145629074958046875 * P**4
    - 6583704417821122560000 * P**2
    - 26093447590576484415176704
) / (23154427662890625 * P**3)
BETA10 = -8 * (
    39878843693546535205078125 * P**6
    - 1832393315703878491921920000000 * P**4
    - 12681601488957070429831550730240000 * P**2
    - 54290188724439791780028059265687093248
) / (7690757619746410400390625 * P**5)

Q = sp.Poly(
    373891487235896675830078125 * Y**3
    - 15097287707154073014589440000000 * Y**2
    - 101452811911656563438652405841920000 * Y
    - 434321509795518334240224474125496745984,
    Y,
    domain=sp.QQ,
)
QP = sp.Poly(Q.as_expr().subs(Y, P**2), P, domain=sp.QQ)

DISC = sp.Poly(
    -818189072042144733182993060302734375 * Y**4
    - 15281178725265792138311029689600000000000 * Y**3
    + 35745234591073505392643546713864273920000000 * Y**2
    + 343583092356524653010661045438183132081684480000 * Y
    + 680868007162161739848062587106828346429167544303616,
    Y,
    domain=sp.QQ,
)

R11 = sp.Poly(
    6846329377771290182382546697998046875 * Y**4
    - 713835723041306505264998768716800000000000 * Y**3
    - 2754991513504883058403611855707575418880000000 * Y**2
    - 31916203206707002973657986739896008646412206080000 * Y
    - 156854967149983010817497418504735580308619473018945536,
    Y,
    domain=sp.QQ,
)


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def exact_zero(expression: sp.Expr) -> bool:
    return sp.cancel(sp.together(sp.expand(expression))) == 0


def check_zero(expression: sp.Expr, label: str) -> None:
    check(exact_zero(expression), label)


def reduce_mod(expression: sp.Expr, variable: sp.Symbol, modulus: sp.Poly) -> sp.Expr:
    """Reduce a rational expression in Frac(other variables)[variable]/(modulus)."""

    value = sp.cancel(sp.together(sp.expand(expression)))
    numerator, denominator = sp.fraction(value)
    mod = sp.Poly(modulus.as_expr(), variable, domain=sp.EX)
    numerator_poly = sp.Poly(numerator, variable, domain=sp.EX)
    denominator_poly = sp.Poly(denominator, variable, domain=sp.EX)
    inverse = sp.invert(denominator_poly, mod)
    return sp.cancel((numerator_poly * inverse).rem(mod).as_expr())


def selected_solution(
    matrix: sp.Matrix,
    rhs: sp.Matrix,
    columns: tuple[int, ...],
) -> tuple[sp.Expr, ...]:
    """Solve the lexicographically inherited full-row-rank square block."""

    row_indices = tuple(matrix.T.rref()[1])
    square = matrix.extract(row_indices, columns)
    check(square.rows == square.cols, "selected square block")
    check(square.det() != 0, "selected block invertible")
    result = square.inv() * rhs.extract(row_indices, (0,))
    return tuple(sp.cancel(value) for value in result)


def hierarchy_value(
    rows: list[sp.Expr],
    m: int,
    ell: int,
    q: int,
    substitutions: dict[sp.Symbol, sp.Expr],
) -> sp.Expr:
    s = (ell + 1) // 2
    value = sum(
        (-1) ** (n - s)
        * sp.binomial(m + q - n, q)
        * R8.xcoeff(rows[n], 2 * n - ell)
        for n in range(s, m + 1)
    )
    return sp.factor(sp.cancel(value.subs(substitutions)))


def even_numerator_polynomial(expression: sp.Expr) -> sp.Poly:
    """Clear a rational expression and replace P^(2j) by Y^j primitively."""

    numerator, _ = sp.fraction(sp.cancel(sp.together(expression)))
    polynomial = sp.Poly(numerator, P, domain=sp.QQ)
    result = sp.Integer(0)
    for (degree,), coefficient in polynomial.terms():
        check(degree % 2 == 0, "even P numerator")
        result += coefficient * Y ** (degree // 2)
    return sp.Poly(result, Y, domain=sp.QQ).primitive()[1]


def add_row_nine_depth(
    fixed_a: list[sp.Expr],
    fixed_c: list[sp.Expr],
) -> tuple[
    list[sp.Expr],
    list[sp.Expr],
    tuple[sp.Symbol, ...],
    list[sp.Expr],
    sp.Matrix,
    sp.Matrix,
]:
    """Append row nine and solve the exact projected P_2/P_3 terminal equations."""

    b9 = R8.B_row(9, fixed_a, fixed_c)
    a9base, c9base = R8.particular_row(9, b9)
    theta9_symbols = tuple(sp.symbols("theta9_0:10"))
    theta9 = sum(theta9_symbols[j] * x**j for j in range(10))
    a9 = sp.expand(a9base + theta9 * sp.diff(R8.A0, x))
    c9 = sp.expand(c9base + theta9 * sp.diff(R8.C0, x))
    row9_a = fixed_a + [a9]
    row9_c = fixed_c + [c9]

    acoords9, amatrix9 = R9.depth_matrix(2, 9)
    ccoords9, cmatrix9 = R9.depth_matrix(3, 9)
    check(
        (len(acoords9), amatrix9.cols, amatrix9.rank(), len(amatrix9.T.nullspace()))
        == (75, 160, 59, 16),
        "row-nine P2 universe",
    )
    check(
        (len(ccoords9), cmatrix9.cols, cmatrix9.rank(), len(cmatrix9.T.nullspace()))
        == (85, 251, 73, 12),
        "row-nine P3 universe",
    )
    avec9 = sp.Matrix([R8.xcoeff(row9_a[n], r) for n, r in acoords9])
    cvec9 = sp.Matrix([R8.xcoeff(row9_c[n], r) for n, r in ccoords9])
    residuals9 = [sp.expand((row.T * avec9)[0]) for row in amatrix9.T.nullspace()]
    residuals9 += [sp.expand((row.T * cvec9)[0]) for row in cmatrix9.T.nullspace()]
    depth9_matrix, depth9_rhs = sp.linear_eq_to_matrix(residuals9, theta9_symbols)
    check(depth9_matrix.rank() == 3, "row-nine terminal rank")
    check(depth9_matrix.rref()[1] == (7, 8, 9), "row-nine terminal pivots")
    check(depth9_matrix[:, :7] == sp.zeros(depth9_matrix.rows, 7), "row-nine free columns")
    depth9_solution = selected_solution(depth9_matrix, depth9_rhs, (7, 8, 9))
    depth9_subs = {
        theta9_symbols[7 + j]: depth9_solution[j]
        for j in range(3)
    }
    return row9_a, row9_c, theta9_symbols, residuals9, depth9_matrix, sp.Matrix(
        [depth9_subs[theta9_symbols[j]] for j in range(7, 10)]
    )


def boundary_phi_zero(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    bracket_subs: dict[sp.Symbol, sp.Expr],
    theta8_symbols: list[sp.Symbol],
    terminal8: dict[sp.Symbol, sp.Expr],
) -> None:
    """The Phi=0 boundary has row-nine points but no row-ten bracket lift."""

    source_subs = {P: 0, R8.xi10: XI_U}
    check_zero(R9.E9_GATE.subs(source_subs) + sp.Rational(11, 243) * F0.as_expr(), "boundary E9")
    check(sp.gcd(F0, F0.diff()) == sp.Poly(1, eta, domain=sp.QQ), "boundary F0 squarefree")
    check(F0.degree() == 2 and F0.TC() != 0, "boundary two nonzero roots")

    response = {
        symbol: sp.factor(value.subs(terminal8).subs(source_subs))
        for symbol, value in bracket_subs.items()
    }
    check(response[R8.U] == 0, "boundary U=0")
    check(response[R8.W] == -sp.Rational(42434609152, 2460375), "boundary W")
    check(response[R8.Z] == sp.Rational(5200877686784, 66430125), "boundary Z")
    check(response[R8.upsilon5] != 0, "boundary upsilon nonzero")
    check(-4 * response[R8.W] * response[R8.upsilon5] != 0, "boundary strict discriminant control")

    g9 = sp.expand(R8.tcoeff(R8.H, 9))
    fixed_a, fixed_c, old_tangent, _ = R9.solve_row_eight_for_g9(
        arows,
        crows,
        bracket_subs,
        terminal8,
        theta8_symbols,
        g9,
    )
    check((old_tangent.rows, old_tangent.cols, old_tangent.rank()) == (10, 7, 7), "boundary G9 old selection")
    fixed_a = [
        sp.expand(row.subs(bracket_subs).subs(terminal8).subs(source_subs))
        for row in fixed_a
    ]
    fixed_c = [
        sp.expand(row.subs(bracket_subs).subs(terminal8).subs(source_subs))
        for row in fixed_c
    ]
    bracket9_defect = R8.predicted_G(9, fixed_a, fixed_c) - g9.subs(bracket_subs).subs(terminal8).subs(source_subs)
    check(reduce_mod(bracket9_defect, eta, F0) == 0, "boundary G9 in F0 quotient")

    row9_a, row9_c, theta9, residuals9, _, solution9 = add_row_nine_depth(fixed_a, fixed_c)
    depth9_subs = {theta9[7 + j]: solution9[j] for j in range(3)}
    check(
        all(reduce_mod(value.subs(depth9_subs), eta, F0) == 0 for value in residuals9),
        "boundary row-nine depth",
    )
    row9_a_depth = fixed_a + [sp.expand(row9_a[9].subs(depth9_subs))]
    row9_c_depth = fixed_c + [sp.expand(row9_c[9].subs(depth9_subs))]

    g10 = sp.expand(R8.tcoeff(R8.H, 10)).subs(bracket_subs).subs(terminal8).subs(source_subs)
    equations10 = [
        R8.xcoeff(g10 - R8.predicted_G(10, row9_a_depth, row9_c_depth), j)
        for j in range(11)
    ]
    matrix10, rhs10 = sp.linear_eq_to_matrix(equations10, theta9[:7])
    check(matrix10.rank() == 7, "boundary G10 free rank")
    check(matrix10.T.rref()[1] == (0, 1, 2, 3, 4, 5, 7), "boundary G10 pivot rows")
    solution10 = selected_solution(matrix10, rhs10, tuple(range(7)))
    substitutions10 = {theta9[j]: solution10[j] for j in range(7)}
    reduced10 = [
        sp.factor(reduce_mod(value.subs(substitutions10), eta, F0))
        for value in equations10
    ]
    expected6 = (
        321853955625 * alpha * eta - 34327712203669504
    ) / 215233605000
    expected8 = -sp.Rational(11388581281792, 332150625)
    check_zero(reduced10[6] - expected6, "boundary G10 residual six")
    check_zero(reduced10[8] - expected8, "boundary G10 fixed residual eight")
    check(
        all(value == 0 or index in (6, 8) for index, value in enumerate(reduced10)),
        "boundary only two G10 residuals",
    )
    check(expected8 != 0, "boundary row-ten extinction")


def nonzero_branch(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    bracket_subs: dict[sp.Symbol, sp.Expr],
    theta8_symbols: list[sp.Symbol],
    terminal8: dict[sp.Symbol, sp.Expr],
) -> None:
    """Continue the Phi!=0 branch through row eleven."""

    alpha9 = sp.factor(sp.solve(R9.E9_GATE.subs(R8.xi10, XI_U), alpha)[0])
    expected_alpha9 = (
        13553385750 * P**2
        - 69677820000 * P * eta
        - 5646560625 * eta**2
        - 379697122115584
    ) / (11293121250 * P)
    check_zero(alpha9 - expected_alpha9, "nonzero E9 alpha graph")
    source9_subs = {R8.xi10: XI_U, alpha: alpha9}
    check_zero(R9.E9_GATE.subs(source9_subs), "nonzero E9 solved")

    response9 = {
        symbol: sp.factor(value.subs(terminal8).subs(source9_subs))
        for symbol, value in bracket_subs.items()
    }
    expected_w = -13 * (820125 * P**2 + 13056802816) / 9841500
    expected_z = -13 * (
        32805000 * P**2 + 36905625 * P * eta - 1600270057472
    ) / 265720500
    check(response9[R8.U] == 0, "nonzero U=0")
    check_zero(response9[R8.W] - expected_w, "nonzero W response")
    check_zero(response9[R8.Z] - expected_z, "nonzero Z response")

    g9 = sp.expand(R8.tcoeff(R8.H, 9))
    fixed_a, fixed_c, old_tangent, _ = R9.solve_row_eight_for_g9(
        arows,
        crows,
        bracket_subs,
        terminal8,
        theta8_symbols,
        g9,
    )
    check((old_tangent.rows, old_tangent.cols, old_tangent.rank()) == (10, 7, 7), "nonzero G9 old selection")
    fixed_a = [
        sp.cancel(row.subs(bracket_subs).subs(terminal8).subs(source9_subs))
        for row in fixed_a
    ]
    fixed_c = [
        sp.cancel(row.subs(bracket_subs).subs(terminal8).subs(source9_subs))
        for row in fixed_c
    ]
    g9_source = g9.subs(bracket_subs).subs(terminal8).subs(source9_subs)
    check_zero(R8.predicted_G(9, fixed_a, fixed_c) - g9_source, "nonzero G9 bracket")

    row9_a, row9_c, theta9, residuals9, _, solution9 = add_row_nine_depth(fixed_a, fixed_c)
    depth9_subs = {theta9[7 + j]: solution9[j] for j in range(3)}
    check(all(exact_zero(value.subs(depth9_subs)) for value in residuals9), "nonzero row-nine depth")
    row9_a_depth = fixed_a + [sp.cancel(row9_a[9].subs(depth9_subs))]
    row9_c_depth = fixed_c + [sp.cancel(row9_c[9].subs(depth9_subs))]

    g10_literal = sp.expand(R8.tcoeff(R8.H, 10))
    check_zero(
        g10_literal
        - (
            (15 * R8.U + 10 * R8.W + 6 * R8.Z) * x**8
            + (5 * alpha + 4 * beta) * x**9
            + (R8.upsilon5 + R8.xi10) * x**10
        ),
        "literal G10",
    )
    g10 = g10_literal.subs(bracket_subs).subs(terminal8).subs(source9_subs)
    equations10 = [
        R8.xcoeff(g10 - R8.predicted_G(10, row9_a_depth, row9_c_depth), j)
        for j in range(11)
    ]
    matrix10, rhs10 = sp.linear_eq_to_matrix(equations10, theta9[:7])
    check(matrix10.rank() == 7, "nonzero G10 free rank")
    check(matrix10.T.rref()[1] == (0, 1, 2, 3, 4, 5, 7), "nonzero G10 pivot rows")
    solution10 = selected_solution(matrix10, rhs10, tuple(range(7)))
    g10_subs = {theta9[j]: solution10[j] for j in range(7)}
    residuals10 = [sp.factor(sp.cancel(value.subs(g10_subs))) for value in equations10]
    expected6 = (
        1772782200000 * P**3
        + 609828547500 * P**2 * beta
        - 29740924426500 * P**2 * eta
        - 3971635740000 * P * eta**2
        - 68655424407339008 * P
        - 321853955625 * eta**3
        - 21642735960588288 * eta
    ) / (430467210000 * P)
    expected8 = 13 * (
        32805000 * P**2 + 36905625 * P * eta - 1752089427968
    ) / 664301250
    check_zero(residuals10[6] - expected6, "nonzero G10 residual six")
    check_zero(residuals10[8] - expected8, "nonzero G10 residual eight")
    check(
        all(value == 0 or index in (6, 8) for index, value in enumerate(residuals10)),
        "nonzero only two G10 residuals",
    )
    solved_eta = sp.factor(sp.solve(residuals10[8], eta)[0])
    solved_beta = sp.factor(sp.solve(residuals10[6].subs(eta, solved_eta), beta)[0])
    check_zero(solved_eta - ETA10, "row-ten eta graph")
    check_zero(solved_beta - BETA10, "row-ten beta graph")
    locus10_subs = {eta: ETA10, beta: BETA10}
    check(all(exact_zero(value.subs(locus10_subs)) for value in residuals10), "row-ten bracket locus")
    check_zero(alpha9.subs(eta, ETA10) - ALPHA10, "row-ten alpha graph")
    check_zero(response9[R8.Z].subs(eta, ETA10) + sp.Rational(225611776, 30375), "row-ten Z constant")

    selected_a9 = [sp.cancel(row.subs(g10_subs).subs(locus10_subs)) for row in row9_a_depth]
    selected_c9 = [sp.cancel(row.subs(g10_subs).subs(locus10_subs)) for row in row9_c_depth]
    b10 = R8.B_row(10, selected_a9, selected_c9)
    a10base, c10base = R8.particular_row(10, b10)
    theta10_symbols = tuple(sp.symbols("theta10_0:11"))
    theta10 = sum(theta10_symbols[j] * x**j for j in range(11))
    a10 = sp.expand(a10base + theta10 * sp.diff(R8.A0, x))
    c10 = sp.expand(c10base + theta10 * sp.diff(R8.C0, x))
    row10_a = selected_a9 + [a10]
    row10_c = selected_c9 + [c10]

    acoords10, amatrix10 = R9.depth_matrix(2, 10)
    ccoords10, cmatrix10 = R9.depth_matrix(3, 10)
    anull10 = amatrix10.T.nullspace()
    cnull10 = cmatrix10.T.nullspace()
    check(
        (len(acoords10), amatrix10.cols, amatrix10.rank(), len(anull10))
        == (88, 193, 68, 20),
        "row-ten P2 universe",
    )
    check(
        (len(ccoords10), cmatrix10.cols, cmatrix10.rank(), len(cnull10))
        == (99, 304, 83, 16),
        "row-ten P3 universe",
    )
    avec10 = sp.Matrix([R8.xcoeff(row10_a[n], r) for n, r in acoords10])
    cvec10 = sp.Matrix([R8.xcoeff(row10_c[n], r) for n, r in ccoords10])
    aresiduals10 = [sp.expand((row.T * avec10)[0]) for row in anull10]
    cresiduals10 = [sp.expand((row.T * cvec10)[0]) for row in cnull10]
    amatrix_terminal, arhs_terminal = sp.linear_eq_to_matrix(aresiduals10, theta10_symbols)
    cmatrix_terminal, crhs_terminal = sp.linear_eq_to_matrix(cresiduals10, theta10_symbols)
    check(amatrix_terminal.rank() == 3, "row-ten standalone P2 rank")
    check(cmatrix_terminal.rank() == 3, "row-ten standalone P3 rank")
    check(amatrix_terminal.rref()[1] == (8, 9, 10), "row-ten P2 pivots")
    check(cmatrix_terminal.rref()[1] == (8, 9, 10), "row-ten P3 pivots")
    check(amatrix_terminal[:, :8] == sp.zeros(amatrix_terminal.rows, 8), "row-ten P2 free columns")
    check(cmatrix_terminal[:, :8] == sp.zeros(cmatrix_terminal.rows, 8), "row-ten P3 free columns")
    joint_matrix = amatrix_terminal.col_join(cmatrix_terminal)
    check(joint_matrix.rank() == 3, "row-ten joint coefficient rank")

    asolution10 = selected_solution(amatrix_terminal, arhs_terminal, (8, 9, 10))
    csolution10 = selected_solution(cmatrix_terminal, crhs_terminal, (8, 9, 10))
    asubs10 = {theta10_symbols[8 + j]: asolution10[j] for j in range(3)}
    csubs10 = {theta10_symbols[8 + j]: csolution10[j] for j in range(3)}
    check(all(exact_zero(value.subs(asubs10)) for value in aresiduals10), "P2 standalone compatible")
    check(all(exact_zero(value.subs(csubs10)) for value in cresiduals10), "P3 standalone compatible")
    expected_delta = 4 * QP.as_expr() / (23072272859239231201171875 * P**5)
    check_zero(asolution10[0] - csolution10[0] - expected_delta, "joint theta10_8 difference")
    check_zero(asolution10[1] - csolution10[1], "joint theta10_9 agreement")
    check_zero(asolution10[2] - csolution10[2], "joint theta10_10 agreement")

    # The new opposite-order hierarchy consumer: select P3, then read P2.
    a_terms = {(6, 1): 35, (7, 3): -20, (8, 5): 10, (9, 7): -4, (10, 9): 1}
    a_vector = sp.zeros(len(acoords10), 1)
    for coordinate, coefficient in a_terms.items():
        a_vector[acoords10.index(coordinate)] = coefficient
    check(a_vector.T * amatrix10 == sp.zeros(1, amatrix10.cols), "L(10,11,3) annihilates P2")
    l_a_after_c = hierarchy_value(row10_a, 10, 11, 3, csubs10)
    expected_l_a = -2 * QP.as_expr() / (23072272859239231201171875 * P**5)
    check_zero(l_a_after_c - expected_l_a, "new A hierarchy identity")

    # The old opposite ordering gives the same source ideal, not an extra cut.
    c_terms = {(5, 0): 56, (6, 2): -35, (7, 4): 20, (8, 6): -10, (9, 8): 4, (10, 10): -1}
    c_vector = sp.zeros(len(ccoords10), 1)
    for coordinate, coefficient in c_terms.items():
        c_vector[ccoords10.index(coordinate)] = coefficient
    check(c_vector.T * cmatrix10 == sp.zeros(1, cmatrix10.cols), "L(10,10,3) annihilates P3")
    l_c_after_a = hierarchy_value(row10_c, 10, 10, 3, asubs10)
    expected_l_c = QP.as_expr() / (15381515239492820800781250 * P**5)
    check_zero(l_c_after_a - expected_l_c, "old C hierarchy identity")
    check_zero(l_a_after_c + sp.Rational(4, 3) * l_c_after_a, "hierarchy consumers same ideal")

    nonzero_a = []
    nonzero_c = []
    for ell in range(2, 21):
        rho = (ell + 2) // 3
        for q in range(rho):
            if 10 + q >= ell and 2 <= 10 + q - ell:
                value = hierarchy_value(row10_a, 10, ell, q, csubs10)
                if value != 0:
                    nonzero_a.append((ell, q))
            if 10 + q >= ell and 3 <= 10 + q - ell:
                value = hierarchy_value(row10_c, 10, ell, q, asubs10)
                if value != 0:
                    nonzero_c.append((ell, q))
    check(nonzero_a == [(11, 3)], "unique nonzero admissible A hierarchy consumer")
    check(nonzero_c == [(10, 3)], "unique nonzero admissible C hierarchy consumer")

    check(Q.primitive()[0] == 1 and Q.degree() == 3, "Q primitive cubic")
    check(Q.TC() != 0, "Q nonzero roots")
    check(sp.gcd(Q, Q.diff()) == sp.Poly(1, Y, domain=sp.QQ), "Q squarefree")
    check(sp.gcd(QP, QP.diff()) == sp.Poly(1, P, domain=sp.QQ), "six Phi sheets distinct")
    q7 = sp.Poly(Q.as_expr(), Y, modulus=7)
    check(q7 == sp.Poly(-Y**3 + 3 * Y**2 + 3 * Y + 2, Y, modulus=7), "Q mod7")
    q7_squarefree = sp.Poly(
        (Y - 1) * q7.as_expr() + (2 * Y**2 + 3 * Y + 1) * q7.diff().as_expr() - 1,
        Y,
        modulus=7,
    )
    check(q7_squarefree.is_zero, "Q squarefree mod7 Bezout")

    w_factor = sp.Poly(820125 * Y + 13056802816, Y, domain=sp.QQ)
    check(sp.gcd(Q, w_factor) == sp.Poly(1, Y, domain=sp.QQ), "Q avoids W=0")
    w7 = sp.Poly(w_factor.as_expr(), Y, modulus=7)
    check(w7 == sp.Poly(1 - 2 * Y, Y, modulus=7), "W factor mod7")
    check(
        sp.Poly(3 * q7.as_expr() + (2 * Y**2 + 2 * Y + 2) * w7.as_expr() - 1, Y, modulus=7).is_zero,
        "Q/W mod7 Bezout",
    )

    alpha10_actual = sp.factor(alpha9.subs(eta, ETA10))
    w10 = sp.factor(response9[R8.W].subs(eta, ETA10))
    z10 = sp.factor(response9[R8.Z].subs(eta, ETA10))
    discriminant = sp.factor(alpha10_actual**2 - 4 * w10 * response9[R8.upsilon5])
    check(even_numerator_polynomial(discriminant) == DISC, "strict discriminant polynomial")
    check(sp.gcd(Q, DISC) == sp.Poly(1, Y, domain=sp.QQ), "Q avoids strict discriminant wall")
    d7 = sp.Poly(DISC.as_expr(), Y, modulus=7)
    check(d7 == sp.Poly(-Y**4 + 2 * Y**3 - Y + 1, Y, modulus=7), "discriminant mod7")
    check(
        sp.Poly((-Y**3 - Y**2 + Y - 3) * q7.as_expr() + Y**2 * d7.as_expr() - 1, Y, modulus=7).is_zero,
        "Q/discriminant mod7 Bezout",
    )
    check(z10 == -sp.Rational(225611776, 30375), "six-point Z strict gate")
    check(response9[R8.upsilon5] != 0 and R8.Delta.subs(terminal8) != 0, "remaining strict gates")

    # P=1 is a rational row-ten bracket/standalone-depth hostile to any
    # premature extinction statement.  It fails only the joint equation Q=0.
    check(Q.eval(1) != 0, "P=1 joint-depth hostile")
    check(w10.subs(P, 1) != 0 and z10 != 0 and discriminant.subs(P, 1) != 0, "P=1 strict hostile")

    # On Q=0 the A-selected terminal also satisfies P3.  Continue its eight
    # free terminal directions to the literal row-eleven source.
    row10_a_joint = selected_a9 + [sp.cancel(a10.subs(asubs10))]
    row10_c_joint = selected_c9 + [sp.cancel(c10.subs(asubs10))]
    check(
        all(reduce_mod(value.subs(asubs10), P, QP) == 0 for value in cresiduals10),
        "joint row-ten depth on Q",
    )
    g11_literal = sp.expand(R8.tcoeff(R8.H, 11))
    check_zero(
        g11_literal - ((6 * R8.U + 5 * R8.W + 4 * R8.Z) * x**10 + (alpha + beta) * x**11),
        "literal G11",
    )
    g11 = g11_literal.subs(bracket_subs).subs(terminal8).subs(source9_subs).subs(locus10_subs)
    equations11 = [
        R8.xcoeff(g11 - R8.predicted_G(11, row10_a_joint, row10_c_joint), j)
        for j in range(12)
    ]
    matrix11, rhs11 = sp.linear_eq_to_matrix(equations11, theta10_symbols[:8])
    check(matrix11.rank() == 8, "G11 free terminal rank")
    check(matrix11.T.rref()[1] == tuple(range(8)), "G11 pivot rows")
    solution11 = selected_solution(matrix11, rhs11, tuple(range(8)))
    g11_subs = {theta10_symbols[j]: solution11[j] for j in range(8)}
    residuals11 = [sp.factor(sp.cancel(value.subs(g11_subs))) for value in equations11]
    expected_r8 = -R11.as_expr().subs(Y, P**2) / (
        71525718603423911567894897460937500 * P**6
    )
    expected_r9 = 8 * QP.as_expr() / (84598333817210514404296875 * P**5)
    check_zero(residuals11[8] - expected_r8, "G11 residual eight")
    check_zero(residuals11[9] - expected_r9, "G11 residual nine")
    check(
        all(value == 0 or index in (8, 9) for index, value in enumerate(residuals11)),
        "only two G11 residuals",
    )
    check(R11.primitive()[0] == 1 and R11.degree() == 4, "R11 primitive quartic")
    check(sp.gcd(Q, R11) == sp.Poly(1, Y, domain=sp.QQ), "Q/R11 coprime")
    r7 = sp.Poly(R11.as_expr(), Y, modulus=7)
    check(
        r7 == sp.Poly(-2 * Y**4 + 3 * Y**3 + 3 * Y**2 + 3 * Y - 2, Y, modulus=7),
        "R11 mod7",
    )
    check(
        sp.Poly(
            (3 * Y**3 - 2 * Y**2 - 3 * Y) * q7.as_expr()
            + (2 * Y**2 - 2 * Y + 3) * r7.as_expr()
            - 1,
            Y,
            modulus=7,
        ).is_zero,
        "Q/R11 mod7 Bezout",
    )

    print(f"ALPHA9={sp.sstr(alpha9)}")
    print(f"ROW10_ETA={sp.sstr(ETA10)}")
    print(f"ROW10_ALPHA={sp.sstr(ALPHA10)}")
    print(f"ROW10_BETA={sp.sstr(BETA10)}")
    print(f"ROW10_W={sp.sstr(w10)} ROW10_Z={sp.sstr(z10)}")
    print(f"ROW10_Q={sp.sstr(Q.as_expr())}")
    print(f"STRICT_DISCRIMINANT_NUMERATOR={sp.sstr(DISC.as_expr())}")
    print(f"ROW11_R={sp.sstr(R11.as_expr())}")


def main() -> None:
    arows, crows, bracket_subs = R8.build_bracket_rows()
    amatrix8, cmatrix8, theta8_symbols, terminal8 = R8.hasse_audit(
        arows,
        crows,
        bracket_subs,
    )
    check(
        (amatrix8.rows, amatrix8.cols, amatrix8.rank()) == (63, 131, 51),
        "inherited row-eight P2",
    )
    check(
        (cmatrix8.rows, cmatrix8.cols, cmatrix8.rank()) == (72, 204, 63),
        "inherited row-eight P3",
    )
    check_zero(bracket_subs[R8.U].subs(terminal8).subs(R8.xi10, XI_U), "XI_U really U=0")

    boundary_phi_zero(arows, crows, bracket_subs, theta8_symbols, terminal8)
    nonzero_branch(arows, crows, bracket_subs, theta8_symbols, terminal8)

    print("THM-4366 PRIMARY: U-ZERO ROW-ELEVEN HIERARCHY-SELECTED EXTINCTION PASS")
    print(f"CHECKS={CHECKS}")
    print(f"U_ZERO xi_10={XI_U}; Phi,eta,alpha_11,beta_11 initially free")
    print(f"PHI_ZERO E9=(-11/243)*({sp.sstr(F0.as_expr())}); two row9 components with fibre A7")
    print("PHI_ZERO G10 selected residual[8]=-11388581281792/332150625 !=0; dies before row10")
    print("PHI_NONZERO row9 source=G_mxA2; row10 bracket solves eta,beta and leaves Phi in G_m")
    print("ROW10_DEPTH P2=88x193/rank68/null20 terminal=3/aug3; P3=99x304/rank83/null16 terminal=3/aug3")
    print("JOINT_DEPTH compatible iff Q(Phi^2)=0; Q squarefree cubic, Q(0)!=0, six sign sheets, fibre A8")
    print("NEW_A_L(10,11,3)=35a_(6,1)-20a_(7,3)+10a_(8,5)-4a_(9,7)+a_(10,9)")
    print("AFTER_P3_SELECTION new_A_L=-2Q(Phi^2)/(23072272859239231201171875*Phi^5)")
    print("AFTER_P2_SELECTION old_C_L=Q(Phi^2)/(15381515239492820800781250*Phi^5)")
    print("REDUNDANCY source ideals agree; new_A_L=(-4/3)old_C_L, so it is an alternative opposite-order certificate, not an extra cut")
    print("HIERARCHY_SCAN at m=10: only nonzero admissible opposite-order consumers are A(ell,q)=(11,3), C(ell,q)=(10,3)")
    print("STRICT Q avoids Phi=0, W=0, Z=0, and alpha_11^2-4W*upsilon_5=0")
    print("ROW11 terminal rank=8; residual[9]=8Q/(84598333817210514404296875*Phi^5)")
    print("ROW11 residual[8]=-R(Phi^2)/(71525718603423911567894897460937500*Phi^6); gcd(Q,R)=1")
    print("MOD7 Q=-Y^3+3Y^2+3Y+2; R=-2Y^4+3Y^3+3Y^2+3Y-2")
    print("MOD7_BEZOUT (3Y^3-2Y^2-3Y)Q+(2Y^2-2Y+3)R=1")
    print("HOSTILES Phi=0 has genuine row9 points; Phi=1 passes row10 bracket and each depth module separately but fails joint Q")
    print("SCOPE exact finite source-normal residual_weight<=12 projected-depth extinction through row11; no all-row lift, termination, seam entry, Keller pair, JC2, or DC2")


if __name__ == "__main__":
    main()
