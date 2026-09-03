#!/usr/bin/env python3
"""Primary exact certificate for THM-4361.

Restrict THM-4308's source-normal row-eight family to THM-4357's
three-dimensional beta-zero endpoint pullback of THM-4334.  The certificate
derives the row-nine bracket graph and projected-depth fibre, proves that the
row-ten bracket plus inherited row-nine depth leave six sign-sheeted source
points, and then proves that the *joint* row-ten P_2/P_3 projected-depth
system is inconsistent there.

The decisive P_3 row is compatible when P_3 is tested alone.  It becomes an
obstruction only after P_2 has selected the three visible terminal tangent
coordinates.  This distinction is part of the theorem.

This is a finite projected-depth result in the fixed residual-weight-at-most-
twelve universe.  It proves no all-row B_2 lift, termination, seam entry,
Keller pair, JC(2), or DC(2) conclusion.
"""

from __future__ import annotations

from math import gcd
from pathlib import Path
import sys

import sympy as sp


# Keep frozen output byte-identical on Windows, normal mode, and -O mode.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_student_stein_row9_thm4315 as R9  # noqa: E402


R8 = R9.R8
x = R8.x
P = R8.Phi
X = R8.xi10
eta = R8.eta
alpha = R8.alpha11
beta = R8.beta11
Y = sp.symbols("Y")

X0 = -sp.Rational(1928704, 10125)

ETA_GRAPH = sp.cancel(
    (
        sp.Integer(12506118074368)
        - sp.Integer(173745000) * P**2
        - sp.Integer(926883000) * X
    )
    / (sp.Integer(195463125) * P)
)

ALPHA_GRAPH = sp.cancel(
    (
        sp.Integer(4085005423617421875) * P**4
        + sp.Integer(24824465812575000000) * P**2 * X
        - sp.Integer(278518552828793671680000) * P**2
        - sp.Integer(7302452813356500000) * X**2
        + sp.Integer(197059040065115394048000) * X
        - sp.Integer(1329425408965288765218095104)
    )
    / (sp.Integer(649499164991015625) * P**3)
)

Q = sp.Poly(
    sp.Integer(2779225183740234375) * Y**3
    - sp.Integer(194721282033880320000000) * Y**2
    - sp.Integer(1868800030080493839974400000) * Y
    - sp.Integer(9659395340042262184105231777792),
    Y,
    domain=sp.QQ,
)

D = sp.Poly(
    sp.Integer(1482253431328125) * Y**2
    - sp.Integer(103851350418069504000) * Y
    - sp.Integer(1468662667625265243357184),
    Y,
    domain=sp.QQ,
)

QP = sp.Poly(Q.as_expr().subs(Y, P**2), P, domain=sp.QQ)
DP = sp.Poly(D.as_expr().subs(Y, P**2), P, domain=sp.QQ)

H_DENOMINATOR = sp.Integer(132327846238037905244160)
Q_CORRECTION_DENOMINATOR = sp.Integer(248114711696321072332800000)

CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def exact_zero(expression: sp.Expr) -> bool:
    return sp.cancel(sp.together(sp.expand(expression))) == 0


def check_zero(expression: sp.Expr, label: str) -> None:
    check(exact_zero(expression), label)


def selected_solution(
    matrix: sp.Matrix,
    rhs: sp.Matrix,
    row_indices: tuple[int, ...],
    column_indices: tuple[int, ...],
    free_symbols: tuple[sp.Symbol, ...] = (),
    free_columns: tuple[int, ...] = (),
) -> sp.Matrix:
    """Solve a declared square pivot block, retaining declared free terms."""

    square = matrix.extract(row_indices, column_indices)
    check(square.rows == square.cols, "selected block is square")
    check(square.det() != 0, "selected block is invertible")
    effective_rhs = rhs.extract(row_indices, (0,))
    if free_symbols:
        free_matrix = matrix.extract(row_indices, free_columns)
        effective_rhs -= free_matrix * sp.Matrix(free_symbols)
    return sp.simplify(square.inv() * effective_rhs)


def reduce_q(expression: sp.Expr) -> sp.Expr:
    """Reduce in QQ[P,P^-1]/(q(P^2)); q(0)!=0 makes P a unit."""

    value = sp.cancel(sp.together(sp.expand(expression)))
    numerator, denominator = sp.fraction(value)
    numerator_poly = sp.Poly(numerator, P, domain=sp.EX)
    denominator_poly = sp.Poly(denominator, P, domain=sp.EX)
    inverse = sp.invert(denominator_poly, QP)
    return sp.cancel((numerator_poly * inverse).rem(QP).as_expr())


def primitive_integer_functional(vector: sp.Matrix) -> tuple[int, ...]:
    denominator = 1
    for value in vector:
        denominator = sp.ilcm(denominator, int(sp.denom(value)))
    entries = [int(value * denominator) for value in vector]
    common = 0
    for value in entries:
        common = gcd(common, abs(value))
    entries = [value // common for value in entries]
    first = next(value for value in entries if value)
    if first < 0:
        entries = [-value for value in entries]
    return tuple(entries)


def residue(expression: sp.Expr, prime: int, substitutions: dict[sp.Symbol, int]) -> int:
    value = sp.cancel(sp.together(expression.subs(substitutions)))
    numerator, denominator = sp.fraction(value)
    numerator_value = int(numerator) % prime
    denominator_value = int(denominator) % prime
    check(denominator_value != 0, f"denominator unit modulo {prime}")
    return numerator_value * pow(denominator_value, -1, prime) % prime


def main() -> None:
    # I. Reconstruct the inherited exact row-eight family.
    arows, crows, bracket_subs = R8.build_bracket_rows()
    amatrix8, cmatrix8, theta8_symbols, terminal_subs = R8.hasse_audit(
        arows, crows, bracket_subs
    )
    check(
        (amatrix8.rows, amatrix8.cols, amatrix8.rank()) == (63, 131, 51),
        "inherited row-eight P2 universe",
    )
    check(
        (cmatrix8.rows, cmatrix8.cols, cmatrix8.rank()) == (72, 204, 63),
        "inherited row-eight P3 universe",
    )

    graph_subs = {eta: ETA_GRAPH, alpha: ALPHA_GRAPH, beta: 0}
    response = {
        symbol: sp.cancel(value.subs(terminal_subs).subs(graph_subs))
        for symbol, value in bracket_subs.items()
    }
    expected_u = sp.cancel(
        (sp.Integer(475515904) - sp.Integer(109350) * X)
        / sp.Integer(200475)
    )
    expected_w = sp.cancel(
        -(
            sp.Integer(4343625) * P**2
            - sp.Integer(17172000) * X
            + sp.Integer(143826305024)
        )
        / sp.Integer(4009500)
    )
    check_zero(response[R8.U] - expected_u, "beta-zero U response")
    check_zero(response[R8.W] - expected_w, "beta-zero W response")
    check_zero(response[R8.Z], "beta-zero Z graph")
    check_zero(terminal_subs[R8.zeta3] + 3 * P / 2, "zeta graph")

    # II. The row-nine bracket fixes alpha and no other source coordinate.
    check_zero(R9.E9_GATE.subs(graph_subs), "displayed alpha solves E9")
    check(
        sp.diff(R9.E9_GATE, alpha) == -sp.Integer(511211250) * P,
        "E9 alpha coefficient",
    )
    check(QP.degree() == 6 and Q.degree() == 3, "declared quotient degrees")

    g9 = sp.expand(R8.tcoeff(R8.H, 9))
    expected_g9 = (
        (20 * R8.U + 10 * R8.W + 4 * R8.Z) * x**6
        + (10 * alpha + 6 * beta) * x**7
        + (5 * R8.upsilon5 + 4 * X) * x**8
        + (eta + R8.zeta3) * x**9
    )
    check_zero(g9 - expected_g9, "literal G9")
    fixed_a, fixed_c, old_tangent, _ = R9.solve_row_eight_for_g9(
        arows,
        crows,
        bracket_subs,
        terminal_subs,
        theta8_symbols,
        g9,
    )
    check((old_tangent.rows, old_tangent.cols, old_tangent.rank()) == (10, 7, 7), "G9 selects old terminal A7")
    fixed_a = [
        sp.cancel(row.subs(bracket_subs).subs(terminal_subs).subs(graph_subs))
        for row in fixed_a
    ]
    fixed_c = [
        sp.cancel(row.subs(bracket_subs).subs(terminal_subs).subs(graph_subs))
        for row in fixed_c
    ]
    g9_graph = sp.cancel(
        g9.subs(bracket_subs).subs(terminal_subs).subs(graph_subs)
    )
    check_zero(R8.predicted_G(9, fixed_a, fixed_c) - g9_graph, "G9 bracket graph")

    # Add row nine and impose both exact row-nine projected modules.
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
    anull9 = amatrix9.T.nullspace()
    cnull9 = cmatrix9.T.nullspace()
    check(
        (len(acoords9), amatrix9.cols, amatrix9.rank(), len(anull9))
        == (75, 160, 59, 16),
        "row-nine P2 universe",
    )
    check(
        (len(ccoords9), cmatrix9.cols, cmatrix9.rank(), len(cnull9))
        == (85, 251, 73, 12),
        "row-nine P3 universe",
    )
    avec9 = sp.Matrix([R8.xcoeff(row9_a[n], r) for n, r in acoords9])
    cvec9 = sp.Matrix([R8.xcoeff(row9_c[n], r) for n, r in ccoords9])
    residuals9 = [sp.expand((row.T * avec9)[0]) for row in anull9]
    residuals9 += [sp.expand((row.T * cvec9)[0]) for row in cnull9]
    depth9_matrix, depth9_rhs = sp.linear_eq_to_matrix(residuals9, theta9_symbols)
    check(depth9_matrix.rank() == 3, "row-nine terminal rank three")
    check(depth9_matrix.rref()[1] == (7, 8, 9), "row-nine terminal pivots")
    depth9_rows = tuple(depth9_matrix.T.rref()[1])
    depth9_pivots = (7, 8, 9)
    depth9_free_columns = tuple(range(7))
    depth9_solution = selected_solution(
        depth9_matrix,
        depth9_rhs,
        depth9_rows,
        depth9_pivots,
        theta9_symbols[:7],
        depth9_free_columns,
    )
    depth9_subs = {
        theta9_symbols[depth9_pivots[j]]: sp.cancel(depth9_solution[j])
        for j in range(3)
    }
    check(
        all(exact_zero(value.subs(depth9_subs)) for value in residuals9),
        "all row-nine depth equations",
    )
    check(depth9_matrix.row_join(depth9_rhs).rank() == 3, "row-nine augmented rank three")
    check(10 - depth9_matrix.rank() == 7, "row-nine affine fibre A7")

    # A rational row-nine control survives every strict source condition.
    hostile_subs = {P: 1, X: 0}
    check(
        ETA_GRAPH.subs(hostile_subs) == sp.Rational(12505944329368, 195463125),
        "rational row-nine eta",
    )
    check(
        ALPHA_GRAPH.subs(hostile_subs)
        == -sp.Rational(1329703923433112135272353229, 649499164991015625),
        "rational row-nine alpha",
    )
    check(expected_u.subs(hostile_subs) != 0, "rational row-nine U gate")
    check(expected_w.subs(hostile_subs) != 0, "rational row-nine W gate")
    check((-3 * P / 2).subs(hostile_subs) != 0, "rational row-nine zeta gate")

    # III. Impose row-nine depth first; seven free directions meet G10.
    row9_a_depth = fixed_a + [sp.cancel(a9.subs(depth9_subs))]
    row9_c_depth = fixed_c + [sp.cancel(c9.subs(depth9_subs))]
    free9 = theta9_symbols[:7]
    g10 = sp.expand(R8.tcoeff(R8.H, 10))
    expected_g10 = (
        (15 * R8.U + 10 * R8.W + 6 * R8.Z) * x**8
        + (5 * alpha + 4 * beta) * x**9
        + (R8.upsilon5 + X) * x**10
    )
    check_zero(g10 - expected_g10, "literal G10")
    g10_graph = sp.cancel(
        g10.subs(bracket_subs).subs(terminal_subs).subs(graph_subs)
    )
    difference10 = sp.expand(
        g10_graph - R8.predicted_G(10, row9_a_depth, row9_c_depth)
    )
    equations10 = [R8.xcoeff(difference10, j) for j in range(11)]
    g10_matrix, g10_rhs = sp.linear_eq_to_matrix(equations10, free9)
    check(g10_matrix.rank() == 7, "G10 free-row-nine rank seven")
    g10_rows = (0, 1, 2, 3, 4, 5, 7)
    g10_columns = tuple(range(7))
    g10_solution = selected_solution(
        g10_matrix,
        g10_rhs,
        g10_rows,
        g10_columns,
    )
    g10_subs = {free9[j]: sp.cancel(g10_solution[j]) for j in range(7)}
    g10_residuals = [sp.factor(sp.cancel(value.subs(g10_subs))) for value in equations10]
    check(
        exact_zero(g10_residuals[8] + sp.Rational(4, 61875) * (10125 * X + 1928704)),
        "G10 xi residual",
    )
    check_zero(
        g10_residuals[6].subs(X, X0)
        - sp.Rational(2, 94585080322265625) * QP.as_expr() / P**4,
        "G10 q residual",
    )
    check(
        all(value == 0 or index in (6, 8) for index, value in enumerate(g10_residuals)),
        "only two G10 residuals",
    )
    check(
        all(
            reduce_q(value.subs(X, X0)) == 0
            for value in g10_residuals
        ),
        "G10 exact on six-point quotient",
    )

    # Select the seven row-nine directions, set X=X0, and retain q off shell.
    row9_a_selected = [
        sp.cancel(row.subs(g10_subs).subs(X, X0)) for row in row9_a_depth
    ]
    row9_c_selected = [
        sp.cancel(row.subs(g10_subs).subs(X, X0)) for row in row9_c_depth
    ]
    check(
        all(
            reduce_q(value.subs(depth9_subs).subs(g10_subs).subs(X, X0)) == 0
            for value in residuals9
        ),
        "row-nine depth retained on q",
    )

    # Six-sheet/open-gate audit in the quotient coordinate Y=P^2.
    check(Q.primitive()[0] == 1, "q primitive")
    check(sp.gcd(Q, Q.diff()) == sp.Poly(1, Y, domain=sp.QQ), "q squarefree")
    check(Q.eval(0) != 0, "q has no zero Y root")
    check(sp.gcd(QP, QP.diff()) == sp.Poly(1, P, domain=sp.QQ), "six Phi sheets distinct")
    u_x0 = sp.factor(expected_u.subs(X, X0))
    w_x0 = sp.factor(expected_w.subs(X, X0))
    lambda_x0 = sp.factor((expected_u + expected_w).subs(X, X0))
    check(u_x0 == sp.Rational(225611776, 91125), "six-sheet U gate")
    check_zero(w_x0 + sp.Rational(13, 40500) * (3375 * P**2 + 114294784), "six-sheet W formula")
    check_zero(lambda_x0 + sp.Rational(13, 364500) * (30375 * P**2 + 959234048), "six-sheet Lambda formula")
    check(
        sp.gcd(Q, sp.Poly(3375 * Y + 114294784, Y, domain=sp.QQ))
        == sp.Poly(1, Y, domain=sp.QQ),
        "q avoids W zero",
    )
    check(
        sp.gcd(Q, sp.Poly(30375 * Y + 959234048, Y, domain=sp.QQ))
        == sp.Poly(1, Y, domain=sp.QQ),
        "q avoids Lambda zero",
    )
    eta_x0 = sp.factor(ETA_GRAPH.subs(X, X0))
    alpha_x0 = sp.factor(ALPHA_GRAPH.subs(X, X0))
    check_zero(eta_x0.subs(P, -P) + eta_x0, "eta sign sheet")
    check_zero(alpha_x0.subs(P, -P) + alpha_x0, "alpha sign sheet")
    check_zero((-3 * P / 2).subs(P, -P) - 3 * P / 2, "zeta sign sheet")
    check_zero(w_x0.subs(P, -P) - w_x0, "W sign invariant")

    # IV. Add row ten and rebuild the two projected depth modules separately.
    b10 = R8.B_row(10, row9_a_selected, row9_c_selected)
    a10base, c10base = R8.particular_row(10, b10)
    theta10_symbols = tuple(sp.symbols("theta10_0:11"))
    theta10 = sum(theta10_symbols[j] * x**j for j in range(11))
    a10 = sp.expand(a10base + theta10 * sp.diff(R8.A0, x))
    c10 = sp.expand(c10base + theta10 * sp.diff(R8.C0, x))
    row10_a = row9_a_selected + [a10]
    row10_c = row9_c_selected + [c10]

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
    amatrix_terminal, arhs_terminal = sp.linear_eq_to_matrix(
        aresiduals10, theta10_symbols
    )
    cmatrix_terminal, crhs_terminal = sp.linear_eq_to_matrix(
        cresiduals10, theta10_symbols
    )
    joint_terminal = amatrix_terminal.col_join(cmatrix_terminal)
    check(amatrix_terminal.rank() == 3, "standalone P2 terminal rank three")
    check(cmatrix_terminal.rank() == 3, "standalone P3 terminal rank three")
    check(joint_terminal.rank() == 3, "joint terminal coefficient rank three")
    check(amatrix_terminal.rref()[1] == (8, 9, 10), "P2 terminal pivots")
    check(cmatrix_terminal.rref()[1] == (8, 9, 10), "P3 terminal pivots")

    terminal_pivots = (8, 9, 10)
    arows_terminal = tuple(amatrix_terminal.T.rref()[1])
    crows_terminal = tuple(cmatrix_terminal.T.rref()[1])
    asolution10 = selected_solution(
        amatrix_terminal,
        arhs_terminal,
        arows_terminal,
        terminal_pivots,
    )
    csolution10 = selected_solution(
        cmatrix_terminal,
        crhs_terminal,
        crows_terminal,
        terminal_pivots,
    )
    free10_zero = {symbol: 0 for symbol in theta10_symbols[:8]}
    asubs10 = {
        **free10_zero,
        **{
            theta10_symbols[terminal_pivots[j]]: sp.cancel(asolution10[j])
            for j in range(3)
        },
    }
    csubs10 = {
        **free10_zero,
        **{
            theta10_symbols[terminal_pivots[j]]: sp.cancel(csolution10[j])
            for j in range(3)
        },
    }
    p2_after_p2 = [reduce_q(value.subs(asubs10)) for value in aresiduals10]
    p3_after_p3 = [reduce_q(value.subs(csubs10)) for value in cresiduals10]
    check(all(value == 0 for value in p2_after_p2), "P2 standalone augmented rank three")
    check(all(value == 0 for value in p3_after_p3), "P3 standalone augmented rank three")

    # The primitive tetrahedral P3 row, evaluated after P2's selection.
    tetrahedral_terms = {
        (5, 0): 56,
        (6, 2): -35,
        (7, 4): 20,
        (8, 6): -10,
        (9, 8): 4,
        (10, 10): -1,
    }
    tetrahedral_vector = sp.zeros(len(ccoords10), 1)
    for coordinate, coefficient in tetrahedral_terms.items():
        tetrahedral_vector[ccoords10.index(coordinate)] = coefficient
    check(
        (tetrahedral_vector.T * cmatrix10) == sp.zeros(1, cmatrix10.cols),
        "tetrahedral row annihilates pi10(P3)",
    )
    check(
        tuple(abs(value) for value in tetrahedral_terms.values())
        == (56, 35, 20, 10, 4, 1),
        "tetrahedral coefficient sequence",
    )
    check(
        tuple(sp.binomial(m, 3) for m in range(8, 2, -1))
        == (56, 35, 20, 10, 4, 1),
        "third-simplex address",
    )
    check(
        primitive_integer_functional(tetrahedral_vector)
        == tuple(int(value) for value in tetrahedral_vector),
        "tetrahedral row primitive",
    )

    row10_c_after_p2 = [row.subs(asubs10) for row in row10_c]
    h_c = sp.expand(
        sum(
            coefficient * R8.xcoeff(row10_c_after_p2[n], r)
            for (n, r), coefficient in tetrahedral_terms.items()
        )
    )
    h_claimed = -P * DP.as_expr() / H_DENOMINATOR
    h_correction = QP.as_expr() / (Q_CORRECTION_DENOMINATOR * P)
    check_zero(h_c - h_claimed - h_correction, "off-q tetrahedral identity")
    check_zero(reduce_q(h_c) - h_claimed, "on-q tetrahedral value")

    # Standalone ranks are 3=augmented rank; joint augmentation is rank four.
    h_mod_q = reduce_q(h_c)
    check(h_mod_q != 0, "joint P2/P3 compatibility fails")
    p3_after_p2 = [reduce_q(value.subs(asubs10)) for value in cresiduals10]
    check(any(value != 0 for value in p3_after_p2), "P2 selection violates P3")
    p2_rank = amatrix_terminal.rank()
    p3_rank = cmatrix_terminal.rank()
    joint_rank = joint_terminal.rank()
    p2_augmented_rank = p2_rank if all(value == 0 for value in p2_after_p2) else p2_rank + 1
    p3_augmented_rank = p3_rank if all(value == 0 for value in p3_after_p3) else p3_rank + 1
    joint_augmented_rank = joint_rank + 1 if h_mod_q != 0 else joint_rank
    check((p2_rank, p2_augmented_rank) == (3, 3), "P2 rank equals augmented rank")
    check((p3_rank, p3_augmented_rank) == (3, 3), "P3 rank equals augmented rank")
    check((joint_rank, joint_augmented_rank) == (3, 4), "joint rank/augmented-rank jump")

    # Exact coprimality and Euclidean stopping witnesses.
    check(sp.gcd(Q, D) == sp.Poly(1, Y, domain=sp.QQ), "gcd(q,d)=1")
    q_remainder = sp.factor(Q.rem(D).as_expr())
    expected_q_remainder = sp.Integer(490103134214955204608) * (
        sp.Integer(1805625) * Y - sp.Integer(19708903424)
    )
    check_zero(q_remainder - expected_q_remainder, "q mod d")
    linear = sp.Poly(1805625 * Y - 19708903424, Y, domain=sp.QQ)
    check(
        D.rem(linear).as_expr()
        == -sp.Rational(138855112274442930681014124544, 57245),
        "d mod linear is a unit",
    )

    # V. Good-prime control: bracket survives while joint depth fails.
    prime = 29
    f29_subs = {P: 1, X: 7}
    check(residue(X0, prime, {}) == 7, "F29 xi condition")
    check(residue(ETA_GRAPH, prime, f29_subs) == 8, "F29 eta graph")
    check(residue(ALPHA_GRAPH, prime, f29_subs) == 8, "F29 alpha graph")
    check(residue(expected_u, prime, f29_subs) == 10, "F29 U gate")
    check(residue(expected_w, prime, f29_subs) == 24, "F29 W gate")
    check(residue(-3 * P / 2, prime, f29_subs) == 13, "F29 zeta gate")
    check(int(Q.eval(1)) % prime == 0, "F29 q survivor")
    check(int(D.eval(1)) % prime == 23, "F29 depth hostile")
    check(residue(h_claimed.subs(P, 1), prime, {}) != 0, "F29 tetrahedral obstruction")

    print("THM-4361 PRIMARY: BETA-ZERO ROW-TEN JOINT DEPTH EXTINCTION PASS")
    print(f"CHECKS={CHECKS}")
    print(
        "ROW9_DEPTH P2=75x160/rank59/null16 "
        "P3=85x251/rank73/null12 terminal_rank=3 fibre=A7"
    )
    print(f"ROW10_BRACKET xi={X0} q_degree={Q.degree()} six_sign_sheets=6")
    print(
        "ROW10_DEPTH P2=88x193/rank68/null20 terminal=3/aug3; "
        "P3=99x304/rank83/null16 terminal=3/aug3; joint=3/aug4"
    )
    print("TETRAHEDRAL=(56,-35,20,-10,4,-1)=alternating_binomials_m_choose_3")
    print(
        "H_C_OFF_Q=-Phi*d(Phi^2)/132327846238037905244160"
        "+q(Phi^2)/(248114711696321072332800000*Phi)"
    )
    print("H_C_ON_Q=-Phi*d(Phi^2)/132327846238037905244160")
    print("COPRIME gcd(q,qprime)=gcd(q,d)=gcd(q,W)=1; q(0)!=0")
    print("SIGN_ADDRESS=three_Phi_squared_roots_times_two_sign_sheets")
    print("F29_CONTROL=(Phi,xi,eta,alpha,beta)=(1,7,8,8,0); q=0; d=23")
    print(
        "SCOPE finite row10 joint projected P2/P3 extinction only; no all-row "
        "B2 lift, termination, seam entry, Keller pair, JC2, or DC2"
    )


if __name__ == "__main__":
    main()
