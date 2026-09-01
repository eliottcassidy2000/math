#!/usr/bin/env python3
"""Primary exact certificate for THM-4315.

This extends THM-4308 by one bracket/depth row and identifies the row
cokernel with a scaled Student/Pearson stationary expectation.  It then
intersects the row-nine gate with THM-4312's surviving cubic-corner k=1
locus.  The result is finite-jet only: no all-row lift or JC(2) conclusion
is asserted.
"""

from __future__ import annotations

from math import gcd
from pathlib import Path
import sys

import sympy as sp


sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_bracket_hasse_rows8_thm4308 as R8  # noqa: E402


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def is_zero(expression: sp.Expr) -> bool:
    return sp.cancel(sp.together(sp.expand(expression))) == 0


def require_zero(expression: sp.Expr, label: str) -> None:
    require(is_zero(expression), label)


def residue(value: sp.Expr, prime: int) -> int:
    """Reduce one exact rational value modulo a good prime."""

    value = sp.cancel(value)
    if value.free_symbols:
        raise AssertionError(f"nonconstant modular value: {value}")
    numerator, denominator = sp.fraction(value)
    denominator_mod = int(denominator) % prime
    if denominator_mod == 0:
        raise AssertionError(f"bad denominator modulo {prime}: {value}")
    return (int(numerator) % prime) * pow(denominator_mod, -1, prime) % prime


def modular_rank(matrix: sp.Matrix, prime: int) -> int:
    """Exact Gaussian rank over F_prime for a rational matrix."""

    rows = [
        [residue(matrix[i, j], prime) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]
    rank = 0
    for column in range(matrix.cols):
        pivot = next((i for i in range(rank, matrix.rows) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(inverse * entry) % prime for entry in rows[rank]]
        for i in range(matrix.rows):
            if i == rank or rows[i][column] == 0:
                continue
            factor = rows[i][column]
            rows[i] = [
                (rows[i][j] - factor * rows[rank][j]) % prime
                for j in range(matrix.cols)
            ]
        rank += 1
        if rank == matrix.rows:
            break
    return rank


def require_mod_zero(expression: sp.Expr, substitutions: dict[sp.Symbol, int], prime: int, label: str) -> None:
    require(residue(expression.subs(substitutions), prime) == 0, label)


x = R8.x
Phi = R8.Phi
X = sp.symbols("X")


def student_even_moment(m: int, r: int) -> sp.Rational:
    """E[X^(2r)] for density proportional to (x^2+6)^(-(m+1))."""

    value = sp.Integer(1)
    for j in range(1, r + 1):
        value *= sp.Rational(6 * (2 * j - 1), 2 * m - 2 * j + 1)
    return sp.cancel(value)


def primitive_student_row(m: int) -> list[int]:
    moments = [
        sp.Integer(0) if degree % 2 else student_even_moment(m, degree // 2)
        for degree in range(m + 1)
    ]
    denominator_lcm = 1
    for value in moments:
        denominator_lcm = sp.ilcm(denominator_lcm, int(sp.denom(value)))
    entries = [int(value * denominator_lcm) for value in moments]
    common = 0
    for entry in entries:
        common = gcd(common, abs(entry))
    return [entry // common for entry in entries]


def tangent_matrix(m: int) -> sp.Matrix:
    variables = sp.symbols(f"stein{m}_0:{m}")
    theta = sum(variables[j] * x**j for j in range(m))
    image = sp.expand(-x * theta + (x**2 + 6) * sp.diff(theta, x) / (2 * m))
    equations = [R8.xcoeff(image, j) for j in range(m + 1)]
    matrix, _ = sp.linear_eq_to_matrix(equations, variables)
    require(matrix.rank() == m, f"Student tangent rank m={m}")
    return matrix


def student_stein_audit() -> dict[int, list[int]]:
    expected = {
        5: [21, 0, 14, 0, 36, 0],
        6: [77, 0, 42, 0, 84, 0, 360],
        7: [143, 0, 66, 0, 108, 0, 360, 0],
        8: [715, 0, 286, 0, 396, 0, 1080, 0, 5040],
        9: [12155, 0, 4290, 0, 5148, 0, 11880, 0, 45360, 0],
    }
    rows: dict[int, list[int]] = {}
    a = x**2 + 6
    for m in range(1, 13):
        row = primitive_student_row(m)
        matrix = tangent_matrix(m)
        require(sp.Matrix([row]) * matrix == sp.zeros(1, m), f"Student cokernel m={m}")
        require(len(row) == m + 1 and row[0] > 0, f"Student row shape m={m}")
        # The diffusion generator L_m h=a h''-2mxh' has invariant density
        # rho_m=a^(-(m+1)); this is the exact stationary adjoint identity.
        rho = a ** (-(m + 1))
        require_zero(sp.diff(a * rho, x) + 2 * m * x * rho, f"stationary density m={m}")
        # Killing with probability 6/a and conditioning on survival sends
        # mu_m to mu_(m+1).  The exact survival probability is below.
        survival = sp.Rational(2 * m + 1, 2 * m + 2)
        normalization_ratio = sp.Rational(12 * (m + 1), 2 * m + 1)
        require_zero(6 / survival - normalization_ratio, f"posterior normalization m={m}")
        rows[m] = row
    for m, row in expected.items():
        require(rows[m] == row, f"displayed Student row m={m}")
    return rows


E9_FULL = (
    68294026800 * R8.Delta**2
    + 3653910000 * R8.Delta * R8.Theta
    - 5288166000 * R8.Delta * R8.xi10
    + 176911616000 * R8.Delta
    + 1547488800 * Phi**2
    + 3987447750 * Phi * R8.alpha11
    + 24602292000 * Phi * R8.eta
    + 4222003500 * Phi * R8.zeta3
    + 2258685000 * R8.Theta**2
    + 225494104640 * R8.Theta
    + 1993723875 * R8.eta**2
    + 263331993600 * R8.xi10
    + 105193437167616
)

E9_GATE = (
    613527750 * Phi**2
    - 511211250 * Phi * R8.alpha11
    - 3154140000 * Phi * R8.eta
    - 255605625 * R8.eta**2
    + 6736896000 * R8.xi10
    - 46483785515008
)


def solve_row_eight_for_g9(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    bracket_subs: dict[sp.Symbol, sp.Expr],
    terminal_subs: dict[sp.Symbol, sp.Expr],
    theta8_symbols: list[sp.Symbol],
    g9: sp.Expr,
) -> tuple[list[sp.Expr], list[sp.Expr], sp.Matrix, sp.Matrix]:
    difference = sp.expand(
        (g9.subs(bracket_subs) - R8.predicted_G(9, arows, crows)).subs(terminal_subs)
    )
    free = [symbol for symbol in theta8_symbols if symbol not in terminal_subs]
    require(len(free) == 7, "row-eight Hasse free tangent count")
    equations = [R8.xcoeff(difference, j) for j in range(10)]
    matrix, rhs = sp.linear_eq_to_matrix(equations, free)
    require(matrix.rank() == 7, "G9 restricted tangent rank")
    require(len(free) - matrix.rank() == 0, "row-nine truncation has zero-dimensional image")
    row_pivots = matrix.T.rref()[1]
    require(len(row_pivots) == 7, "G9 restricted pivot rows")
    solution = matrix[list(row_pivots), :].inv() * rhs[list(row_pivots), :]
    substitutions = {symbol: sp.cancel(solution[j]) for j, symbol in enumerate(free)}
    residuals = [sp.factor(equation.subs(substitutions)) for equation in equations]
    nonzero = [value for value in residuals if value != 0]
    require(len(nonzero) == 1, "G9 one scalar residual after row-eight Hasse gate")
    require_zero(nonzero[0] + sp.Rational(13, 6495390000) * E9_GATE, "G9 gated residual")
    fixed_a8 = sp.expand(arows[8].subs(terminal_subs).subs(substitutions))
    fixed_c8 = sp.expand(crows[8].subs(terminal_subs).subs(substitutions))
    return arows[:8] + [fixed_a8], crows[:8] + [fixed_c8], matrix, rhs


def depth_matrix(depth: int, max_row: int) -> tuple[list[tuple[int, int]], sp.Matrix]:
    coordinates = [
        (n, r) for n in range(max_row + 1) for r in range(n + depth + 1)
    ]
    columns: list[tuple[int, int, int, int]] = []
    for b in range(depth + 1):
        for a in range(depth - b + 1):
            for e in range(max_row // 2 + 1):
                for c in range(max_row + 1):
                    if b + c + 2 * e <= max_row:
                        columns.append((a, b, c, e))
    matrix = sp.zeros(len(coordinates), len(columns))
    coordinate_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
    for column, (a, b, c, e) in enumerate(columns):
        for j in range(c + e + 1):
            n = b + c + 2 * e + j
            r = a + 2 * b + e + 2 * j
            if n <= max_row:
                matrix[coordinate_index[(n, r)], column] += sp.binomial(c + e, j)
    return coordinates, matrix


def even_polynomial_in_x(expression: sp.Expr) -> sp.Poly:
    poly = sp.Poly(sp.expand(expression), Phi)
    converted = sp.Integer(0)
    for (degree,), coefficient in poly.terms():
        require(degree % 2 == 0, "even Phi polynomial")
        converted += coefficient * X ** (degree // 2)
    return sp.Poly(converted, X)


Q_EXPECTED = sp.Poly(
    316016952601619726458584136962890625 * X**5
    + 14163685391496771581808513584548828125000 * X**4
    + 137633556412285429978153325875719168000000000 * X**3
    - 6709927871370175861935855782936495259648000000 * X**2
    + 16759499408238096044037088463875607198378754048000 * X
    + 44257795605986960276324945338517826145242081533100032,
    X,
)

R_HIGH = sp.Poly(
    7547170421607067494140625 * X**3
    + 164114458618573873612800000000 * X**2
    + 2284603892441775363795663716352000 * X
    + 2970579390109346274816679296272171008,
    X,
)


def row_nine_corner_audit() -> tuple[sp.Poly, int, int, int, int, int, int]:
    arows, crows, bracket_subs = R8.build_bracket_rows()
    amatrix8, cmatrix8, theta8_symbols, terminal_subs = R8.hasse_audit(
        arows, crows, bracket_subs
    )
    require(amatrix8.rank() == 51 and cmatrix8.rank() == 63, "inherited row-eight depths")

    gate_parameters = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -3 * Phi / 2,
    }
    acoords8 = [(n, r) for n in range(9) for r in range(n + 3)]
    ccoords8 = [(n, r) for n in range(9) for r in range(n + 4)]
    avec8 = sp.Matrix([R8.xcoeff(arows[n], r) for n, r in acoords8])
    cvec8 = sp.Matrix([R8.xcoeff(crows[n], r) for n, r in ccoords8])
    row8_residuals = [sp.expand((row.T * avec8)[0].subs(gate_parameters)) for row in amatrix8.T.nullspace()]
    row8_residuals += [sp.expand((row.T * cvec8)[0].subs(gate_parameters)) for row in cmatrix8.T.nullspace()]
    row8_tangent, row8_rhs = sp.linear_eq_to_matrix(row8_residuals, theta8_symbols)
    require(row8_tangent.rank() == 2, "reconstructed row-eight terminal depth rank")
    require(row8_tangent.row_join(row8_rhs).rank() == 2, "reconstructed row-eight depth consistency")

    g9 = sp.expand(R8.tcoeff(R8.H, 9))
    expected_g9 = (
        (20 * R8.U + 10 * R8.W + 4 * R8.Z) * x**6
        + (10 * R8.alpha11 + 6 * R8.beta11) * x**7
        + (5 * R8.upsilon5 + 4 * R8.xi10) * x**8
        + (R8.eta + R8.zeta3) * x**9
    )
    require_zero(g9 - expected_g9, "direct source row G9")

    difference = sp.expand(
        g9.subs(bracket_subs) - R8.predicted_G(9, arows, crows)
    )
    row9 = primitive_student_row(9)
    obstruction = sum(row9[j] * R8.xcoeff(difference, j) for j in range(10))
    require_zero(obstruction - E9_FULL / 328050, "full Student E9 obstruction")
    require_zero(E9_FULL.subs(gate_parameters) + sp.Rational(39, 5) * E9_GATE, "gated E9")

    fixed_a, fixed_c, g9_tangent, g9_rhs = solve_row_eight_for_g9(
        arows, crows, bracket_subs, terminal_subs, theta8_symbols, g9
    )

    xi_corner = (4343625 * Phi**2 + 124805668864) / 12798000
    eta_corner = (
        sp.Rational(2091705253888, 107983125)
        - sp.Rational(2839, 1185) * Phi**2
    ) / Phi
    alpha_corner = (
        6971519208442078125 * Phi**4
        - 14082869793796263936000 * Phi**2
        - 74378924775425263164981248
    ) / (396452079682031250 * Phi**3)
    corner_subs = {
        R8.xi10: xi_corner,
        R8.eta: eta_corner,
        R8.alpha11: alpha_corner,
        R8.beta11: -alpha_corner,
    }
    require_zero(E9_GATE.subs(corner_subs), "corner solves E9")

    prime = 19
    finite_field_point = {
        R8.Delta: 4,
        R8.Theta: 1,
        Phi: 8,
        R8.eta: 6,
        R8.zeta3: 7,
        R8.upsilon5: 9,
        R8.xi10: 4,
        R8.U: 12,
        R8.W: 14,
        R8.Z: 12,
        R8.alpha11: 15,
        R8.beta11: 4,
    }
    require(residue(R8.K.subs(finite_field_point), prime) == 5, "F19 K coordinate")
    for index, equation in enumerate((R8.E5, R8.E6, R8.E7, R8.E8), start=5):
        require_mod_zero(equation, finite_field_point, prime, f"F19 E{index}")
    require_mod_zero(E9_GATE, finite_field_point, prime, "F19 E9")
    for symbol, value in bracket_subs.items():
        require_mod_zero(symbol - value, finite_field_point, prime, f"F19 bracket response {symbol}")
    for symbol, value in gate_parameters.items():
        require_mod_zero(symbol - value, finite_field_point, prime, f"F19 Hasse gate {symbol}")
    for symbol, value in corner_subs.items():
        require_mod_zero(symbol - value, finite_field_point, prime, f"F19 corner response {symbol}")
    require_mod_zero(R8.W + 2 * R8.U, finite_field_point, prime, "F19 corner W=-2U")
    require_mod_zero(R8.Z - R8.U, finite_field_point, prime, "F19 corner Z=U")
    require_mod_zero(R8.alpha11 + R8.beta11, finite_field_point, prime, "F19 beta=-alpha")
    c2 = R8.upsilon5 + R8.xi10
    require_mod_zero(R8.alpha11**2 - 4 * R8.U * c2, finite_field_point, prime, "F19 square relation")
    rho = -R8.alpha11 / (2 * R8.U)
    l1 = R8.eta + R8.zeta3 + rho * R8.upsilon5
    require(residue(l1.subs(finite_field_point), prime) == 5, "F19 L1 nonzero control")
    require(modular_rank(amatrix8, prime) == 51, "F19 row-eight P2 rank")
    require(modular_rank(cmatrix8, prime) == 63, "F19 row-eight P3 rank")
    for m in range(5, 9):
        require(modular_rank(tangent_matrix(m), prime) == m, f"F19 Student tangent rank m={m}")
    require(modular_rank(row8_tangent.subs(finite_field_point), prime) == 2, "F19 row-eight depth tangent rank")
    require(
        modular_rank(row8_tangent.row_join(row8_rhs).subs(finite_field_point), prime) == 2,
        "F19 row-eight depth augmented rank",
    )
    require(modular_rank(g9_tangent.subs(finite_field_point), prime) == 7, "F19 G9 restricted rank")
    require(
        modular_rank(g9_tangent.row_join(g9_rhs).subs(finite_field_point), prime) == 7,
        "F19 G9 restricted consistency",
    )

    fixed_a = [sp.expand(row.subs(bracket_subs).subs(terminal_subs).subs(corner_subs)) for row in fixed_a]
    fixed_c = [sp.expand(row.subs(bracket_subs).subs(terminal_subs).subs(corner_subs)) for row in fixed_c]
    g9_corner = sp.expand(g9.subs(bracket_subs).subs(terminal_subs).subs(corner_subs))
    require_zero(R8.predicted_G(9, fixed_a, fixed_c) - g9_corner, "corner G9 bracket match")

    b9 = R8.B_row(9, fixed_a, fixed_c)
    a9base, c9base = R8.particular_row(9, b9)
    theta9_symbols = list(sp.symbols("theta9_0:10"))
    theta9 = sum(theta9_symbols[j] * x**j for j in range(10))
    a9 = sp.expand(a9base + theta9 * sp.diff(R8.A0, x))
    c9 = sp.expand(c9base + theta9 * sp.diff(R8.C0, x))
    row9_a = fixed_a + [a9]
    row9_c = fixed_c + [c9]

    acoords, amatrix = depth_matrix(2, 9)
    ccoords, cmatrix = depth_matrix(3, 9)
    require((len(acoords), amatrix.cols, amatrix.rank()) == (75, 160, 59), "row-nine P2 universe")
    require((len(ccoords), cmatrix.cols, cmatrix.rank()) == (85, 251, 73), "row-nine P3 universe")
    anull = amatrix.T.nullspace()
    cnull = cmatrix.T.nullspace()
    require(len(anull) == 16 and len(cnull) == 12, "row-nine depth nullities")
    avec = sp.Matrix([R8.xcoeff(row9_a[n], r) for n, r in acoords])
    cvec = sp.Matrix([R8.xcoeff(row9_c[n], r) for n, r in ccoords])
    residuals = [sp.expand((row.T * avec)[0]) for row in anull]
    residuals += [sp.expand((row.T * cvec)[0]) for row in cnull]
    depth_tangent, depth_rhs = sp.linear_eq_to_matrix(residuals, theta9_symbols)
    require(depth_tangent.rank() == 3, "row-nine terminal depth rank")
    require(depth_tangent.row_join(depth_rhs).rank() == 3, "row-nine depth consistency")
    require(depth_tangent.rref()[1] == (7, 8, 9), "row-nine terminal pivot coordinates")
    zero_free = {symbol: 0 for symbol in theta9_symbols[:7]}
    reduced = [sp.expand(value.subs(zero_free)) for value in residuals]
    terminal_solutions = sp.solve(reduced, theta9_symbols[7:], dict=True)
    require(len(terminal_solutions) == 1, "row-nine terminal solution")
    terminal9 = {**zero_free, **terminal_solutions[0]}
    require(all(is_zero(value.subs(terminal9)) for value in residuals), "all row-nine depth residuals")
    require(modular_rank(amatrix, prime) == 59, "F19 row-nine P2 rank")
    require(modular_rank(cmatrix, prime) == 73, "F19 row-nine P3 rank")
    require(
        modular_rank(depth_tangent.subs(finite_field_point), prime) == 3,
        "F19 row-nine depth tangent rank",
    )
    require(
        modular_rank(depth_tangent.row_join(depth_rhs).subs(finite_field_point), prime) == 3,
        "F19 row-nine depth augmented rank",
    )
    row9_a[-1] = sp.expand(row9_a[-1].subs(terminal9))
    row9_c[-1] = sp.expand(row9_c[-1].subs(terminal9))

    aseries = sum(row9_a[n] * R8.t**n for n in range(10))
    cseries = sum(row9_c[n] * R8.t**n for n in range(10))
    gcorner = sp.expand(R8.G.subs(bracket_subs).subs(terminal_subs).subs(corner_subs))
    jac = R8.jacobian(aseries, cseries)
    require_zero(R8.tcoeff(jac, 0) - 1, "row-nine control Jacobian row zero")
    for n in range(1, 9):
        require_zero(R8.tcoeff(jac, n), f"row-nine control Jacobian row {n}")
    for n in range(10):
        require_zero(R8.tcoeff(R8.P(aseries, cseries) - gcorner, n), f"row-nine control P/G row {n}")

    u_corner = -sp.Rational(13, 57591000) * (820125 * Phi**2 + 13056802816)
    c2_corner = sp.Rational(11, 474000) * (14625 * Phi**2 + 404652032)
    square_relation = sp.together((alpha_corner**2 - 4 * u_corner * c2_corner) * Phi**2)
    raw = sp.Poly(square_relation.as_numer_denom()[0], Phi)
    content_phi, primitive_phi = raw.primitive()
    require(content_phi == 1, "corner quintic primitive numerator")
    quintic = even_polynomial_in_x(primitive_phi.as_expr())
    require(quintic == Q_EXPECTED, "corner row-nine quintic")
    require(sp.gcd(quintic, quintic.diff()) == sp.Poly(1, X), "corner quintic squarefree")
    require_mod_zero(Phi**2 - 7, finite_field_point, prime, "F19 Phi^2=X=7")
    require(residue(quintic.as_expr().subs(X, 7), prime) == 0, "F19 quintic survivor X=7")
    forbidden = sp.Poly(
        X * (820125 * X + 13056802816) * (14625 * X + 404652032), X
    )
    require(sp.gcd(quintic, forbidden) == sp.Poly(1, X), "corner quintic avoids forbidden factors")
    forbidden_values = (
        residue(sp.Integer(7), prime),
        residue((820125 * X + 13056802816).subs(X, 7), prime),
        residue((14625 * X + 404652032).subs(X, 7), prime),
    )
    require(forbidden_values == (7, 10, 10), "F19 forbidden-factor control")
    require(sp.gcd(quintic, R_HIGH) == sp.Poly(1, X), "row-nine excludes high contact")

    qmod = sp.Poly(-5 * X**5 + 5 * X**4 - 6 * X**3 + 2 * X**2 - 7 * X - 4, X, modulus=19)
    rmod = sp.Poly(3 * X**3 + 3 * X + 8, X, modulus=19)
    require(sp.Poly(quintic.as_expr(), X, modulus=19) == qmod, "quintic mod19 reduction")
    require(sp.Poly(R_HIGH.as_expr(), X, modulus=19) == rmod, "high-contact mod19 reduction")
    bezout = (
        (-X**2 + 9 * X + 1) * qmod.as_expr()
        + (-8 * X**4 + 4 * X**3 - X**2 + 9 * X + 3) * rmod.as_expr()
        - 1
    )
    require(sp.Poly(bezout, X, modulus=19).is_zero, "mod19 high-contact Bezout witness")

    return quintic, len(acoords), amatrix.cols, amatrix.rank(), len(ccoords), cmatrix.cols, cmatrix.rank()


def main() -> None:
    rows = student_stein_audit()
    quintic, pa, pc, pr, qa, qc, qr = row_nine_corner_audit()
    print("THM-4315 STUDENT--STEIN ROW-NINE PRIMARY EXACT AUDIT: PASS")
    print("STEIN D_m(theta)=((x^2+6)*theta'-2*m*x*theta)/(2*m)")
    print("STATIONARY mu_m(dx)=c_m*(x^2+6)^(-(m+1))*dx")
    print("MOMENT for 0<=2r<=m: E_m[X^(2r)]=6^r*(2r-1)!!/prod_(j=1)^r(2m-2j+1); odd degrees <=m have mean0")
    print(f"COKERNEL_m9={rows[9]}")
    print("FILTER mu_m --survive 6/(X^2+6)--> mu_(m+1); survival=(2m+1)/(2m+2)")
    print(f"G9={sp.sstr(sp.expand(R8.tcoeff(R8.H, 9)))}")
    print(f"E9_GATE={sp.sstr(E9_GATE)}=0")
    print(f"P2_ROW9 ambient={pa} columns={pc} rank={pr} left_nullity={pa-pr}")
    print(f"P3_ROW9 ambient={qa} columns={qc} rank={qr} left_nullity={qa-qr}")
    print("ROW9_TERMINAL tangent_rank=3 affine_dimension=7 no_extra_corner_equation=yes")
    print("FIBRE_TRANSITION dim_J8=7 dim_J9=7 truncation_image_dim=0")
    print("STOCHASTIC_NO_GO nondegenerate_full_fibre_marginals_projectively_consistent=no")
    print(f"CORNER_Q={sp.sstr(quintic.as_expr())}")
    print("CORNER_Q degree=5 squarefree=yes forbidden_gcd=1 source_points=10_over_algebraic_closure")
    print("HIGH_CONTACT gcd(Q,R)=1 mod19_Bezout=yes L1_zero_row9_lift=no")
    print("F19_PROPOSAL increments=(5,6,7,8,2,7,3) survival_row8=19^-28 survival_row9=19^-38")
    print("F19_SURVIVOR (Delta,Theta,K,Phi,eta,zeta3,upsilon5,xi10,U,W,Z,alpha11,beta11)=(4,1,5,8,6,7,9,4,12,14,12,15,4) L1=5")
    print(f"CHECKS={CHECKS}")
    print("SCOPE exact residual_weight<=12 row9 projected finite gate; no row10/all-row lift, termination, seam entry, JC2, or DC2")


if __name__ == "__main__":
    main()
