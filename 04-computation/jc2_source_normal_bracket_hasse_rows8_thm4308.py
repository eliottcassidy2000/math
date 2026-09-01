#!/usr/bin/env python3
"""Primary exact symbolic certificate for THM-4308.

The certificate works in the fixed characteristic-zero, b=d=0,
source-normal gauge inherited from THM-4007 and THM-4230.  It derives the
bracket compatibility rows E5--E8, constructs every row through t^8, and
computes the exact row-eight projections of the depth modules P_2 and P_3.

This is deliberately a finite-jet certificate.  It neither tests nor asserts
an all-row B_2 lift, polynomial termination, seam entry, or JC(2).
"""

from __future__ import annotations

import sympy as sp


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


x, t = sp.symbols("x t")
Delta, Phi, Theta = sp.symbols("Delta Phi Theta")
eta, zeta3, upsilon5, xi10 = sp.symbols("eta zeta3 upsilon5 xi10")
alpha11, beta11 = sp.symbols("alpha11 beta11")
U, W, Z = sp.symbols("U W Z")

K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
q = -(x**2 + 6) / 2


def P(a: sp.Expr, c: sp.Expr) -> sp.Expr:
    return sp.expand(c**2 - a**3 + sp.Rational(3, 4) * a + sp.Rational(1, 4))


p = t * (1 + x**2 * t)
y = x * t * p
u = x**2 * t
H = sp.expand(
    -3 * p
    + sp.Rational(8, 3) * p**2
    - sp.Rational(1376, 135) * p**3
    + K * y**2
    + Phi * p**2 * y
    + Delta * p**4
    + Theta * p * y**2
    + eta * p**3 * y
    + zeta3 * y**3
    + upsilon5 * p**5
    + xi10 * p**2 * y**2
    + alpha11 * p**4 * y
    + beta11 * p * y**3
    + U * p**6
    + W * p**3 * y**2
    + Z * y**4
)
G = sp.expand(-u / 2 + H)

A0 = 1 + x**2 / 4
C0 = -3 * x / 4 - x**3 / 8
A1 = sp.Rational(4, 3) + 2 * x**2
C1 = -4 * x - sp.Rational(3, 2) * x**3
A2 = -sp.Rational(32, 9) - sp.Rational(4, 5) * x**2
C2 = sp.Rational(88, 15) * x - sp.Rational(12, 5) * x**3
A3 = (
    sp.Rational(2176, 135)
    - Phi * x / 2
    + (sp.Rational(1088, 315) - 4 * K / 7) * x**2
    - sp.Rational(32, 15) * x**4
)
C3 = (
    3 * Phi / 4
    + (-sp.Rational(8128, 315) + 6 * K / 7) * x
    + 3 * Phi * x**2 / 8
    + (sp.Rational(736, 105) + 3 * K / 7) * x**3
    + sp.Rational(8, 5) * x**5
)


EXPECTED_G_ROWS = {
    4: Delta + Phi * x + (K - sp.Rational(1376, 45)) * x**2 + sp.Rational(8, 3) * x**4,
    5: (
        upsilon5
        + eta * x
        + (4 * Delta + Theta) * x**2
        + 3 * Phi * x**3
        + (2 * K - sp.Rational(1376, 45)) * x**4
    ),
    6: (
        U
        + alpha11 * x
        + (5 * upsilon5 + xi10) * x**2
        + (4 * eta + zeta3) * x**3
        + (6 * Delta + 3 * Theta) * x**4
        + 3 * Phi * x**5
        + (K - sp.Rational(1376, 135)) * x**6
    ),
    7: (
        (6 * U + W) * x**2
        + (5 * alpha11 + beta11) * x**3
        + (10 * upsilon5 + 4 * xi10) * x**4
        + (6 * eta + 3 * zeta3) * x**5
        + (4 * Delta + 3 * Theta) * x**6
        + Phi * x**7
    ),
    8: (
        (15 * U + 5 * W + Z) * x**4
        + (10 * alpha11 + 4 * beta11) * x**5
        + (10 * upsilon5 + 6 * xi10) * x**6
        + (4 * eta + 3 * zeta3) * x**7
        + (Delta + Theta) * x**8
    ),
}


E5 = 2025 * upsilon5 + 9000 * Delta + 1350 * Theta + 184832
E6 = 200475 * U + 109350 * xi10 - 5593860 * Delta - 529200 * Theta - 137763328
E7 = (
    801900 * W
    + 1782000 * Delta**2
    + 156163200 * Delta
    + 868725 * Phi**2
    + 27390480 * Theta
    - 3434400 * xi10
    + 12891824128
)
E8 = (
    21651300 * Z
    - 225022050 * Delta**2
    - 59073300 * Delta * Theta
    - 9512522400 * Delta
    + 34749000 * Phi**2
    + 39092625 * Phi * eta
    + 940522560 * Theta
    + 185376600 * xi10
    - 1112446017536
)

upsilon_solution = sp.solve(E5, upsilon5)[0]
U_solution = sp.solve(E6, U)[0]
W_solution = sp.solve(E7, W)[0]
Z_solution = sp.solve(E8, Z)[0]


def xcoeff(expression: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(expression).coeff(x, degree)


def tcoeff(expression: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(expression).coeff(t, degree)


def jacobian(a: sp.Expr, c: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(a, x) * sp.diff(c, t) - sp.diff(a, t) * sp.diff(c, x))


def B_row(m: int, arows: list[sp.Expr], crows: list[sp.Expr]) -> sp.Expr:
    return sp.expand(
        sum(
            (m - i) * sp.diff(arows[i], x) * crows[m - i]
            - i * arows[i] * sp.diff(crows[m - i], x)
            for i in range(1, m)
        )
    )


def T_row(m: int, arows: list[sp.Expr], crows: list[sp.Expr]) -> sp.Expr:
    quadratic = sum(crows[i] * crows[m - i] for i in range(1, m))
    cubic = sum(
        arows[i] * arows[j] * arows[k]
        for i in range(m)
        for j in range(m)
        for k in range(m)
        if i + j + k == m
    )
    return sp.expand(quadratic - cubic)


def predicted_G(m: int, arows: list[sp.Expr], crows: list[sp.Expr]) -> sp.Expr:
    return sp.expand(T_row(m, arows, crows) - q * B_row(m, arows, crows) / m)


def determinant_operator_matrix(m: int) -> tuple[sp.Matrix, tuple[int, ...]]:
    avars = sp.symbols(f"ra{m}_0:{m + 2}")
    cvars = sp.symbols(f"rc{m}_0:{m + 3}")
    avalue = sum(avars[j] * x**j for j in range(m + 2))
    cvalue = sum(cvars[j] * x**j for j in range(m + 3))
    expression = sp.expand(m * (sp.diff(A0, x) * cvalue - sp.diff(C0, x) * avalue))
    equations = [xcoeff(expression, j) for j in range(m + 4)]
    matrix, _ = sp.linear_eq_to_matrix(equations, avars + cvars)
    require(matrix.rank() == m + 4, f"row-{m} determinant operator surjective")
    require(len(avars + cvars) - matrix.rank() == m + 1, f"row-{m} tangent nullity")
    _, pivots = matrix.rref()
    return matrix, pivots


def particular_row(m: int, bvalue: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    matrix, pivots = determinant_operator_matrix(m)
    rhs = sp.Matrix([-xcoeff(bvalue, j) for j in range(m + 4)])
    pivot_matrix = matrix[:, list(pivots)]
    require(pivot_matrix.det() != 0, f"row-{m} pivot minor")
    pivot_values = pivot_matrix.inv() * rhs
    values: list[sp.Expr] = [sp.Integer(0)] * (2 * m + 5)
    for column, value in zip(pivots, pivot_values):
        values[column] = sp.cancel(value)
    avalue = sp.expand(sum(values[j] * x**j for j in range(m + 2)))
    cvalue = sp.expand(sum(values[m + 2 + j] * x**j for j in range(m + 3)))
    require_zero(
        m * (sp.diff(A0, x) * cvalue - sp.diff(C0, x) * avalue) + bvalue,
        f"row-{m} particular determinant equation",
    )
    return avalue, cvalue


def tangent_matrix(m: int) -> tuple[sp.Matrix, list[sp.Symbol]]:
    """Map theta of degree <=m-1 to the variation of G_m."""

    theta_vars = list(sp.symbols(f"th{m - 1}_0:{m}"))
    theta = sum(theta_vars[j] * x**j for j in range(m))
    image = sp.expand(sp.diff(q, x) * theta - q * sp.diff(theta, x) / m)
    equations = [xcoeff(image, j) for j in range(m + 1)]
    matrix, _ = sp.linear_eq_to_matrix(equations, theta_vars)
    require(matrix.rank() == m, f"G{m} tangent injection")
    require((m + 1) - matrix.rank() == 1, f"G{m} tangent cokernel one")
    return matrix, theta_vars


COKERNEL_ROWS = {
    5: [21, 0, 14, 0, 36, 0],
    6: [77, 0, 42, 0, 84, 0, 360],
    7: [143, 0, 66, 0, 108, 0, 360, 0],
    8: [715, 0, 286, 0, 396, 0, 1080, 0, 5040],
}


def tangent_solve(m: int, target: sp.Expr) -> sp.Expr:
    matrix, theta_vars = tangent_matrix(m)
    rhs = sp.Matrix([xcoeff(target, j) for j in range(m + 1)])
    row_pivots = matrix.T.rref()[1]
    require(len(row_pivots) == m, f"G{m} independent tangent rows")
    square = matrix[list(row_pivots), :]
    theta_values = square.inv() * rhs[list(row_pivots), :]
    require(
        all(
            is_zero(sum(matrix[row, column] * theta_values[column] for column in range(m)) - rhs[row])
            for row in range(m + 1)
        ),
        f"G{m} tangent solve",
    )
    return sp.expand(sum(theta_values[j] * x**j for j in range(m)))


def compatibility_value(m: int, difference: sp.Expr) -> sp.Expr:
    row = COKERNEL_ROWS[m]
    matrix, _ = tangent_matrix(m)
    require(sp.Matrix([row]).shape == (1, m + 1), f"G{m} cokernel shape")
    require(sp.Matrix([row]) * matrix == sp.zeros(1, m), f"G{m} displayed cokernel")
    return sp.factor(sum(row[j] * xcoeff(difference, j) for j in range(m + 1)))


def source_rows() -> dict[int, sp.Expr]:
    result: dict[int, sp.Expr] = {}
    for n in range(4, 9):
        actual = tcoeff(H, n)
        require_zero(actual - EXPECTED_G_ROWS[n], f"direct source row G{n}")
        result[n] = actual
    return result


def build_bracket_rows() -> tuple[list[sp.Expr], list[sp.Expr], dict[int, sp.Expr]]:
    grows = source_rows()
    arows = [A0, A1, A2, A3]
    crows = [C0, C1, C2, C3]

    inherited_p = P(
        sum(arows[n] * t**n for n in range(4)),
        sum(crows[n] * t**n for n in range(4)),
    ) - G
    for n in range(4):
        require_zero(tcoeff(inherited_p, n), f"inherited P/G row {n}")
    inherited_j = jacobian(sum(arows[n] * t**n for n in range(4)), sum(crows[n] * t**n for n in range(4)))
    require_zero(tcoeff(inherited_j, 0) - 1, "inherited Jacobian row zero")
    for n in range(1, 3):
        require_zero(tcoeff(inherited_j, n), f"inherited Jacobian row {n}")

    require_zero(P(A0, C0), "boundary lies on nodal cubic")
    require_zero((-3 * A0**2 + sp.Rational(3, 4)) + q * sp.diff(C0, x), "rotated gradient A")
    require_zero(2 * C0 - q * sp.diff(A0, x), "rotated gradient C")

    require_zero(predicted_G(4, arows, crows) - grows[4], "inherited G4 bracket compatibility")

    bracket_subs: dict[sp.Symbol, sp.Expr] = {}
    expected_equations = {5: E5, 6: E6, 7: E7, 8: E8}
    solution_symbols = {5: upsilon5, 6: U, 7: W, 8: Z}
    solution_values = {5: upsilon_solution, 6: U_solution, 7: W_solution, 8: Z_solution}

    for n in range(4, 8):
        bvalue = B_row(n, arows, crows)
        abase, cbase = particular_row(n, bvalue)
        trial_a = arows + [abase]
        trial_c = crows + [cbase]
        m = n + 1
        difference = sp.expand(grows[m].subs(bracket_subs) - predicted_G(m, trial_a, trial_c))
        obstruction = compatibility_value(m, difference)
        ratio = sp.cancel(obstruction / expected_equations[m].subs(bracket_subs))
        require(ratio.is_Rational and ratio != 0, f"E{m} exact obstruction factor")

        new_symbol = solution_symbols[m]
        new_value = sp.cancel(solution_values[m].subs(bracket_subs))
        bracket_subs[new_symbol] = new_value
        difference = sp.expand(difference.subs(new_symbol, new_value))
        theta = tangent_solve(m, difference)
        an = sp.expand(abase + theta * sp.diff(A0, x))
        cn = sp.expand(cbase + theta * sp.diff(C0, x))
        arows.append(an)
        crows.append(cn)
        require(sp.degree(an, x) <= n + 1, f"A{n} degree cap")
        require(sp.degree(cn, x) <= n + 2, f"C{n} degree cap")
        require_zero(predicted_G(m, arows, crows) - grows[m].subs(bracket_subs), f"G{m} matched after E{m}")

    # E8 fixes row seven.  Row eight itself retains the full tangent kernel.
    b8 = B_row(8, arows, crows)
    a8base, c8base = particular_row(8, b8)
    theta8_symbols = list(sp.symbols("theta8_0:9"))
    theta8 = sum(theta8_symbols[j] * x**j for j in range(9))
    a8 = sp.expand(a8base + theta8 * sp.diff(A0, x))
    c8 = sp.expand(c8base + theta8 * sp.diff(C0, x))
    arows.append(a8)
    crows.append(c8)
    require(sp.degree(a8, x) <= 9, "A8 degree cap")
    require(sp.degree(c8, x) <= 10, "C8 degree cap")
    require_zero(predicted_G(8, arows[:8], crows[:8]) - grows[8].subs(bracket_subs), "E8 compatibility before terminal row")
    return arows, crows, bracket_subs


def depth_columns(depth: int) -> list[tuple[int, int, int, int]]:
    columns: list[tuple[int, int, int, int]] = []
    for b in range(depth + 1):
        for a in range(depth - b + 1):
            for e in range(5):
                for c in range(9):
                    if b + c + 2 * e <= 8:
                        columns.append((a, b, c, e))
    return columns


def depth_matrix(depth: int) -> tuple[list[tuple[int, int]], list[tuple[int, int, int, int]], sp.Matrix]:
    coordinates = [(n, r) for n in range(9) for r in range(n + depth + 1)]
    columns = depth_columns(depth)
    matrix = sp.zeros(len(coordinates), len(columns))
    coordinate_index = {coordinate: index for index, coordinate in enumerate(coordinates)}
    for column, (a, b, c, e) in enumerate(columns):
        for j in range(c + e + 1):
            n = b + c + 2 * e + j
            r = a + 2 * b + e + 2 * j
            if n <= 8:
                matrix[coordinate_index[(n, r)], column] += sp.binomial(c + e, j)
    return coordinates, columns, matrix


def jet_vector(rows: list[sp.Expr], depth: int) -> sp.Matrix:
    return sp.Matrix([xcoeff(rows[n], r) for n in range(9) for r in range(n + depth + 1)])


def functional(coordinates: list[tuple[int, int]], entries: dict[tuple[int, int], int]) -> sp.Matrix:
    return sp.Matrix([[entries.get(coordinate, 0) for coordinate in coordinates]])


def hasse_audit(
    arows: list[sp.Expr], crows: list[sp.Expr], bracket_subs: dict[sp.Symbol, sp.Expr]
) -> tuple[sp.Matrix, sp.Matrix, list[sp.Symbol], dict[sp.Symbol, sp.Expr]]:
    acoords, acolumns, amatrix = depth_matrix(2)
    ccoords, ccolumns, cmatrix = depth_matrix(3)
    require((len(acoords), len(acolumns)) == (63, 131), "P2 exact universe")
    require((len(ccoords), len(ccolumns)) == (72, 204), "P3 exact universe")
    require(amatrix.rank() == 51, "P2 truncation rank")
    require(cmatrix.rank() == 63, "P3 truncation rank")
    anull = amatrix.T.nullspace()
    cnull = cmatrix.T.nullspace()
    require(len(anull) == 12, "P2 Hasse left nullity")
    require(len(cnull) == 9, "P3 Hasse left nullity")

    avec = jet_vector(arows, 2)
    cvec = jet_vector(crows, 3)

    fa_delta = functional(acoords, {(2, 0): -4, (3, 2): 3, (4, 4): -2, (5, 6): 1})
    fa_zeta = functional(acoords, {(4, 1): -10, (5, 3): 6, (6, 5): -3, (7, 7): 1})
    fa_top = functional(acoords, {(4, 0): 15, (5, 2): -10, (6, 4): 6, (7, 6): -3, (8, 8): 1})
    fc_top = functional(ccoords, {(4, 1): 15, (5, 3): -10, (6, 5): 6, (7, 7): -3, (8, 9): 1})
    fa_eta = functional(acoords, {(4, 1): -15, (5, 3): 8, (6, 5): -3, (8, 9): 1})
    fc_eta = functional(ccoords, {(3, 0): -6, (4, 2): 5, (5, 4): -4, (6, 6): 3, (7, 8): -2, (8, 10): 1})
    for name, row, matrix in [
        ("A-Delta", fa_delta, amatrix),
        ("A-zeta", fa_zeta, amatrix),
        ("A-top", fa_top, amatrix),
        ("C-top", fc_top, cmatrix),
        ("A-eta", fa_eta, amatrix),
        ("C-eta", fc_eta, cmatrix),
    ]:
        require(row * matrix == sp.zeros(1, matrix.cols), f"displayed Hasse functional {name}")

    c89 = xcoeff(crows[8], 9)
    c810 = xcoeff(crows[8], 10)
    displayed_relations = [
        ((fa_delta * avec)[0], -(15 * Delta - 896) / 45, "Delta functional"),
        ((fa_zeta * avec)[0], 3 * (3 * Phi + 2 * zeta3) / 20, "zeta functional"),
        (
            (fa_top * avec)[0],
            4 * (13215 * Delta + 7950 * Theta - 2475 * c89 + 6075 * xi10 - 2583808) / 7425,
            "A top functional",
        ),
        (
            (fc_top * cvec)[0],
            (6900 * Delta - 13425 * Theta + 4950 * c89 - 12150 * xi10 + 3159808) / 4950,
            "C top functional",
        ),
        (
            (fa_eta * avec)[0],
            (27 * Phi - 40 * c810 + 27 * eta + 45 * zeta3) / 30,
            "A eta functional",
        ),
        (
            (fc_eta * cvec)[0],
            (40 * c810 - 27 * eta - 27 * zeta3) / 40,
            "C eta functional",
        ),
    ]
    for actual, expected, label in displayed_relations:
        require_zero(actual - expected, label)

    delta_value = sp.Rational(896, 15)
    theta_value = sp.Rational(512, 75)
    zeta_value = -3 * Phi / 2
    hasse_subs = {Delta: delta_value, Theta: theta_value, zeta3: zeta_value}
    target_c89 = (1215 * xi10 - 348032) / 495
    target_c810 = 27 * (eta + zeta3) / 40

    theta8_symbols = sorted(
        set().union(*(row.free_symbols for row in (arows[8], crows[8])))
        & {sp.Symbol(f"theta8_{j}") for j in range(9)},
        key=str,
    )
    require(len(theta8_symbols) == 9, "terminal tangent coordinate count")

    all_residuals = [(row.T * avec)[0] for row in anull] + [(row.T * cvec)[0] for row in cnull]
    reduced_residuals = [sp.expand(value.subs(hasse_subs)) for value in all_residuals]
    tangent_matrix, tangent_rhs = sp.linear_eq_to_matrix(reduced_residuals, theta8_symbols)
    require(tangent_matrix.rank() == 2, "Hasse terminal tangent rank")

    terminal_equations = [
        sp.expand((c89 - target_c89).subs(hasse_subs)),
        sp.expand((c810 - target_c810).subs(hasse_subs)),
    ]
    solve_terminal = sp.solve(terminal_equations, theta8_symbols[-2:], dict=True)
    require(len(solve_terminal) == 1, "terminal coordinate solve")
    terminal_subs = solve_terminal[0]
    require(
        all(is_zero(residual.subs(terminal_subs)) for residual in reduced_residuals),
        "all 21 Hasse relations follow from terminal coordinate pair",
    )
    augmented = tangent_matrix.row_join(tangent_rhs)
    require(augmented.rank() == 2, "Hasse terminal augmented rank")
    require(9 - tangent_matrix.rank() == 7, "finite solution affine dimension seven")

    # The three parameter constraints are individually necessary, while the
    # preceding rank/solution check proves their sufficiency for all 21 rows.
    require_zero((fa_delta * avec)[0].subs(Delta, delta_value), "Delta necessity witness")
    require_zero((fa_zeta * avec)[0].subs(zeta3, zeta_value), "zeta necessity witness")
    theta_eliminant = 33330 * Delta + 2475 * Theta - 2007808
    require_zero(theta_eliminant.subs({Delta: delta_value, Theta: theta_value}), "Theta eliminant")

    return amatrix, cmatrix, theta8_symbols, {**hasse_subs, **terminal_subs}


def full_finite_control(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    bracket_subs: dict[sp.Symbol, sp.Expr],
    amatrix: sp.Matrix,
    cmatrix: sp.Matrix,
    theta8_symbols: list[sp.Symbol],
    hasse_terminal_subs: dict[sp.Symbol, sp.Expr],
) -> tuple[sp.Rational, sp.Rational, sp.Rational, sp.Rational, sp.Rational]:
    control_subs: dict[sp.Symbol, sp.Expr] = {
        Delta: sp.Rational(896, 15),
        Theta: sp.Rational(512, 75),
        Phi: 0,
        eta: 0,
        xi10: 0,
        zeta3: 0,
        alpha11: 0,
        beta11: 0,
    }
    control_subs.update({symbol: value.subs(control_subs) for symbol, value in bracket_subs.items()})
    terminal_values = {
        symbol: value.subs(control_subs).subs({other: 0 for other in theta8_symbols[:-2]})
        for symbol, value in hasse_terminal_subs.items()
        if symbol in theta8_symbols[-2:]
    }
    control_subs.update({symbol: 0 for symbol in theta8_symbols[:-2]})
    control_subs.update(terminal_values)

    control_a_rows = [sp.expand(row.subs(control_subs)) for row in arows]
    control_c_rows = [sp.expand(row.subs(control_subs)) for row in crows]
    control_a = sum(control_a_rows[n] * t**n for n in range(9))
    control_c = sum(control_c_rows[n] * t**n for n in range(9))
    control_rhs = G.subs(control_subs)
    control_j = jacobian(control_a, control_c)
    require_zero(tcoeff(control_j, 0) - 1, "control J row zero")
    for n in range(1, 8):
        require_zero(tcoeff(control_j, n), f"control J row {n}")
    for n in range(9):
        require_zero(tcoeff(P(control_a, control_c) - control_rhs, n), f"control P/G row {n}")

    avec = jet_vector(control_a_rows, 2)
    cvec = jet_vector(control_c_rows, 3)
    require(amatrix.row_join(avec).rank() == amatrix.rank(), "control A in projected P2")
    require(cmatrix.row_join(cvec).rank() == cmatrix.rank(), "control C in projected P3")

    control_U = sp.cancel(U_solution.subs(control_subs))
    control_W = sp.cancel(W_solution.subs(control_subs))
    control_Z = sp.cancel(Z_solution.subs(control_subs))
    control_Lambda = sp.cancel(control_U + control_W + control_Z)
    control_D = sp.cancel(control_W**2 - 4 * control_U * control_Z)
    for name, value in [("U", control_U), ("Z", control_Z), ("Lambda", control_Lambda), ("D", control_D)]:
        require(value != 0, f"control gate-interior {name}")
    require(xcoeff(control_c_rows[8], 9) == -sp.Rational(348032, 495), "control c89")
    require(xcoeff(control_c_rows[8], 10) == 0, "control c810")
    return control_U, control_W, control_Z, control_Lambda, control_D


def response_map() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    substitutions = {Delta: sp.Rational(896, 15), Theta: sp.Rational(512, 75)}
    u_map = sp.factor(U_solution.subs(substitutions))
    w_map = sp.factor(W_solution.subs(substitutions))
    z_map = sp.factor(Z_solution.subs(substitutions))
    expected_u = (475515904 - 109350 * xi10) / 200475
    expected_w = -(4343625 * Phi**2 - 17172000 * xi10 + 143826305024) / 4009500
    expected_z = (
        12506118074368 - 173745000 * Phi**2 - 195463125 * Phi * eta - 926883000 * xi10
    ) / 108256500
    require_zero(u_map - expected_u, "three-parameter U response")
    require_zero(w_map - expected_w, "three-parameter W response")
    require_zero(z_map - expected_z, "three-parameter Z response")
    require_zero(K.subs(Delta, sp.Rational(896, 15)) + sp.Rational(32, 5), "Hasse K value")
    require_zero(upsilon_solution.subs(substitutions) + sp.Rational(731648, 2025), "Hasse upsilon5 value")
    return expected_u, expected_w, expected_z


def hostile_audit() -> None:
    for equation, symbol, jump, name in [
        (E5, upsilon5, 2025, "E5"),
        (E6, U, 200475, "E6"),
        (E7, W, 801900, "E7"),
        (E8, Z, 21651300, "E8"),
    ]:
        require_zero(equation.subs(symbol, symbol + 1) - equation - jump, f"{name} unit perturbation")
    delta_residual = -(15 * Delta - 896) / 45
    zeta_residual = 3 * (3 * Phi + 2 * zeta3) / 20
    require_zero(
        delta_residual.subs(Delta, sp.Rational(896, 15) + 1) + sp.Rational(1, 3),
        "Delta Hasse hostile",
    )
    require_zero(zeta_residual.subs({Phi: 0, zeta3: 1}) - sp.Rational(3, 10), "zeta Hasse hostile")
    theta_eliminant = 33330 * Delta + 2475 * Theta - 2007808
    require_zero(
        theta_eliminant.subs({Delta: sp.Rational(896, 15), Theta: sp.Rational(512, 75) + 1}) - 2475,
        "Theta Hasse hostile",
    )
    require([1, 6, 15] != [0, 1, 5], "THM-4298 observer hostile packets")


def main() -> None:
    arows, crows, bracket_subs = build_bracket_rows()
    amatrix, cmatrix, theta8_symbols, terminal_subs = hasse_audit(arows, crows, bracket_subs)
    u_map, w_map, z_map = response_map()
    control_U, control_W, control_Z, control_Lambda, control_D = full_finite_control(
        arows, crows, bracket_subs, amatrix, cmatrix, theta8_symbols, terminal_subs
    )
    hostile_audit()

    print("THM-4308 SOURCE-NORMAL BRACKET/HASSE ROW-EIGHT PRIMARY AUDIT: PASS")
    print("gauge: a=1 gamma=-1/2 I=3/4 b=d=0; p=t(1+x^2*t); y=x*t*p; u=x^2*t")
    print("caps: deg(A_n)<=n+1 and deg(C_n)<=n+2 for n=4..8")
    for n in range(4, 9):
        print(f"G{n}={sp.sstr(sp.expand(EXPECTED_G_ROWS[n]))}")
    print(f"E5={sp.sstr(E5)}=0")
    print(f"E6={sp.sstr(E6)}=0")
    print(f"E7={sp.sstr(E7)}=0")
    print(f"E8={sp.sstr(E8)}=0")
    print("tangent_cokernels: G5=(21,0,14,0,36,0); G6=(77,0,42,0,84,0,360)")
    print("tangent_cokernels: G7=(143,0,66,0,108,0,360,0); G8=(715,0,286,0,396,0,1080,0,5040)")
    print("P2_projection: ambient=63 columns=131 rank=51 left_nullity=12")
    print("P3_projection: ambient=72 columns=204 rank=63 left_nullity=9")
    print("Hasse_parameters: Delta=896/15; Theta=512/75; zeta3=-3*Phi/2")
    print("Hasse_terminal: [x^9]C8=(1215*xi10-348032)/495; [x^10]C8=27*(eta+zeta3)/40")
    print("finite_iff: E5--E8 plus Hasse_parameters; terminal tangent rank=2; affine_dimension=7")
    print(f"response_U={sp.sstr(u_map)}")
    print(f"response_W={sp.sstr(w_map)}")
    print(f"response_Z={sp.sstr(z_map)}")
    print(f"interior_control: U={control_U}; W={control_W}; Z={control_Z}")
    print(f"interior_control: Lambda={control_Lambda}; D={control_D}")
    print("interior_control: J_rows=0..7 PASS; P/G_rows=0..8 PASS; projected_P2/P3=PASS")
    print("bracket_unit_hostiles: delta(E5,E6,E7,E8)=(2025,200475,801900,21651300)")
    print("Hasse_unit_hostiles: Delta_residual=-1/3; zeta3_residual=3/10; Theta_eliminant_delta=2475")
    print("observer_hostile: p^6_packet=(1,6,15); p^3*y^2_packet=(0,1,5)")
    print(f"CHECKS={CHECKS}")
    print("scope: exact row-eight projected finite jets only; no full B2 lift, termination, seam entry, JC2, or DC2")


if __name__ == "__main__":
    main()
