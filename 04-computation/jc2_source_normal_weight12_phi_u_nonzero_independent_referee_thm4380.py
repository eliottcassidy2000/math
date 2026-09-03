#!/usr/bin/env python3
"""Standalone exact audit of the THM-4380 Phi*U != 0 scratch candidate.

This file intentionally imports no repository-local program.  It starts from
the canonical source polynomial and boundary rows in THM-4308, reconstructs
the bracket recursion and the complete projected P_2/P_3 matrices through
row eleven, and then checks the row-twelve Student obstruction.
"""

from __future__ import annotations

from math import gcd

import sympy as sp


CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def zz(expr: sp.Expr) -> sp.Expr:
    return sp.cancel(sp.together(sp.expand(expr)))


def check0(expr: sp.Expr, label: str) -> None:
    check(zz(expr) == 0, label)


x, t = sp.symbols("x t")
Delta, Phi, Theta = sp.symbols("Delta Phi Theta")
eta, zeta3, upsilon5, xi10 = sp.symbols("eta zeta3 upsilon5 xi10")
alpha11, beta11 = sp.symbols("alpha11 beta11")
U, W, Z = sp.symbols("U W Z")

K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
q = -(x**2 + 6) / 2
p = t * (1 + x**2 * t)
y = x * t * p
u = x**2 * t

Hsrc = sp.expand(
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
Gsrc = sp.expand(-u / 2 + Hsrc)

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


def xcoef(expr: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(expr).coeff(x, degree)


def tcoef(expr: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(expr).coeff(t, degree)


def cubic(a: sp.Expr, c: sp.Expr) -> sp.Expr:
    return sp.expand(c**2 - a**3 + sp.Rational(3, 4) * a + sp.Rational(1, 4))


def jacobian(a: sp.Expr, c: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(a, x) * sp.diff(c, t) - sp.diff(a, t) * sp.diff(c, x))


def brow(m: int, aa: list[sp.Expr], cc: list[sp.Expr]) -> sp.Expr:
    return sp.expand(sum(
        (m - i) * sp.diff(aa[i], x) * cc[m - i]
        - i * aa[i] * sp.diff(cc[m - i], x)
        for i in range(1, m)
    ))


def trow(m: int, aa: list[sp.Expr], cc: list[sp.Expr]) -> sp.Expr:
    quadratic = sum(cc[i] * cc[m - i] for i in range(1, m))
    cub = sum(
        aa[i] * aa[j] * aa[k]
        for i in range(m)
        for j in range(m)
        for k in range(m)
        if i + j + k == m
    )
    return sp.expand(quadratic - cub)


def predicted_g(m: int, aa: list[sp.Expr], cc: list[sp.Expr]) -> sp.Expr:
    return sp.expand(trow(m, aa, cc) - q * brow(m, aa, cc) / m)


def determinant_particular(m: int, bvalue: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    avars = list(sp.symbols(f"a{m}b_0:{m + 2}"))
    cvars = list(sp.symbols(f"c{m}b_0:{m + 3}"))
    aval = sum(avars[j] * x**j for j in range(m + 2))
    cval = sum(cvars[j] * x**j for j in range(m + 3))
    op = sp.expand(m * (sp.diff(A0, x) * cval - sp.diff(C0, x) * aval))
    mat, _ = sp.linear_eq_to_matrix([xcoef(op, j) for j in range(m + 4)], avars + cvars)
    pivots = mat.rref()[1]
    check(mat.rank() == m + 4 and len(avars + cvars) - mat.rank() == m + 1,
          f"determinant row {m} exact rank")
    rhs = sp.Matrix([-xcoef(bvalue, j) for j in range(m + 4)])
    vals = [sp.Integer(0)] * (2 * m + 5)
    sol = mat[:, list(pivots)].inv() * rhs
    for j, value in zip(pivots, sol):
        vals[j] = zz(value)
    a = sp.expand(sum(vals[j] * x**j for j in range(m + 2)))
    c = sp.expand(sum(vals[m + 2 + j] * x**j for j in range(m + 3)))
    check0(m * (sp.diff(A0, x) * c - sp.diff(C0, x) * a) + bvalue,
           f"determinant row {m} particular")
    return a, c


def tangent_image_matrix(m: int) -> tuple[sp.Matrix, list[sp.Symbol]]:
    variables = list(sp.symbols(f"d{m - 1}_0:{m}"))
    theta = sum(variables[j] * x**j for j in range(m))
    image = sp.expand(sp.diff(q, x) * theta - q * sp.diff(theta, x) / m)
    mat, _ = sp.linear_eq_to_matrix([xcoef(image, j) for j in range(m + 1)], variables)
    check(mat.rank() == m, f"Student tangent rank {m}")
    return mat, variables


def student_row(m: int) -> list[int]:
    values = []
    for degree in range(m + 1):
        if degree % 2:
            values.append(sp.Rational(0))
        else:
            r = degree // 2
            num = 6**r * sp.factorial2(2 * r - 1)
            den = sp.prod(2 * m - 2 * j + 1 for j in range(1, r + 1))
            values.append(sp.Rational(num, den))
    common = sp.ilcm(*[int(sp.denom(v)) for v in values])
    ints = [int(common * v) for v in values]
    divisor = 0
    for value in ints:
        divisor = gcd(divisor, abs(value))
    ints = [value // divisor for value in ints]
    if ints[0] < 0:
        ints = [-value for value in ints]
    mat, _ = tangent_image_matrix(m)
    check(sp.Matrix([ints]) * mat == sp.zeros(1, m), f"Student cokernel {m}")
    return ints


def student_value(m: int, expr: sp.Expr) -> sp.Expr:
    row = student_row(m)
    return zz(sum(row[j] * xcoef(expr, j) for j in range(m + 1)))


def solve_tangent(m: int, target: sp.Expr) -> sp.Expr:
    mat, variables = tangent_image_matrix(m)
    rhs = sp.Matrix([xcoef(target, j) for j in range(m + 1)])
    rows = mat.T.rref()[1]
    solution = mat[list(rows), :].inv() * rhs[list(rows), :]
    check(all(zz((mat * solution - rhs)[j]) == 0 for j in range(m + 1)),
          f"tangent solve {m}")
    return sp.expand(sum(solution[j] * x**j for j in range(m)))


def depth_projection(depth: int, maxrow: int) -> tuple[list[tuple[int, int]], sp.Matrix]:
    coords = [(n, r) for n in range(maxrow + 1) for r in range(n + depth + 1)]
    index = {coord: j for j, coord in enumerate(coords)}
    columns: list[tuple[int, int, int, int]] = []
    for b in range(depth + 1):
        for a in range(depth - b + 1):
            for e in range(maxrow // 2 + 1):
                for c in range(maxrow + 1):
                    if b + c + 2 * e <= maxrow:
                        columns.append((a, b, c, e))
    matrix = sp.zeros(len(coords), len(columns))
    for col, (a, b, c, e) in enumerate(columns):
        for j in range(c + e + 1):
            n = b + c + 2 * e + j
            if n <= maxrow:
                r = a + 2 * b + e + 2 * j
                matrix[index[(n, r)], col] += sp.binomial(c + e, j)
    return coords, matrix


def depth_equations(rows: list[sp.Expr], depth: int, maxrow: int) -> tuple[list[sp.Expr], sp.Matrix]:
    coords, matrix = depth_projection(depth, maxrow)
    jet = sp.Matrix([xcoef(rows[n], r) for n, r in coords])
    equations = [zz((null.T * jet)[0]) for null in matrix.T.nullspace()]
    return equations, matrix


def solve_affine(equations: list[sp.Expr], variables: list[sp.Symbol]) -> tuple[dict[sp.Symbol, sp.Expr], sp.Matrix, sp.Matrix]:
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    pivot_columns = list(matrix.rref()[1])
    free_columns = [j for j in range(len(variables)) if j not in pivot_columns]
    # In every depth solve below the free terminal directions are genuinely
    # silent columns, so retain them instead of choosing the zero section.
    check(matrix[:, free_columns] == sp.zeros(matrix.rows, len(free_columns)),
          "free affine directions are silent")
    reduced = equations
    pivars = [variables[j] for j in pivot_columns]
    cmat, crhs = sp.linear_eq_to_matrix(reduced, pivars)
    rows = cmat.T.rref()[1]
    sol = cmat[list(rows), :].inv() * crhs[list(rows), :]
    result = {symbol: zz(sol[j]) for j, symbol in enumerate(pivars)}
    return result, matrix, rhs


def primitive_expr(expr: sp.Expr, variables: tuple[sp.Symbol, ...]) -> sp.Expr:
    numerator = sp.fraction(zz(expr))[0]
    poly = sp.Poly(numerator, *variables, domain=sp.QQ)
    den = sp.ilcm(*[int(sp.denom(c)) for c in poly.coeffs()])
    ipoly = sp.Poly(sp.expand(den * poly.as_expr()), *variables, domain=sp.ZZ)
    _, primitive = ipoly.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return sp.factor(primitive.as_expr())


def unique_scalars(expressions: list[sp.Expr]) -> list[sp.Expr]:
    answer: list[sp.Expr] = []
    for expr in expressions:
        if expr == 0:
            continue
        if not any(zz(expr / old).is_Rational for old in answer):
            answer.append(expr)
    return answer


def even_phi_to_y(expr: sp.Expr, yy: sp.Symbol) -> sp.Expr:
    source = sp.Poly(sp.expand(expr), Phi)
    result = 0
    for (degree,), coefficient in source.terms():
        check(degree % 2 == 0, "even Phi parity")
        result += coefficient * yy ** (degree // 2)
    return sp.expand(result)


def integer_primitive(poly: sp.Poly) -> sp.Poly:
    rat = sp.Poly(poly.as_expr(), poly.gens[0], domain=sp.QQ)
    den = sp.ilcm(*[int(sp.denom(c)) for c in rat.all_coeffs()])
    integral = sp.Poly(sp.expand(den * rat.as_expr()), poly.gens[0], domain=sp.ZZ)
    _, primitive = integral.primitive()
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive


def main() -> None:
    # Direct source rows 4--12, reconstructed from H(p,y).
    grows = {m: sp.expand(tcoef(Hsrc, m)) for m in range(4, 13)}
    expected_sparse = {
        9: (20 * U + 10 * W + 4 * Z) * x**6
           + (10 * alpha11 + 6 * beta11) * x**7
           + (5 * upsilon5 + 4 * xi10) * x**8 + (eta + zeta3) * x**9,
        10: (15 * U + 10 * W + 6 * Z) * x**8
            + (5 * alpha11 + 4 * beta11) * x**9 + (upsilon5 + xi10) * x**10,
        11: (6 * U + 5 * W + 4 * Z) * x**10 + (alpha11 + beta11) * x**11,
        12: (U + W + Z) * x**12,
    }
    for m, expected in expected_sparse.items():
        check0(grows[m] - expected, f"direct source row G{m}")

    aa = [A0, A1, A2, A3]
    cc = [C0, C1, C2, C3]
    atrunc = sum(aa[n] * t**n for n in range(4))
    ctrunc = sum(cc[n] * t**n for n in range(4))
    for n in range(4):
        check0(tcoef(cubic(atrunc, ctrunc) - Gsrc, n), f"inherited cubic row {n}")
    jtrunc = jacobian(atrunc, ctrunc)
    check0(tcoef(jtrunc, 0) - 1, "inherited Jacobian row 0")
    for n in (1, 2):
        check0(tcoef(jtrunc, n), f"inherited Jacobian row {n}")
    check0(predicted_g(4, aa, cc) - grows[4], "row 4 compatibility")

    solved_source: dict[sp.Symbol, sp.Expr] = {}
    next_symbols = {5: upsilon5, 6: U, 7: W, 8: Z}
    for n in range(4, 8):
        abase, cbase = determinant_particular(n, brow(n, aa, cc))
        m = n + 1
        trial_a, trial_c = aa + [abase], cc + [cbase]
        difference = zz(grows[m].subs(solved_source) - predicted_g(m, trial_a, trial_c))
        obstruction = student_value(m, difference)
        symbol = next_symbols[m]
        value = zz(sp.solve(obstruction, symbol)[0])
        solved_source[symbol] = value
        difference = zz(difference.subs(symbol, value))
        theta = solve_tangent(m, difference)
        aa.append(sp.expand(abase + theta * sp.diff(A0, x)))
        cc.append(sp.expand(cbase + theta * sp.diff(C0, x)))
        check0(predicted_g(m, aa, cc) - grows[m].subs(solved_source), f"row {m} matched")

    # Append determinant row 8 and impose the exact row-eight depth gate.
    a8base, c8base = determinant_particular(8, brow(8, aa, cc))
    th8 = list(sp.symbols("theta8_0:9"))
    theta8 = sum(th8[j] * x**j for j in range(9))
    aa.append(sp.expand(a8base + theta8 * sp.diff(A0, x)))
    cc.append(sp.expand(c8base + theta8 * sp.diff(C0, x)))
    gate = {Delta: sp.Rational(896, 15), Theta: sp.Rational(512, 75), zeta3: -3 * Phi / 2}
    aeq8, amat8 = depth_equations(aa, 2, 8)
    ceq8, cmat8 = depth_equations(cc, 3, 8)
    check((amat8.shape, amat8.rank(), cmat8.shape, cmat8.rank()) ==
          ((63, 131), 51, (72, 204), 63), "row 8 complete depth universes")
    eq8 = [zz(eq.subs(gate)) for eq in aeq8 + ceq8]
    terminal8, dmat8, drhs8 = solve_affine(eq8, th8)
    check(dmat8.rank() == dmat8.row_join(drhs8).rank() == 2, "row 8 depth rank")
    check(all(zz(eq.subs(terminal8)) == 0 for eq in eq8), "row 8 all depth residuals")
    base_subs = {**solved_source, **gate, **terminal8}
    depth8_a = [zz(row.subs(solved_source).subs(gate).subs(terminal8)) for row in aa]
    depth8_c = [zz(row.subs(solved_source).subs(gate).subs(terminal8)) for row in cc]

    # Row 9 Student equation, then its seven-direction bracket selector.
    g9 = zz(grows[9].subs(solved_source).subs(gate).subs(terminal8))
    diff9 = zz(g9 - predicted_g(9, depth8_a, depth8_c))
    student9 = student_value(9, diff9)
    alpha_solution = zz(sp.solve(student9, alpha11)[0])
    diff9 = zz(diff9.subs(alpha11, alpha_solution))
    eq9br = [xcoef(diff9, j) for j in range(10)]
    section9, bmat9, brhs9 = solve_affine(eq9br, th8[:7])
    if bmat9.row_join(brhs9).rank() != 7:
        print("DEBUG_student9=", sp.sstr(student9))
        print("DEBUG_alpha=", sp.sstr(alpha_solution))
        print("DEBUG_rem9=", [sp.sstr(zz(eq.subs(section9))) for eq in eq9br])
    check(bmat9.rank() == bmat9.row_join(brhs9).rank() == 7,
          f"row 9 bracket rank {bmat9.rank()}/{bmat9.row_join(brhs9).rank()}")
    check(all(zz(eq.subs(section9)) == 0 for eq in eq9br), "row 9 bracket all rows")
    selected8_a = [zz(row.subs(alpha11, alpha_solution).subs(section9)) for row in depth8_a]
    selected8_c = [zz(row.subs(alpha11, alpha_solution).subs(section9)) for row in depth8_c]

    # Append determinant row 9, impose complete depth, retain seven free directions.
    a9base, c9base = determinant_particular(9, brow(9, selected8_a, selected8_c))
    th9 = list(sp.symbols("theta9_0:10"))
    theta9 = sum(th9[j] * x**j for j in range(10))
    row9_a = selected8_a + [sp.expand(a9base + theta9 * sp.diff(A0, x))]
    row9_c = selected8_c + [sp.expand(c9base + theta9 * sp.diff(C0, x))]
    aeq9, amat9 = depth_equations(row9_a, 2, 9)
    ceq9, cmat9 = depth_equations(row9_c, 3, 9)
    check((amat9.shape, amat9.rank(), cmat9.shape, cmat9.rank()) ==
          ((75, 160), 59, (85, 251), 73), "row 9 complete depth universes")
    terminal9, dmat9, drhs9 = solve_affine(aeq9 + ceq9, th9)
    check(dmat9.rank() == dmat9.row_join(drhs9).rank() == 3, "row 9 depth rank")
    check(tuple(dmat9.rref()[1]) == (7, 8, 9), "row 9 depth pivots")
    check(all(zz(eq.subs(terminal9)) == 0 for eq in aeq9 + ceq9), "row 9 all depth residuals")
    top9 = {th9[j]: terminal9[th9[j]] for j in (7, 8, 9)}
    depth9_a = [zz(row.subs(top9)) for row in row9_a]
    depth9_c = [zz(row.subs(top9)) for row in row9_c]

    # Row 10 bracket gives two source cuts; solve them as the Phi,eta graph.
    g10 = zz(grows[10].subs(solved_source).subs(gate).subs(terminal8).subs(alpha11, alpha_solution))
    diff10 = zz(g10 - predicted_g(10, depth9_a, depth9_c))
    student10 = student_value(10, diff10)
    eq10br = [xcoef(diff10, j) for j in range(11)]
    section10, bmat10, brhs10 = solve_affine(eq10br, th9[:7])
    check(bmat10.rank() == 7, "row 10 restricted bracket tangent rank")
    rem10 = [zz(eq.subs(section10)) for eq in eq10br]
    nums10 = unique_scalars([
        primitive_expr(value, (Phi, eta, xi10, beta11)) for value in rem10 if value != 0
    ])
    check(len(nums10) == 2, "row 10 has two source numerators")
    no_beta = next(poly for poly in nums10 if beta11 not in poly.free_symbols)
    with_beta = next(poly for poly in nums10 if beta11 in poly.free_symbols)
    xi_surface = zz(sp.solve(no_beta, xi10)[0])
    beta_surface = zz(sp.solve(with_beta.subs(xi10, xi_surface), beta11)[0])
    surface = {xi10: xi_surface, beta11: beta_surface}
    alpha_surface = zz(alpha_solution.subs(surface))
    u_response = zz(solved_source[U].subs(gate).subs(surface))
    check0(student9.subs({alpha11: alpha_solution}), "row 9 Student source graph")
    check0(student10.subs(surface), "row 10 Student consequence")
    check(u_response != 0, "row 10 U nonzero rational-function control")
    row10_positive = zz(u_response.subs({Phi: 1, eta: 0}))
    check(row10_positive != 0, "row 10 explicit Phi=1 eta=0 U control")
    sec10_surface = {symbol: zz(value.subs(surface)) for symbol, value in section10.items()}
    selected9_a = [zz(row.subs(surface).subs(sec10_surface)) for row in depth9_a]
    selected9_c = [zz(row.subs(surface).subs(sec10_surface)) for row in depth9_c]
    g10_surface = zz(g10.subs(surface))
    check0(predicted_g(10, selected9_a, selected9_c) - g10_surface, "row 10 bracket match")

    # Append row 10 and impose complete joint depth.  The unique remnant is H.
    a10base, c10base = determinant_particular(10, brow(10, selected9_a, selected9_c))
    th10 = list(sp.symbols("theta10_0:11"))
    theta10 = sum(th10[j] * x**j for j in range(11))
    row10_a = selected9_a + [sp.expand(a10base + theta10 * sp.diff(A0, x))]
    row10_c = selected9_c + [sp.expand(c10base + theta10 * sp.diff(C0, x))]
    aeq10, amat10 = depth_equations(row10_a, 2, 10)
    ceq10, cmat10 = depth_equations(row10_c, 3, 10)
    check((amat10.shape, amat10.rank(), cmat10.shape, cmat10.rank()) ==
          ((88, 193), 68, (99, 304), 83), "row 10 complete depth universes")
    terminal10, dmat10, drhs10 = solve_affine(aeq10 + ceq10, th10)
    check(dmat10.rank() == 3 and dmat10.row_join(drhs10).rank() == 4,
          "row 10 joint depth generic rank jump")
    check(tuple(dmat10.rref()[1]) == (8, 9, 10), "row 10 depth pivots")
    rem10d = [zz(eq.subs(terminal10)) for eq in aeq10 + ceq10]
    h_candidates = unique_scalars([
        primitive_expr(value, (Phi, eta)) for value in rem10d if value != 0
    ])
    check(len(h_candidates) == 1, "row 10 depth has one gluing polynomial")
    h10 = h_candidates[0]
    upoly = sp.Poly(sp.fraction(u_response)[0], Phi, eta).primitive()[1]
    check(sp.gcd(sp.Poly(h10, Phi, eta), upoly) == 1, "H coprime to U")
    top10 = {th10[j]: terminal10[th10[j]] for j in (8, 9, 10)}
    depth10_a = [zz(row.subs(top10)) for row in row10_a]
    depth10_c = [zz(row.subs(top10)) for row in row10_c]

    # Row 11 bracket: eight tangent directions and one new source cut modulo H.
    g11 = zz(grows[11].subs(solved_source).subs(gate).subs(terminal8).subs(alpha11, alpha_solution).subs(surface))
    diff11 = zz(g11 - predicted_g(11, depth10_a, depth10_c))
    student11 = student_value(11, diff11)
    eq11br = [xcoef(diff11, j) for j in range(12)]
    section11, bmat11, brhs11 = solve_affine(eq11br, th10[:8])
    check(bmat11.rank() == 8, "row 11 restricted bracket tangent rank")
    rem11 = [zz(eq.subs(section11)) for eq in eq11br]
    hgb = sp.groebner([h10], eta, Phi, order="lex")
    s_candidates = []
    for value in rem11:
        if value == 0:
            continue
        primitive = primitive_expr(value, (eta, Phi))
        remainder = sp.factor(hgb.reduce(primitive)[1])
        if remainder != 0:
            s_candidates.append(remainder)
    s_candidates = unique_scalars(s_candidates)
    check(len(s_candidates) == 1, "row 11 has one new source cut modulo H")
    s11 = s_candidates[0]
    student11_primitive = primitive_expr(student11, (eta, Phi))
    student11_mod_h = sp.factor(hgb.reduce(student11_primitive)[1])
    check(zz(student11_mod_h / s11).is_Rational, "row 11 Student generates source cut")

    # Sign quotient lambda=eta/Phi, Y=Phi^2 and exact degree-seven eliminant.
    lam, yy = sp.symbols("lambda Y")
    hly = even_phi_to_y(zz(h10.subs(eta, lam * Phi) / Phi), yy)
    sly = even_phi_to_y(sp.expand(s11.subs(eta, lam * Phi)), yy)
    hpoly = sp.Poly(hly, yy)
    check(hpoly.degree() == 1, "H is linear in Y")
    alam = sp.Poly(hpoly.coeff_monomial(yy), lam)
    blam = sp.Poly(-hpoly.coeff_monomial(1), lam)
    check(sp.gcd(alam, blam).degree() == 0, "A and B coprime")
    ylam = zz(blam.as_expr() / alam.as_expr())
    eraw = sp.Poly(sp.fraction(zz(sly.subs(yy, ylam)))[0], lam)
    eliminant = integer_primitive(eraw)
    check(eliminant.degree() == 7, "row 11 eliminant degree seven")
    check(sp.gcd(eliminant, eliminant.diff()).degree() == 0, "row 11 eliminant squarefree")
    uly = even_phi_to_y(sp.expand(upoly.as_expr().subs(eta, lam * Phi)), yy)
    ulam = integer_primitive(sp.Poly(sp.fraction(zz(uly.subs(yy, ylam)))[0], lam))
    check(sp.gcd(eliminant, alam * blam * ulam).degree() == 0,
          "row 11 eliminant avoids A B U")

    expected_e = sp.Poly(
        21252176198679866250754006276839556755825 * lam**7
        + 799311827675117522149997435401077574131600 * lam**6
        + 14384863896403857958176347858723924433398460 * lam**5
        + 142730433788981669223548142320603830110956220 * lam**4
        + 786944231209420107856657052701244375708027892 * lam**3
        + 1852564642916723756803328543705267257790149632 * lam**2
        - 147733098192443646925107791876239203619548432 * lam
        + 3714269896529642422852685702214695613036368,
        lam,
    )
    check(eliminant == expected_e, "row 11 eliminant exact displayed coefficients")

    # The actual row-eleven bracket vanishes on (H,S); append determinant row 11.
    jointgb = sp.groebner([h10, s11], eta, Phi, order="lex")
    selected10_a = [zz(row.subs(section11)) for row in depth10_a]
    selected10_c = [zz(row.subs(section11)) for row in depth10_c]
    bracket11_error = zz(predicted_g(11, selected10_a, selected10_c) - g11)
    for degree in range(12):
        num = sp.fraction(zz(xcoef(bracket11_error, degree)))[0]
        check(jointgb.reduce(num)[1] == 0, f"row 11 bracket coefficient {degree} modulo H,S")
    a11base, c11base = determinant_particular(11, brow(11, selected10_a, selected10_c))
    th11 = list(sp.symbols("theta11_0:12"))
    theta11 = sum(th11[j] * x**j for j in range(12))
    row11_a = selected10_a + [sp.expand(a11base + theta11 * sp.diff(A0, x))]
    row11_c = selected10_c + [sp.expand(c11base + theta11 * sp.diff(C0, x))]
    aeq11, amat11 = depth_equations(row11_a, 2, 11)
    ceq11, cmat11 = depth_equations(row11_c, 3, 11)
    check((amat11.shape, amat11.rank(), cmat11.shape, cmat11.rank()) ==
          ((102, 228), 77, (114, 361), 94), "row 11 complete depth universes")
    acoef11, _ = sp.linear_eq_to_matrix(aeq11, th11)
    ccoef11, _ = sp.linear_eq_to_matrix(ceq11, th11)
    check((acoef11.rank(), ccoef11.rank()) == (4, 3), "row 11 separate terminal ranks")
    terminal11, dmat11, drhs11 = solve_affine(aeq11 + ceq11, th11)
    check(dmat11.rank() == 4 and tuple(dmat11.rref()[1]) == (8, 9, 10, 11),
          "row 11 joint depth rank and pivots")
    rem11d = [zz(eq.subs(terminal11)) for eq in aeq11 + ceq11]
    nonzero_depth = []
    for value in rem11d:
        if value == 0:
            continue
        num = sp.fraction(value)[0]
        remainder = sp.factor(jointgb.reduce(num)[1])
        if remainder != 0:
            nonzero_depth.append(remainder)
    check(not nonzero_depth, "all 45 row 11 depth residuals vanish modulo H,S")

    # Row 12 Student residual and its sign-quotient pullback.
    top11 = {th11[j]: terminal11[th11[j]] for j in (8, 9, 10, 11)}
    depth11_a = [zz(row.subs(top11)) for row in row11_a]
    depth11_c = [zz(row.subs(top11)) for row in row11_c]
    g12 = zz(grows[12].subs(solved_source).subs(gate).subs(terminal8).subs(alpha11, alpha_solution).subs(surface))
    diff12 = zz(g12 - predicted_g(12, depth11_a, depth11_c))
    student12 = student_value(12, diff12)
    check(not (set(th11[:8]) & student12.free_symbols), "row 12 Student independent of free tangent")

    def signed_pull(expr: sp.Expr) -> tuple[int, sp.Poly]:
        primitive = primitive_expr(expr, (Phi, eta))
        source = sp.Poly(sp.expand(primitive.subs(eta, lam * Phi)), Phi)
        parities = {degree % 2 for (degree,), _ in source.terms()}
        check(len(parities) == 1, "signed Phi parity")
        parity = next(iter(parities))
        ly = sum(coef * yy ** ((degree - parity) // 2)
                 for (degree,), coef in source.terms())
        pulled = sp.Poly(sp.fraction(zz(ly.subs(yy, ylam)))[0], lam)
        return parity, integer_primitive(pulled)

    parity12, t12 = signed_pull(student12)
    check(t12.degree() == 7, "row 12 Student pull degree seven")
    check(sp.gcd(eliminant, t12).degree() == 0, "row 12 Student coprime over Q")

    # Explicit compact F_29 certificate, checked both against the derived
    # reductions and by direct multiplication.
    E29 = sp.Poly(8*lam**7 + 9*lam**6 - 10*lam**5 + 12*lam**4 + 2*lam**3 + 8*lam**2 - 12*lam + 3, lam, modulus=29)
    T29 = sp.Poly(8*lam**7 - 13*lam**6 + 9*lam**5 - 2*lam**4 + 2*lam**3 + 12*lam**2 - 2*lam + 8, lam, modulus=29)
    a29 = sp.Poly(-6*lam**6 - 10*lam**5 - 13*lam**4 - 14*lam**3 + lam**2 - 11*lam - 8, lam, modulus=29)
    b29 = sp.Poly(6*lam**6 + 12*lam**5 - 14*lam**4 - 6*lam**3 + 8*lam**2 - 8*lam + 14, lam, modulus=29)
    derived_e29 = sp.Poly(eliminant.as_expr(), lam, modulus=29)
    derived_t29 = sp.Poly(t12.as_expr(), lam, modulus=29)
    check(derived_e29 == E29, "derived E reduction equals E29")
    check(derived_t29 == T29, "derived T reduction equals T29")
    check(a29 * E29 + b29 * T29 == sp.Poly(1, lam, modulus=29), "F29 Bezout identity")
    check(sp.gcd(E29, E29.diff()) == sp.Poly(1, lam, modulus=29), "E29 squarefree")
    check(sp.gcd(E29, T29) == sp.Poly(1, lam, modulus=29), "E29 T29 coprime")
    A29 = sp.Poly(-9*lam**3 - 4*lam**2 - 11*lam - 5, lam, modulus=29)
    B29 = sp.Poly(13 - 8*lam, lam, modulus=29)
    U29 = sp.Poly(5*lam**3 - 2*lam**2 + 1, lam, modulus=29)
    check(sp.Poly(alam.as_expr(), lam, modulus=29) == A29, "derived A reduction equals A29")
    check(sp.Poly(blam.as_expr(), lam, modulus=29) == B29, "derived B reduction equals B29")
    derived_u29 = sp.Poly(ulam.as_expr(), lam, modulus=29)
    check(sp.monic(derived_u29) == sp.monic(U29),
          f"derived U reduction equals U29 up to unit: {derived_u29.as_expr()}")
    for label, forbidden in (("A29", A29), ("B29", B29), ("U29", U29)):
        check(sp.gcd(E29, forbidden) == sp.Poly(1, lam, modulus=29), f"E29 avoids {label}")
    # The promoted theorem uses a smaller row-twelve residual J29.  Verify
    # its published Bezout identity and show that our independently selected
    # Student residual differs from it by a unit on the E29 carrier.
    J29 = sp.Poly(-11*lam**2 + 8*lam + 5, lam, modulus=29)
    jbez_e = sp.Poly(11*lam + 11, lam, modulus=29)
    jbez_j = sp.Poly(8*lam**6 + 7*lam**5 + 13*lam**4 + 12*lam**3
                     - 3*lam**2 + 8*lam + 11, lam, modulus=29)
    check(jbez_e * E29 + jbez_j * J29 == sp.Poly(1, lam, modulus=29),
          "canonical E29 J29 Bezout identity")
    R29 = sp.Poly(9*lam**6 - 4*lam**5 - 5*lam**4 - 5*lam**3
                  + 5*lam**2 - 14*lam + 6, lam, modulus=29)
    check(sp.rem(T29 - R29 * J29, E29) == sp.Poly(0, lam, modulus=29),
          "Student T29 equals unit R29 times J29 on E29")
    check(sp.gcd(E29, R29) == sp.Poly(1, lam, modulus=29),
          "R29 is a unit on E29")

    # Verify full row-twelve selector rank and that Student is one of its
    # source residuals; extinction already follows from gcd(E,T)=1.
    eq12br = [xcoef(diff12, j) for j in range(13)]
    section12, bmat12, brhs12 = solve_affine(eq12br, th11[:8])
    check(bmat12.rank() == 8, "row 12 restricted bracket tangent rank")
    rem12 = [zz(eq.subs(section12)) for eq in eq12br]
    nonzero12 = sum(value != 0 for value in rem12)
    check(nonzero12 > 0, "row 12 selector has nonzero source residual")

    print("THM4380 CLEANROOM PASS")
    print(f"checks={CHECKS}")
    print("source_rows=G9,G10,G11,G12 direct from H(p,y)")
    print(f"depth_universes=row8:{amat8.shape}/{amat8.rank()},{cmat8.shape}/{cmat8.rank()} row9:{amat9.shape}/{amat9.rank()},{cmat9.shape}/{cmat9.rank()} row10:{amat10.shape}/{amat10.rank()},{cmat10.shape}/{cmat10.rank()} row11:{amat11.shape}/{amat11.rank()},{cmat11.shape}/{cmat11.rank()}")
    print(f"row10_depth_H={sp.sstr(h10)}")
    print(f"row11_source_S={sp.sstr(s11)}")
    print(f"row11_eliminant={sp.sstr(eliminant.as_expr())}")
    print(f"row11_forbidden_gcds=squarefree:{sp.gcd(eliminant, eliminant.diff()).degree()} A:{sp.gcd(eliminant, alam).degree()} B:{sp.gcd(eliminant, blam).degree()} U:{sp.gcd(eliminant, ulam).degree()}")
    print("row11_geometric_points=14 (7 simple lambda roots, nonzero Y, two Phi sign sheets)")
    print(f"row11_depth=annihilator_rows:{len(aeq11)+len(ceq11)} ranks_P2_P3_joint:{acoef11.rank()}/{ccoef11.rank()}/{dmat11.rank()} pivots:{dmat11.rref()[1]} residuals_mod_HS:{len(nonzero_depth)} affine_fibre:8")
    print(f"row12_student_parity={parity12} degree={t12.degree()} gcd_degree={sp.gcd(eliminant,t12).degree()}")
    print(f"E29={sp.sstr(E29.as_expr())}")
    print(f"T29={sp.sstr(T29.as_expr())}")
    print(f"a29={sp.sstr(a29.as_expr())}")
    print(f"b29={sp.sstr(b29.as_expr())}")
    print("bezout29=a29*E29+b29*T29=1")
    print(f"A29={sp.sstr(A29.as_expr())}")
    print(f"B29={sp.sstr(B29.as_expr())}")
    print(f"U29={sp.sstr(U29.as_expr())}")
    print(f"J29={sp.sstr(J29.as_expr())}")
    print("canonical_bezout29=(11*lambda+11)*E29+(8*lambda**6+7*lambda**5+13*lambda**4+12*lambda**3-3*lambda**2+8*lambda+11)*J29=1")
    print(f"R29={sp.sstr(R29.as_expr())}")
    print("student_vs_canonical=T29=R29*J29 mod E29; gcd(E29,R29)=1")
    print(f"row12_selector=rank:{bmat12.rank()} nonzero_residuals:{nonzero12}")
    print(f"u_positive_phi1_eta0={row10_positive}")
    print("verdict=all 14 genuine row-11 bracket+joint-depth points die at row 12")
    print("scope=FINITE-EXACT; residual_weight<=12; source-normal D(Phi*U); no chart/seam entry, all-row membership, termination, Keller pair, JC2, or DC2")


if __name__ == "__main__":
    main()
