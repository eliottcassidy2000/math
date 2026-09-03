#!/usr/bin/env python3
"""Clean-room exact audit for the proposed THM-4360.

This file deliberately imports no repository computation.  It rebuilds the
source-normal rows, determinant recursion, projected depth matrices, and the
row-nine/ten tests from the displayed formulas in THM-4308 and THM-4315.
"""

from __future__ import annotations

from math import gcd
import sys

import sympy as s


sys.stdout.reconfigure(newline="\n")

checks = 0


def require(ok: bool, label: str) -> None:
    global checks
    if not ok:
        raise AssertionError(label)
    checks += 1


def z(expr: s.Expr) -> bool:
    return s.cancel(s.together(s.expand(expr))) == 0


def require_zero(expr: s.Expr, label: str) -> None:
    require(z(expr), label)


x, t = s.symbols("x t")
eta, alpha, beta = s.symbols("eta alpha11 beta11")

Delta = s.Rational(896, 15)
Theta = s.Rational(512, 75)
K = -s.Rational(32, 5)
upsilon = -s.Rational(731648, 2025)
Phi = s.Integer(0)
zeta = s.Integer(0)
xi = s.Rational(1563264759296, 115860375)

U = s.cancel((s.Integer(475515904) - 109350 * xi) / 200475)
W = s.cancel(-(s.Integer(4343625) * Phi**2 - 17172000 * xi + 143826305024) / 4009500)
Z = s.cancel(
    (
        s.Integer(12506118074368)
        - 173745000 * Phi**2
        - 195463125 * Phi * eta
        - 926883000 * xi
    )
    / 108256500
)

p = t * (1 + x**2 * t)
y = x * t * p
u = x**2 * t
H = s.expand(
    -3 * p
    + s.Rational(8, 3) * p**2
    - s.Rational(1376, 135) * p**3
    + K * y**2
    + Phi * p**2 * y
    + Delta * p**4
    + Theta * p * y**2
    + eta * p**3 * y
    + zeta * y**3
    + upsilon * p**5
    + xi * p**2 * y**2
    + alpha * p**4 * y
    + beta * p * y**3
    + U * p**6
    + W * p**3 * y**2
    + Z * y**4
)

A0 = 1 + x**2 / 4
C0 = -3 * x / 4 - x**3 / 8
A1 = s.Rational(4, 3) + 2 * x**2
C1 = -4 * x - s.Rational(3, 2) * x**3
A2 = -s.Rational(32, 9) - s.Rational(4, 5) * x**2
C2 = s.Rational(88, 15) * x - s.Rational(12, 5) * x**3
A3 = s.Rational(2176, 135) + (s.Rational(1088, 315) - 4 * K / 7) * x**2 - s.Rational(32, 15) * x**4
C3 = (-s.Rational(8128, 315) + 6 * K / 7) * x + (s.Rational(736, 105) + 3 * K / 7) * x**3 + s.Rational(8, 5) * x**5

qboundary = -(x**2 + 6) / 2


def xc(expr: s.Expr, degree: int) -> s.Expr:
    return s.expand(expr).coeff(x, degree)


def tc(expr: s.Expr, degree: int) -> s.Expr:
    return s.expand(expr).coeff(t, degree)


def brow(m: int, aa: list[s.Expr], cc: list[s.Expr]) -> s.Expr:
    return s.expand(
        sum(
            (m - i) * s.diff(aa[i], x) * cc[m - i]
            - i * aa[i] * s.diff(cc[m - i], x)
            for i in range(1, m)
        )
    )


def trow(m: int, aa: list[s.Expr], cc: list[s.Expr]) -> s.Expr:
    quadratic = sum(cc[i] * cc[m - i] for i in range(1, m))
    cubic = sum(
        aa[i] * aa[j] * aa[k]
        for i in range(m)
        for j in range(m)
        for k in range(m)
        if i + j + k == m
    )
    return s.expand(quadratic - cubic)


def predicted(m: int, aa: list[s.Expr], cc: list[s.Expr]) -> s.Expr:
    return s.expand(trow(m, aa, cc) - qboundary * brow(m, aa, cc) / m)


def determinant_particular(m: int, b: s.Expr) -> tuple[s.Expr, s.Expr]:
    av = s.symbols(f"a{m}_0:{m + 2}")
    cv = s.symbols(f"c{m}_0:{m + 3}")
    ap = sum(av[j] * x**j for j in range(m + 2))
    cp = sum(cv[j] * x**j for j in range(m + 3))
    expr = s.expand(m * (s.diff(A0, x) * cp - s.diff(C0, x) * ap))
    equations = [xc(expr, j) for j in range(m + 4)]
    matrix, _ = s.linear_eq_to_matrix(equations, av + cv)
    require(matrix.rank() == m + 4, f"determinant rank row {m}")
    pivots = matrix.rref()[1]
    pm = matrix[:, list(pivots)]
    rhs = s.Matrix([-xc(b, j) for j in range(m + 4)])
    vals = pm.inv() * rhs
    allvals = [s.Integer(0)] * (2 * m + 5)
    for col, value in zip(pivots, vals):
        allvals[col] = s.cancel(value)
    anew = s.expand(sum(allvals[j] * x**j for j in range(m + 2)))
    cnew = s.expand(sum(allvals[m + 2 + j] * x**j for j in range(m + 3)))
    require_zero(m * (s.diff(A0, x) * cnew - s.diff(C0, x) * anew) + b, f"determinant solution row {m}")
    return anew, cnew


def tangent_solve(m: int, target: s.Expr) -> s.Expr:
    tv = s.symbols(f"r{m}_0:{m}")
    theta = sum(tv[j] * x**j for j in range(m))
    image = s.expand(s.diff(qboundary, x) * theta - qboundary * s.diff(theta, x) / m)
    eqs = [xc(image, j) for j in range(m + 1)]
    matrix, _ = s.linear_eq_to_matrix(eqs, tv)
    require(matrix.rank() == m, f"tangent rank row {m}")
    rhs = s.Matrix([xc(target, j) for j in range(m + 1)])
    row_pivots = matrix.T.rref()[1]
    vals = matrix[list(row_pivots), :].inv() * rhs[list(row_pivots), :]
    require(all(z((matrix * vals - rhs)[j]) for j in range(m + 1)), f"tangent consistency row {m}")
    return s.expand(sum(vals[j] * x**j for j in range(m)))


def build_rows_to_eight() -> tuple[list[s.Expr], list[s.Expr], tuple[s.Symbol, ...]]:
    aa = [A0, A1, A2, A3]
    cc = [C0, C1, C2, C3]
    for n in range(4, 8):
        abase, cbase = determinant_particular(n, brow(n, aa, cc))
        trial_a = aa + [abase]
        trial_c = cc + [cbase]
        target = s.expand(tc(H, n + 1) - predicted(n + 1, trial_a, trial_c))
        theta = tangent_solve(n + 1, target)
        aa.append(s.expand(abase + theta * s.diff(A0, x)))
        cc.append(s.expand(cbase + theta * s.diff(C0, x)))
        require_zero(predicted(n + 1, aa, cc) - tc(H, n + 1), f"source match G{n + 1}")
    abase, cbase = determinant_particular(8, brow(8, aa, cc))
    r8 = s.symbols("r8_0:9")
    theta8 = sum(r8[j] * x**j for j in range(9))
    aa.append(s.expand(abase + theta8 * s.diff(A0, x)))
    cc.append(s.expand(cbase + theta8 * s.diff(C0, x)))
    require_zero(predicted(8, aa[:8], cc[:8]) - tc(H, 8), "source match G8")
    return aa, cc, r8


def depth_matrix(depth: int, maxrow: int) -> tuple[list[tuple[int, int]], s.Matrix]:
    coords = [(n, r) for n in range(maxrow + 1) for r in range(n + depth + 1)]
    columns: list[tuple[int, int, int, int]] = []
    for b in range(depth + 1):
        for a in range(depth - b + 1):
            for e in range(maxrow // 2 + 1):
                for c in range(maxrow + 1):
                    if b + c + 2 * e <= maxrow:
                        columns.append((a, b, c, e))
    mat = s.zeros(len(coords), len(columns))
    location = {coord: i for i, coord in enumerate(coords)}
    for col, (a, b, c, e) in enumerate(columns):
        for j in range(c + e + 1):
            n = b + c + 2 * e + j
            r = a + 2 * b + e + 2 * j
            if n <= maxrow:
                mat[location[(n, r)], col] += s.binomial(c + e, j)
    return coords, mat


def jet(rows: list[s.Expr], coords: list[tuple[int, int]]) -> s.Matrix:
    return s.Matrix([xc(rows[n], r) for n, r in coords])


F = s.expand(s.Integer(2393096045625) * eta**2 - s.Integer(415832184456871936))
NE = s.expand(s.Integer(1010418330375) * alpha * eta + s.Integer(619241095293435904))
NH = s.expand(s.Integer(5052091651875) * alpha * eta + s.Integer(3193963923683016704))


def mod_f(expr: s.Expr) -> s.Expr:
    val = s.cancel(s.together(expr))
    num, den = s.fraction(val)
    require(s.degree(den, eta) <= 0, "F reduction denominator independent of eta")
    remainder = s.Poly(num, eta, domain="EX").rem(s.Poly(F, eta, domain="EX")).as_expr()
    return s.cancel(remainder / den)


def solve_linear_subset(equations: list[s.Expr], variables: tuple[s.Symbol, ...]) -> tuple[dict[s.Symbol, s.Expr], s.Matrix, s.Matrix, tuple[int, ...]]:
    matrix, rhs = s.linear_eq_to_matrix(equations, variables)
    row_pivots = matrix.T.rref()[1]
    require(len(row_pivots) == matrix.rank(), "independent row count")
    square = matrix[list(row_pivots), :]
    require(square.rows == square.cols == len(variables), "square selected system")
    vals = square.inv() * rhs[list(row_pivots), :]
    return {var: s.cancel(vals[i]) for i, var in enumerate(variables)}, matrix, rhs, row_pivots


def student_row(m: int) -> list[int]:
    moments: list[s.Rational] = []
    for degree in range(m + 1):
        if degree % 2:
            moments.append(s.Integer(0))
            continue
        r = degree // 2
        value = s.Integer(1)
        for j in range(1, r + 1):
            value *= s.Rational(6 * (2 * j - 1), 2 * m - 2 * j + 1)
        moments.append(s.cancel(value))
    common_den = 1
    for value in moments:
        common_den = s.ilcm(common_den, int(s.denom(value)))
    entries = [int(value * common_den) for value in moments]
    common_num = 0
    for entry in entries:
        common_num = gcd(common_num, abs(entry))
    return [entry // common_num for entry in entries]


def main() -> None:
    aa, cc, r8 = build_rows_to_eight()

    a8coords, p2_8 = depth_matrix(2, 8)
    c8coords, p3_8 = depth_matrix(3, 8)
    require((len(a8coords), p2_8.cols, p2_8.rank()) == (63, 131, 51), "P2 row8 universe")
    require((len(c8coords), p3_8.cols, p3_8.rank()) == (72, 204, 63), "P3 row8 universe")
    old_residuals = [s.expand((v.T * jet(aa, a8coords))[0]) for v in p2_8.T.nullspace()]
    old_residuals += [s.expand((v.T * jet(cc, c8coords))[0]) for v in p3_8.T.nullspace()]
    m8, b8 = s.linear_eq_to_matrix(old_residuals, r8)
    require((len(old_residuals), m8.rank(), m8.row_join(b8).rank()) == (21, 2, 2), "row8 depth affine fibre")
    solved8 = s.solve(old_residuals, r8, dict=True)
    require(len(solved8) == 1 and len(solved8[0]) == 2, "row8 two pivot coordinates")
    pivot8 = solved8[0]
    free8 = tuple(v for v in r8 if v not in pivot8)
    require(len(free8) == 7, "row8 seven free tangent coordinates")
    aa8depth = [s.expand(row.subs(pivot8)) for row in aa]
    cc8depth = [s.expand(row.subs(pivot8)) for row in cc]
    require(all(z(r.subs(pivot8)) for r in old_residuals), "all row8 depth residuals")

    expected_g9 = (
        (20 * U + 10 * W + 4 * Z) * x**6
        + (10 * alpha + 6 * beta) * x**7
        + (5 * upsilon + 4 * xi) * x**8
        + eta * x**9
    )
    g9 = tc(H, 9)
    require_zero(g9 - expected_g9, "literal G9")
    defect9 = s.expand(g9 - predicted(9, aa8depth, cc8depth))
    equations9 = [xc(defect9, j) for j in range(10)]
    m9, b9lin = s.linear_eq_to_matrix(equations9, free8)
    require(m9.rank() == 7, "old G9 tangent rank seven")
    require(m9.row_join(b9lin).rank() == 8, "old G9 generic augmented rank eight")
    sr9 = student_row(9)
    require(sr9 == [12155, 0, 4290, 0, 5148, 0, 11880, 0, 45360, 0], "row9 Student row")
    e9_student = s.factor(sum(sr9[j] * equations9[j] for j in range(10)))
    require_zero(e9_student - s.Rational(143, 56308142250) * F, "literal row9 Student normalization")
    e9 = s.factor(-s.Rational(5, 39) * 328050 * e9_student)
    require_zero(e9 + s.Rational(11, 102987) * F, "row9 F obstruction")
    sol8, _, _, rows9 = solve_linear_subset(equations9, free8)
    residual9 = [s.cancel(eq.subs(sol8)) for eq in equations9]
    require(all(z(mod_f(value)) for value in residual9), "G9 consistency modulo F")
    require(s.Matrix([mod_f(v) for v in residual9]).rank() == 0, "no extra G9 source equation")
    aa8fixed = [s.expand(row.subs(sol8)) for row in aa8depth]
    cc8fixed = [s.expand(row.subs(sol8)) for row in cc8depth]

    b9det = brow(9, aa8fixed, cc8fixed)
    a9base, c9base = determinant_particular(9, b9det)
    qv = s.symbols("q0:10")
    theta9 = sum(qv[j] * x**j for j in range(10))
    a9 = s.expand(a9base + theta9 * s.diff(A0, x))
    c9 = s.expand(c9base + theta9 * s.diff(C0, x))
    aa9 = aa8fixed + [a9]
    cc9 = cc8fixed + [c9]

    a9coords, p2_9 = depth_matrix(2, 9)
    c9coords, p3_9 = depth_matrix(3, 9)
    require((len(a9coords), p2_9.cols, p2_9.rank()) == (75, 160, 59), "P2 row9 universe")
    require((len(c9coords), p3_9.cols, p3_9.rank()) == (85, 251, 73), "P3 row9 universe")
    depth9 = [s.expand((v.T * jet(aa9, a9coords))[0]) for v in p2_9.T.nullspace()]
    depth9 += [s.expand((v.T * jet(cc9, c9coords))[0]) for v in p3_9.T.nullspace()]
    require(len(depth9) == 28, "row9 28 depth residuals")
    depth9mod = [mod_f(v) for v in depth9]
    md, bd = s.linear_eq_to_matrix(depth9mod, qv)
    require(md.rank() == 3, "row9 depth tangent rank three")
    require(md.row_join(bd).rank() == 3, "row9 depth augmented rank three")
    require(md.rref()[1] == (7, 8, 9), "row9 depth pivots q7 q8 q9")
    depth_particular = {
        **{qv[j]: s.Integer(0) for j in range(7)},
        qv[7]: -s.Rational(317075357581312, 347581125),
        qv[8]: -(125 * alpha + 100 * beta + 198 * eta) / 15,
        qv[9]: -s.Rational(553237643264, 23172075),
    }
    require(all(z(mod_f(v.subs(depth_particular))) for v in depth9), "all row9 depth residuals modulo F")

    expected_g10 = (
        (15 * U + 10 * W + 6 * Z) * x**8
        + (5 * alpha + 4 * beta) * x**9
        + (upsilon + xi) * x**10
    )
    g10 = tc(H, 10)
    require_zero(g10 - expected_g10, "literal G10")
    defect10 = s.expand(g10 - predicted(10, aa9, cc9))
    equations10 = [xc(defect10, j) for j in range(11)]
    sr10 = student_row(10)
    require(sr10 == [46189, 0, 14586, 0, 15444, 0, 30888, 0, 99792, 0, 489888], "row10 Student row")
    e10raw = s.factor(sum(sr10[j] * equations10[j] for j in range(11)))
    expected_raw_numerator = (
        27281294920125 * alpha * eta
        + 241702700608125 * eta**2
        - 25279541057221296128
    )
    require_zero(e10raw - s.Rational(143, 84462213375) * expected_raw_numerator, "raw E10")
    e10 = mod_f(e10raw)
    require_zero(e10 - s.Rational(143, 3128230125) * NE, "E10 modulo F")
    require_zero(s.diff(e10raw, beta), "raw beta cancellation")

    m10, b10 = s.linear_eq_to_matrix(equations10, qv)
    require(m10.rank() == 10, "full row10 tangent rank ten")
    require(m10.row_join(b10).rank() == 11, "full row10 generic augmented rank eleven")
    qselected, _, _, rows10 = solve_linear_subset(equations10, qv)
    selected_defect = [mod_f(eq.subs(qselected)) for eq in equations10]
    alpha_on_ne = -s.Rational(619241095293435904, 1010418330375) / eta
    require(all(z(mod_f(value.subs(alpha, alpha_on_ne))) for value in selected_defect), "row10 consistency on F and NE")

    dq_dbeta = tuple(s.cancel(s.diff(qselected[q], beta)) for q in qv)
    expected_translation = (
        s.Rational(2384, 945),
        0,
        -s.Rational(504928, 8505),
        0,
        s.Rational(9484, 105),
        0,
        -s.Rational(568, 7),
        0,
        -s.Rational(20, 3),
        0,
    )
    require(dq_dbeta == expected_translation, "beta translation derivative")

    selected_depth_raw = [s.cancel(value.subs(qselected)) for value in depth9]
    selected_depth = [mod_f(value) for value in selected_depth_raw]
    nonzero_indices = [i for i, value in enumerate(selected_depth) if not z(value)]
    require(nonzero_indices == [12, 14, 26], "only three selected depth residuals")
    require_zero(selected_depth[12] - NH / 14189651847000, "selected residual index 12")
    require_zero(selected_depth[14] - s.Rational(13, 153248239947600) * NE, "selected residual index 14")
    require_zero(selected_depth[26] + s.Rational(13, 204330986596800) * NE, "selected residual index 26")
    require(all(z(s.diff(value, beta)) for value in selected_depth), "selected depth residuals beta-invariant")
    require(all(z(s.diff(value, beta)) for value in selected_depth_raw), "selected raw depth residuals beta-invariant")

    ha_coefficients = {(5, 0): 35, (6, 2): -20, (7, 4): 10, (8, 6): -4, (9, 8): 1}
    ha_row = s.Matrix([[ha_coefficients.get(coord, 0) for coord in a9coords]])
    require(ha_row * p2_9 == s.zeros(1, p2_9.cols), "H_A annihilates all P2 columns")
    ha_content = 0
    for coefficient in ha_coefficients.values():
        ha_content = gcd(ha_content, abs(coefficient))
    require(ha_content == 1, "H_A primitive")
    ha = 35 * xc(aa9[5], 0) - 20 * xc(aa9[6], 2) + 10 * xc(aa9[7], 4) - 4 * xc(aa9[8], 6) + xc(aa9[9], 8)
    ha_selected = mod_f(ha.subs(qselected))
    require_zero(ha_selected - NH / 14189651847000, "primitive H_A evaluation")
    require_zero(NH - 5 * NE - s.Integer(97758447215837184), "unit witness")
    ha_on_scalar_gate = s.cancel(s.Integer(97758447215837184) / s.Integer(14189651847000))
    require(ha_on_scalar_gate == s.Rational(9854451712, 1430375), "forced H_A nonzero value")

    fpoly = s.Poly(F, eta, domain=s.QQ)
    require(fpoly.is_irreducible, "F irreducible over Q")
    require(s.gcd(fpoly, fpoly.diff()) == s.Poly(1, eta, domain=s.QQ), "F squarefree")
    require(fpoly.eval(0) != 0, "F roots nonzero")
    require(s.degree(F, eta) == 2, "F has two geometric roots")

    print("THM-4360 CLEAN-ROOM EXACT REFEREE: PASS")
    print(f"U={U}; W={W}; Z={Z}")
    print(f"G9_OLD_TANGENT rank={m9.rank()} generic_augmented_rank={m9.row_join(b9lin).rank()} selected_rows={rows9}")
    print(f"F={F}")
    print(f"E9_NAMED={e9}")
    print(f"E9_LITERAL_STUDENT={e9_student}")
    print("E9_NORMALIZATION named and literal are nonzero rational multiples with identical zero locus")
    print(f"ROW9_DEPTH P2=(75,160,59,16) P3=(85,251,73,12) tangent_rank={md.rank()} pivots={md.rref()[1]}")
    print(f"ROW9_DEPTH_PARTICULAR q7={depth_particular[qv[7]]}; q8={depth_particular[qv[8]]}; q9={depth_particular[qv[9]]}")
    print(f"E10_RAW={e10raw}")
    print(f"E10_MOD_F={e10}")
    print(f"ROW10 tangent_rank={m10.rank()} generic_augmented_rank={m10.row_join(b10).rank()} selected_rows={rows10}")
    print(f"NONZERO_DEPTH_INDICES={nonzero_indices}")
    print(f"H_A={ha_selected}")
    print(f"UNIT={s.expand(NH - 5 * NE)}")
    print(f"BETA_TRANSLATION={dq_dbeta}")
    print("GEOMETRY row9_closure=two_A2 row9_Dbeta=two_(A1xGm); scalar_row10_closure=two_A1 scalar_row10_Dbeta=two_Gm; depth_consumer=empty")
    print(f"CHECKS={checks}")


if __name__ == "__main__":
    main()

