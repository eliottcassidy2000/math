#!/usr/bin/env python3
"""Clean-room exact referee for the proposed S_4339 row-9/10 result.

This script deliberately does not import the exploratory certificate.  It
reconstructs the source rows, bracket rows, depth-projection matrices, and
Student cokernels directly from the displayed definitions in THM-4308 and
THM-4315.
"""

from __future__ import annotations

from math import gcd, isqrt
import sys

import sympy as sp

sys.stdout.reconfigure(newline="\n")


x, t = sp.symbols("x t")
eta, alpha = sp.symbols("eta alpha11")

DELTA = sp.Rational(896, 15)
THETA = sp.Rational(512, 75)
K = -sp.Rational(32, 5)
UPSILON = -sp.Rational(731648, 2025)
XI = sp.Rational(1563264759296, 115860375)
U = -sp.Rational(5200877686784, 1042743375)
W = sp.Rational(2539122786304, 115860375)
Z = sp.Integer(0)

F = 2393096045625 * eta**2 - 415832184456871936
E10_NUM = 1010418330375 * alpha * eta + 619241095293435904
H_NUM = 5052091651875 * alpha * eta + 3193963923683016704

checks = 0


def require(condition: bool, label: str) -> None:
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1


def zero(value: sp.Expr) -> bool:
    return sp.cancel(sp.together(sp.expand(value))) == 0


def require_zero(value: sp.Expr, label: str) -> None:
    require(zero(value), label)


def xc(poly: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(poly).coeff(x, degree)


def tc(poly: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(poly).coeff(t, degree)


def quotient_reduce(value: sp.Expr) -> sp.Expr:
    """Canonical eta-degree < 2 representative modulo F."""
    numerator, denominator = sp.together(sp.expand(value)).as_numer_denom()
    remainder = sp.rem(numerator, F, eta)
    return sp.cancel(remainder / denominator)


def reduce_x(poly: sp.Expr) -> sp.Expr:
    expanded = sp.Poly(sp.expand(poly), x)
    return sp.expand(
        sum(quotient_reduce(coefficient) * x**degree[0] for degree, coefficient in expanded.terms())
    )


def curve_zero(value: sp.Expr) -> bool:
    return zero(quotient_reduce(value))


def P(a: sp.Expr, c: sp.Expr) -> sp.Expr:
    return sp.expand(c**2 - a**3 + sp.Rational(3, 4) * a + sp.Rational(1, 4))


p = t * (1 + x**2 * t)
y = x * t * p
H = sp.expand(
    -3 * p
    + sp.Rational(8, 3) * p**2
    - sp.Rational(1376, 135) * p**3
    + K * y**2
    + DELTA * p**4
    + THETA * p * y**2
    + eta * p**3 * y
    + UPSILON * p**5
    + XI * p**2 * y**2
    + alpha * p**4 * y
    + U * p**6
    + W * p**3 * y**2
)

G = {m: tc(H, m) for m in range(4, 11)}

A0 = 1 + x**2 / 4
C0 = -3 * x / 4 - x**3 / 8
A1 = sp.Rational(4, 3) + 2 * x**2
C1 = -4 * x - sp.Rational(3, 2) * x**3
A2 = -sp.Rational(32, 9) - sp.Rational(4, 5) * x**2
C2 = sp.Rational(88, 15) * x - sp.Rational(12, 5) * x**3
A3 = sp.Rational(2176, 135) + (sp.Rational(1088, 315) - 4 * K / 7) * x**2 - sp.Rational(32, 15) * x**4
C3 = (-sp.Rational(8128, 315) + 6 * K / 7) * x + (sp.Rational(736, 105) + 3 * K / 7) * x**3 + sp.Rational(8, 5) * x**5

q = -(x**2 + 6) / 2


def brow(m: int, aa: list[sp.Expr], cc: list[sp.Expr]) -> sp.Expr:
    return sp.expand(
        sum(
            (m - i) * sp.diff(aa[i], x) * cc[m - i]
            - i * aa[i] * sp.diff(cc[m - i], x)
            for i in range(1, m)
        )
    )


def trow(m: int, aa: list[sp.Expr], cc: list[sp.Expr]) -> sp.Expr:
    quadratic = sum(cc[i] * cc[m - i] for i in range(1, m))
    cubic = sum(
        aa[i] * aa[j] * aa[k]
        for i in range(m)
        for j in range(m)
        for k in range(m)
        if i + j + k == m
    )
    return sp.expand(quadratic - cubic)


def predicted(m: int, aa: list[sp.Expr], cc: list[sp.Expr]) -> sp.Expr:
    return sp.expand(trow(m, aa, cc) - q * brow(m, aa, cc) / m)


def determinant_particular(m: int, forcing: sp.Expr, *, on_curve: bool = False) -> tuple[sp.Expr, sp.Expr]:
    av = list(sp.symbols(f"a{m}_0:{m + 2}"))
    cv = list(sp.symbols(f"c{m}_0:{m + 3}"))
    ap = sum(av[j] * x**j for j in range(m + 2))
    cp = sum(cv[j] * x**j for j in range(m + 3))
    expression = sp.expand(m * (sp.diff(A0, x) * cp - sp.diff(C0, x) * ap) + forcing)
    equations = [xc(expression, j) for j in range(m + 4)]
    matrix, rhs = sp.linear_eq_to_matrix(equations, av + cv)
    require(matrix.rank() == m + 4, f"determinant operator rank m={m}")
    pivots = matrix.rref()[1]
    square = matrix[:, list(pivots)]
    values = square.inv() * rhs
    all_values = [sp.Integer(0)] * (2 * m + 5)
    for column, value in zip(pivots, values):
        all_values[column] = quotient_reduce(value) if on_curve else sp.cancel(value)
    aresult = sp.expand(sum(all_values[j] * x**j for j in range(m + 2)))
    cresult = sp.expand(sum(all_values[m + 2 + j] * x**j for j in range(m + 3)))
    residual = sp.expand(m * (sp.diff(A0, x) * cresult - sp.diff(C0, x) * aresult) + forcing)
    if on_curve:
        require(all(curve_zero(xc(residual, j)) for j in range(m + 4)), f"determinant solve mod F m={m}")
    else:
        require_zero(residual, f"determinant solve m={m}")
    return aresult, cresult


def tangent_matrix(m: int) -> tuple[list[sp.Symbol], sp.Matrix]:
    variables = list(sp.symbols(f"s{m}_0:{m}"))
    theta = sum(variables[j] * x**j for j in range(m))
    image = sp.expand(-x * theta + (x**2 + 6) * sp.diff(theta, x) / (2 * m))
    matrix, _ = sp.linear_eq_to_matrix([xc(image, j) for j in range(m + 1)], variables)
    require(matrix.rank() == m, f"Student tangent rank m={m}")
    return variables, matrix


def primitive_student(m: int) -> list[int]:
    moments: list[sp.Rational] = []
    for degree in range(m + 1):
        if degree % 2:
            moments.append(sp.Integer(0))
            continue
        value = sp.Integer(1)
        for j in range(1, degree // 2 + 1):
            value *= sp.Rational(6 * (2 * j - 1), 2 * m - 2 * j + 1)
        moments.append(sp.cancel(value))
    lcm = 1
    for value in moments:
        lcm = sp.ilcm(lcm, int(sp.denom(value)))
    integers = [int(value * lcm) for value in moments]
    common = 0
    for value in integers:
        common = gcd(common, abs(value))
    result = [value // common for value in integers]
    _, matrix = tangent_matrix(m)
    require(sp.Matrix([result]) * matrix == sp.zeros(1, m), f"Student row annihilates m={m}")
    return result


def solve_injective(matrix: sp.Matrix, rhs: sp.Matrix, label: str, *, on_curve: bool = False) -> tuple[sp.Matrix, list[sp.Expr]]:
    rank = matrix.rank()
    require(rank == matrix.cols, f"{label} full column rank")
    pivot_rows = matrix.T.rref()[1]
    square = matrix[list(pivot_rows), :]
    solution = square.inv() * rhs[list(pivot_rows), :]
    if on_curve:
        solution = solution.applyfunc(quotient_reduce)
    residuals = [sp.expand((matrix * solution - rhs)[j]) for j in range(matrix.rows)]
    return solution, residuals


def depth_projection(depth: int, max_row: int) -> tuple[list[tuple[int, int]], list[tuple[int, int, int, int]], sp.Matrix]:
    coordinates = [(n, r) for n in range(max_row + 1) for r in range(n + depth + 1)]
    columns: list[tuple[int, int, int, int]] = []
    for b in range(depth + 1):
        for a in range(depth - b + 1):
            for e in range(max_row // 2 + 1):
                for c in range(max_row + 1):
                    if b + c + 2 * e <= max_row:
                        columns.append((a, b, c, e))
    matrix = sp.zeros(len(coordinates), len(columns))
    index = {coordinate: i for i, coordinate in enumerate(coordinates)}
    for column, (a, b, c, e) in enumerate(columns):
        for j in range(c + e + 1):
            n = b + c + 2 * e + j
            r = a + 2 * b + e + 2 * j
            if n <= max_row:
                matrix[index[(n, r)], column] += sp.binomial(c + e, j)
    return coordinates, columns, matrix


def jet(rows: list[sp.Expr], coordinates: list[tuple[int, int]]) -> sp.Matrix:
    return sp.Matrix([xc(rows[n], r) for n, r in coordinates])


def parameterize_linear_system(
    equations: list[sp.Expr], variables: list[sp.Symbol], label: str, expected_rank: int, *, on_curve: bool = False
) -> tuple[dict[sp.Symbol, sp.Expr], list[sp.Symbol], sp.Matrix, sp.Matrix]:
    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    require(matrix.rank() == expected_rank, f"{label} rank")
    pivot_columns = list(matrix.rref()[1])
    require(len(pivot_columns) == expected_rank, f"{label} pivot count")
    free_columns = [j for j in range(len(variables)) if j not in pivot_columns]
    pivot_rows = matrix[:, pivot_columns].T.rref()[1]
    square = matrix[list(pivot_rows), pivot_columns]
    adjusted = rhs[list(pivot_rows), :] - matrix[list(pivot_rows), free_columns] * sp.Matrix(
        [variables[j] for j in free_columns]
    )
    values = square.inv() * adjusted
    if on_curve:
        values = values.applyfunc(quotient_reduce)
    substitutions = {variables[column]: values[i] for i, column in enumerate(pivot_columns)}
    return substitutions, [variables[j] for j in free_columns], matrix, rhs


def main() -> None:
    # Literal source rows and fixed response constants.
    expected_g9 = (
        sp.Rational(24900699406336, 208548675) * x**6
        + 10 * alpha * x**7
        + sp.Rational(2014584278528, 38620125) * x**8
        + eta * x**9
    )
    expected_g10 = (
        sp.Rational(2006771806208, 13903245) * x**8
        + 5 * alpha * x**9
        + sp.Rational(1521403518976, 115860375) * x**10
    )
    require_zero(G[9] - expected_g9, "literal G9")
    require_zero(G[10] - expected_g10, "literal G10")

    # Rebuild rows four through seven without importing either canonical script.
    aa = [A0, A1, A2, A3]
    cc = [C0, C1, C2, C3]
    require_zero(predicted(4, aa, cc) - G[4], "inherited G4")
    for n in range(4, 8):
        abase, cbase = determinant_particular(n, brow(n, aa, cc))
        trial_a = aa + [abase]
        trial_c = cc + [cbase]
        m = n + 1
        target = sp.Matrix([xc(G[m] - predicted(m, trial_a, trial_c), j) for j in range(m + 1)])
        variables, matrix = tangent_matrix(m)
        solution, residuals = solve_injective(matrix, target, f"row {m}")
        require(all(zero(value) for value in residuals), f"row {m} bracket compatibility")
        theta = sum(solution[j] * x**j for j in range(m))
        aa.append(sp.expand(abase + theta * sp.diff(A0, x)))
        cc.append(sp.expand(cbase + theta * sp.diff(C0, x)))
        require_zero(predicted(m, aa, cc) - G[m], f"row {m} source match")

    # Rebuild the row-eight terminal affine fibre and impose pi_8(P2/P3).
    a8base, c8base = determinant_particular(8, brow(8, aa, cc))
    theta8_vars = list(sp.symbols("theta8_0:9"))
    theta8 = sum(theta8_vars[j] * x**j for j in range(9))
    aa8 = aa + [sp.expand(a8base + theta8 * sp.diff(A0, x))]
    cc8 = cc + [sp.expand(c8base + theta8 * sp.diff(C0, x))]

    acoords8, acols8, amat8 = depth_projection(2, 8)
    ccoords8, ccols8, cmat8 = depth_projection(3, 8)
    require((len(acoords8), len(acols8), amat8.rank()) == (63, 131, 51), "pi8 P2 dimensions")
    require((len(ccoords8), len(ccols8), cmat8.rank()) == (72, 204, 63), "pi8 P3 dimensions")
    residual8 = [sp.expand((left.T * jet(aa8, acoords8))[0]) for left in amat8.T.nullspace()]
    residual8 += [sp.expand((left.T * jet(cc8, ccoords8))[0]) for left in cmat8.T.nullspace()]
    sub8, free8, depth8_matrix, depth8_rhs = parameterize_linear_system(
        residual8, theta8_vars, "row8 depth", 2
    )
    require(len(free8) == 7, "row8 depth affine dimension seven")
    require(all(zero(value.subs(sub8)) for value in residual8), "row8 depth consistency")
    aa8[-1] = sp.expand(aa8[-1].subs(sub8))
    cc8[-1] = sp.expand(cc8[-1].subs(sub8))

    # Independently derive the m=9 Student factor and the rank-seven old-fibre selection.
    student9 = primitive_student(9)
    require(student9 == [12155, 0, 4290, 0, 5148, 0, 11880, 0, 45360, 0], "primitive m9 row")
    diff9 = sp.expand(G[9] - predicted(9, aa8, cc8))
    obstruction9 = sp.factor(sum(student9[j] * xc(diff9, j) for j in range(10)))
    expected_gate9 = -sp.Rational(11, 102987) * F
    # E9_GATE itself has this normalization.  The actual Student scalar is a
    # separate nonzero rational multiple, checked by a constant quotient.
    e9_gate_direct = -255605625 * eta**2 + 6736896000 * XI - 46483785515008
    require_zero(e9_gate_direct - expected_gate9, "restricted E9_GATE factor")
    quotient9 = sp.factor(obstruction9 / F)
    require(quotient9.is_Rational and quotient9 != 0, "Student m9 obstruction is exactly F up to unit")

    equations9 = [xc(diff9, j) for j in range(10)]
    matrix9, rhs9 = sp.linear_eq_to_matrix(equations9, free8)
    require(matrix9.rank() == 7, "G9 old-tangent rank seven")
    solution9, remaining9 = solve_injective(matrix9, rhs9, "G9 old tangent", on_curve=True)
    require(all(curve_zero(value) for value in remaining9), "G9 consistency modulo F")
    selected8 = {variable: solution9[j] for j, variable in enumerate(free8)}
    aa8[-1] = reduce_x(aa8[-1].subs(selected8))
    cc8[-1] = reduce_x(cc8[-1].subs(selected8))
    require(all(curve_zero(xc(predicted(9, aa8, cc8) - G[9], j)) for j in range(10)), "selected row8 gives G9")

    # The rational-square claim is audited on the reduced fraction.
    eta_squared = sp.Rational(415832184456871936, 2393096045625)
    num, den = map(int, (sp.numer(eta_squared), sp.denom(eta_squared)))
    require(gcd(num, den) == 1, "eta squared fraction reduced")
    require(isqrt(num) ** 2 != num and isqrt(den) ** 2 != den, "eta squared is not rational square")
    require(num != 0, "two distinct geometric eta roots")

    # Construct row nine over Q[alpha,eta]/(F), then rebuild the two pi_9 matrices.
    b9 = reduce_x(brow(9, aa8, cc8))
    a9base, c9base = determinant_particular(9, b9, on_curve=True)
    theta9_vars = list(sp.symbols("theta9_0:10"))
    theta9 = sum(theta9_vars[j] * x**j for j in range(10))
    aa9 = aa8 + [sp.expand(a9base + theta9 * sp.diff(A0, x))]
    cc9 = cc8 + [sp.expand(c9base + theta9 * sp.diff(C0, x))]

    acoords9, acols9, amat9 = depth_projection(2, 9)
    ccoords9, ccols9, cmat9 = depth_projection(3, 9)
    require((len(acoords9), len(acols9), amat9.rank()) == (75, 160, 59), "pi9 P2 dimensions")
    require((len(ccoords9), len(ccols9), cmat9.rank()) == (85, 251, 73), "pi9 P3 dimensions")
    anull9 = amat9.T.nullspace()
    cnull9 = cmat9.T.nullspace()
    require((len(anull9), len(cnull9)) == (16, 12), "pi9 left nullities")
    residual9 = [quotient_reduce((left.T * jet(aa9, acoords9))[0]) for left in anull9]
    residual9 += [quotient_reduce((left.T * jet(cc9, ccoords9))[0]) for left in cnull9]
    depth9_matrix, depth9_rhs = sp.linear_eq_to_matrix(residual9, theta9_vars)
    require(depth9_matrix.rank() == 3, "row9 depth tangent rank")
    require(depth9_matrix.rref()[1] == (7, 8, 9), "row9 depth pivots")
    solution_depth9, depth9_remaining = solve_injective(
        depth9_matrix[:, [7, 8, 9]],
        depth9_rhs - depth9_matrix[:, :7] * sp.Matrix(theta9_vars[:7]),
        "row9 depth pivot system",
        on_curve=True,
    )
    depth9_subs = {theta9_vars[7 + j]: solution_depth9[j] for j in range(3)}
    reduced_depth9 = [quotient_reduce(value.subs(depth9_subs)) for value in residual9]
    require(all(zero(value) for value in reduced_depth9), "all 28 depth residuals vanish modulo F")
    require(10 - depth9_matrix.rank() == 7, "row9 terminal affine dimension seven")

    # Student m=10 obstruction.  This is invariant under the row-nine tangent.
    student10 = primitive_student(10)
    require(student10 == [46189, 0, 14586, 0, 15444, 0, 30888, 0, 99792, 0, 489888], "primitive m10 row")
    base_aa9 = aa8 + [a9base]
    base_cc9 = cc8 + [c9base]
    diff10_base = reduce_x(G[10] - predicted(10, base_aa9, base_cc9))
    obstruction10 = quotient_reduce(sum(student10[j] * xc(diff10_base, j) for j in range(11)))
    expected_obstruction10 = sp.Rational(143, 3128230125) * E10_NUM
    require_zero(obstruction10 - expected_obstruction10, "row10 scalar obstruction modulo F")

    # Full G10 selects all ten row-nine tangent coordinates.  Solve ten
    # independent coefficient equations; the eleventh is exactly E10.
    _, matrix10 = tangent_matrix(10)
    require(matrix10.rank() == 10, "G10 full tangent rank")
    rhs10 = sp.Matrix([xc(diff10_base, j) for j in range(11)])
    selected10, residual10 = solve_injective(matrix10, rhs10, "G10 tangent", on_curve=True)
    nonzero_residual10 = [quotient_reduce(value) for value in residual10 if not curve_zero(value)]
    require(len(nonzero_residual10) == 1, "one G10 scalar residual")
    require(
        sp.factor(nonzero_residual10[0] / E10_NUM).is_Rational
        and sp.factor(nonzero_residual10[0] / E10_NUM) != 0,
        "G10 residual is E10 up to unit",
    )

    selected_a9 = reduce_x(a9base + sum(selected10[j] * x**j for j in range(10)) * sp.diff(A0, x))
    selected_rows_a = aa8 + [selected_a9]

    # Rebuild and test the displayed primitive P2 left functional.
    functional_entries = {(5, 0): 35, (6, 2): -20, (7, 4): 10, (8, 6): -4, (9, 8): 1}
    require(gcd(*[abs(value) for value in functional_entries.values()]) == 1, "H_A is primitive")
    hrow = sp.Matrix([[functional_entries.get(coordinate, 0) for coordinate in acoords9]])
    require(hrow * amat9 == sp.zeros(1, amat9.cols), "displayed H_A annihilates all 160 P2 columns")
    hvalue = quotient_reduce((hrow * jet(selected_rows_a, acoords9))[0])
    expected_hvalue = H_NUM / sp.Integer(14189651847000)
    require_zero(hvalue - expected_hvalue, "selected G10 row H_A value modulo F")

    incompatibility = sp.expand(H_NUM - 5 * E10_NUM)
    require(incompatibility == 97758447215837184, "integer incompatibility witness")
    imposed_value = sp.cancel(
        expected_hvalue.subs(alpha * eta, -sp.Rational(619241095293435904, 1010418330375))
    )
    # SymPy substitution on a product is not always structural, so check by
    # the numerator identity as well and compute directly.
    imposed_value = sp.cancel(
        (-5 * 619241095293435904 + 3193963923683016704)
        / sp.Integer(14189651847000)
    )
    require(imposed_value == sp.Rational(9854451712, 1430375), "normalized delayed-depth residual")
    require(imposed_value != 0, "delayed-depth residual nonzero")
    alpha_coefficient = sp.cancel(
        -sp.Rational(619241095293435904, 1010418330375) / eta_squared
    )
    require(
        alpha_coefficient == -sp.Rational(212599558168065, 60278408086247),
        "row10 alpha as a multiple of eta modulo F",
    )

    # Logical geometry: F gives two nonzero conjugate eta values; E10 then
    # fixes alpha uniquely on each line, while H_A excludes both points.
    require(sp.degree(F, eta) == 2 and sp.discriminant(F, eta) != 0, "row9 curve is two geometric lines")
    require(sp.Poly(F, eta, domain=sp.QQ).is_irreducible, "row9 curve has no Q eta root")

    print("THM4358 CLEANROOM REFEREE: PASS")
    print(f"E9_GATE_RESTRICTED={sp.factor(expected_gate9)}")
    print(f"STUDENT9_UNIT={quotient9}")
    print("E9_NORMALIZATION_CAUTION=named_gate_and_primitive_Student_differ_by_nonzero_rational_unit")
    print("ROW8_DEPTH_DIM=7; G9_OLD_TANGENT_RANK=7; IMAGE_DIM=0")
    print("PI9_P2=75x160 rank59 nullity16")
    print("PI9_P3=85x251 rank73 nullity12")
    print("ROW9_DEPTH_RANK=3 pivots=(7,8,9) affine_dim=7; EXTRA_SOURCE=NONE_MOD_F")
    print(f"E10_MOD_F={sp.factor(obstruction10)}")
    print(f"H_A_MOD_F={sp.factor(hvalue)}")
    print(f"INCOMPATIBILITY={incompatibility}")
    print(f"NORMALIZED_RESIDUAL={imposed_value}")
    print(f"CHECKS={checks}")


if __name__ == "__main__":
    main()
