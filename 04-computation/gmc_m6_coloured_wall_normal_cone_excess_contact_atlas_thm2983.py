#!/usr/bin/env python3
"""Exact companion for THM-2983's coloured M=6 wall/contact atlas.

The computation uses the fixed ORIGINAL_F width-six Macaulay chart imported
from the immutable THM-2969 companion.  It verifies:

* the full-chart and pure-resultant top-n quadratic sections;
* exact QQ/fmpz wall ranks, with three modular ranks as cross-checks;
* the half-wall and six integer-wall first normal maps;
* exact block-Toeplitz liftable-kernel barcodes on the stated arcs; and
* the root-six projective first-death arrangement and its non-SNC boundary.

All mathematical guards use ``require`` and remain active under ``python -O``.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import sys
from fractions import Fraction as F
from pathlib import Path

import sympy as sp
from flint import fmpq_mat, fmpz_mat


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.py"
BASE_LF_SHA256 = "5be5df0fd058f436593339f78cd6240a47975082d22929c135ba811458753bf5"
WIDTH = 6
N = 36
COLOURS = {2: 0, 3: 1, 4: 2}
PRIMES = (1_000_000_007, 1_000_000_009, 1_000_000_033)

# Root-six linear-factor coefficients.
B6 = 2_517_055_682_762_701_080
C6 = 871_966_204_362_527_579
A6_COEFF = 7_748_852_908_937_866_554


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n")


def load_base():
    actual = hashlib.sha256(lf_bytes(BASE_PATH)).hexdigest()
    require(actual == BASE_LF_SHA256, ("THM-2969 dependency hash", actual))
    spec = importlib.util.spec_from_file_location("thm2969_for_thm2983", BASE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-2969")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def coefficient(poly, degree: int) -> int:
    if hasattr(poly, "__getitem__"):
        return int(poly[degree])
    return int(poly) if degree == 0 else 0


def compact(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def as_fraction(value) -> F:
    if hasattr(value, "p") and hasattr(value, "q"):
        return F(int(value.p), int(value.q))
    return F(int(value.numerator), int(value.denominator))


def trace(matrix) -> F:
    return sum((as_fraction(matrix[i, i]) for i in range(matrix.nrows())), F(0))


def matmul_trace(left, right) -> F:
    return trace(left * right)


def homogeneous_values(poly, degree: int, n: int, colour: int) -> tuple[int, int, int]:
    value = dn = dc = 0
    for exponent in range(degree + 1):
        a = coefficient(poly, exponent)
        value += a * n**exponent * colour ** (degree - exponent)
        if exponent:
            dn += exponent * a * n ** (exponent - 1) * colour ** (degree - exponent)
        if exponent < degree:
            dc += (degree - exponent) * a * n**exponent * colour ** (degree - exponent - 1)
    return value, dn, dc


def modular_rank(matrix: list[list[int]], prime: int) -> int:
    rows = [[value % prime for value in row] for row in matrix]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inverse = pow(rows[rank][column], -1, prime)
        rows[rank] = [(value * inverse) % prime for value in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][column]:
                scale = rows[i][column]
                rows[i] = [
                    (rows[i][j] - scale * rows[rank][j]) % prime
                    for j in range(len(rows[0]))
                ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def first_two_colour_jets(poly, colour: int) -> dict[tuple[int, int, int], int]:
    degree = poly.degree()
    answer = {}
    for deficit in range(3):
        monomial = [0, 0, 0]
        monomial[colour] = deficit
        answer[tuple(monomial)] = int(poly[degree - deficit])
    return answer


def add_jet(left, right, scale: int = 1):
    answer = dict(left)
    for monomial, value in right.items():
        answer[monomial] = answer.get(monomial, 0) + scale * value
    return {monomial: value for monomial, value in answer.items() if value}


def mul_jet(left, right):
    answer = {}
    for a, ca in left.items():
        for b, cb in right.items():
            monomial = tuple(a[i] + b[i] for i in range(3))
            if sum(monomial) <= 2:
                answer[monomial] = answer.get(monomial, 0) + ca * cb
    return {monomial: value for monomial, value in answer.items() if value}


def power_jet(jet, exponent: int):
    answer = {(0, 0, 0): 1}
    for _ in range(exponent):
        answer = mul_jet(answer, jet)
    return answer


def linear_square_coefficient(linear: list[F], monomial: tuple[int, int, int]) -> F:
    support = [i for i, exponent in enumerate(monomial) for _ in range(exponent)]
    require(len(support) == 2, monomial)
    if support[0] == support[1]:
        return linear[support[0]] ** 2
    return 2 * linear[support[0]] * linear[support[1]]


def normalized_log_quadratic(linear, ordinary):
    return {
        monomial: ordinary[monomial] - F(1, 2) * linear_square_coefficient(linear, monomial)
        for monomial in ordinary
    }


def restore_ordinary_quadratic(linear, logarithmic):
    return {
        monomial: logarithmic[monomial] + F(1, 2) * linear_square_coefficient(linear, monomial)
        for monomial in logarithmic
    }


def quadratic_hessian(degree: int, linear: list[F], quadratic):
    rows = [[F(0) for _ in range(4)] for _ in range(4)]
    rows[0][0] = F(degree * (degree - 1))
    for colour in range(3):
        rows[0][colour + 1] = rows[colour + 1][0] = F(degree - 1) * linear[colour]
        mono = tuple(2 if i == colour else 0 for i in range(3))
        rows[colour + 1][colour + 1] = 2 * quadratic[mono]
    for left in range(3):
        for right in range(left + 1, 3):
            mono = tuple(int(i in (left, right)) for i in range(3))
            rows[left + 1][right + 1] = rows[right + 1][left + 1] = quadratic[mono]
    return rows


def sympy_matrix(rows):
    return sp.Matrix([[sp.Rational(value.numerator, value.denominator) for value in row] for row in rows])


def exact_inertia(rows) -> tuple[int, int]:
    matrix = sympy_matrix(rows)
    minors = [sp.Rational(1)]
    for size in range(1, matrix.rows + 1):
        minors.append(matrix[:size, :size].det())
    require(all(value != 0 for value in minors), minors)
    pivots = [minors[i] / minors[i - 1] for i in range(1, len(minors))]
    return sum(bool(value > 0) for value in pivots), sum(bool(value < 0) for value in pivots)


def topjet_checks(base, t, forms, rows, metadata, selected):
    a0_rows: list[list[int]] = []
    a1_rows = [[[] for _ in range(N)] for _ in range(3)]
    a2_rows = [[[] for _ in range(N)] for _ in range(3)]
    for local_row, global_row in enumerate(selected):
        order = metadata[global_row][0]
        colour = COLOURS[order]
        degree = (order - 1) * WIDTH - 1
        row0, row1, row2 = [], [], []
        for poly in rows[global_row]:
            row0.append(coefficient(poly, degree))
            row1.append(coefficient(poly, degree - 1))
            row2.append(coefficient(poly, degree - 2))
        a0_rows.append(row0)
        for which in range(3):
            a1_rows[which][local_row] = row1 if which == colour else [0] * N
            a2_rows[which][local_row] = row2 if which == colour else [0] * N

    leading = fmpz_mat(a0_rows)
    det0 = int(leading.det())
    require(det0 > 0, "leading chart determinant is not positive")
    inverse = fmpq_mat(a0_rows).inv()
    x = [inverse * fmpq_mat(matrix) for matrix in a1_rows]
    y = [inverse * fmpq_mat(matrix) for matrix in a2_rows]
    linear_p = [trace(matrix) for matrix in x]
    quadratic_p = {}
    for colour in range(3):
        mono = tuple(2 if i == colour else 0 for i in range(3))
        quadratic_p[mono] = trace(y[colour]) + F(1, 2) * (
            linear_p[colour] ** 2 - matmul_trace(x[colour], x[colour])
        )
    for left in range(3):
        for right in range(left + 1, 3):
            mono = tuple(int(i in (left, right)) for i in range(3))
            quadratic_p[mono] = linear_p[left] * linear_p[right] - matmul_trace(x[left], x[right])

    quadratic, cubic, _quartic = forms
    q200, q110, q020 = (quadratic[key] for key in ((2, 0, 0), (1, 1, 0), (0, 2, 0)))
    c300, c210, c120 = (cubic[key] for key in ((3, 0, 0), (2, 1, 0), (1, 2, 0)))
    q200j, q110j, q020j = (first_two_colour_jets(poly, 0) for poly in (q200, q110, q020))
    c300j, c210j, c120j = (first_two_colour_jets(poly, 1) for poly in (c300, c210, c120))
    curvature = {}
    curvature = add_jet(curvature, mul_jet(c120j, power_jet(q200j, 2)))
    curvature = add_jet(curvature, mul_jet(mul_jet(c210j, q110j), q200j), -1)
    curvature = add_jet(curvature, mul_jet(mul_jet(c300j, q020j), q200j), -1)
    curvature = add_jet(curvature, mul_jet(c300j, power_jet(q110j, 2)))
    flag = mul_jet(mul_jet(power_jet(q200j, 6), c300j), curvature)
    flag0 = flag[(0, 0, 0)]
    require(flag0 > 0, flag0)
    units = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    quadratics = ((2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1))
    linear_f = [F(flag.get(monomial, 0), flag0) for monomial in units]
    quadratic_f = {monomial: F(flag.get(monomial, 0), flag0) for monomial in quadratics}
    log_p = normalized_log_quadratic(linear_p, quadratic_p)
    log_f = normalized_log_quadratic(linear_f, quadratic_f)
    linear_r = [linear_p[i] - linear_f[i] for i in range(3)]
    log_r = {monomial: log_p[monomial] - log_f[monomial] for monomial in quadratics}
    quadratic_r = restore_ordinary_quadratic(linear_r, log_r)

    degree_p = 58 * WIDTH - 36
    degree_r = 46 * WIDTH - 26
    original, _mutated = t.interpolate_pair(rows, degree_p)
    qfull, cfull, kfull, _alternate = t.flag_polynomials(forms)
    resultant = original // (qfull**6 * cfull * kfull)
    require(original.degree() == degree_p and resultant.degree() == degree_r, "degree")
    require(int(original[degree_p]) == det0, "leading determinant mismatch")
    require(sum(linear_p, F(0)) == F(int(original[degree_p - 1]), det0), "chart linear")
    require(sum(quadratic_p.values(), F(0)) == F(int(original[degree_p - 2]), det0), "chart quadratic")
    r0 = int(resultant[degree_r])
    require(sum(linear_r, F(0)) == F(int(resultant[degree_r - 1]), r0), "resultant linear")
    require(sum(quadratic_r.values(), F(0)) == F(int(resultant[degree_r - 2]), r0), "resultant quadratic")
    hp = quadratic_hessian(degree_p, linear_p, quadratic_p)
    hr = quadratic_hessian(degree_r, linear_r, quadratic_r)
    inertia_p = exact_inertia(hp)
    inertia_r = exact_inertia(hr)
    require(inertia_p == inertia_r == (1, 3), (inertia_p, inertia_r))
    colour_r = sympy_matrix(hr)[1:, 1:]
    transverse_basis = sp.Matrix([[1, 1], [-1, 0], [0, -1]])
    transverse = transverse_basis.T * colour_r * transverse_basis
    expected_trace = sp.Rational(244051230785617, 5627343750)
    expected_det = sp.Rational(-10882890806368883653776584, 376988067626953125)
    require(transverse.trace() == expected_trace and transverse.det() == expected_det, "transverse")
    return degree_p, degree_r, expected_trace, expected_det


def point_matrices(rows, metadata, selected, n0: int, colour0: int):
    matrices = [[[0] * N for _ in range(N)] for _ in range(5)]
    for local_row, global_row in enumerate(selected):
        order = metadata[global_row][0]
        colour_index = COLOURS[order]
        degree = (order - 1) * WIDTH - 1
        for column, poly in enumerate(rows[global_row]):
            value, dn, dc = homogeneous_values(poly, degree, n0, colour0)
            matrices[0][local_row][column] = value
            matrices[1][local_row][column] = dn
            matrices[colour_index + 2][local_row][column] = dc
    return matrices


def null_basis(matrix: fmpz_mat) -> tuple[fmpz_mat, int]:
    full, nullity = matrix.nullspace()
    basis = fmpz_mat([[full[i, j] for j in range(nullity)] for i in range(full.nrows())])
    return basis, nullity


def normal_pieces(matrices):
    a0 = fmpz_mat(matrices[0])
    right, nr = null_basis(a0)
    left, nl = null_basis(a0.transpose())
    require(nr == nl, (nr, nl))
    pieces = [left.transpose() * fmpz_mat(matrix) * right for matrix in matrices[1:]]
    return a0, pieces, nr


def common_kernel_dimensions(pieces) -> tuple[int, int]:
    r = pieces[0].nrows()
    stacked_right = fmpz_mat(
        [[pieces[p][i, j] for j in range(r)] for p in range(4) for i in range(r)]
    )
    stacked_left = fmpz_mat(
        [[pieces[p][j, i] for j in range(r)] for p in range(4) for i in range(r)]
    )
    _right, kr = null_basis(stacked_right)
    _left, kl = null_basis(stacked_left)
    return kr, kl


def symbolic_normal_determinant(pieces):
    z, q, c, f = sp.symbols("z q c f")
    variables = (z, q, c, f)
    size = pieces[0].nrows()
    matrix = sp.zeros(size)
    for variable, piece in zip(variables, pieces):
        matrix += variable * sp.Matrix([[int(piece[i, j]) for j in range(size)] for i in range(size)])
    return sp.factor(matrix.det(method="domain-ge")), variables


def verify_first_normal_maps(base, t, forms, rows, metadata, selected):
    q0, c0, curvature, _alternate = t.flag_polynomials(forms)
    flag = q0**6 * c0 * curvature
    w6 = (base.smith_factor(WIDTH) * base.seam_factor(WIDTH)) // base.expected_flag(WIDTH)

    half_matrices = point_matrices(rows, metadata, selected, -1, 2)
    half_a0, half_pieces, half_nullity = normal_pieces(half_matrices)
    require(int(half_a0.rank()) == 30 and half_nullity == 6, "half rank")
    require(all(modular_rank(half_matrices[0], prime) == 30 for prime in PRIMES), "half modular")
    half_det, variables = symbolic_normal_determinant(half_pieces)
    z, q, c, f = variables
    half_ratio = sp.cancel(half_det / (f + 2 * z) ** 6)
    require(half_ratio.is_number and half_ratio != 0, "half first normal")

    expected_nullities = (2, 10, 15, 10, 9, 15)
    expected_flag_orders = (0, 1, 1, 0, 0, 1)
    expected_w_orders = (6, 24, 26, 24, 23, 20)
    expected_kernels = ((0, 0), (5, 2), (8, 10), (9, 8), (8, 8), (0, 0))
    rank_records = []
    for root in range(1, 7):
        matrices = point_matrices(rows, metadata, selected, -root, 1)
        a0, pieces, nullity = normal_pieces(matrices)
        exact_rank = int(a0.rank())
        mod_ranks = tuple(modular_rank(matrices[0], prime) for prime in PRIMES)
        require(nullity == expected_nullities[root - 1], (root, nullity))
        require(exact_rank == N - nullity and mod_ranks == (exact_rank,) * 3, (root, exact_rank, mod_ranks))
        factor = base.X + root
        flag_order = base.valuation(flag, factor)
        w_order = base.valuation(w6, factor)
        require(flag_order == expected_flag_orders[root - 1], (root, flag_order))
        require(w_order == expected_w_orders[root - 1], (root, w_order))
        common = common_kernel_dimensions(pieces)
        require(common == expected_kernels[root - 1], (root, common))
        if root == 1:
            determinant, variables = symbolic_normal_determinant(pieces)
            z, q, c, f = variables
            ratio = sp.cancel(determinant / (q + f - 2 * c) ** 2)
            require(ratio.is_number and ratio != 0, "root1 first normal")
        elif root == 6:
            determinant, variables = symbolic_normal_determinant(pieces)
            z, q, c, f = variables
            plane = z + 6 * c
            linear = A6_COEFF * c - B6 * q + C6 * z
            quadric = -6 * c * f - 6 * c * q - 2 * c * z + 12 * f * q + f * z + q * z
            ratio = sp.cancel(determinant / (plane**2 * linear * quadric**6))
            require(ratio.is_number and ratio != 0, "root6 first normal")
            require(A6_COEFF == B6 + 6 * C6, "root6 hidden incidence")
            require(sp.expand(linear - (C6 * plane + B6 * (c - q))) == 0, "linear rewrite")
            require(
                sp.expand(quadric - (plane * (f - 2 * c + q) + 12 * (c - q) * (c - f))) == 0,
                "quadric rewrite",
            )
            gradient = [sp.diff(quadric, variable) for variable in variables]
            singular = sp.linsolve(gradient, variables)
            require(singular == sp.FiniteSet((-6 * f, f, f, f)), singular)
        rank_records.append((root, exact_rank, nullity, flag_order, w_order, common))
    return rank_records


def line_coefficient(poly, degree: int, n0: int, c0: int, zn: int, zc: int, order: int) -> int:
    answer = 0
    for exponent in range(degree + 1):
        value = coefficient(poly, exponent)
        if not value:
            continue
        for left_order in range(max(0, order - (degree - exponent)), min(order, exponent) + 1):
            right_order = order - left_order
            answer += (
                value
                * math.comb(exponent, left_order)
                * n0 ** (exponent - left_order)
                * zn**left_order
                * math.comb(degree - exponent, right_order)
                * c0 ** (degree - exponent - right_order)
                * zc**right_order
            )
    return answer


def coefficient_matrices(rows, metadata, selected, root: int, direction: tuple[int, int, int, int]):
    zn, zq, zc, zf = direction
    colour_directions = (zq, zc, zf)
    maximum_degree = max((metadata[index][0] - 1) * WIDTH - 1 for index in selected)
    matrices = [[[0 for _ in range(N)] for _ in range(N)] for _ in range(maximum_degree + 1)]
    for local_row, global_row in enumerate(selected):
        degree = (metadata[global_row][0] - 1) * WIDTH - 1
        colour = COLOURS[metadata[global_row][0]]
        for column, poly in enumerate(rows[global_row]):
            for order in range(degree + 1):
                matrices[order][local_row][column] = line_coefficient(
                    poly, degree, -root, 1, zn, colour_directions[colour], order
                )
    return matrices


def block_toeplitz(matrices, depth: int) -> fmpz_mat:
    size = N * depth
    block = [[0 for _ in range(size)] for _ in range(size)]
    zero = [[0] * N for _ in range(N)]
    for equation in range(depth):
        for unknown in range(equation + 1):
            order = equation - unknown
            matrix = matrices[order] if order < len(matrices) else zero
            ro = equation * N
            co = unknown * N
            for row in range(N):
                block[ro + row][co : co + N] = matrix[row]
    return fmpz_mat(block)


def barcode(matrices, maximum_depth: int = 30):
    previous_rank = 0
    lifetimes = []
    for depth in range(1, maximum_depth + 1):
        rank = int(block_toeplitz(matrices, depth).rank())
        liftable = N - rank + previous_rank
        require(liftable >= 0, (depth, rank, previous_rank, liftable))
        lifetimes.append(liftable)
        previous_rank = rank
        if liftable == 0:
            return tuple(lifetimes), True
    return tuple(lifetimes), False


def verify_arc_barcodes(base, t, forms, rows, metadata, selected):
    q0, c0, curvature, _alternate = t.flag_polynomials(forms)
    flag = q0**6 * c0 * curvature
    w6 = (base.smith_factor(WIDTH) * base.seam_factor(WIDTH)) // base.expected_flag(WIDTH)
    expected_n = {
        1: (2, 2, 1, 1, 0), 2: (10, 9, 4, 2, 0), 3: (15, 10, 2, 0),
        4: (10, 9, 4, 1, 0), 5: (9, 8, 5, 1, 0), 6: (15, 6, 0),
    }
    expected_generic = {
        1: (2, 0), 2: (10, 5, 2, 0), 3: (15, 10, 2, 0),
        4: (10, 9, 4, 1, 0), 5: (9, 8, 5, 1, 0), 6: (15, 0),
    }
    records = []
    directions = (
        ("n", (1, 0, 0, 0)),
        ("1235", (1, 2, 3, 5)),
        ("2357", (2, 3, 5, 7)),
        ("signed", (1, -1, 2, -3)),
    )
    for root in range(1, 7):
        diagonal_order = base.valuation(flag, base.X + root) + base.valuation(w6, base.X + root)
        for name, direction in directions:
            lifetimes, terminated = barcode(coefficient_matrices(rows, metadata, selected, root, direction))
            require(terminated, (root, name, lifetimes))
            if name == "n":
                require(lifetimes == expected_n[root] and sum(lifetimes) == diagonal_order, (root, lifetimes))
            elif name == "2357" and root == 1:
                require(lifetimes == (2, 1, 0), lifetimes)
            else:
                require(lifetimes == expected_generic[root], (root, name, lifetimes))
            records.append((root, name, lifetimes, sum(lifetimes)))
    return records


def factor_values(direction: tuple[int, int, int, int]) -> tuple[int, int, int]:
    z, q, c, f = direction
    plane = z + 6 * c
    linear = A6_COEFF * c - B6 * q + C6 * z
    quadric = -6 * c * f - 6 * c * q - 2 * c * z + 12 * f * q + f * z + q * z
    return plane, linear, quadric


def homogeneous_form_at(form, n_value, colour_value, variables):
    answer = sp.Integer(0)
    for monomial, poly in form.items():
        degree = poly.degree()
        scalar = sum(
            coefficient(poly, exponent)
            * n_value**exponent
            * colour_value ** (degree - exponent)
            for exponent in range(degree + 1)
        )
        answer += scalar * sp.prod(variables[i] ** monomial[i] for i in range(3))
    return sp.expand(answer)


def verify_null_arc_factorizations(forms):
    x0, x1, x2, s = sp.symbols("x0 x1 x2 s")
    variables = (x0, x1, x2)
    ell = x0 + x1 + x2
    h = (
        63299 * x0**2
        + 121066 * x0 * x1
        + 114169 * x0 * x2
        + 57772 * x1**2
        + 108669 * x1 * x2
        + 50919 * x2**2
    )
    q0 = homogeneous_form_at(forms[0], -6, 1, variables)
    c0 = homogeneous_form_at(forms[1], -6, 1, variables)
    f0 = homogeneous_form_at(forms[2], -6, 1, variables)
    require(sp.expand(q0 + 24 * h) == 0, "root6 quadratic factor")
    require(sp.expand(c0 - 641571840 * ell * h) == 0, "root6 cubic factor")
    require(
        sp.expand(f0 + 93261243646771200 * ell**2 * h) == 0,
        "root6 quartic factor",
    )

    # On the A cap Q, off-L line direction (-6,0,1,1), put s=1+t.
    # The cubic and quartic colours scale with n, while the quadratic colour
    # stays one.  Hence C_s and F_s retain the fixed conic H for every s.
    c_arc = homogeneous_form_at(forms[1], -6 * s, s, variables)
    f_arc = homogeneous_form_at(forms[2], -6 * s, s, variables)
    require(sp.expand(c_arc - 641571840 * s**11 * ell * h) == 0, "A_Q cubic arc")
    require(
        sp.expand(f_arc + 93261243646771200 * s**17 * ell**2 * h) == 0,
        "A_Q quartic arc",
    )
    return (63299, 121066, 114169, 57772, 108669, 50919)


def verify_death_arrangement(rows, metadata, selected):
    d = B6 + 12 * C6
    directions = {
        "generic": (1, 2, 3, 5),
        "A": (-6, 2, 1, 3),
        "L": (B6, C6, 0, 1),
        "Q": (1, 1, 1, 1),
        "A_Q": (-6, 0, 1, 1),
        "L_Q": (B6 * d, C6 * d, 0, -B6 * C6),
        "triple_generic": (-6, 1, 1, 2),
        "triple_singular": (-6, 1, 1, 1),
    }
    expected_zero = {
        "generic": (), "A": ("A",), "L": ("L",), "Q": ("Q",),
        "A_Q": ("A", "Q"), "L_Q": ("L", "Q"),
        "triple_generic": ("A", "L", "Q"),
        "triple_singular": ("A", "L", "Q"),
    }
    expected = {
        "generic": ((15, 0), True),
        "A": ((15, 2, 0), True),
        "L": ((15, 1, 0), True),
        "Q": ((15, 6, 0), True),
        "A_Q": ((15, 6, 6, 6, 6, 6), False),
        "L_Q": ((15, 6, 1, 0), True),
        "triple_generic": ((15, 9, 9, 9, 9, 9), False),
        "triple_singular": ((15, 15, 15, 15, 15, 15), False),
    }
    records = []
    labels = ("A", "L", "Q")
    for name, direction in directions.items():
        zero_set = tuple(label for label, value in zip(labels, factor_values(direction)) if value == 0)
        require(zero_set == expected_zero[name], (name, zero_set))
        maximum = 6 if not expected[name][1] else 30
        lifetimes, terminated = barcode(
            coefficient_matrices(rows, metadata, selected, 6, direction), maximum
        )
        require((lifetimes, terminated) == expected[name], (name, lifetimes, terminated))
        records.append((name, zero_set, lifetimes, terminated))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    base = load_base()
    t = base.thm2943
    forms = t.polynomial_forms(WIDTH, (0, 1, 2))
    rows, metadata = t.thm2942.macaulay_rows(forms)
    selected = t.selection(t.ORIGINAL_F)
    require(len(selected) == N, len(selected))

    degree_p, degree_r, transverse_trace, transverse_det = topjet_checks(
        base, t, forms, rows, metadata, selected
    )
    ranks = verify_first_normal_maps(base, t, forms, rows, metadata, selected)
    arcs = verify_arc_barcodes(base, t, forms, rows, metadata, selected)
    strata = verify_death_arrangement(rows, metadata, selected)
    h_coefficients = verify_null_arc_factorizations(forms)

    lines = [
        "M6 COLOURED WALL NORMAL-CONE AND EXCESS-CONTACT ATLAS",
        f"dependency_thm2969_lf_sha256={BASE_LF_SHA256}",
        f"chart=ORIGINAL_F;degree={degree_p};resultant_degree={degree_r};topjet_inertia=1+,3-",
        f"transverse_trace={transverse_trace};transverse_det={transverse_det}",
        "halfwall=rank30;nullity6;first_normal=(f+2z)^6",
    ]
    for root, rank, nullity, flag_order, w_order, common in ranks:
        lines.append(
            f"root={root};exact_rank={rank};nullity={nullity};modular_ranks={rank},{rank},{rank};"
            f"flag_order={flag_order};W6_order={w_order};common_kernel=R{common[0]}L{common[1]}"
        )
    for root, name, lifetimes, contact in arcs:
        lines.append(f"arc_root={root};direction={name};liftable={lifetimes};contact={contact}")
    lines.append("root6_factor_geometry=A=z+6c;L=C*A+B*(c-q);Q=A*(f-2c+q)+12*(c-q)*(c-f)")
    lines.append("root6_incidence=A_cap_L_subset_Q;quadric_singular=(-6,1,1,1)")
    lines.append(
        "root6_AQ_null_conic_coefficients=" + ",".join(str(value) for value in h_coefficients)
    )
    lines.append("root6_AQ_null_arc=C_s=s^11*ell*H;F_s=s^17*ell^2*H;resultant=IDENTICALLY_ZERO")
    for name, zero_set, lifetimes, terminated in strata:
        if name in {"A_Q", "triple_generic", "triple_singular"}:
            contact = "NULL_ARC_PREFIX"
        else:
            contact = str(sum(lifetimes)) if terminated else "OPEN_PREFIX"
        lines.append(
            f"root6_stratum={name};zero_set={'+'.join(zero_set) if zero_set else 'none'};"
            f"liftable={lifetimes};contact={contact}"
        )
    digest = hashlib.sha256("\n".join(lines).encode()).hexdigest()
    lines.extend(
        [
            f"record_digest={digest}",
            "exact_qq_fmpz_ranks=PRIMARY;three_prime_ranks=CROSS_CHECK",
            "all_exact_checks=PASS",
        ]
    )
    transcript = "\n".join(lines) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
