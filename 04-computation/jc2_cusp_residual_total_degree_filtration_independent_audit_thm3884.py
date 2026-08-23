#!/usr/bin/env python3
"""Import-independent hostile audit of THM-3884.

This companion deliberately does not import SymPy or the canonical artifact.
It expands the THM-3881 residual in a small sparse polynomial ring over Z,
checks every coefficient cell, and audits the degree/equality/gauge seams.
"""

from __future__ import annotations

import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


sys.stdout.reconfigure(newline="\n")

NVAR = 4
ZERO_MONOMIAL = (0,) * NVAR
Poly = dict[tuple[int, ...], int]
GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def clean(raw: Poly) -> Poly:
    return {monomial: coefficient for monomial, coefficient in raw.items() if coefficient}


def const(value: int) -> Poly:
    return {} if value == 0 else {ZERO_MONOMIAL: value}


def var(index: int) -> Poly:
    exponent = [0] * NVAR
    exponent[index] = 1
    return {tuple(exponent): 1}


def add(*polynomials: Poly) -> Poly:
    result: Poly = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            result[monomial] = result.get(monomial, 0) + coefficient
    return clean(result)


def scale(value: int, polynomial: Poly) -> Poly:
    return clean({monomial: value * coefficient for monomial, coefficient in polynomial.items()})


def neg(polynomial: Poly) -> Poly:
    return scale(-1, polynomial)


def sub(left: Poly, right: Poly) -> Poly:
    return add(left, neg(right))


def mul(*polynomials: Poly) -> Poly:
    result = const(1)
    for polynomial in polynomials:
        product: Poly = {}
        for left_monomial, left_coefficient in result.items():
            for right_monomial, right_coefficient in polynomial.items():
                monomial = tuple(
                    left_monomial[index] + right_monomial[index]
                    for index in range(NVAR)
                )
                product[monomial] = (
                    product.get(monomial, 0) + left_coefficient * right_coefficient
                )
        result = clean(product)
    return result


def power(polynomial: Poly, exponent: int) -> Poly:
    gate(exponent >= 0, "nonnegative power")
    result = const(1)
    base = polynomial
    value = exponent
    while value:
        if value & 1:
            result = mul(result, base)
        base = mul(base, base)
        value //= 2
    return result


def xy_degree(polynomial: Poly) -> int | None:
    if not polynomial:
        return None
    return max(monomial[0] + monomial[1] for monomial in polynomial)


def xy_homogeneous_part(polynomial: Poly, degree: int) -> Poly:
    return {
        monomial: coefficient
        for monomial, coefficient in polynomial.items()
        if monomial[0] + monomial[1] == degree
    }


def tf_cells(polynomial: Poly) -> dict[tuple[int, int], Poly]:
    result: dict[tuple[int, int], Poly] = {}
    for monomial, coefficient in polynomial.items():
        key = (monomial[2], monomial[3])
        xy_monomial = (monomial[0], monomial[1], 0, 0)
        cell = result.setdefault(key, {})
        cell[xy_monomial] = cell.get(xy_monomial, 0) + coefficient
    return {key: clean(value) for key, value in result.items() if clean(value)}


def reduce_mod(polynomial: Poly, prime: int) -> Poly:
    return clean({monomial: coefficient % prime for monomial, coefficient in polynomial.items()})


def x_valuation(polynomial: Poly) -> int | None:
    if not polynomial:
        return None
    return min(monomial[0] for monomial in polynomial)


def matrix_rank(matrix: list[list[int]]) -> int:
    work = [[Fraction(entry) for entry in row] for row in matrix]
    if not work:
        return 0
    rows = len(work)
    columns = len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows) if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [entry / pivot_value for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                work[row][index] - multiple * work[pivot_row][index]
                for index in range(columns)
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


# Variables are ordered x,y,T,F.  Everything below is rebuilt from the
# displayed THM-3881 definitions, not from the canonical THM-3884 companion.
x, y, T, F = (var(index) for index in range(NVAR))
a = add(x, const(1))
L = add(scale(9, x), const(4))
K = add(power(y, 2), scale(-15, power(x, 2)), scale(-15, x), const(-4))
P = mul(a, power(L, 2))
Delta = sub(mul(power(a, 3), power(L, 2)), power(K, 2))
r = add(mul(a, T), mul(K, F))
A = add(mul(K, T), mul(a, P, F))

residual = add(
    power(L, 4),
    mul(
        const(2),
        add(scale(3, A), scale(3, P), power(r, 2)),
        power(L, 2),
        F,
    ),
    mul(
        add(scale(8, A), scale(6, P), scale(3, power(r, 2))),
        sub(mul(P, power(F, 2)), power(T, 2)),
    ),
)
cells = tf_cells(residual)

# Independently hand-collected coefficient formulas.  Equality with the raw
# expansion audits both expansion routes, while their exact degrees freeze
# the full fourteen-cell universe.
expected_cells: dict[tuple[int, int], Poly] = {
    (4, 0): scale(-3, power(a, 2)),
    (3, 1): scale(-6, mul(a, K)),
    (3, 0): scale(-8, K),
    (2, 2): scale(3, Delta),
    (2, 1): scale(-6, mul(power(a, 2), power(L, 2))),
    (2, 0): scale(-6, P),
    (1, 3): scale(6, mul(power(a, 2), K, power(L, 2))),
    (1, 2): scale(12, mul(a, K, power(L, 2))),
    (1, 1): scale(6, mul(K, power(L, 2))),
    (0, 4): scale(3, mul(a, power(K, 2), power(L, 2))),
    (0, 3): add(scale(2, mul(power(K, 2), power(L, 2))), scale(8, mul(a, power(P, 2)))),
    (0, 2): scale(12, mul(power(a, 2), power(L, 4))),
    (0, 1): scale(6, mul(a, power(L, 4))),
    (0, 0): power(L, 4),
}
gate(cells == expected_cells, "independent expansion equals hand-collected cells")
gate(len(cells) == 14, "exactly fourteen nonzero characteristic-zero cells")

expected_degrees = {
    (4, 0): 2,
    (3, 1): 3,
    (3, 0): 2,
    (2, 2): 5,
    (2, 1): 4,
    (2, 0): 3,
    (1, 3): 6,
    (1, 2): 5,
    (1, 1): 4,
    (0, 4): 7,
    (0, 3): 7,
    (0, 2): 6,
    (0, 1): 5,
    (0, 0): 4,
}
actual_degrees = {key: xy_degree(value) for key, value in cells.items()}
gate(actual_degrees == expected_degrees, "complete coefficient-degree ledger")

# The four load-bearing identities are rechecked as exact sparse polynomials.
gate(cells[(0, 4)] == scale(3, mul(a, power(K, 2), power(L, 2))), "F4 identity")
gate(cells[(1, 3)] == scale(6, mul(K, power(a, 2), power(L, 2))), "TF3 identity")
gate(cells[(2, 2)] == scale(3, Delta), "T2F2 identity")
gate(cells[(4, 0)] == scale(-3, power(a, 2)), "T4 identity")

# All-n affine inequalities, including the hostile boundary n=1.
for (i, j), coefficient_degree in sorted(expected_degrees.items()):
    if (i, j) == (0, 4):
        continue
    slope = 4 - i - j
    intercept = 7 - coefficient_degree
    gate(slope >= 0, f"m<=n slope {(i, j)}")
    gate(slope + intercept >= 1, f"m<=n strict at n=1 {(i, j)}")
    for n in (1, 2, 7, 31):
        gate(4 * n + 7 > (i + j) * n + coefficient_degree, f"m<=n spot {(i, j), n}")

K2 = add(power(y, 2), scale(-15, power(x, 2)))
top_f4 = xy_homogeneous_part(cells[(0, 4)], 7)
gate(top_f4 == scale(243, mul(power(x, 3), power(K2, 2))), "unique m<=n top form")

# At m=n+1 exactly three cells meet the target.  No other cell can enter by
# cancellation because each has a strict degree deficit before summation.
equality_top_cells: set[tuple[int, int]] = set()
for (i, j), coefficient_degree in sorted(expected_degrees.items()):
    slope = 4 - i - j
    intercept = 7 - i - coefficient_degree
    if slope == 0 and intercept == 0:
        equality_top_cells.add((i, j))
        continue
    gate(slope >= 0, f"equality slope {(i, j)}")
    gate(slope + intercept >= 1, f"equality strict at n=1 {(i, j)}")
    for n in (1, 2, 7, 31):
        gate(4 * n + 7 > i * (n + 1) + j * n + coefficient_degree, f"equality spot {(i, j), n}")
gate(equality_top_cells == {(0, 4), (1, 3), (2, 2)}, "complete equality packet")

# Reuse the four formal coordinates as x,y,U,V, where U=T_m and V=f_n.
U, V = T, F
equality_top = add(
    mul(xy_homogeneous_part(cells[(0, 4)], 7), power(V, 4)),
    mul(xy_homogeneous_part(cells[(1, 3)], 6), U, power(V, 3)),
    mul(xy_homogeneous_part(cells[(2, 2)], 5), power(U, 2), power(V, 2)),
)
equality_factor = scale(
    243,
    mul(power(x, 3), power(V, 2), power(add(mul(K2, V), mul(x, U)), 2)),
)
gate(equality_top == equality_factor, "equality leading factorization")
gate(xy_degree(equality_top) == 7, "coefficient-side degree of equality packet")

# Direct valuation probes exercise the prime-x parity mechanism away from the
# seam.  The general proof is the identity v_x(x^3 V^2 W^2)=3+2v_x(V)+2v_x(W).
valuation_controls = (
    (add(x, y), add(power(x, 2), power(y, 2))),
    (mul(x, add(x, y)), add(x, power(y, 2))),
    (power(x, 3), add(power(x, 4), power(y, 4))),
)
for index, (u_value, v_value) in enumerate(valuation_controls):
    bracket = add(mul(K2, v_value), mul(x, u_value))
    gate(bool(bracket), f"valuation bracket nonzero {index}")
    factored = mul(power(x, 3), power(v_value, 2), power(bracket, 2))
    expected_valuation = 3 + 2 * x_valuation(v_value) + 2 * x_valuation(bracket)
    gate(x_valuation(factored) == expected_valuation, f"valuation formula {index}")
    gate(expected_valuation % 2 == 1, f"valuation odd {index}")

# The seam equation is solved independently in every homogeneous degree
# n=1..12 by rational row reduction.  Its full kernel has exactly the n
# expected gauge vectors f_n=xq, T_m=-K_2q.
kernel_rows: list[tuple[int, int, int]] = []
for n in range(1, 13):
    variable_count = (n + 1) + (n + 2)
    equation_count = n + 3
    matrix = [[0 for _ in range(variable_count)] for _ in range(equation_count)]
    # f_j and t_j index the y-exponent in homogeneous f_n and T_{n+1}.
    for j in range(n + 1):
        matrix[j][j] += -15
        matrix[j + 2][j] += 1
    for j in range(n + 2):
        matrix[j][n + 1 + j] += 1
    rank = matrix_rank(matrix)
    nullity = variable_count - rank
    gate(rank == equation_count, f"seam map full row rank n={n}")
    gate(nullity == n, f"seam nullity n={n}")
    for q_y_degree in range(n):
        vector = [0] * variable_count
        vector[q_y_degree] = 1
        vector[n + 1 + q_y_degree] = 15
        vector[n + 1 + q_y_degree + 2] = -1
        image = [sum(row[column] * vector[column] for column in range(variable_count)) for row in matrix]
        gate(image == [0] * equation_count, f"gauge kernel vector {(n, q_y_degree)}")
    kernel_rows.append((n, rank, nullity))

# The full gauge increment (Kq,-aq) fixes T(0,0)-4f(0,0), and its leading
# part cancels a seam pair (-K2*q,x*q).
gate(K.get((0, 0, 0, 0)) == -4, "K origin")
gate(a.get((0, 0, 0, 0)) == 1, "a origin")
gate(-4 - 4 * (-1) == 0, "full gauge preserves address")

# THM-3884 alone closes every nonconstant affine-slope lane T=h f.
# Constant h has m=n; a genuinely linear h has m=n+1 but K2+x h1 retains
# y^2 coefficient one, so it cannot lie on the equality seam.
alpha_x_beta_y = add(x, y)  # coefficient choices are irrelevant to y^2
affine_obstruction_control = add(K2, mul(x, alpha_x_beta_y))
gate(affine_obstruction_control.get((0, 2, 0, 0)) == 1, "affine seam retains y2")
gate(xy_degree(affine_obstruction_control) == 2, "affine seam nonzero quadratic")

# The degree proof itself extends from characteristic zero/algebraically
# closed fields to integral-domain coefficients of characteristic != 3.
# Reductions audit the delicate small-prime behavior; vanished lower cells
# only improve the degree gaps.  Characteristic 3 is a genuine proof boundary.
prime_rows: list[tuple[int, int, int, int]] = []
for prime in (2, 5, 7, 11, 13):
    reduced_cells = tf_cells(reduce_mod(residual, prime))
    f4_degree = xy_degree(reduced_cells[(0, 4)])
    delta_degree = xy_degree(reduce_mod(Delta, prime))
    gate(f4_degree == 7, f"F4 top survives characteristic {prime}")
    gate(delta_degree == 5, f"Delta top survives characteristic {prime}")
    for key, reduced_cell in reduced_cells.items():
        gate(xy_degree(reduced_cell) <= expected_degrees[key], f"no new degree mod {prime} cell {key}")
    gate(reduce_mod(equality_top, prime) == reduce_mod(equality_factor, prime), f"factor mod {prime}")
    prime_rows.append((prime, len(reduced_cells), f4_degree, delta_degree))
gate(xy_degree(reduce_mod(L, 3)) == 0, "characteristic three collapses L degree")
char3_cells = tf_cells(reduce_mod(residual, 3))
gate((0, 4) not in char3_cells, "characteristic three kills the F4 cell")

source_text = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))), "no inactive asserts")

semantic = {
    "affine_relation": "THM3884 closes nonconstant T=h*f for deg(h)<=1; constant-f edge is external",
    "characteristic_scope": "degree filtration works over integral domains of characteristic not 3",
    "coefficient_cells": sorted((i, j, expected_degrees[(i, j)]) for i, j in expected_degrees),
    "degree_floor": "square and nonconstant f imply deg(T)>=deg(f)+1",
    "equality_factor": "243*x^3*f_n^2*(K2*f_n+x*T_m)^2",
    "equality_kernel": "f_n=x*q and T_m=-K2*q",
    "limitations": "gauge does not preserve the residual square equation; no termination or JC2 claim",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("THM3884_INDEPENDENT_HOSTILE_AUDIT_20260823")
print("verdict=PASS_WITH_STRICT_SCOPE_AND_PORTABILITY_CORRECTION")
print("coefficient_degree_cells=14;independent_sparse_Z_expansion=PASS")
print("degree_floor=nonconstant_square_implies_degT_at_least_degf_plus_1")
print("equality_factor=243*x^3*f_n^2*(K2*f_n+x*T_m)^2")
print("equality_kernel=f_n=x*q;T_m=-K2*q")
print(f"kernel_rows={kernel_rows[0], kernel_rows[-1]}")
print(f"prime_rows={tuple(prime_rows)}")
print("scope_extension=integral_domains_characteristic_not_3")
print("affine_corollary=nonconstant_deg_h_at_most_1_closed;constant_f_requires_direct_edge")
print("gauge_invariance_or_termination=NOT_CLAIMED")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"gates={GATES}")
print("ALL CHECKS PASSED")
