#!/usr/bin/env python3
"""Exact referee for THM-3392.

The theorem compares the synchronized quadratic sign norm

    Q(A) = max_z |z^T A z|

with its bipartite bilinear lift

    B(A) = max_{x,y} |x^T A y|

for symmetric zero-diagonal matrices.  The analytic proof is in the theorem
file.  This companion exhausts every switching class of signed complete
graphs through order six, checks the first strict synchronization loss, and
checks the Sylvester-Hadamard data used in the sharpness construction.

Every truth-bearing check uses ``require`` and remains active under
``python -O``.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations, product


Vector = tuple[int, ...]
Matrix = tuple[tuple[int, ...], ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dot_form(matrix: Matrix, left: Vector, right: Vector) -> int:
    n = len(matrix)
    return sum(
        left[i] * matrix[i][j] * right[j]
        for i in range(n)
        for j in range(n)
    )


def sign_norms(matrix: Matrix) -> tuple[int, int, Vector, Vector, Vector]:
    """Return Q, B and exact maximizers, using B=max_x ||Ax||_1."""

    n = len(matrix)
    best_q = -1
    best_b = -1
    best_z: Vector = ()
    best_x: Vector = ()
    best_y: Vector = ()
    for x in product((-1, 1), repeat=n):
        image = tuple(
            sum(matrix[i][j] * x[j] for j in range(n)) for i in range(n)
        )
        q_value = abs(sum(x[i] * image[i] for i in range(n)))
        b_value = sum(abs(entry) for entry in image)
        if q_value > best_q:
            best_q = q_value
            best_z = x
        if b_value > best_b:
            best_b = b_value
            best_x = x
            best_y = tuple(1 if entry >= 0 else -1 for entry in image)
    require(
        abs(dot_form(matrix, best_x, best_y)) == best_b,
        "bilinear maximizer reconstruction failed",
    )
    return best_q, best_b, best_z, best_x, best_y


def normalized_sign_matrix(n: int, tail: Vector) -> Matrix:
    """Unique first-row-positive representative of a switching class."""

    matrix = [[0] * n for _ in range(n)]
    for j in range(1, n):
        matrix[0][j] = 1
        matrix[j][0] = 1
    tail_edges = tuple(combinations(range(1, n), 2))
    require(len(tail) == len(tail_edges), "wrong switching-tail length")
    for sign, (i, j) in zip(tail, tail_edges, strict=True):
        matrix[i][j] = sign
        matrix[j][i] = sign
    return tuple(tuple(row) for row in matrix)


def sylvester_hadamard(order: int) -> Matrix:
    require(order >= 1 and order & (order - 1) == 0, "order is not a power of 2")
    matrix: Matrix = ((1,),)
    while len(matrix) < order:
        size = len(matrix)
        matrix = tuple(
            tuple(
                matrix[i % size][j % size]
                * (-1 if i >= size and j >= size else 1)
                for j in range(2 * size)
            )
            for i in range(2 * size)
        )
    return matrix


def check_hadamard(matrix: Matrix) -> None:
    order = len(matrix)
    require(all(len(row) == order for row in matrix), "Hadamard is not square")
    require(
        all(entry in (-1, 1) for row in matrix for entry in row),
        "Hadamard entry is not a sign",
    )
    for i in range(order):
        for j in range(order):
            inner = sum(matrix[i][k] * matrix[j][k] for k in range(order))
            require(inner == (order if i == j else 0), "Hadamard Gram identity failed")


ratio_histograms: dict[int, Counter[Fraction]] = {}
switching_classes = 0
sign_cube_evaluations = 0
first_strict_order = None
for n in range(2, 7):
    tail_size = (n - 1) * (n - 2) // 2
    histogram: Counter[Fraction] = Counter()
    for tail in product((-1, 1), repeat=tail_size):
        matrix = normalized_sign_matrix(n, tail)
        q_norm, b_norm, _, _, _ = sign_norms(matrix)
        require(q_norm <= b_norm <= 2 * q_norm, "factor-two comparison failed")
        histogram[Fraction(b_norm, q_norm)] += 1
        switching_classes += 1
        sign_cube_evaluations += 1 << n
        if b_norm > q_norm and first_strict_order is None:
            first_strict_order = n
    ratio_histograms[n] = histogram

require(first_strict_order == 6, "unexpected first strict order")
for n in range(2, 6):
    require(ratio_histograms[n] == Counter({Fraction(1): 1 << ((n - 1) * (n - 2) // 2)}),
            "unexpected pre-six ratio histogram")
require(
    ratio_histograms[6] == Counter({Fraction(1): 1012, Fraction(6, 5): 12}),
    "unexpected order-six ratio histogram",
)

# A canonical first-row-positive strict witness.
witness: Matrix = (
    (0, 1, 1, 1, 1, 1),
    (1, 0, -1, -1, 1, 1),
    (1, -1, 0, 1, -1, 1),
    (1, -1, 1, 0, 1, -1),
    (1, 1, -1, 1, 0, -1),
    (1, 1, 1, -1, -1, 0),
)
q_witness, b_witness, _, x_witness, y_witness = sign_norms(witness)
require((q_witness, b_witness) == (10, 12), "strict witness norms changed")
u_witness = tuple((x + y) // 2 for x, y in zip(x_witness, y_witness, strict=True))
v_witness = tuple((x - y) // 2 for x, y in zip(x_witness, y_witness, strict=True))
require(
    all(u * v == 0 for u, v in zip(u_witness, v_witness, strict=True)),
    "ternary supports are not disjoint",
)
require(
    dot_form(witness, x_witness, y_witness)
    == dot_form(witness, u_witness, u_witness)
    - dot_form(witness, v_witness, v_witness),
    "ternary polarization identity failed",
)

# Assumption hostiles: dropping either zero diagonal or symmetry permits Q=0<B.
diagonal_hostile: Matrix = ((1, 0), (0, -1))
skew_hostile: Matrix = ((0, 1), (-1, 0))
require(sign_norms(diagonal_hostile)[:2] == (0, 2), "diagonal hostile failed")
require(sign_norms(skew_hostile)[:2] == (0, 2), "skew hostile failed")

# Exact algebraic input for the asymptotically sharp sign-matrix family.
# If C=H_m is Sylvester-Hadamard, ||C||_2=sqrt(m), hence
# max_{s,t in {+-1}^m}|s^T C t| <= m^(3/2).  The theorem file embeds C
# between a negative clique and a positive clique.
hadamard_orders = (2, 4, 8, 16, 32, 64)
for order in hadamard_orders:
    check_hadamard(sylvester_hadamard(order))

print("THM-3392 bipartite/synchronized sign referee")
print(f"switching_normalized_classes={switching_classes}")
print(f"sign_cube_evaluations={sign_cube_evaluations}")
for n in range(2, 7):
    encoded = ",".join(
        f"{ratio.numerator}/{ratio.denominator}:{count}"
        for ratio, count in sorted(ratio_histograms[n].items())
    )
    print(f"n={n}_B_over_Q_histogram={encoded}")
print("first_strict_order=6_FINITE_EXACT")
print("strict_witness_Q=10_B=12_ratio=6/5")
print("ternary_disjoint_support_identity=PASS")
print("zero_diagonal_and_symmetry_hostiles_Q=0_B=2_PASS")
print("sylvester_hadamard_orders=" + ",".join(map(str, hadamard_orders)) + "_PASS")
print("sharpness_lower_bound=2m(m-1)/(m^2+2m^(3/2))_TENDS_TO_2")
print("ALL CHECKS PASSED")
