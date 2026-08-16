#!/usr/bin/env python3
"""Exact bridge from 2x2 mixed Haar data to tetrahedral K4 cycles.

The four cells of F_2^2 are the vertices of K4.  Each nontrivial Walsh
character selects a perfect matching by equal signs; its opposite-face image
is twice the oriented complementary four-cycle.  This identifies the local
transportation-kernel direction with one marked H1 coordinate without
claiming any physical LRC realization.
"""

from __future__ import annotations

import hashlib
from itertools import combinations, permutations, product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


VERTICES = tuple(range(4))
BITS = ((0, 0), (0, 1), (1, 0), (1, 1))
EDGES = tuple(combinations(VERTICES, 2))
EDGE_INDEX = {edge: index for index, edge in enumerate(EDGES)}
EXPECTED_TOURNAMENT_LEDGER_SHA256 = "8da9354e69971ab3e27b54564567ce3f431721c77412480b667f655d93c55e85"


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def matrix_vector(matrix, vector):
    return tuple(
        sum(row[index] * vector[index] for index in range(len(vector)))
        for row in matrix
    )


def determinant_integer(matrix):
    work = [row[:] for row in matrix]
    size = len(work)
    require(all(len(row) == size for row in work), "determinant matrix is not square")
    sign, previous = 1, 1
    for column in range(size - 1):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        pivot_value = work[column][column]
        for row in range(column + 1, size):
            for target in range(column + 1, size):
                numerator = (
                    work[row][target] * pivot_value
                    - work[row][column] * work[column][target]
                )
                require(numerator % previous == 0, "Bareiss division ceased to be integral")
                work[row][target] = numerator // previous
            work[row][column] = 0
        previous = pivot_value
    return sign * work[-1][-1]


def rank_mod(matrix, prime):
    work = [[value % prime for value in row] for row in matrix]
    row = 0
    for column in range(len(work[0]) if work else 0):
        pivot = next((index for index in range(row, len(work)) if work[index][column]), None)
        if pivot is None:
            continue
        work[row], work[pivot] = work[pivot], work[row]
        inverse = pow(work[row][column], prime - 2, prime)
        work[row] = [(value * inverse) % prime for value in work[row]]
        for index in range(len(work)):
            if index == row or work[index][column] == 0:
                continue
            scale = work[index][column]
            work[index] = [
                (work[index][target] - scale * work[row][target]) % prime
                for target in range(len(work[index]))
            ]
        row += 1
    return row


def incidence():
    matrix = [[0] * 4 for _ in EDGES]
    for row, (tail, head) in enumerate(EDGES):
        matrix[row][tail] = -1
        matrix[row][head] = 1
    return matrix


def triangle_boundary(face):
    a, b, c = face
    vector = [0] * 6
    for tail, head, coefficient in ((b, c, 1), (a, c, -1), (a, b, 1)):
        edge = (min(tail, head), max(tail, head))
        vector[EDGE_INDEX[edge]] += coefficient * (1 if tail < head else -1)
    return vector


def opposite_face_matrix():
    matrix = [[0] * 4 for _ in EDGES]
    for omitted in VERTICES:
        face = tuple(vertex for vertex in VERTICES if vertex != omitted)
        boundary = triangle_boundary(face)
        sign = -1 if omitted % 2 else 1
        for edge, coefficient in enumerate(boundary):
            matrix[edge][omitted] = sign * coefficient
    return matrix


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def canonical_line(vector):
    negative = tuple(-value for value in vector)
    return min(tuple(vector), negative)


B = incidence()
BOUNDARY = transpose(B)
OMEGA = opposite_face_matrix()

# Omega is a scaled isometry on the sum-zero vertex space.
gram = [
    [dot(tuple(row[i] for row in OMEGA), tuple(row[j] for row in OMEGA)) for j in VERTICES]
    for i in VERTICES
]
require(
    gram == [[3 if i == j else -1 for j in VERTICES] for i in VERTICES],
    "opposite-face Gram identity changed",
)

characters = {}
for frequency in BITS[1:]:
    character = tuple(
        -1 if (frequency[0] * point[0] + frequency[1] * point[1]) % 2 else 1
        for point in BITS
    )
    image = matrix_vector(OMEGA, character)
    require(all(value % 2 == 0 for value in image), "Walsh image ceased to be even")
    cycle = tuple(value // 2 for value in image)
    require(matrix_vector(BOUNDARY, cycle) == (0, 0, 0, 0), "Walsh image is not a cycle")
    matching = tuple(
        edge for edge in EDGES if character[edge[0]] == character[edge[1]]
    )
    support = tuple(edge for edge, value in zip(EDGES, cycle) if value)
    require(len(matching) == 2 and len(support) == 4, "matching/cycle support changed")
    require(set(matching).isdisjoint(support) and set(matching) | set(support) == set(EDGES), "support is not the complementary square")
    characters[frequency] = (character, matching, cycle)

character_vectors = [characters[frequency][0] for frequency in BITS[1:]]
cycle_vectors = [characters[frequency][2] for frequency in BITS[1:]]
require(all(dot(left, right) == 0 for left, right in combinations(character_vectors, 2)), "Walsh characters lost orthogonality")
require(all(dot(left, right) == 0 for left, right in combinations(cycle_vectors, 2)), "square cycles lost orthogonality")
require(all(dot(vector, vector) == 4 for vector in character_vectors), "Walsh norm changed")
require(all(dot(vector, vector) == 4 for vector in cycle_vectors), "cycle norm changed")

# The ordinary 2x2 row/column marginal has the XOR character as its exact
# one-dimensional kernel over every odd characteristic.
marginal = [
    [int(point[0] == bit) for point in BITS] for bit in (0, 1)
] + [
    [int(point[1] == bit) for point in BITS] for bit in (0, 1)
]
xor_character = characters[(1, 1)][0]
require(rank_mod(marginal, 13) == 3, "2x2 marginal rank changed")
require(matrix_vector(marginal, xor_character) == (0, 0, 0, 0), "XOR left the transportation kernel")

# The three Walsh lines are exactly the three balanced 2+2 cuts.  S4 acts on
# them through S3, with the Klein four group as kernel.
walsh_lines = {canonical_line(vector) for vector in character_vectors}
line_actions = set()
line_kernel = []
for permutation in permutations(VERTICES):
    action = []
    for vector in character_vectors:
        moved = [0] * 4
        for source, target in enumerate(permutation):
            moved[target] = vector[source]
        action.append(character_vectors.index(next(v for v in character_vectors if canonical_line(v) == canonical_line(moved))))
    image = tuple(action)
    line_actions.add(image)
    if image == (0, 1, 2):
        line_kernel.append(permutation)
require(len(walsh_lines) == 3 and len(line_actions) == 6, "S4 did not induce S3 on Haar lines")
require(len(line_kernel) == 4, "Haar-line kernel is not V4")

# Scores (three independent cuts) plus the three primitive Haar cycles recover
# all six edge coordinates over Q; the integral lattice index is 32.
split_matrix = [B[row][:3] + [cycle[row] for cycle in cycle_vectors] for row in range(6)]
split_determinant = determinant_integer(split_matrix)
require(abs(split_determinant) == 32, "cut/Haar-cycle lattice index changed")

tournament_ledger = []
for signs in product((-1, 1), repeat=6):
    score = matrix_vector(transpose(B), signs)
    holonomy = tuple(dot(cycle, signs) for cycle in cycle_vectors)
    tournament_ledger.append((signs, score, holonomy))
require(len({(row[1], row[2]) for row in tournament_ledger}) == 64, "score plus Haar holonomy lost a tournament")
tournament_digest = hashlib.sha256(repr(tournament_ledger).encode("ascii")).hexdigest()
require(
    tournament_digest == EXPECTED_TOURNAMENT_LEDGER_SHA256,
    "tournament score/Haar ledger changed",
)

print("== tetrahedral K4 Haar/XOR cycle bridge ==")
print(f"vertices={BITS}; edges={EDGES}")
print("Omega^T Omega=4I-J; Walsh and primitive square-cycle bases are orthogonal")
for frequency in BITS[1:]:
    character, matching, cycle = characters[frequency]
    print(f"frequency={frequency}; character={character}; equal-sign matching={matching}; complementary cycle={cycle}")
print("2x2 row/column marginal rank over F13=3; kernel=<character_(1,1)>")
print(f"S4 -> S3 on three Haar/matching lines; image={len(line_actions)}; kernel={len(line_kernel)}=V4")
print(f"cut plus primitive-Haar-cycle determinant={split_determinant}; lattice index={abs(split_determinant)}")
print(f"score plus three Haar holonomies reconstructs tournaments=64/64")
print(f"tournament score/Haar ledger sha256={tournament_digest}")
print("scope: exact representation bridge only; no actual U_full address map, physical current, row exclusion, or LRC(14)")
print("all exact checks passed")
