#!/usr/bin/env python3
"""Exact tetrahedral atlas for T4, XOR, H1, and Berggren calibration.

The four vertices, six oriented edges, three perfect matchings, and four
opposite faces of K4 are different S4-sets/representations.  This companion
checks their exact maps, enumerates all 64 labelled tournaments, and records
which information the variable-translation Berggren language discards.
"""

from __future__ import annotations

import hashlib
from itertools import combinations, permutations, product


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


VERTICES = tuple(range(4))
EDGES = tuple(combinations(VERTICES, 2))
EDGE_INDEX = {edge: index for index, edge in enumerate(EDGES)}
EXPECTED_TOURNAMENT_LEDGER_SHA256 = "0e3878081e5cfd5ee63a24e85554605e1f06c5dcb54d389cb8c29634e2a4a2df"


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


def matrix_vector(matrix, vector):
    return tuple(sum(row[index] * vector[index] for index in range(len(vector))) for row in matrix)


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
                numerator = work[row][target] * pivot_value - work[row][column] * work[column][target]
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
        sign = 1 if tail < head else -1
        vector[EDGE_INDEX[edge]] += coefficient * sign
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


def permutation_parity(permutation):
    return sum(
        permutation[i] > permutation[j]
        for i in range(4)
        for j in range(i + 1, 4)
    ) % 2


def cycle_type(permutation):
    seen, lengths = set(), []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current, length = start, 0
        while current not in seen:
            seen.add(current)
            current = permutation[current]
            length += 1
        lengths.append(length)
    return tuple(sorted(lengths))


def edge_action_matrix(permutation):
    matrix = [[0] * 6 for _ in EDGES]
    for source, (tail, head) in enumerate(EDGES):
        moved_tail, moved_head = permutation[tail], permutation[head]
        target = (min(moved_tail, moved_head), max(moved_tail, moved_head))
        sign = 1 if moved_tail < moved_head else -1
        matrix[EDGE_INDEX[target]][source] = sign
    return matrix


MATCHINGS = (
    frozenset((frozenset((0, 1)), frozenset((2, 3)))),
    frozenset((frozenset((0, 2)), frozenset((1, 3)))),
    frozenset((frozenset((0, 3)), frozenset((1, 2)))),
)


def matching_action(permutation):
    lookup = {matching: index for index, matching in enumerate(MATCHINGS)}
    result = []
    for matching in MATCHINGS:
        moved = frozenset(frozenset((permutation[a], permutation[b])) for a, b in matching)
        result.append(lookup[moved])
    return tuple(result)


def affine_permutation(linear_bit, translation):
    def swap_bits(value):
        return ((value & 1) << 1) | ((value & 2) >> 1)

    return tuple((swap_bits(value) if linear_bit else value) ^ translation for value in range(4))


def relabel_tournament(signs, permutation):
    result = [0] * 6
    for source, (tail, head) in enumerate(EDGES):
        moved_tail, moved_head = permutation[tail], permutation[head]
        target = (min(moved_tail, moved_head), max(moved_tail, moved_head))
        sign = 1 if moved_tail < moved_head else -1
        result[EDGE_INDEX[target]] = sign * signs[source]
    return tuple(result)


def outdegrees(signs):
    degrees = []
    for vertex in VERTICES:
        wins = 0
        for other in VERTICES:
            if other == vertex:
                continue
            if vertex < other:
                wins += signs[EDGE_INDEX[(vertex, other)]] == -1
            else:
                wins += signs[EDGE_INDEX[(other, vertex)]] == 1
        degrees.append(wins)
    return tuple(degrees)


B = incidence()
OMEGA = opposite_face_matrix()
BOUNDARY = transpose(B)

# B is the cut map and Omega consists of oriented opposite-face boundaries.
for column in transpose(OMEGA):
    require(matrix_vector(BOUNDARY, column) == (0, 0, 0, 0), "Omega left the cycle space")
require(matrix_vector(OMEGA, (1, 1, 1, 1)) == (0,) * 6, "Omega did not kill constants")
split_matrix = [B[row][:3] + OMEGA[row][:3] for row in range(6)]
require(determinant_integer(split_matrix) == -16, "tetrahedral cut/cycle index changed")

# Sign-twisted S4 covariance and faithfulness on graph H1.
h1_matrices = set()
all_permutations = tuple(permutations(VERTICES))
for permutation in all_permutations:
    action = edge_action_matrix(permutation)
    orientation_sign = -1 if permutation_parity(permutation) else 1
    for vertex in VERTICES:
        moved = matrix_vector(action, tuple(row[vertex] for row in OMEGA))
        expected = tuple(orientation_sign * row[permutation[vertex]] for row in OMEGA)
        require(moved == expected, "opposite-face covariance failed")
    induced_columns = []
    for vertex in range(3):
        target = permutation[vertex]
        coordinates = [0, 0, 0]
        if target < 3:
            coordinates[target] = orientation_sign
        else:
            coordinates = [-orientation_sign] * 3
        induced_columns.append(tuple(coordinates))
    induced = tuple(tuple(induced_columns[column][row] for column in range(3)) for row in range(3))
    expected_trace = (-1 if permutation_parity(permutation) else 1) * (
        sum(permutation[index] == index for index in VERTICES) - 1
    )
    require(sum(induced[index][index] for index in range(3)) == expected_trace, "H1 character changed")
    h1_matrices.add(induced)
require(len(h1_matrices) == 24, "the odd-characteristic H1 action lost S4 information")

# The action on perfect matchings is S4 -> S3 with kernel V4.
matching_images = {matching_action(permutation) for permutation in all_permutations}
matching_kernel = [permutation for permutation in all_permutations if matching_action(permutation) == (0, 1, 2)]
require(len(matching_images) == 6, "matching action is not onto S3")
require(
    sorted(cycle_type(permutation) for permutation in matching_kernel)
    == [(1, 1, 1, 1), (2, 2), (2, 2), (2, 2)],
    "matching kernel is not the normal V4",
)

# Affine lifts with true linear part I have source classes 1 and 2^2;
# those with J have classes 2*1^2 and 4.  These are exactly the two fibres
# used by THM-3497's variable-translation acceptance predicate.
affine_types = {
    bit: {cycle_type(affine_permutation(bit, translation)) for translation in range(4)}
    for bit in (0, 1)
}
require(affine_types[0] == {(1, 1, 1, 1), (2, 2)}, "I-affine classes changed")
require(affine_types[1] == {(1, 1, 2), (4,)}, "J-affine classes changed")

class_order = ((1, 1, 1, 1), (1, 1, 2), (2, 2), (1, 3), (4,))
class_table = []
for kind in class_order:
    representatives = [permutation for permutation in all_permutations if cycle_type(permutation) == kind]
    representative = representatives[0]
    matching_kind = cycle_type(matching_action(representative))
    h1_trace = (-1 if permutation_parity(representative) else 1) * (
        sum(representative[index] == index for index in VERTICES) - 1
    )
    accepted_bits = tuple(bit for bit in (0, 1) if kind in affine_types[bit])
    class_table.append((kind, len(representatives), matching_kind, h1_trace, accepted_bits))

require(
    [(row[0], row[2], row[3], row[4]) for row in class_table]
    == [
        ((1, 1, 1, 1), (1, 1, 1), 3, (0,)),
        ((1, 1, 2), (1, 2), -1, (1,)),
        ((2, 2), (1, 1, 1), -1, (0,)),
        ((1, 3), (3,), 0, ()),
        ((4,), (1, 2), 1, (1,)),
    ],
    "S4/matching/H1/affine class table changed",
)

# THM-3497's projective generators: A and C rotate an opposite face, while B
# is a four-cycle.  The matching quotient sees A,C as 3-cycles and B as a
# reflection.
berggren = {
    "A": (1, 3, 2, 0),
    "B": (1, 3, 0, 2),
    "C": (3, 1, 0, 2),
}
require(cycle_type(berggren["A"]) == (1, 3) and berggren["A"][2] == 2, "A chart changed")
require(cycle_type(berggren["C"]) == (1, 3) and berggren["C"][1] == 1, "C chart changed")
require(cycle_type(berggren["B"]) == (4,), "B chart changed")
require(cycle_type(matching_action(berggren["A"])) == (3,), "q(A) changed")
require(cycle_type(matching_action(berggren["C"])) == (3,), "q(C) changed")
require(cycle_type(matching_action(berggren["B"])) == (1, 2), "q(B) changed")

# Enumerate all labelled tournaments.  Edge sign -1 means the lower-labelled
# endpoint beats the higher-labelled endpoint.  B^T*s is 2*outdegree-3,
# while Omega^T*s records oriented circulation on each opposite triangle.
all_tournaments = set(product((-1, 1), repeat=6))
orbits = []
while all_tournaments:
    representative = min(all_tournaments)
    orbit = {relabel_tournament(representative, permutation) for permutation in all_permutations}
    all_tournaments -= orbit
    orbits.append(tuple(sorted(orbit)))
require(len(orbits) == 4 and sorted(len(orbit) for orbit in orbits) == [8, 8, 24, 24], "T4 orbit census changed")

orbit_rows = []
ledger_rows = []
for orbit in sorted(orbits, key=lambda values: (len(values), values[0])):
    representative = orbit[0]
    score = matrix_vector(BOUNDARY, representative)
    face = matrix_vector(transpose(OMEGA), representative)
    degree_sequence = tuple(sorted(outdegrees(representative), reverse=True))
    cyclic_faces = sum(abs(value) == 3 for value in face)
    require(all(abs(value) in (1, 3) for value in face), "triangle circulation left {1,3}")
    decoded_degrees = tuple(sorted(((value + 3) // 2 for value in score), reverse=True))
    require(decoded_degrees == degree_sequence, "score decoding failed")
    orbit_rows.append((len(orbit), degree_sequence, cyclic_faces, score, face, representative))
    for signs in orbit:
        scores = matrix_vector(BOUNDARY, signs)
        faces = matrix_vector(transpose(OMEGA), signs)
        ledger_rows.append(f"{signs}:{scores}:{faces}")

expected_orbit_profiles = {
    ((3, 1, 1, 1), 1, 8),
    ((2, 2, 2, 0), 1, 8),
    ((3, 2, 1, 0), 0, 24),
    ((2, 2, 1, 1), 2, 24),
}
require(
    {(degrees, cyclic, size) for size, degrees, cyclic, _score, _face, _rep in orbit_rows}
    == expected_orbit_profiles,
    "unlabelled T4 profiles changed",
)

ledger = "\n".join(sorted(ledger_rows))
digest = hashlib.sha256(ledger.encode("ascii")).hexdigest()
require(digest == EXPECTED_TOURNAMENT_LEDGER_SHA256, "tournament atlas ledger changed")

print("== tetrahedral K4 representation atlas ==")
print("sets: vertices=4 edges=6 perfect_matchings=3 opposite_faces=4 tournaments=64")
print("cut_plus_opposite_face_cycle determinant=-16; odd-characteristic H1 action faithful of order 24")
print("matching action: S4 -> S3 onto; kernel cycle types=1^4 plus 3*(2^2)")
print("S4 class table: class,size,matching_class,H1_trace,affine_bits")
for row in class_table:
    print(row)
print("Berggren: A,C are fixed-vertex 3-cycles with matching 3-cycles; B is a 4-cycle with matching reflection")
print("T4 orbit table: size,score_sequence,cyclic_faces,score_vector,face_vector,representative")
for row in orbit_rows:
    print(row)
print(f"labelled tournament score/face ledger sha256={digest}")
print("score plus opposite-face circulation reconstructs every labelled tournament; lattice index=16")
print("scope: finite representation atlas only; no physical current, LRC row exclusion, or JC consequence")
print("all exact checks passed")
