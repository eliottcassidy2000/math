#!/usr/bin/env python3
"""Exact K4 atlas for missing, single, and double directed edges.

Each unordered pair carries two arc bits.  Their centered symmetric sum and
antisymmetric difference split the four states into a cross.  On K4, vertex
degrees plus matching contrasts recover the symmetric channel, while net
scores plus the three primitive Haar cycles recover the skew channel.
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
ARC_STATES = ((0, 0), (1, 0), (0, 1), (1, 1))
STATE_NAMES = {(0, 0): "missing", (1, 0): "forward", (0, 1): "reverse", (1, 1): "double"}
MATCHINGS = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)
HAAR_CYCLES = (
    (1, 0, -1, 1, 0, 1),
    (0, -1, 1, 1, -1, 0),
    (-1, 1, 0, 0, -1, 1),
)
EXPECTED_SEMANTIC_SHA256 = "5e6b7523f0da756849e870e7180b457761430c7b62618fca3c488d41e4a5e08a"


def matrix_vector(matrix, vector):
    return tuple(
        sum(row[index] * vector[index] for index in range(len(vector)))
        for row in matrix
    )


def transpose(matrix):
    return [list(column) for column in zip(*matrix)]


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


def incidence(signed):
    matrix = [[0] * 4 for _ in EDGES]
    for row, (left, right) in enumerate(EDGES):
        matrix[row][left] = -1 if signed else 1
        matrix[row][right] = 1
    return matrix


def relabel(states, permutation):
    result = [None] * 6
    for source, ((left, right), (forward, reverse)) in enumerate(zip(EDGES, states)):
        moved_left, moved_right = permutation[left], permutation[right]
        target = (min(moved_left, moved_right), max(moved_left, moved_right))
        if moved_left < moved_right:
            result[EDGE_INDEX[target]] = (forward, reverse)
        else:
            result[EDGE_INDEX[target]] = (reverse, forward)
    require(all(state is not None for state in result), "relabel lost an edge")
    return tuple(result)


SIGNED_INCIDENCE = incidence(True)
UNSIGNED_INCIDENCE = incidence(False)
SIGNED_BOUNDARY = transpose(SIGNED_INCIDENCE)
UNSIGNED_DEGREES = transpose(UNSIGNED_INCIDENCE)
MATCHING_ROWS = [
    [int(edge in matching) for edge in EDGES]
    for matching in MATCHINGS
]

for cycle in HAAR_CYCLES:
    require(matrix_vector(SIGNED_BOUNDARY, cycle) == (0, 0, 0, 0), "Haar vector is not a cycle")

symmetry_matrix = UNSIGNED_DEGREES + MATCHING_ROWS[:2]
skew_matrix = [row[:] for row in transpose(SIGNED_INCIDENCE)[:3]] + [list(cycle) for cycle in HAAR_CYCLES]
symmetry_determinant = determinant_integer(symmetry_matrix)
skew_determinant = determinant_integer(skew_matrix)
require(abs(symmetry_determinant) == 8, "symmetric degree/matching determinant changed")
require(abs(skew_determinant) == 32, "skew score/Haar determinant changed")
require(abs(symmetry_determinant * skew_determinant) == 256, "combined channel index changed")

records = []
minimal_signatures = set()
symmetric_signatures = {}
skew_signatures = {}
score_only = {}
degree_only = {}
for states in product(ARC_STATES, repeat=6):
    symmetric = tuple(forward + reverse - 1 for forward, reverse in states)
    skew = tuple(forward - reverse for forward, reverse in states)
    require(all(abs(symmetric[i]) + abs(skew[i]) == 1 for i in range(6)), "four-state cross law failed")
    xor_bits = tuple(forward ^ reverse for forward, reverse in states)
    require(xor_bits == tuple(abs(value) for value in skew), "XOR/skew-support identity failed")

    undirected_degrees = matrix_vector(UNSIGNED_DEGREES, symmetric)
    matching_sums = matrix_vector(MATCHING_ROWS, symmetric)
    net_scores = matrix_vector(transpose(SIGNED_INCIDENCE), skew)
    holonomies = tuple(sum(cycle[i] * skew[i] for i in range(6)) for cycle in HAAR_CYCLES)
    minimal = (
        undirected_degrees,
        matching_sums[:2],
        net_scores[:3],
        holonomies,
    )
    minimal_signatures.add(minimal)
    symmetric_signatures.setdefault((undirected_degrees, matching_sums), []).append(states)
    skew_signatures.setdefault((net_scores, holonomies), []).append(states)
    score_only.setdefault(net_scores, []).append(states)
    degree_only.setdefault(undirected_degrees, []).append(states)
    records.append((states, symmetric, skew, undirected_degrees, matching_sums, net_scores, holonomies))

require(len(records) == 4**6 == 4096, "four-state edge census changed")
require(len(minimal_signatures) == 4096, "combined symmetric/skew observables lost a digraph")
require(max(map(len, symmetric_signatures.values())) > 1, "symmetric channel unexpectedly recovered arrows")
require(max(map(len, skew_signatures.values())) > 1, "skew channel unexpectedly distinguished missing from double")
require(max(map(len, score_only.values())) > 1 and max(map(len, degree_only.values())) > 1, "hostile scalar quotient became injective")

# Partial orientations and semicomplete digraphs are two global ternary charts
# meeting on the tournament equator.
partial_orientations = sum(all(state != (1, 1) for state in states) for states, *_ in records)
semicomplete = sum(all(state != (0, 0) for state in states) for states, *_ in records)
tournaments = sum(all(state in ((1, 0), (0, 1)) for state in states) for states, *_ in records)
require((partial_orientations, semicomplete, tournaments) == (3**6, 3**6, 2**6), "ternary-chart census changed")

# Full S4 orbit census of loopless directed graphs on four labelled vertices.
unseen = {record[0] for record in records}
orbit_sizes = []
all_permutations = tuple(permutations(VERTICES))
while unseen:
    representative = min(unseen)
    orbit = {relabel(representative, permutation) for permutation in all_permutations}
    unseen -= orbit
    orbit_sizes.append(len(orbit))
require(len(orbit_sizes) == 218, "unlabelled four-vertex digraph census changed")

semantic_ledger = (
    tuple(sorted(orbit_sizes)),
    tuple(sorted((len(values) for values in symmetric_signatures.values()))),
    tuple(sorted((len(values) for values in skew_signatures.values()))),
    tuple(records),
)
semantic_digest = hashlib.sha256(repr(semantic_ledger).encode("ascii")).hexdigest()
require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "directed-graph semantic ledger changed")

missing_all = ((0, 0),) * 6
double_all = ((1, 1),) * 6
tournament_forward = ((1, 0),) * 6
tournament_reverse = ((0, 1),) * 6
record_by_states = {record[0]: record for record in records}
require(
    record_by_states[missing_all][2] == record_by_states[double_all][2] == (0,) * 6,
    "missing/double skew hostile changed",
)
require(
    record_by_states[tournament_forward][1]
    == record_by_states[tournament_reverse][1]
    == (0,) * 6,
    "converse symmetric hostile changed",
)

print("== K4 missing/single/double directed-edge channel split ==")
print("pair state -> (centered symmetric, skew, XOR):")
for state in ARC_STATES:
    forward, reverse = state
    print(f"  {STATE_NAMES[state]}={state} -> ({forward+reverse-1},{forward-reverse},{forward^reverse})")
print(f"labelled loopless digraphs={len(records)}; S4 orbits={len(orbit_sizes)}")
print(f"partial-orientation/semicomplete/tournament charts={partial_orientations}/{semicomplete}/{tournaments}")
print(f"symmetric degree+two-matching determinant={symmetry_determinant}; index={abs(symmetry_determinant)}")
print(f"skew score+three-Haar determinant={skew_determinant}; index={abs(skew_determinant)}")
print(f"combined continuous channel index={abs(symmetry_determinant*skew_determinant)}; discrete reconstruction=4096/4096")
print(f"hostile skew collision: all-missing vs all-double={missing_all != double_all}; both have skew zero")
print(f"hostile symmetric collision: forward vs converse tournament={tournament_forward != tournament_reverse}; both have symmetric zero")
print(f"semantic ledger sha256={semantic_digest}")
print("scope: complete K4 representation atlas for loopless digraphs; no physical LRC/JC transport")
print("all exact checks passed")
