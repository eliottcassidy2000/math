#!/usr/bin/env python3
"""Exact finite referee for THM-3072.

The regular A4 flag torsor has three conjugate C2 quotients.  This companion
checks their edge/oriented-cycle geometry and the exact characteristic-zero
projector identity which reconstructs a flag signal from all three quotient
tables while leaving a three-dimensional blind sector for any two.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import permutations


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


V = tuple(range(4))
DIRECTIONS = (1, 2, 3)
RHO = (0, 2, 3, 1)


def rho_power(x, exponent):
    for _ in range(exponent % 3):
        x = RHO[x]
    return x


FLAGS = tuple((x, direction) for x in V for direction in DIRECTIONS)
FLAG_INDEX = {flag: index for index, flag in enumerate(FLAGS)}


def tau(index, flag):
    x, direction = flag
    return (x ^ rho_power(direction, index), direction)


def rotate(flag):
    x, direction = flag
    return (x, RHO[direction])


def quotient_class(index, flag):
    return tuple(sorted((flag, tau(index, flag))))


QUOTIENTS = tuple(tuple(sorted({quotient_class(j, flag) for flag in FLAGS}))
                  for j in range(3))
require(tuple(len(q) for q in QUOTIENTS) == (6, 6, 6), "quotient census")

# The three nonzero translations form V4 and are cycled by the ternary move.
for flag in FLAGS:
    require(tau(0, tau(0, flag)) == flag, ("tau square", 0, flag))
    require(tau(1, tau(1, flag)) == flag, ("tau square", 1, flag))
    require(tau(2, tau(2, flag)) == flag, ("tau square", 2, flag))
    require(tau(0, tau(1, flag)) == tau(2, flag), ("V4 product", flag))
    require(rotate(tau(0, flag)) == tau(2, rotate(flag)),
            ("C3 conjugacy", flag))


def edge_of(flag):
    x, direction = flag
    return tuple(sorted((x, x ^ direction)))


def canonical_oriented_cycle(sequence):
    sequence = tuple(sequence)
    rotations = tuple(sequence[j:] + sequence[:j] for j in range(4))
    return min(rotations)


def oriented_cycle_of(flag):
    x, direction = flag
    next_direction = RHO[direction]
    return canonical_oriented_cycle((
        x,
        x ^ direction,
        x ^ direction ^ next_direction,
        x ^ next_direction,
    ))


EDGES = tuple(sorted({edge_of(flag) for flag in FLAGS}))
ALL_ORIENTED_CYCLES = tuple(sorted({
    canonical_oriented_cycle(sequence) for sequence in permutations(V)
}))
FLAG_CYCLES = tuple(sorted({oriented_cycle_of(flag) for flag in FLAGS}))
require(len(EDGES) == 6, ("edge count", EDGES))
require(len(ALL_ORIENTED_CYCLES) == 6, ("oriented cycle count", ALL_ORIENTED_CYCLES))
require(FLAG_CYCLES == ALL_ORIENTED_CYCLES, "flag cycles are incomplete")

# X0 is the ordinary edge set.  X2 is the oriented Hamiltonian-cycle set.
for first in FLAGS:
    for second in FLAGS:
        require((quotient_class(0, first) == quotient_class(0, second)) ==
                (edge_of(first) == edge_of(second)),
                ("edge deck", first, second))
        require((quotient_class(2, first) == quotient_class(2, second)) ==
                (oriented_cycle_of(first) == oriented_cycle_of(second)),
                ("cycle deck", first, second))

# X1 is the third conjugate edge chart.
def skew_edge_of(flag):
    x, direction = flag
    return tuple(sorted((x, x ^ RHO[direction])))


require({skew_edge_of(flag) for flag in FLAGS} == set(EDGES), "skew edge census")
for first in FLAGS:
    for second in FLAGS:
        require((quotient_class(1, first) == quotient_class(1, second)) ==
                (skew_edge_of(first) == skew_edge_of(second)),
                ("skew deck", first, second))

# Left A4 equivariance.
AFFINE = tuple((translation, power) for translation in V for power in range(3))


def left_affine(element, flag):
    translation, power = element
    x, direction = flag
    return (translation ^ rho_power(x, power), rho_power(direction, power))


def act_on_edge(element, edge):
    translation, power = element
    return tuple(sorted(translation ^ rho_power(x, power) for x in edge))


def act_on_cycle(element, cycle):
    translation, power = element
    return canonical_oriented_cycle(tuple(
        translation ^ rho_power(x, power) for x in cycle
    ))


for element in AFFINE:
    for flag in FLAGS:
        require(edge_of(left_affine(element, flag)) ==
                act_on_edge(element, edge_of(flag)), ("edge equivariance", element, flag))
        require(oriented_cycle_of(left_affine(element, flag)) ==
                act_on_cycle(element, oriented_cycle_of(flag)),
                ("cycle equivariance", element, flag))

# The joint edge/cycle label recovers a flag, but its image is only three
# disjoint K2,2 components, one for each perfect matching direction.
JOINT = {(edge_of(flag), oriented_cycle_of(flag)) for flag in FLAGS}
require(len(JOINT) == 12, ("joint map is not injective", len(JOINT)))
component_records = []
for direction in DIRECTIONS:
    matching_edges = {edge_of((x, direction)) for x in V}
    direction_cycles = {oriented_cycle_of((x, direction)) for x in V}
    pairs = {(edge, cycle) for edge, cycle in JOINT
             if edge in matching_edges and cycle in direction_cycles}
    require(len(matching_edges) == 2 and len(direction_cycles) == 2,
            ("component sides", direction, matching_edges, direction_cycles))
    require(pairs == {(edge, cycle) for edge in matching_edges
                      for cycle in direction_cycles},
            ("component is not K2,2", direction, pairs))
    component_records.append((direction, tuple(sorted(matching_edges)),
                              tuple(sorted(direction_cycles))))
require(sum(4 for _ in component_records) == len(JOINT), "component total")

# Exact rational matrix helpers.
N = len(FLAGS)


def zero_matrix(rows, columns):
    return [[Fraction(0) for _ in range(columns)] for _ in range(rows)]


def identity_matrix(size):
    result = zero_matrix(size, size)
    for j in range(size):
        result[j][j] = Fraction(1)
    return result


def permutation_matrix(action):
    result = zero_matrix(N, N)
    for row, flag in enumerate(FLAGS):
        result[row][FLAG_INDEX[action(flag)]] = Fraction(1)
    return result


def matrix_add(*matrices):
    return [[sum(matrix[row][column] for matrix in matrices)
             for column in range(len(matrices[0][0]))]
            for row in range(len(matrices[0]))]


def matrix_scale(scalar, matrix):
    scalar = Fraction(scalar)
    return [[scalar * value for value in row] for row in matrix]


def matrix_multiply(left, right):
    return [[sum(left[row][middle] * right[middle][column]
                 for middle in range(len(right)))
             for column in range(len(right[0]))]
            for row in range(len(left))]


def matrix_rank(matrix):
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column] != 0), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        pivot_value = work[pivot_row][column]
        work[pivot_row] = [value / pivot_value for value in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or work[row][column] == 0:
                continue
            factor = work[row][column]
            work[row] = [work[row][j] - factor * work[pivot_row][j]
                         for j in range(columns)]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def horizontal(left, right):
    return [left[row] + right[row] for row in range(len(left))]


def vertical(*matrices):
    return [row[:] for matrix in matrices for row in matrix]


IDENTITY = identity_matrix(N)
TAU = tuple(permutation_matrix(lambda flag, j=j: tau(j, flag)) for j in range(3))
PROJECTORS = tuple(matrix_scale(Fraction(1, 2), matrix_add(IDENTITY, TAU[j]))
                   for j in range(3))
P_V4 = matrix_scale(Fraction(1, 4), matrix_add(IDENTITY, *TAU))

for projector in PROJECTORS:
    require(matrix_multiply(projector, projector) == projector, "Pj not idempotent")
require(matrix_multiply(P_V4, P_V4) == P_V4, "PV4 not idempotent")
for i in range(3):
    for j in range(i + 1, 3):
        require(matrix_multiply(PROJECTORS[i], PROJECTORS[j]) == P_V4,
                ("pair intersection projector", i, j))

projector_sum = matrix_add(*PROJECTORS)
reconstruction_right = matrix_add(IDENTITY, matrix_scale(2, P_V4))
require(projector_sum == reconstruction_right, "projector reconstruction identity")

projector_ranks = tuple(matrix_rank(projector) for projector in PROJECTORS)
pair_sum_ranks = tuple(matrix_rank(horizontal(PROJECTORS[i], PROJECTORS[j]))
                       for i in range(3) for j in range(i + 1, 3))
pair_observation_ranks = tuple(matrix_rank(vertical(PROJECTORS[i], PROJECTORS[j]))
                               for i in range(3) for j in range(i + 1, 3))
triple_observation_rank = matrix_rank(vertical(*PROJECTORS))
require(projector_ranks == (6, 6, 6), ("projector ranks", projector_ranks))
require(matrix_rank(P_V4) == 3, ("common rank", matrix_rank(P_V4)))
require(pair_sum_ranks == (9, 9, 9), ("pair sum ranks", pair_sum_ranks))
require(pair_observation_ranks == (9, 9, 9),
        ("pair observation ranks", pair_observation_ranks))
require(triple_observation_rank == 12,
        ("triple observation rank", triple_observation_rank))


def apply(matrix, vector):
    return tuple(sum(matrix[row][column] * vector[column]
                     for column in range(len(vector)))
                 for row in range(len(matrix)))


# For the pair of ordinary edge and oriented-cycle tables, X1 is the exact
# three-dimensional blind sector.  Its basis is supported one direction at a
# time and uses the V4 character whose kernel is <tau_1>.
blind_vectors = []
for supported_direction in DIRECTIONS:
    values = []
    for x, direction in FLAGS:
        if direction != supported_direction:
            values.append(Fraction(0))
        else:
            kernel_vector = RHO[direction]
            values.append(Fraction(1 if x in (0, kernel_vector) else -1))
    vector = tuple(values)
    require(apply(PROJECTORS[0], vector) == (Fraction(0),) * N,
            ("edge hostile not killed", supported_direction))
    require(apply(PROJECTORS[2], vector) == (Fraction(0),) * N,
            ("cycle hostile not killed", supported_direction))
    require(apply(PROJECTORS[1], vector) == vector,
            ("missing chart does not retain hostile", supported_direction))
    blind_vectors.append(vector)
require(matrix_rank([list(vector) for vector in blind_vectors]) == 3,
        "blind hostile basis rank")

# The odd reflection fixes the ordinary edge quotient and swaps the other two.
REFLECTION = (0, 2, 1, 3)


def reflected(flag):
    x, direction = flag
    return (REFLECTION[x], REFLECTION[direction])


for flag in FLAGS:
    require(reflected(tau(0, reflected(flag))) == tau(0, flag),
            ("reflection H0", flag))
    require(reflected(tau(1, reflected(flag))) == tau(2, flag),
            ("reflection H1/H2", flag))

semantic = repr((
    FLAGS,
    QUOTIENTS,
    EDGES,
    ALL_ORIENTED_CYCLES,
    sorted(JOINT),
    component_records,
    projector_ranks,
    pair_sum_ranks,
    pair_observation_ranks,
    triple_observation_rank,
    blind_vectors,
)).encode("ascii")
semantic_sha256 = sha256(semantic).hexdigest()

print("THM3072_A4_FLAG_THREE_C2_TOMOGRAPHY")
print("flags=12;C2_quotients=3;quotient_sizes=6,6,6;right_kernel=V4")
print("geometry=X0:K4_edges,X1:skew_edges,X2:oriented_Hamiltonian_cycles")
print("edge_cycle_joint_map=injective;image=3_disjoint_K2,2;size=12")
print("projector_identity=P0+P1+P2=I+2PV4;PiPj=PV4")
print("ranks=Pi:6,PV4:3,pair_sum:9,pair_observation:9,triple_observation:12")
print("edge_plus_cycle_blind_dimension=3;blind_sector=X1_nontrivial_character")
print("C3_cycles_the_three_quotients=1;odd_reflection_fixes_X0_and_swaps_X1_X2=1")
print(f"semantic_sha256={semantic_sha256}")
print("scope=finite_A4_permutation_module;common_atom_not_physical_quartic_LRC_or_JC_carrier")
