#!/usr/bin/env python3
"""Assertion-independent exact audit of the THM-3280 lattice/type claim.

This deliberately does not import the theorem companion or SymPy.  Smith
claims are certified by Bareiss determinants and determinantal divisors.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


N = 12
SELECTOR_SUPPORT = (1, 2, 4, 5, 10)
GEODESIC_SUPPORT = (0, 1, 4, 8, 11)
ORDINARY_HOP_DISTANCES = (1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1)

# The canonical nonlinear rank ordering from THM-3269 and the genuine J_12
# phase chart from THM-3277, both normalized by 17 -> 0 and 16 -> 1.
ELL_ORDER = (17, 16, 18, 13, 11, 22, 19, 10, 7, 21, 2, 3)
ELL = {vertex: residue for residue, vertex in enumerate(ELL_ORDER)}
J12 = {
    2: 10, 3: 10, 7: 7, 10: 10, 11: 4, 13: 4,
    16: 1, 17: 0, 18: 1, 19: 4, 21: 8, 22: 1,
}


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def cyclic_shift(vector, amount):
    return tuple(vector[(index - amount) % N] for index in range(N))


def cyclic_difference(vector):
    shifted = cyclic_shift(vector, 1)
    return tuple(shifted[index] - vector[index] for index in range(N))


def matrix_from_columns(columns):
    return [[columns[column][row] for column in range(len(columns))]
            for row in range(len(columns[0]))]


def augmentation_coordinates(vector):
    require(sum(vector) == 0, ("nonzero augmentation", vector))
    return tuple(vector[index] for index in range(1, N))


def augmentation_multiplication(vector):
    columns = []
    for power in range(1, N):
        shifted = cyclic_shift(vector, power)
        columns.append(augmentation_coordinates(tuple(
            shifted[index] - vector[index] for index in range(N)
        )))
    return matrix_from_columns(tuple(columns))


def difference_orbit(vector):
    difference = cyclic_difference(vector)
    return matrix_from_columns(tuple(
        augmentation_coordinates(cyclic_shift(difference, power))
        for power in range(N)
    ))


def bareiss_determinant(matrix):
    require(len(matrix) == len(matrix[0]), "square determinant")
    size = len(matrix)
    if size == 0:
        return 1
    work = [list(map(int, row)) for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        pivot_row = next((row for row in range(pivot_index, size)
                          if work[row][pivot_index] != 0), None)
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (work[row][column] * pivot
                             - work[row][pivot_index] * work[pivot_index][column])
                require(numerator % previous == 0,
                        ("nonintegral Bareiss step", numerator, previous))
                work[row][column] = numerator // previous
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def minor(matrix, rows, columns):
    return [[matrix[row][column] for column in columns] for row in rows]


def determinantal_divisor(matrix, size):
    value = 0
    row_sets = tuple(combinations(range(len(matrix)), size))
    column_sets = tuple(combinations(range(len(matrix[0])), size))
    for rows in row_sets:
        for columns in column_sets:
            value = gcd(value, abs(bareiss_determinant(minor(matrix, rows, columns))))
            if value == 1:
                return 1
    return value


def first_minor_with_abs(matrix, size, target):
    for rows in combinations(range(len(matrix)), size):
        for columns in combinations(range(len(matrix[0])), size):
            determinant = bareiss_determinant(minor(matrix, rows, columns))
            if abs(determinant) == target:
                return rows, columns, determinant
    return None


def rank_over_q(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    columns = len(work[0])
    rank = 0
    for column in range(columns):
        pivot = next((row for row in range(rank, rows)
                      if work[row][column] != 0), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][column]
        work[rank] = [entry / scale for entry in work[rank]]
        for row in range(rows):
            if row != rank and work[row][column] != 0:
                factor = work[row][column]
                work[row] = [left - factor * right
                             for left, right in zip(work[row], work[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def column(matrix, index):
    return tuple(row[index] for row in matrix)


selector_word = tuple(int(index in SELECTOR_SUPPORT) for index in range(N))
geodesic_word = tuple(int(index in GEODESIC_SUPPORT) for index in range(N))
selector_aug = augmentation_multiplication(selector_word)
geodesic_aug = augmentation_multiplication(geodesic_word)
selector_orbit = difference_orbit(selector_word)
geodesic_orbit = difference_orbit(geodesic_word)

selector_det = bareiss_determinant(selector_aug)
geodesic_det = bareiss_determinant(geodesic_aug)
selector_delta10 = determinantal_divisor(selector_aug, 10)
geodesic_delta10 = determinantal_divisor(geodesic_aug, 10)
selector_unit10 = first_minor_with_abs(selector_aug, 10, 1)
geodesic_unit9 = first_minor_with_abs(geodesic_aug, 9, 1)

require(abs(selector_det) == 35, ("selector determinant", selector_det))
require(abs(geodesic_det) == 81, ("geodesic determinant", geodesic_det))
require(selector_delta10 == 1 and selector_unit10 is not None,
        ("selector delta10", selector_delta10, selector_unit10))
require(geodesic_delta10 == 9 and geodesic_unit9 is not None,
        ("geodesic determinantal divisors", geodesic_delta10, geodesic_unit9))

# A direct unimodular minor of the combined 11 x 24 difference-orbit bank.
sparse_selector_shifts = (0, 1, 2, 3, 5, 6, 7, 8, 10, 11)
combined_columns = tuple(column(selector_orbit, shift)
                         for shift in sparse_selector_shifts) + (
                             column(geodesic_orbit, 0),
                         )
combined_certificate = bareiss_determinant(matrix_from_columns(combined_columns))
require(abs(combined_certificate) == 1,
        ("combined unimodular certificate", combined_certificate))

# Hostile type test: if the target-residue parity marker is incorrectly pulled
# back as an absolute vertex marker through J_12 and then written in ell order,
# the resulting word has mass eight, not five.  This is deliberately *not* the
# typed word used in THM-3280.
physical_pullback_word = tuple(
    int(J12[vertex] in GEODESIC_SUPPORT) for vertex in ELL_ORDER
)
physical_pullback_aug = augmentation_multiplication(physical_pullback_word)
physical_pullback_det = bareiss_determinant(physical_pullback_aug)
physical_pullback_rank = rank_over_q(physical_pullback_aug)
physical_pullback_delta8 = determinantal_divisor(physical_pullback_aug, 8)
physical_pullback_unit7 = first_minor_with_abs(physical_pullback_aug, 7, 1)

# A separate semantic hostile: THM-3277's transcript distinguishes ordinary
# hop distance from the number of edges in the clutch-product minimizer.  The
# defect 81 belongs to the latter.  Literal ordinary-hop parity has support
# {4,8} and a different 2-primary Smith signature.
ordinary_hop_support = tuple(
    residue for residue, distance in enumerate(ORDINARY_HOP_DISTANCES)
    if distance % 2 == 0
)
ordinary_hop_word = tuple(int(index in ordinary_hop_support) for index in range(N))
ordinary_hop_aug = augmentation_multiplication(ordinary_hop_word)
ordinary_hop_det = bareiss_determinant(ordinary_hop_aug)
ordinary_delta8 = determinantal_divisor(ordinary_hop_aug, 8)
ordinary_delta9 = determinantal_divisor(ordinary_hop_aug, 9)
ordinary_delta10 = determinantal_divisor(ordinary_hop_aug, 10)

require(ELL[17] == J12[17] == 0 and ELL[16] == J12[16] == 1,
        "root/generator normalization")
require(J12[16] == J12[18] == 1 and ELL[16] != ELL[18],
        "minimal vertex-fibre hostile")
ell_edge_increment = (ELL[16] - ELL[3]) % N
j12_edge_increment = (J12[16] - J12[3]) % N
require((ell_edge_increment, j12_edge_increment) == (2, 3),
        "physical edge mismatch")
require(sum(physical_pullback_word) == 8,
        ("physical pullback mass", sum(physical_pullback_word)))
require((physical_pullback_rank, physical_pullback_delta8,
         physical_pullback_unit7 is not None) == (8, 2, True),
        ("physical pullback Smith data", physical_pullback_rank,
         physical_pullback_delta8, physical_pullback_unit7))
require((ordinary_hop_support, abs(ordinary_hop_det), ordinary_delta8,
         ordinary_delta9, ordinary_delta10)
        == ((4, 8), 8, 1, 2, 4),
        ("ordinary-hop parity signature", ordinary_hop_support,
         ordinary_hop_det, ordinary_delta8, ordinary_delta9, ordinary_delta10))

print("INDEPENDENT THM-3280 LATTICE/TYPE AUDIT (NO SYMPY)")
print("selector_det=%d,delta10=%d,SNF=(1^10,35),unit10=%s" %
      (selector_det, selector_delta10, selector_unit10))
print("geodesic_det=%d,delta10=%d,unit9=%s,SNF=(1^9,9,9)" %
      (geodesic_det, geodesic_delta10, geodesic_unit9))
print("combined_unimodular_minor=%d,selector_shifts=%s,geodesic_shift=0,SNF=(1^11)" %
      (combined_certificate, sparse_selector_shifts))
print("bezout=1=(-37)*35+16*81")
print("vertex_fibre_hostile=(j12(16)=j12(18)=1;ell(16)=1,ell(18)=2)")
print("edge_hostile=(3->16;ell_increment=%d,j12_increment=%d)" %
      (ell_edge_increment, j12_edge_increment))
print("wrong_physical_pullback_support=%s,mass=%d,rank=%d,SNF=(1^7,2,0,0,0),det=%d" %
      (tuple(index for index, value in enumerate(physical_pullback_word) if value),
       sum(physical_pullback_word), physical_pullback_rank, physical_pullback_det))
print("ordinary_hop_parity_support=%s,index=%d,SNF=(1^8,2,2,2)" %
      (ordinary_hop_support, abs(ordinary_hop_det)))
print("metric_scope=defect81_is_edge_count_parity_of_clutch_product_minimizers_not_ordinary_hop_distance")
print("verdict=abstract_gauge_identity_valid;physical_vertex_identification_false")
