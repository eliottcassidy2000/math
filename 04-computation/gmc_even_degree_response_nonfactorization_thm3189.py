#!/usr/bin/env python3
"""Exact even-degree response nonfactorization controls for THM-3189."""

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations_with_replacement
from math import gcd
from pathlib import Path

import numpy as np


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
UPSTREAM = HERE / "gmc_pole_prefix_hasse_current_scout.py"
UPSTREAM_BANK = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
THM3184_SCRIPT = HERE / "gmc_depth_seven_degree_fourteen_farkas_thm3184.py"
THM3184_OUTPUT = ROOT / "05-knowledge/results/gmc_depth_seven_degree_fourteen_farkas_thm3184.out"
DEPENDENCIES = (
    (UPSTREAM,
     "151edb9b8cee4807d3f08dc17af32e45420021ba30dfd116c04d9fcaf8bbd5b7"),
    (UPSTREAM_BANK,
     "15b94691d53afbcdc6aefda89fc7cd5497534ca70fb780a686892dabb5646d6f"),
    (THM3184_SCRIPT,
     "9d560fc04fd7772a0b92b5973cd85b7ed61a2ec76328badfe5420ae2fcb7d129"),
    (THM3184_OUTPUT,
     "d599057d44b12cd222a0d6ecf3755cb5d2ae59d1faa1f7bd6731ad295ab667b3"),
)


def lf_hash(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for path, expected in DEPENDENCIES:
    require(lf_hash(path) == expected, ("dependency hash drift", str(path)))


def load_upstream_prefix(maximum_degree):
    tree = ast.parse(UPSTREAM.read_text(encoding="utf-8"))
    prefix = []
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "MAXIMUM_DEGREE"
                        for target in node.targets)):
            node.value = ast.Constant(maximum_degree)
        prefix.append(node)
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "UNIVERSE"
                        for target in node.targets)):
            break
    module = ast.fix_missing_locations(ast.Module(body=prefix, type_ignores=[]))
    namespace = {"__file__": str(UPSTREAM)}
    exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


UP = load_upstream_prefix(14)
UP["residual_roots"] = lru_cache(maxsize=None)(UP["residual_roots"])
UP["all_monomial_values"] = lru_cache(maxsize=None)(UP["all_monomial_values"])
coefficient_vectors = UP["coefficient_vectors"]
reduced_poles = UP["reduced_poles"]
BANK = UP["BANKS"][1]


@lru_cache(maxsize=None)
def partitions(total, maximum=None):
    if total == 0:
        return ((),)
    if maximum is None or maximum > total:
        maximum = total
    answer = []
    for first in range(maximum, 0, -1):
        for tail in partitions(total - first, first):
            answer.append((first,) + tail)
    return tuple(answer)


@lru_cache(maxsize=None)
def coarsens(fine, coarse):
    if sum(fine) != sum(coarse) or len(fine) < len(coarse):
        return False
    pieces = tuple(sorted(fine, reverse=True))
    bins = tuple(sorted(coarse, reverse=True))

    @lru_cache(maxsize=None)
    def place(index, capacities):
        if index == len(pieces):
            return not any(capacities)
        piece = pieces[index]
        tried = set()
        for slot, capacity in enumerate(capacities):
            if capacity < piece or capacity in tried:
                continue
            tried.add(capacity)
            changed = list(capacities)
            changed[slot] -= piece
            changed.sort(reverse=True)
            if place(index + 1, tuple(changed)):
                return True
        return False

    return place(0, bins)


def ones(count):
    return (1,) * count


# (source degree, target degree, expected upset size, minimal generators).
TARGETS = (
    (6, 8, 21, ((2,) + ones(6),)),
    (8, 10, 40, ((3,) + ones(7), (2, 2) + ones(6))),
    (10, 12, 76, ((2,) + ones(10),)),
    (10, 12, 74, ((2, 2) + ones(8),)),
    (12, 14, 130, ((3, 2) + ones(9), (2, 2, 2) + ones(8))),
    (12, 14, 128, ((4,) + ones(10), (3, 2) + ones(9),
                   (2, 2, 2, 2, 2, 2, 1, 1))),
    (12, 14, 121, ((7,) + ones(7), (3, 3, 2) + ones(6),
                   (2, 2, 2, 2) + ones(6))),
    (12, 14, 132, ((2, 2) + ones(10),)),
    (12, 14, 129, ((5,) + ones(9), (3, 3) + ones(8),
                   (2, 2, 2) + ones(8))),
    (12, 14, 127, ((5,) + ones(9), (4, 2) + ones(8),
                   (3, 3) + ones(8), (2, 2, 2, 2) + ones(6))),
)


UPSETS = []
for _, degree, expected_size, generators in TARGETS:
    bank = partitions(degree)
    upset = frozenset(
        shape for shape in bank
        if any(coarsens(generator, shape) for generator in generators)
    )
    require(len(upset) == expected_size, ("upset size drift", degree))
    minimal = tuple(
        shape for shape in bank
        if shape in upset and not any(
            other != shape and other in upset and coarsens(other, shape)
            for other in bank
        )
    )
    require(set(minimal) == set(generators),
            ("minimal generator drift", degree, generators, minimal))
    UPSETS.append(upset)
UPSETS = tuple(UPSETS)


POLES, _ = reduced_poles(1, BANK, 1, 3)
MULTIPLICITY = Counter(POLES)
VALUES = tuple(sorted(MULTIPLICITY))
BY_DEPTH = tuple(
    tuple(
        state
        for state in combinations_with_replacement(VALUES, depth)
        if all(Counter(state)[value] <= MULTIPLICITY[value]
               for value in set(state))
    )
    for depth in range(1, 8)
)
STATES = tuple(state for layer in BY_DEPTH for state in layer)
require(tuple(map(len, BY_DEPTH)) == (8, 33, 93, 200, 348, 507, 631)
        and len(STATES) == 1820, "state census drift")


SOURCE_DEGREES = (6, 8, 10, 12)
SOURCE_SHAPES = {degree: partitions(degree) for degree in SOURCE_DEGREES}
SOURCE_ROWS = {
    degree: [[] for _ in SOURCE_SHAPES[degree]]
    for degree in SOURCE_DEGREES
}
TARGET_ROWS = [[] for _ in TARGETS]
BALANCE_CHECKS = 0
for state in STATES:
    vectors = coefficient_vectors(1, BANK, 1, 3, state)
    for degree in SOURCE_DEGREES:
        vector = vectors[degree]
        require(sum(vector.values()) == 0,
                ("source zero-mass drift", state, degree))
        BALANCE_CHECKS += 1
        for row, shape in zip(SOURCE_ROWS[degree], SOURCE_SHAPES[degree]):
            value = vector[shape]
            require(value.denominator == 1, "source response lost integrality")
            row.append(value.numerator)
    for row, (_, degree, _, _), upset in zip(TARGET_ROWS, TARGETS, UPSETS):
        value = sum(vectors[degree][shape] for shape in upset)
        require(value.denominator == 1, "target response lost integrality")
        row.append(value.numerator)


SOURCE_ROWS = {degree: tuple(tuple(row) for row in rows)
               for degree, rows in SOURCE_ROWS.items()}
TARGET_ROWS = tuple(tuple(row) for row in TARGET_ROWS)
PRIME = 1000003


def rank_mod(rows):
    matrix = np.array(rows, dtype=object)
    matrix = np.vectorize(lambda value: int(value) % PRIME,
                          otypes=[np.int64])(matrix)
    rank = 0
    for column in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        if pivot != rank:
            matrix[[rank, pivot]] = matrix[[pivot, rank]]
        inverse = pow(int(matrix[rank, column]), -1, PRIME)
        matrix[rank] = matrix[rank] * inverse % PRIME
        active = np.flatnonzero(matrix[:, column])
        active = active[active != rank]
        for row in active:
            matrix[row] = (
                matrix[row] - int(matrix[row, column]) * matrix[rank]
            ) % PRIME
        rank += 1
        if rank == matrix.shape[0]:
            break
    return rank


SOURCE_RANKS = []
AFFINE_SOURCE_RANKS = []
INJECTIVE_COUNTS = []
for degree in SOURCE_DEGREES:
    rows = SOURCE_ROWS[degree]
    expected = len(SOURCE_SHAPES[degree]) - 1
    source_rank = rank_mod(rows)
    affine_rank = rank_mod(rows + ((1,) * len(STATES),))
    require(source_rank == expected and affine_rank == expected + 1,
            ("source rank drift", degree, source_rank, affine_rank))
    SOURCE_RANKS.append((degree, len(rows), source_rank))
    AFFINE_SOURCE_RANKS.append((degree, affine_rank))
    columns = {
        tuple(row[index] for row in rows)
        for index in range(len(STATES))
    }
    require(len(columns) == len(STATES),
            ("individual state responses collided", degree))
    INJECTIVE_COUNTS.append((degree, len(columns)))


RANK_JUMPS = []
AFFINE_RANK_JUMPS = []
for target, row in zip(TARGETS, TARGET_ROWS):
    source_degree, target_degree, upset_size, generators = target
    source = SOURCE_ROWS[source_degree]
    source_rank = len(SOURCE_SHAPES[source_degree]) - 1
    augmented_rank = rank_mod(source + (row,))
    affine = source + ((1,) * len(STATES),)
    affine_rank = source_rank + 1
    affine_augmented_rank = rank_mod(affine + (row,))
    require(augmented_rank == source_rank + 1,
            ("linear rank failed to jump", target))
    require(affine_augmented_rank == affine_rank + 1,
            ("affine rank failed to jump", target))
    RANK_JUMPS.append((source_degree, target_degree, upset_size,
                       source_rank, augmented_rank, generators))
    AFFINE_RANK_JUMPS.append((source_degree, target_degree,
                              affine_rank, affine_augmented_rank))


# Circuit-minimal affine hostile for the first transition.  The positive and
# negative parts give two rational probability laws with the same complete
# degree-six response and different selected degree-eight response.
CIRCUIT_INDICES = tuple(range(12))
CIRCUIT_STATES = tuple(STATES[index] for index in CIRCUIT_INDICES)
require(CIRCUIT_STATES == (
    (1,), (2,), (3,), (4,), (5,), (6,), (7,), (8,),
    (1, 1), (1, 2), (1, 3), (1, 4),
), "affine circuit state order drift")
CIRCUIT = (
    -1047277306414574923871547479704943667621575439892091,
    7504006088796364741516689429373777425450644438777969,
    -23076608092958877009276603275374222862197557645092271,
    39446840155066170582433740489847512927510928685260765,
    -40481016835146433236610937331785735797691389653391585,
    24940417689554928500780094388795198827192237041614851,
    -8542034917490033954128648539285677801015113585173749,
    1254725791744637777888128377808698770041845615340271,
    4534007100750039738438417212526753888724041742080,
    -8125859585348642127289323848561206793955708255360,
    6660162235070956799665346800833144498192277589760,
    -2120882902654833141730499839406513262980068520640,
)
require(reduce(gcd, (abs(value) for value in CIRCUIT), 0) == 1
        and all(CIRCUIT) and sum(CIRCUIT) == 0,
        "affine circuit lost primitivity or normalization")
degree_six = SOURCE_ROWS[6]
require(all(sum(coefficient * row[index]
                for coefficient, index in zip(CIRCUIT, CIRCUIT_INDICES)) == 0
            for row in degree_six),
        "affine circuit does not cancel complete degree-six response")
common_mass = sum(value for value in CIRCUIT if value > 0)
require(common_mass == -sum(value for value in CIRCUIT if value < 0),
        "affine circuit masses differ")
target_pairing = sum(
    coefficient * TARGET_ROWS[0][index]
    for coefficient, index in zip(CIRCUIT, CIRCUIT_INDICES)
)
require(target_pairing
        == 145762976976221426428610795350527471841760941217242422769712640,
        "affine circuit target pairing drift")
restricted_affine = tuple(
    tuple(row[index] for index in CIRCUIT_INDICES)
    for row in degree_six
) + ((1,) * len(CIRCUIT_INDICES),)
PROPER_SUBSET_CHECKS = 0
for omitted in range(len(CIRCUIT_INDICES)):
    restricted = tuple(
        tuple(value for index, value in enumerate(row) if index != omitted)
        for row in restricted_affine
    )
    require(rank_mod(restricted) == 11,
            ("affine circuit proper subset became dependent", omitted))
    PROPER_SUBSET_CHECKS += 1
target_difference = Fraction(target_pairing, common_mass)


print("THM-3189 even-degree response nonfactorization exact control")
print("dependency_hash_checks=" + repr(len(DEPENDENCIES)))
print("prime=" + repr(PRIME))
print("poles=" + repr(POLES))
print("state_counts_depth1_depth7_total="
      + repr((tuple(map(len, BY_DEPTH)), len(STATES))))
print("source_balance_checks=" + repr(BALANCE_CHECKS))
print("source_degree_partition_modular_ranks=" + repr(tuple(SOURCE_RANKS)))
print("affine_source_modular_ranks=" + repr(tuple(AFFINE_SOURCE_RANKS)))
print("individual_response_injective_counts=" + repr(tuple(INJECTIVE_COUNTS)))
print("linear_target_rank_jumps=" + repr(tuple(RANK_JUMPS)))
print("affine_target_rank_jumps=" + repr(tuple(AFFINE_RANK_JUMPS)))
print("affine_circuit_states=" + repr(CIRCUIT_STATES))
print("affine_circuit_coefficients=" + repr(CIRCUIT))
print("affine_circuit_proper_subset_checks=" + repr(PROPER_SUBSET_CHECKS))
print("affine_circuit_common_probability_mass=" + repr(common_mass))
print("affine_circuit_target_pairing=" + repr(target_pairing))
print("affine_probability_target_difference=" + repr(target_difference))
print("scope=no_linear_or_affine_two_step_factorization_on_fixed_response_bank")
print("all_exact_checks=PASS")
