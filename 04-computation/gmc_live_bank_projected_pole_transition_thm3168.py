#!/usr/bin/env python3
"""Exact controls for THM-3168's live-bank projected transition wall."""

import ast
from collections import Counter
from functools import reduce
from itertools import combinations_with_replacement
from math import gcd
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_pole_prefix_hasse_current_scout.py"


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
    module = ast.fix_missing_locations(
        ast.Module(body=prefix, type_ignores=[]))
    namespace = {"__file__": str(UPSTREAM)}
    exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


UP = load_upstream_prefix(13)
partitions = UP["partitions"]
coefficient_vectors = UP["coefficient_vectors"]
reduced_poles = UP["reduced_poles"]
BANK = UP["BANKS"][1]

POLES, _ = reduced_poles(1, BANK, 1, 3)
MULTIPLICITY = Counter(POLES)
VALUES = tuple(sorted(MULTIPLICITY))
require(POLES == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "support-(1,3) pole multiset drift")

STATES = ((),) + tuple(
    state
    for depth in range(1, 4)
    for state in combinations_with_replacement(VALUES, depth)
    if all(Counter(state)[value] <= MULTIPLICITY[value]
           for value in set(state))
)
require(len(STATES) == 135, "empty-through-depth-three state census drift")

VECTORS = {
    state: coefficient_vectors(1, BANK, 1, 3, state)
    for state in STATES
}
for state, degrees in VECTORS.items():
    for degree, vector in degrees.items():
        require(all(value.denominator == 1 for value in vector.values()),
                f"current lost integrality at {(state, degree)}")
        require(sum(vector.values()) == 0,
                f"current lost total-mass cancellation at {(state, degree)}")


def pairs(pole, cap):
    answer = []
    for state in STATES:
        if len(state) > cap:
            continue
        if Counter(state)[pole] >= MULTIPLICITY[pole]:
            continue
        child = tuple(sorted(state + (pole,)))
        require(child in VECTORS, "legal child left reconstructed state bank")
        answer.append((state, child))
    return tuple(answer)


def matrices(pole, cap, degree):
    pair_bank = pairs(pole, cap)
    shapes = partitions(degree)[:-1]
    parent = sp.Matrix([
        [VECTORS[state][degree][shape].numerator
         for state, _ in pair_bank]
        for shape in shapes
    ])
    child = sp.Matrix([
        [VECTORS[next_state][degree][shape].numerator
         for _, next_state in pair_bank]
        for shape in shapes
    ])
    return pair_bank, parent, child


PRIME = 65521
require(all(PRIME % divisor for divisor in range(2, int(PRIME ** .5) + 1)),
        "finite-field rank modulus lost primality")


def rank_mod_prime(matrix):
    work = [[int(value) % PRIME for value in row]
            for row in matrix.tolist()]
    row_count = len(work)
    column_count = len(work[0]) if work else 0
    rank = 0
    for column in range(column_count):
        pivot = next((row for row in range(rank, row_count)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], PRIME - 2, PRIME)
        work[rank] = [(value * inverse) % PRIME
                      for value in work[rank]]
        for row in range(row_count):
            if row == rank or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [
                (left - scale * right) % PRIME
                for left, right in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == row_count:
            break
    return rank


# ---------------------------------------------------------------------------
# 1. Exact finite-span survivor through parent depth one.
# ---------------------------------------------------------------------------

SHALLOW_PAIR_COUNTS = []
SHALLOW_RANKS = []
SHALLOW_KERNEL_DIMENSIONS_N5_N6 = []
for pole in VALUES:
    pole_ranks = []
    kernel_dimensions = []
    for degree in range(5, 14):
        pair_bank, parent, child = matrices(pole, 1, degree)
        if degree <= 6:
            parent_rank = parent.rank()
            kernel = parent.nullspace()
            require(all(child * relation == sp.zeros(child.rows, 1)
                        for relation in kernel),
                    f"shallow child lost parent relation at {(pole, degree)}")
            stacked_rank = parent.col_join(child).rank()
            require(stacked_rank == parent_rank,
                    f"shallow exact rank equality failed at {(pole, degree)}")
            kernel_dimensions.append(len(kernel))
        else:
            parent_rank = rank_mod_prime(parent)
            require(parent_rank == len(pair_bank),
                    f"shallow parent columns lost independence at {(pole, degree)}")
            stacked_rank = parent_rank
        pole_ranks.append((parent_rank, stacked_rank))
    SHALLOW_PAIR_COUNTS.append(len(pairs(pole, 1)))
    SHALLOW_RANKS.append(tuple(pole_ranks))
    SHALLOW_KERNEL_DIMENSIONS_N5_N6.append(tuple(kernel_dimensions))

require(tuple(SHALLOW_PAIR_COUNTS) == (9, 9, 9, 9, 9, 8, 8, 8),
        "shallow fixed-pole pair census drift")
require(all(profile == ((6, 6), (8, 8)) + ((9, 9),) * 7
            for profile in SHALLOW_RANKS[:5]),
        "repeated-pole shallow rank profile drift")
require(all(profile == ((6, 6),) + ((8, 8),) * 8
            for profile in SHALLOW_RANKS[5:]),
        "singleton-pole shallow rank profile drift")
require(tuple(SHALLOW_KERNEL_DIMENSIONS_N5_N6)
        == ((3, 1),) * 5 + ((2, 0),) * 3,
        "shallow exact kernel-dimension profile drift")


# ---------------------------------------------------------------------------
# 2. Degree-five depth-two rank wall for every physical pole value.
# ---------------------------------------------------------------------------

DEPTH_TWO_ROWS = []
for pole in VALUES:
    pair_bank, parent, child = matrices(pole, 2, 5)
    parent_rank = rank_mod_prime(parent)
    child_rank = rank_mod_prime(child)
    stacked_rank = rank_mod_prime(parent.col_join(child))
    require((parent_rank, child_rank, stacked_rank) == (6, 6, 12),
            f"depth-two rank wall drift at pole {pole}")
    DEPTH_TWO_ROWS.append(
        (pole, len(pair_bank), parent_rank, child_rank, stacked_rank))

require(tuple(row[1] for row in DEPTH_TWO_ROWS)
        == (42, 42, 41, 41, 41, 34, 34, 34),
        "depth-two fixed-pole pair census drift")


# ---------------------------------------------------------------------------
# 3. Direct seven-pair determinant for pole one.
# ---------------------------------------------------------------------------

WITNESS_PAIRS = (
    ((), (1,)),
    ((1,), (1, 1)),
    ((2,), (1, 2)),
    ((3,), (1, 3)),
    ((4,), (1, 4)),
    ((5,), (1, 5)),
    ((1, 1), (1, 1, 1)),
)
SHAPES_FIVE = partitions(5)
NORMALIZED_PARENT_COLUMNS = []
NORMALIZED_CHILD_COLUMNS = []
NORMALIZERS = []
for parent_state, child_state in WITNESS_PAIRS:
    parent_column = [
        VECTORS[parent_state][5][shape].numerator for shape in SHAPES_FIVE
    ]
    child_column = [
        VECTORS[child_state][5][shape].numerator for shape in SHAPES_FIVE
    ]
    normalizer = reduce(gcd, (
        abs(value) for value in parent_column + child_column if value
    ))
    NORMALIZERS.append(normalizer)
    NORMALIZED_PARENT_COLUMNS.append([
        value // normalizer for value in parent_column
    ])
    NORMALIZED_CHILD_COLUMNS.append([
        value // normalizer for value in child_column
    ])

require(tuple(NORMALIZERS) == (1440, 720, 1440, 720, 1440, 720, 720),
        "seven-pair normalization drift")

DETERMINANT_ROWS = tuple(
    tuple(column[row] for column in NORMALIZED_PARENT_COLUMNS)
    for row in range(6)
) + (tuple(column[0] for column in NORMALIZED_CHILD_COLUMNS),)
EXPECTED_MATRIX = (
    (61987, 123972, 61955, 123488, 60963, 117724, 123970),
    (264982, 512220, 246780, 471920, 222210, 407844, 494480),
    (212848, 422652, 205980, 391712, 180480, 320196, 419612),
    (497718, 917412, 419790, 759408, 338238, 591756, 842154),
    (431766, 816756, 375090, 670584, 292542, 500412, 770878),
    (-500833, -986466, -478553, -908966, -420417, -753946, -970113),
    (61986, 123970, 61954, 123486, 60962, 117722, 123968),
)
require(DETERMINANT_ROWS == EXPECTED_MATRIX,
        "seven-pair normalized matrix drift")
DETERMINANT = int(sp.Matrix(DETERMINANT_ROWS).det())
EXPECTED_DETERMINANT = -3821168719915360450560
require(DETERMINANT == EXPECTED_DETERMINANT,
        "seven-pair determinant drift")
require(sp.factorint(abs(DETERMINANT))
        == {2: 11, 3: 6, 5: 1, 7: 1, 11: 1,
            37: 1, 877: 1, 204869207: 1},
        "seven-pair determinant factorization drift")
ROW_DIFFERENCE = tuple(
    child - parent
    for parent, child in zip(DETERMINANT_ROWS[0], DETERMINANT_ROWS[-1])
)
require(ROW_DIFFERENCE == (-1, -2, -1, -2, -1, -2, -2),
        "normalized child/parent top-row difference drift")


print("THM-3168 live-bank projected pole transition wall")
print("poles=" + repr(POLES))
print("empty_through_depth3_states=" + repr(len(STATES)))
print("depth1_pair_counts_M1_M8=" + repr(tuple(SHALLOW_PAIR_COUNTS)))
print("depth1_kernel_dimensions_N5_N6="
      + repr(tuple(SHALLOW_KERNEL_DIMENSIONS_N5_N6)))
print("depth1_rank_profile_repeated_M1_M5=((6,6),(8,8),7x(9,9))")
print("depth1_rank_profile_singleton_M6_M8=((6,6),8x(8,8))")
print("depth2_N5_rows=" + repr(tuple(DEPTH_TWO_ROWS)))
print("witness_pairs=" + repr(WITNESS_PAIRS))
print("witness_normalizers=" + repr(tuple(NORMALIZERS)))
print("witness_determinant=" + repr(DETERMINANT))
print("witness_determinant_factorization="
      + repr(tuple(sorted(sp.factorint(abs(DETERMINANT)).items()))))
print("witness_child_minus_parent_top=" + repr(ROW_DIFFERENCE))
print("scope=finite_projected_linear_span_not_positive_or_Markov_transition")
print("all_exact_checks=PASS")
