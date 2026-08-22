#!/usr/bin/env python3
"""Exact decision-policy compiler and trace obstruction for the THM-3244 atlas."""

import ast
import contextlib
import hashlib
import io
import runpy
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / "04-computation/gmc_unique_reset_rips_nonmorse_thm3244.py":
        "3ff0babc41e35e6a185b0ff442cfb9284d9688360c0b96cd947c1128e16400ba",
    ROOT / "05-knowledge/results/gmc_unique_reset_rips_nonmorse_thm3244.out":
        "27dcd7c68e628465a1f09a564be0be366ded6075ef009a3d85d029f8f18605c9",
}


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


for dependency, expected in DEPENDENCIES.items():
    require(lf_sha256(dependency) == expected, ("dependency drift", dependency.name))

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)),
        "optimization-sensitive assert")
require(not any(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
), "floating literal")

# Loading the exact companion also runs all of its dependency, provenance,
# certificate, and hostile checks.  Its transcript is captured so this
# companion has one stable output surface.
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    atlas = runpy.run_path(str(next(iter(DEPENDENCIES))))
require("all_exact_checks=PASS" in captured.getvalue(), "THM-3244 replay failed")

capacities = atlas["CAPACITIES"]
states = atlas["states"]
count_vectors = atlas["count_vectors"]
state_index = atlas["state_index"]
reset = atlas["RESET"]
reset_index = state_index[reset]
reset_counts = atlas["reset_counts"]
q_distances = atlas["q_distances"]
q_neighbor_index = atlas["q_neighbor_index"]
row_direction_masks = atlas["row_direction_masks"]
row_two = atlas["row_two"]
row_ten = atlas["row_ten"]
nonreset_indices = atlas["nonreset_indices"]

# Only exclusive states constrain a chart classifier.  On an overlap state
# either label is lawful.
row_two_only = tuple(sorted(row_two - row_ten))
row_ten_only = tuple(sorted(row_ten - row_two))
classified_indices = row_two_only + row_ten_only
classified_counts = tuple(count_vectors[index] for index in classified_indices)
classified_size = len(classified_indices)
all_mask = (1 << classified_size) - 1
positive_mask = (1 << len(row_two_only)) - 1
negative_mask = all_mask ^ positive_mask


def is_pure(mask):
    return not (mask & positive_mask) or not (mask & negative_mask)


# Every nontrivial real threshold in one multiplicity coordinate is equivalent
# on this integer box to exactly one test n_j <= t with 0 <= t < capacity_j.
thresholds = []
for coordinate, capacity in enumerate(capacities):
    for threshold in range(capacity):
        mask = sum(
            1 << position
            for position, counts in enumerate(classified_counts)
            if counts[coordinate] <= threshold
        )
        thresholds.append((coordinate, threshold, mask))
thresholds = tuple(thresholds)
require(len(thresholds) == 16, "primitive threshold census changed")

infinity = classified_size + 1
bounded_choice = {}


@lru_cache(None)
def minimum_leaves_at_depth(mask, depth):
    if is_pure(mask):
        return 1
    if depth == 0:
        return infinity
    best = infinity
    best_split = None
    for split_index, (_, _, left_mask) in enumerate(thresholds):
        left = mask & left_mask
        right = mask & (all_mask ^ left_mask)
        if not left or not right:
            continue
        cost = (
            minimum_leaves_at_depth(left, depth - 1)
            + minimum_leaves_at_depth(right, depth - 1)
        )
        if cost < best:
            best = cost
            best_split = split_index
    bounded_choice[(mask, depth)] = best_split
    return best


depth_profile = tuple(
    minimum_leaves_at_depth(all_mask, depth) for depth in range(6)
)
require(depth_profile == (infinity, infinity, infinity, infinity, infinity, 15),
        "minimum-depth profile changed")

unbounded_choice = {}


@lru_cache(None)
def minimum_leaves(mask):
    if is_pure(mask):
        return 1
    best = infinity
    best_split = None
    for split_index, (_, _, left_mask) in enumerate(thresholds):
        left = mask & left_mask
        right = mask & (all_mask ^ left_mask)
        if not left or not right:
            continue
        cost = minimum_leaves(left) + minimum_leaves(right)
        if cost < best:
            best = cost
            best_split = split_index
    unbounded_choice[mask] = best_split
    return best


require(minimum_leaves(all_mask) == 15, "global minimum leaf count changed")


def recover_bounded_tree(mask, depth):
    if is_pure(mask):
        return 2 if mask & positive_mask else 10
    split_index = bounded_choice[(mask, depth)]
    require(split_index is not None, "missing bounded split")
    coordinate, threshold, left_mask = thresholds[split_index]
    left = mask & left_mask
    right = mask & (all_mask ^ left_mask)
    return (
        coordinate + 1,
        threshold,
        recover_bounded_tree(left, depth - 1),
        recover_bounded_tree(right, depth - 1),
    )


tree = recover_bounded_tree(all_mask, 5)
expected_tree = (
    8, 0,
    (7, 0,
     (3, 0,
      10,
      (6, 0, (4, 1, 10, 2), (5, 0, 10, 2))),
     2),
    (5, 1,
     (1, 1, (2, 0, (4, 1, 2, 10), 10), 10),
     (2, 1, (4, 1, 10, 2), (2, 2, (3, 1, 2, 10), 2))),
)
require(tree == expected_tree, "canonical minimum tree changed")


def tree_statistics(node):
    if isinstance(node, int):
        return (0, 1)
    left_depth, left_leaves = tree_statistics(node[2])
    right_depth, right_leaves = tree_statistics(node[3])
    return (1 + max(left_depth, right_depth), left_leaves + right_leaves)


require(tree_statistics(tree) == (5, 15), "tree shape changed")


def select_row(counts, node=tree):
    while not isinstance(node, int):
        coordinate, threshold, left, right = node
        node = left if counts[coordinate - 1] <= threshold else right
    return node


assignments = {2: set(), 10: set()}
for index in nonreset_indices:
    assignments[select_row(count_vectors[index])].add(index)
require(assignments[2] <= row_two, "tree selected unavailable row 2")
require(assignments[10] <= row_ten, "tree selected unavailable row 10")
require(assignments[2] | assignments[10] == nonreset_indices,
        "tree failed to classify the nonreset bank")
require((len(assignments[2]), len(assignments[10])) == (1955, 2363),
        "tree assignment census changed")

# Compile one deterministic operation trace: at each state use the tree's row
# and then the least-index Q-directed coordinate certified by that row.  This
# tie-break is a control, not part of the optimal-tree assertion.
chart_words = []
unsigned_coordinate_words = []
signed_edit_words = []
chart_coordinate_words = []
switch_counts = []
route_lengths = []
unsigned_word_fibres = defaultdict(list)

for start in sorted(nonreset_indices):
    index = start
    chart_word = []
    unsigned_word = []
    signed_word = []
    chart_coordinate_word = []
    while index != reset_index:
        counts = count_vectors[index]
        row = select_row(counts)
        chart = 0 if row == 2 else 1
        direction_mask = row_direction_masks[chart][index]
        require(direction_mask, "selected chart has no certified direction")
        coordinate = (direction_mask & -direction_mask).bit_length() - 1
        sign = 1 if counts[coordinate] < reset_counts[coordinate] else -1
        target = q_neighbor_index(index, coordinate)
        require(target is not None, "canonical direction has no target")
        require(q_distances[target] + 1 == q_distances[index],
                "canonical route lost the Q grading")
        chart_word.append(row)
        unsigned_word.append(coordinate + 1)
        signed_word.append(sign * (coordinate + 1))
        chart_coordinate_word.append((row, coordinate + 1))
        index = target
    require(len(chart_word) == q_distances[start], "route length is not exact")
    chart_word = tuple(chart_word)
    unsigned_word = tuple(unsigned_word)
    signed_word = tuple(signed_word)
    chart_coordinate_word = tuple(chart_coordinate_word)
    chart_words.append(chart_word)
    unsigned_coordinate_words.append(unsigned_word)
    signed_edit_words.append(signed_word)
    chart_coordinate_words.append(chart_coordinate_word)
    switch_counts.append(sum(
        left != right for left, right in zip(chart_word, chart_word[1:])
    ))
    route_lengths.append(len(chart_word))
    unsigned_word_fibres[chart_coordinate_word].append(states[start])

require(Counter(route_lengths) == {
    1: 11, 2: 55, 3: 169, 4: 365, 5: 598, 6: 775, 7: 810,
    8: 685, 9: 467, 10: 250, 11: 101, 12: 28, 13: 4,
}, "route-length histogram changed")
require(Counter(switch_counts) == {0: 553, 1: 3158, 2: 387, 3: 220},
        "canonical switch histogram changed")

empty_word = ()
chart_classes = len(set(chart_words) | {empty_word})
unsigned_coordinate_classes = len(set(unsigned_coordinate_words) | {empty_word})
signed_edit_classes = len(set(signed_edit_words) | {empty_word})
chart_coordinate_classes = len(set(chart_coordinate_words) | {empty_word})
require((chart_classes, unsigned_coordinate_classes, chart_coordinate_classes,
         signed_edit_classes) == (296, 4002, 4318, 4319),
        "continuation-class census changed")

collisions = tuple(
    tuple(fibre) for fibre in unsigned_word_fibres.values() if len(fibre) > 1
)
require(collisions == (((2, 3, 3), (1, 1, 2, 3, 3)),),
        "chart-coordinate collision changed")

# Policy-independent trace lower bounds.  Every Q-monotone route has unsigned
# action Parikh vector |n-Q| and signed action Parikh vector Q-n.  The latter
# reconstructs the full state.  The former has the exact product census below.
absolute_deviations = {
    tuple(abs(value - target) for value, target in zip(counts, reset_counts))
    for counts in count_vectors
}
signed_displacements = {
    tuple(target - value for value, target in zip(counts, reset_counts))
    for counts in count_vectors
}
require(len(absolute_deviations) == 1536, "unsigned Parikh census changed")
require(len(signed_displacements) == len(states) == 4319,
        "signed edit trace ceased to reconstruct the state")


def multiply_polynomials(left, right):
    answer = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            answer[left_degree + right_degree] += left_coefficient * right_coefficient
    return answer


# The commutative distance enumerator is compact even though the signed trace
# continuation quotient is not.  The eight factors are the coordinatewise
# deviation enumerators; the missing zero physical state has Q-distance eight.
distance_factors = (
    (1, 2, 1, 1),
    (1, 1, 1, 1),
    (1, 1, 1),
    (1, 2),
    (1, 2),
    (1, 1),
    (1, 1),
    (1, 1),
)
distance_enumerator = [1]
for factor in distance_factors:
    distance_enumerator = multiply_polynomials(distance_enumerator, factor)
distance_enumerator[8] -= 1
distance_enumerator = tuple(distance_enumerator)
require(distance_enumerator == (
    1, 11, 55, 169, 365, 598, 775, 810, 685, 467, 250, 101, 28, 4,
), "factored distance enumerator changed")
require(tuple(distance_enumerator[1:]) == tuple(
    Counter(route_lengths).get(distance, 0)
    for distance in range(1, len(distance_enumerator))
), "factored enumerator disagrees with canonical routes")

tree_digest = hashlib.sha256(repr(tree).encode("ascii")).hexdigest()
trace_digest = hashlib.sha256(repr(tuple(
    zip(chart_words, unsigned_coordinate_words, signed_edit_words)
)).encode("ascii")).hexdigest()

print("FC3 THM-3244 two-row decision-policy compiler beach audit")
print("dependency_hash_checks=2")
print("assert_nodes=0,float_literals=0")
print("exclusive_classes=(row2_only=304,row10_only=103),primitive_thresholds=16")
print("axis_tree_depth_profile=(0:IMPOSSIBLE,1:IMPOSSIBLE,2:IMPOSSIBLE,3:IMPOSSIBLE,4:IMPOSSIBLE,5:15_leaves)")
print("axis_tree_global_minimum_leaves=15")
print("axis_tree=(coord,threshold,true,false)=" + repr(tree))
print("axis_tree_sha256=" + tree_digest)
print("full_bank_assignment=(row2=1955,row10=2363),unavailable_choices=0")
print("canonical_least_coordinate_routes=4318/4318,exact_Q_distance=PASS")
print("canonical_route_length_histogram=" + repr(tuple(sorted(Counter(route_lengths).items()))))
print("canonical_chart_switch_histogram=" + repr(tuple(sorted(Counter(switch_counts).items()))))
print("continuation_classes=(chart=296,unsigned_coordinate=4002,chart_plus_coordinate=4318,signed_edit=4319)")
print("unique_chart_coordinate_collision=" + repr(collisions[0]))
print("policy_independent_parikh_lower_bounds=(unsigned_coordinate=1536,signed_edit=4319)")
print("distance_enumerator=(1+2z+z^2+z^3)(1+z+z^2+z^3)(1+z+z^2)(1+2z)^2(1+z)^3-z^8")
print("distance_coefficients=" + repr(distance_enumerator))
print("canonical_trace_sha256=" + trace_digest)
print("typed_scope=optimal_axis_threshold_chart_classifier; signed operation trace remains injective; no FC3 closure")
print("all_exact_checks=PASS")
