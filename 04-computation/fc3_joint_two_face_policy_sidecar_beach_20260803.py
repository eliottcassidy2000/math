#!/usr/bin/env python3
"""Exact joint policy/sidecar compiler for the THM-3249 two-face atlas."""

import ast
import contextlib
import hashlib
import importlib.util
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
    ROOT / "04-computation/gmc_cross_support_upset_atlas_holonomy_thm3249.py":
        "bb901e92687544c69d67a55d057f8293ecbf516b80d491a636c4a62af19eebef",
    ROOT / "05-knowledge/results/gmc_cross_support_upset_atlas_holonomy_thm3249.out":
        "76d037d4f0737b37ad48531e23be1a0a37509ce0a3d029d3cefdd42662ec04f2",
    ROOT / "04-computation/gmc_first_shell_pair_clutch_thm3254.py":
        "05efd37eeedeca7e3be581977a894592a7873d94a966f06d9533482cc8498fee",
    ROOT / "05-knowledge/results/gmc_first_shell_pair_clutch_thm3254.out":
        "dd415c8ce6e2e196c115421d3508addabb724305a16843509264a8b3205beee9",
    ROOT / "04-computation/gmc_unique_reset_rips_nonmorse_thm3244.py":
        "3ff0babc41e35e6a185b0ff442cfb9284d9688360c0b96cd947c1128e16400ba",
    ROOT / "05-knowledge/results/gmc_unique_reset_rips_nonmorse_thm3244.out":
        "27dcd7c68e628465a1f09a564be0be366ded6075ef009a3d85d029f8f18605c9",
    ROOT / "04-computation/fc3_two_row_decision_policy_compiler_beach_20260803.py":
        "8910a0bc1a6bc10ac27742781142c28d89209b57bc0f6f5e461da314d0e0e197",
    ROOT / "05-knowledge/results/fc3_two_row_decision_policy_compiler_beach_20260803.out":
        "b4171403580c05a34e33f5cc3b0d4b67402dea81986c0459a7f3b29afc19529a",
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


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


cross = load_module("thm3249_cross_support", next(iter(DEPENDENCIES)))

captured = io.StringIO()
full_script = ROOT / "04-computation/gmc_unique_reset_rips_nonmorse_thm3244.py"
with contextlib.redirect_stdout(captured):
    full = runpy.run_path(str(full_script))
require("all_exact_checks=PASS" in captured.getvalue(), "THM-3244 replay failed")


def padded_counts(state):
    counts = Counter(state)
    return tuple(counts[value] for value in range(1, 9))


def coordinate_of_move(source, target):
    source_counts = padded_counts(source)
    target_counts = padded_counts(target)
    changed = tuple(
        coordinate for coordinate, (left, right) in enumerate(
            zip(source_counts, target_counts)
        ) if left != right
    )
    require(len(changed) == 1, "one-pole move changed multiple coordinates")
    require(abs(source_counts[changed[0]] - target_counts[changed[0]]) == 1,
            "one-pole move was not unit")
    return changed[0]


# Reconstruct the support-(1,2) responses from the pinned product-Gamma
# formula.  Unlike the large face, this face has no embedded direction bank.
small_states = cross.SMALL_STATES
small_state_index = {state: index for index, state in enumerate(small_states)}
small_count_vectors = tuple(padded_counts(state) for state in small_states)
small_reset = cross.SMALL_RESET
small_reset_index = small_state_index[small_reset]
small_reset_counts = small_count_vectors[small_reset_index]
small_rows = tuple(
    cross.row_coordinates(cross.T.coefficient_vectors(
        1, cross.BANK, 1, 2, state
    ))
    for state in small_states
)

small_direction_masks = []
for row_index in (1, 9):
    masks = []
    for state, response in zip(small_states, small_rows):
        mask = 0
        if state != small_reset:
            for target in cross.toward_reset_neighbours(
                    state, cross.SMALL_POLES, small_reset):
                target_index = small_state_index[target]
                if small_rows[target_index][row_index] > response[row_index]:
                    mask |= 1 << coordinate_of_move(state, target)
        masks.append(mask)
    small_direction_masks.append(tuple(masks))
small_direction_masks = tuple(small_direction_masks)
small_row_two = frozenset(
    index for index, mask in enumerate(small_direction_masks[0]) if mask
)
small_row_ten = frozenset(
    index for index, mask in enumerate(small_direction_masks[1]) if mask
)
small_nonreset = frozenset(range(len(small_states))) - {small_reset_index}
require((len(small_row_two), len(small_row_ten),
         len(small_row_two & small_row_ten),
         len(small_row_two - small_row_ten),
         len(small_row_ten - small_row_two),
         len(small_nonreset - (small_row_two | small_row_ten)))
        == (219, 186, 167, 52, 19, 0), "small row cover drift")

full_states = full["states"]
full_count_vectors = full["count_vectors"]
full_state_index = full["state_index"]
full_reset = full["RESET"]
full_reset_index = full_state_index[full_reset]
full_reset_counts = full["reset_counts"]
full_direction_masks = full["row_direction_masks"]
full_row_two = full["row_two"]
full_row_ten = full["row_ten"]
full_nonreset = full["nonreset_indices"]


def exclusive_records(face, counts, row_two, row_ten):
    records = []
    for index in sorted(row_two - row_ten):
        records.append((face, counts[index], 2))
    for index in sorted(row_ten - row_two):
        records.append((face, counts[index], 10))
    return tuple(records)


small_records = exclusive_records(
    0, small_count_vectors, small_row_two, small_row_ten
)
full_records = exclusive_records(
    1, full_count_vectors, full_row_two, full_row_ten
)
joint_records = small_records + full_records
require((len(small_records), len(full_records), len(joint_records)) == (71, 407, 478),
        "exclusive record census drift")


def make_tests(records, include_face, capacities):
    tests = []
    for coordinate, capacity in enumerate(capacities):
        for threshold in range(capacity):
            mask = sum(
                1 << position
                for position, (_, counts, _) in enumerate(records)
                if counts[coordinate] <= threshold
            )
            tests.append(("n", coordinate + 1, threshold, mask))
    if include_face:
        mask = sum(
            1 << position
            for position, (face, _, _) in enumerate(records)
            if face == 0
        )
        tests.append(("face", 0, 0, mask))
    return tuple(tests)


class TreeOptimizer:
    def __init__(self, records, tests):
        self.records = records
        self.tests = tests
        self.size = len(records)
        self.all_mask = (1 << self.size) - 1
        self.row_two_mask = sum(
            1 << position
            for position, (_, _, label) in enumerate(records)
            if label == 2
        )
        self.row_ten_mask = self.all_mask ^ self.row_two_mask
        self.infinity = self.size + 1
        self.bounded_choice = {}
        self.unbounded_choice = {}

    def is_pure(self, mask):
        return not (mask & self.row_two_mask) or not (mask & self.row_ten_mask)

    def minimum_leaves_at_depth(self, mask, depth):
        key = (mask, depth)
        if key in self._bounded_cache:
            return self._bounded_cache[key]
        if self.is_pure(mask):
            answer = 1
        elif depth == 0:
            answer = self.infinity
        else:
            answer = self.infinity
            best_split = None
            for split_index, (_, _, _, left_mask) in enumerate(self.tests):
                left = mask & left_mask
                right = mask & (self.all_mask ^ left_mask)
                if not left or not right:
                    continue
                cost = (
                    self.minimum_leaves_at_depth(left, depth - 1)
                    + self.minimum_leaves_at_depth(right, depth - 1)
                )
                if cost < answer:
                    answer = cost
                    best_split = split_index
            self.bounded_choice[key] = best_split
        self._bounded_cache[key] = answer
        return answer

    def minimum_leaves(self, mask):
        if mask in self._unbounded_cache:
            return self._unbounded_cache[mask]
        if self.is_pure(mask):
            answer = 1
        else:
            answer = self.infinity
            best_split = None
            for split_index, (_, _, _, left_mask) in enumerate(self.tests):
                left = mask & left_mask
                right = mask & (self.all_mask ^ left_mask)
                if not left or not right:
                    continue
                cost = self.minimum_leaves(left) + self.minimum_leaves(right)
                if cost < answer:
                    answer = cost
                    best_split = split_index
            self.unbounded_choice[mask] = best_split
        self._unbounded_cache[mask] = answer
        return answer

    def solve(self, maximum_depth):
        self._bounded_cache = {}
        self._unbounded_cache = {}
        profile = tuple(
            self.minimum_leaves_at_depth(self.all_mask, depth)
            for depth in range(maximum_depth + 1)
        )
        leaves = self.minimum_leaves(self.all_mask)
        return profile, leaves

    def recover(self, mask, depth):
        if self.is_pure(mask):
            return 2 if mask & self.row_two_mask else 10
        split_index = self.bounded_choice[(mask, depth)]
        require(split_index is not None, "missing bounded tree split")
        kind, coordinate, threshold, left_mask = self.tests[split_index]
        left = mask & left_mask
        right = mask & (self.all_mask ^ left_mask)
        return (
            kind, coordinate, threshold,
            self.recover(left, depth - 1),
            self.recover(right, depth - 1),
        )


small_tests = make_tests(small_records, False, (4, 3, 2, 1, 1, 0, 0, 0))
full_tests = make_tests(full_records, False, full["CAPACITIES"])
blind_tests = make_tests(joint_records, False, full["CAPACITIES"])
aware_tests = make_tests(joint_records, True, full["CAPACITIES"])
require((len(small_tests), len(full_tests), len(blind_tests), len(aware_tests))
        == (11, 16, 16, 17), "primitive test census drift")

small_optimizer = TreeOptimizer(small_records, small_tests)
full_optimizer = TreeOptimizer(full_records, full_tests)
blind_optimizer = TreeOptimizer(joint_records, blind_tests)
aware_optimizer = TreeOptimizer(joint_records, aware_tests)
small_profile, small_minimum_leaves = small_optimizer.solve(7)
full_profile, full_minimum_leaves = full_optimizer.solve(7)
blind_profile, blind_minimum_leaves = blind_optimizer.solve(8)
aware_profile, aware_minimum_leaves = aware_optimizer.solve(8)
require(small_profile == (72, 72, 72, 72, 8, 8, 8, 8)
        and small_minimum_leaves == 8, "small policy optimum drift")
require(full_profile == (408, 408, 408, 408, 408, 15, 15, 15)
        and full_minimum_leaves == 15, "full policy optimum drift")
require(blind_profile == (479,) * 9 and blind_minimum_leaves == 479,
        "face-blind impossibility profile drift")
require(aware_profile == (479, 479, 479, 479, 479, 479, 23, 23, 22)
        and aware_minimum_leaves == 22, "face-aware policy frontier drift")


def first_finite(profile, infinity):
    options = tuple(
        (depth, leaves) for depth, leaves in enumerate(profile) if leaves < infinity
    )
    return options[0] if options else None


small_depth_result = first_finite(small_profile, small_optimizer.infinity)
full_depth_result = first_finite(full_profile, full_optimizer.infinity)
aware_depth_result = first_finite(aware_profile, aware_optimizer.infinity)

feature_labels = defaultdict(set)
feature_sources = defaultdict(list)
for face, counts, label in joint_records:
    feature_labels[counts].add(label)
    feature_sources[counts].append((face, label))
blind_conflicts = tuple(sorted(
    (counts, tuple(sorted(feature_sources[counts])))
    for counts, labels in feature_labels.items() if len(labels) > 1
))
require(blind_conflicts, "face-blind policy unexpectedly became feasible")
require(first_finite(blind_profile, blind_optimizer.infinity) is None,
        "face-blind threshold tree unexpectedly exists")
require(blind_minimum_leaves == blind_optimizer.infinity,
        "unbounded face-blind tree unexpectedly exists")
require(len(blind_conflicts) == 19, "cross-face conflict census drift")
require(blind_conflicts[0] == (
    (0, 0, 0, 0, 1, 0, 0, 0), ((0, 2), (1, 10)),
), "minimal state-(5) face conflict drift")
require(all(sources == ((0, 2), (1, 10))
            for _, sources in blind_conflicts),
        "cross-face conflict orientation drift")
raw_count_intersection = set(small_count_vectors) & set(full_count_vectors)
require(len(raw_count_intersection) == 239
        and len(set(small_count_vectors) | set(full_count_vectors)) == 4319,
        "raw count-vector collision tax drift")

small_depth, _ = small_depth_result
full_depth, _ = full_depth_result
aware_depth, _ = aware_depth_result
small_tree = small_optimizer.recover(small_optimizer.all_mask, small_depth)
full_tree = full_optimizer.recover(full_optimizer.all_mask, full_depth)
aware_depth_tree = aware_optimizer.recover(aware_optimizer.all_mask, aware_depth)
aware_leaf_depth = next(
    depth for depth, leaves in enumerate(aware_profile)
    if leaves == aware_minimum_leaves
)
aware_leaf_tree = aware_optimizer.recover(
    aware_optimizer.all_mask, aware_leaf_depth
)


def tree_statistics(tree):
    if isinstance(tree, int):
        return 0, 1, 0
    left = tree_statistics(tree[3])
    right = tree_statistics(tree[4])
    return (
        1 + max(left[0], right[0]),
        left[1] + right[1],
        int(tree[0] == "face") + left[2] + right[2],
    )


def select_row(face, counts, tree):
    while not isinstance(tree, int):
        kind, coordinate, threshold, left, right = tree
        truth = face == 0 if kind == "face" else counts[coordinate - 1] <= threshold
        tree = left if truth else right
    return tree


def consulted_face(face, counts, tree):
    while not isinstance(tree, int):
        kind, coordinate, threshold, left, right = tree
        if kind == "face":
            return True
        tree = left if counts[coordinate - 1] <= threshold else right
    return False


def verify_tree(face, counts, nonreset, row_two, row_ten, tree):
    assignments = Counter()
    unavailable = []
    face_consults = 0
    for index in sorted(nonreset):
        row = select_row(face, counts[index], tree)
        assignments[row] += 1
        face_consults += int(consulted_face(face, counts[index], tree))
        if (row == 2 and index not in row_two) or (row == 10 and index not in row_ten):
            unavailable.append(index)
    require(not unavailable, "tree selected an unavailable row")
    return assignments, face_consults


small_assignment, small_face_consults = verify_tree(
    0, small_count_vectors, small_nonreset, small_row_two, small_row_ten,
    aware_depth_tree,
)
full_assignment, full_face_consults = verify_tree(
    1, full_count_vectors, full_nonreset, full_row_two, full_row_ten,
    aware_depth_tree,
)
small_leaf_assignment, small_leaf_face_consults = verify_tree(
    0, small_count_vectors, small_nonreset, small_row_two, small_row_ten,
    aware_leaf_tree,
)
full_leaf_assignment, full_leaf_face_consults = verify_tree(
    1, full_count_vectors, full_nonreset, full_row_two, full_row_ten,
    aware_leaf_tree,
)
require(tree_statistics(small_tree) == (4, 8, 0), "small tree shape drift")
require(tree_statistics(full_tree) == (5, 15, 0), "full tree shape drift")
require(tree_statistics(aware_depth_tree) == (6, 23, 1),
        "minimum-depth joint tree shape drift")
require(tree_statistics(aware_leaf_tree) == (8, 22, 2),
        "minimum-leaf joint tree shape drift")
require((small_face_consults, full_face_consults) == (238, 2159),
        "minimum-depth face-consult census drift")
require((small_leaf_face_consults, full_leaf_face_consults) == (198, 599),
        "minimum-leaf face-consult census drift")
require((small_assignment[2], small_assignment[10],
         full_assignment[2], full_assignment[10]) == (171, 67, 1955, 2363),
        "minimum-depth assignment census drift")


def q_neighbor_index(counts, state_index, reset_counts, index, coordinate):
    changed = list(counts[index])
    if changed[coordinate] < reset_counts[coordinate]:
        changed[coordinate] += 1
    elif changed[coordinate] > reset_counts[coordinate]:
        changed[coordinate] -= 1
    else:
        return None
    state = tuple(
        value for value, multiplicity in enumerate(changed, 1)
        for _ in range(multiplicity)
    )
    return state_index.get(state)


def compile_face(face, states, counts, state_index, reset_index, reset_counts,
                 nonreset, direction_masks, tree):
    records = []
    for start in sorted(range(len(states))):
        index = start
        chart_word = []
        coordinate_word = []
        signed_word = []
        chart_coordinate_word = []
        while index != reset_index:
            row = select_row(face, counts[index], tree)
            chart = 0 if row == 2 else 1
            mask = direction_masks[chart][index]
            require(mask, "joint controller selected an unavailable direction")
            coordinate = (mask & -mask).bit_length() - 1
            sign = 1 if counts[index][coordinate] < reset_counts[coordinate] else -1
            target = q_neighbor_index(
                counts, state_index, reset_counts, index, coordinate
            )
            require(target is not None, "joint route lost its target")
            chart_word.append(row)
            coordinate_word.append(coordinate + 1)
            signed_word.append(sign * (coordinate + 1))
            chart_coordinate_word.append((row, coordinate + 1))
            index = target
        distance = sum(abs(left - right) for left, right in zip(counts[start], reset_counts))
        require(len(chart_word) == distance, "joint route lost exact Q distance")
        records.append((
            face, states[start], distance,
            tuple(chart_word), tuple(coordinate_word), tuple(signed_word),
            tuple(chart_coordinate_word),
        ))
    return tuple(records)


small_routes = compile_face(
    0, small_states, small_count_vectors, small_state_index, small_reset_index,
    small_reset_counts, small_nonreset, small_direction_masks, aware_depth_tree,
)
full_routes = compile_face(
    1, full_states, full_count_vectors, full_state_index, full_reset_index,
    full_reset_counts, full_nonreset, full_direction_masks, aware_depth_tree,
)
joint_routes = small_routes + full_routes
require(len(joint_routes) == 4558, "joint route census drift")


def class_count(position, tag_face):
    return len({
        ((record[0], record[position]) if tag_face else record[position])
        for record in joint_routes
    })


continuation_classes = tuple(
    (class_count(position, False), class_count(position, True))
    for position in (3, 4, 6, 5)
)
require(continuation_classes == (
    (313, 345), (4125, 4208), (4520, 4548), (4536, 4558),
), "joint continuation-class census drift")

signed_word_fibres = defaultdict(list)
for record in joint_routes:
    signed_word_fibres[record[5]].append((record[0], record[1]))
signed_word_collisions = tuple(sorted(
    tuple(fibre) for fibre in signed_word_fibres.values() if len(fibre) > 1
))
require(len(signed_word_collisions) == 22
        and all(len(fibre) == 2 for fibre in signed_word_collisions),
        "ordered signed-word collision census drift")


def displacement_set(counts, reset_counts, signed):
    if signed:
        return {
            tuple(target - value for value, target in zip(vector, reset_counts))
            for vector in counts
        }
    return {
        tuple(abs(target - value) for value, target in zip(vector, reset_counts))
        for vector in counts
    }


small_signed = displacement_set(small_count_vectors, small_reset_counts, True)
full_signed = displacement_set(full_count_vectors, full_reset_counts, True)
small_unsigned = displacement_set(small_count_vectors, small_reset_counts, False)
full_unsigned = displacement_set(full_count_vectors, full_reset_counts, False)
signed_parikh_intersection = small_signed & full_signed
unsigned_parikh_intersection = small_unsigned & full_unsigned
require((len(small_signed), len(full_signed), len(signed_parikh_intersection),
         len(small_signed | full_signed)) == (239, 4319, 80, 4478),
        "signed displacement collision tax drift")
require((len(small_unsigned), len(full_unsigned),
         len(unsigned_parikh_intersection), len(small_unsigned | full_unsigned))
        == (96, 1536, 96, 1536), "unsigned displacement quotient drift")


def multiply_polynomials(left, right):
    answer = [0] * (len(left) + len(right) - 1)
    for left_degree, left_coefficient in enumerate(left):
        for right_degree, right_coefficient in enumerate(right):
            answer[left_degree + right_degree] += left_coefficient * right_coefficient
    return answer


def distance_data(capacities, reset_counts):
    factors = []
    polynomial = [1]
    for capacity, target in zip(capacities, reset_counts):
        counter = Counter(abs(value - target) for value in range(capacity + 1))
        factor = tuple(counter.get(degree, 0) for degree in range(max(counter) + 1))
        factors.append(factor)
        polynomial = multiply_polynomials(polynomial, factor)
    empty_distance = sum(reset_counts)
    polynomial[empty_distance] -= 1
    return tuple(factors), tuple(polynomial)


small_factors, small_distance = distance_data(
    (4, 3, 2, 1, 1, 0, 0, 0), small_reset_counts
)
full_factors, full_distance = distance_data(full["CAPACITIES"], full_reset_counts)
joint_distance = tuple(
    (small_distance[degree] if degree < len(small_distance) else 0)
    + (full_distance[degree] if degree < len(full_distance) else 0)
    for degree in range(max(len(small_distance), len(full_distance)))
)
require(sum(small_distance) == 239 and sum(full_distance) == 4319
        and sum(joint_distance) == 4558, "distance enumerator mass drift")
require(small_distance == (1, 8, 27, 51, 61, 50, 28, 11, 2),
        "small distance sequence drift")
require(full_distance == (
    1, 11, 55, 169, 365, 598, 775, 810, 685, 467, 250, 101, 28, 4,
), "full distance sequence drift")
require(joint_distance == (
    2, 19, 82, 220, 426, 648, 803, 821, 687, 467, 250, 101, 28, 4,
), "joint distance sequence drift")
require(Counter(record[2] for record in joint_routes) == Counter({
    degree: coefficient for degree, coefficient in enumerate(joint_distance)
}), "joint distance enumerator disagrees with routes")

response_by_face_distance = tuple(
    (face, distance,
     sum(record[0] == face and record[2] == distance
         and (record[3][0] if record[3] else 0) == 2 for record in joint_routes),
     sum(record[0] == face and record[2] == distance
         and (record[3][0] if record[3] else 0) == 10 for record in joint_routes))
    for face, maximum in ((0, len(small_distance) - 1), (1, len(full_distance) - 1))
    for distance in range(1, maximum + 1)
)

switch_distance_enumerator = tuple(sorted(Counter(
    (record[2], sum(left != right for left, right in zip(record[3], record[3][1:])))
    for record in joint_routes
).items()))

depth_tree_digest = hashlib.sha256(
    repr(aware_depth_tree).encode("ascii")
).hexdigest()
leaf_tree_digest = hashlib.sha256(
    repr(aware_leaf_tree).encode("ascii")
).hexdigest()
conflict_digest = hashlib.sha256(
    repr(blind_conflicts).encode("ascii")
).hexdigest()
signed_collision_digest = hashlib.sha256(
    repr(signed_word_collisions).encode("ascii")
).hexdigest()
route_digest = hashlib.sha256(repr(tuple(
    (record[0], record[1], record[3], record[4], record[5])
    for record in joint_routes
)).encode("ascii")).hexdigest()
require(conflict_digest
        == "7bd246a6267514657b32af3cc4993fd2516a8ad577a71f1a1a348e8ac2966d45",
        "cross-face conflict digest drift")
require(depth_tree_digest
        == "61eac5c65580bc41cd2a62a4d129d25ba1990fa7ce0d3479084c9e35362b7c10",
        "minimum-depth tree digest drift")
require(leaf_tree_digest
        == "94e83bfeece6fde6c0c5182479644560caca287de3b7257542a52ad556abd3d0",
        "minimum-leaf tree digest drift")
require(signed_collision_digest
        == "8bac417122ef8df947f856212e263724dcb8989217cad1f332d041d879f63624",
        "ordered signed-word collision digest drift")
require(route_digest
        == "091cd79646808e5eee2dc96c0722d337fcbdc2b1ad728f2d09fdbac51831fa67",
        "joint route digest drift")

print("FC3 joint two-face row-(2,10) policy/sidecar beach audit")
print("dependency_hash_checks=8")
print("faces=(support_(1,2)_I2:239,support_(1,3)_I2:4319),joint=4558")
print("exclusive_records=(small=71,full=407,joint=478)")
print("primitive_tests=(small_axis=11,full_axis=16,joint_axis=16,joint_plus_face=17)")
print("small_axis_depth_profile=" + repr(tuple(
    None if value == small_optimizer.infinity else value for value in small_profile
)))
print("small_axis_minimum_leaves=" + repr(small_minimum_leaves))
print("small_axis_optimal_tree=" + repr(small_tree))
print("full_axis_depth_profile=" + repr(tuple(
    None if value == full_optimizer.infinity else value for value in full_profile
)))
print("full_axis_minimum_leaves=" + repr(full_minimum_leaves))
print("face_blind_conflict_count=" + repr(len(blind_conflicts)))
print("face_blind_conflicts=" + repr(blind_conflicts))
print("face_blind_conflict_sha256=" + conflict_digest)
print("minimal_face_blind_witness=state_(5):small_requires_row2,full_requires_row10")
print("face_blind_axis_policy=IMPOSSIBLE_even_for_arbitrary_count_vector_function")
print("raw_count_vector_collision_tax=4558-4319=239_cross_face_pairs")
print("face_aware_depth_profile=" + repr(tuple(
    None if value == aware_optimizer.infinity else value for value in aware_profile
)))
print("face_aware_minimum_leaves=" + repr(aware_minimum_leaves))
print("face_aware_minimum_depth_tree=" + repr(aware_depth_tree))
print("face_aware_minimum_depth_tree_stats=(depth,leaves,face_nodes)="
      + repr(tree_statistics(aware_depth_tree)))
print("face_aware_minimum_depth_tree_sha256=" + depth_tree_digest)
print("face_aware_minimum_leaf_tree=" + repr(aware_leaf_tree))
print("face_aware_minimum_leaf_tree_stats=(depth,leaves,face_nodes)="
      + repr(tree_statistics(aware_leaf_tree)))
print("face_aware_minimum_leaf_tree_sha256=" + leaf_tree_digest)
print("gate_first_face_bit=(depth=%d,leaves=%d)" % (
    1 + max(small_depth, full_depth), small_minimum_leaves + full_minimum_leaves
))
print("conditional_face_bit_consults=(small=%d/238,full=%d/4318)" % (
    small_face_consults, full_face_consults
))
print("minimum_leaf_tree_face_bit_consults=(small=%d/238,full=%d/4318)" % (
    small_leaf_face_consults, full_leaf_face_consults
))
print("joint_assignments=(small_row2=%d,small_row10=%d,full_row2=%d,full_row10=%d)" % (
    small_assignment[2], small_assignment[10], full_assignment[2], full_assignment[10]
))
print("minimal_selector_sidecar=one_persistent_face_bit; zero_bits_impossible")
print("uniformly_initialized_previous_chart_bit=IMPOSSIBLE_at_state_(5)")
print("face_initialized_memory_bit=EQUIVALENT_to_persistent_face_bit_and_SUFFICIENT")
print("continuation_classes=(chart,unsigned_coordinate,chart_plus_coordinate,signed_edit)")
print("continuation_classes_(untagged,face_tagged)=" + repr(continuation_classes))
print("signed_word_collision_count=" + repr(len(signed_word_collisions)))
print("signed_word_collisions=" + repr(signed_word_collisions))
print("signed_word_collision_sha256=" + signed_collision_digest)
print("parikh_classes=(unsigned_untagged=%d,unsigned_face_tagged=%d,signed_untagged=%d,signed_face_tagged=4558)" % (
    len(small_unsigned | full_unsigned), len(small_unsigned) + len(full_unsigned),
    len(small_signed | full_signed)
))
print("parikh_cross_face_intersections=(unsigned=%d,signed=%d)" % (
    len(unsigned_parikh_intersection), len(signed_parikh_intersection)
))
print("small_distance_factors=" + repr(small_factors))
print("small_distance_coefficients=" + repr(small_distance))
print("full_distance_factors=" + repr(full_factors))
print("full_distance_coefficients=" + repr(full_distance))
print("joint_distance_coefficients=" + repr(joint_distance))
print("joint_distance_closed_form=(1+2z+z^2+z^3)(1+z)^4(1+2z)[1+(1+z^2)(1+z+z^2)(1+2z)]-z^6-z^8")
print("face_catalytic_enumerator=u*D12(z)+v*D13(z)")
print("signed_parikh_collision_tax=4558-4478=80_cross_face_pairs")
print("initial_response_by_face_distance=" + repr(response_by_face_distance))
print("joint_(distance,switch)_enumerator=" + repr(switch_distance_enumerator))
print("joint_route_sha256=" + route_digest)
print("typed_scope=two promoted I2 faces; row selector and commutative route responses only; no FC3 closure")
print("all_exact_checks=PASS")
