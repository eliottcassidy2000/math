#!/usr/bin/env python3
"""Exact companion for THM-3269's weighted bispanning polarization."""

import ast
import contextlib
import io
import runpy
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from pathlib import Path


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
THM3254_SCRIPT = ROOT / "04-computation/gmc_first_shell_pair_clutch_thm3254.py"
DEPENDENCIES = {
    ROOT / "01-canon/theorems/THM-3254-first-shell-two-row-clutch-and-graded-gauge-no-go.md":
        "c7c5948d5181dc845da8f40683517c06e3eecfe837e17cb84d356b91eae1a081",
    THM3254_SCRIPT:
        "05efd37eeedeca7e3be581977a894592a7873d94a966f06d9533482cc8498fee",
    ROOT / "05-knowledge/results/gmc_first_shell_pair_clutch_thm3254.out":
        "dd415c8ce6e2e196c115421d3508addabb724305a16843509264a8b3205beee9",
    ROOT / "01-canon/theorems/THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary.md":
        "754250ba866b2190a51540e9b4f75234079c6cd0a20bc23083d45d2a8dffe0c9",
    ROOT / "04-computation/gmc_bispanning_reset_link_holotopy_thm3260.py":
        "7adf81a692bf477d493860483f50b291c01dd267ed578caa12dd31bca9128616",
    ROOT / "05-knowledge/results/gmc_bispanning_reset_link_holotopy_thm3260.out":
        "eb0b2d2de4808ae1aabdede94640390aeb35cc7a7deec55b193da44dcd043e37",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected in DEPENDENCIES.items():
    require(sha256(lf_bytes(dependency)).hexdigest() == expected,
            ("dependency hash drift", dependency.name))

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0, "optimization-sensitive assert")
require(float_literals == 0, "floating literal")


# Execute the pinned THM-3254 exact companion to recover its directly rebuilt
# reset-link responses and trap-interval function.  Its transcript is hidden;
# none of its selected two-state covers or graph conclusions is trusted below.
sink = io.StringIO()
with contextlib.redirect_stdout(sink):
    clutch = runpy.run_path(str(THM3254_SCRIPT))

trap_interval = clutch["trap_interval"]
link_states = tuple(clutch["link_states"])
covering_pairs = tuple(clutch["ROW_COVERING_PAIRS"])
delayed_pairs = set(clutch["DELAYED_WITNESSES"])
require((len(link_states), len(covering_pairs), len(delayed_pairs)) == (11, 31, 8),
        "pinned clutch bank census")


def strength_and_witnesses(pair):
    """Canonical maximum overlap U/L for [0,U] and [L,infinity] traps."""
    intervals = tuple(
        (state, interval)
        for state in link_states
        if (interval := trap_interval(state, *pair)) is not None
    )
    require(not any(lower == 0 and upper is None
                    for _, (lower, upper) in intervals),
            ("one-state trap", pair))
    left = tuple((state, upper) for state, (lower, upper) in intervals
                 if lower == 0 and upper is not None)
    right = tuple((state, lower) for state, (lower, upper) in intervals
                  if upper is None and lower > 0)
    require(left and right, ("missing endpoint trap", pair))
    candidates = tuple(
        (Fraction(upper, lower), left_state, right_state)
        for left_state, upper in left
        for right_state, lower in right
    )
    strength = max(value for value, _, _ in candidates)
    witnesses = tuple(sorted(
        (left_state, right_state)
        for value, left_state, right_state in candidates if value == strength
    ))
    require(strength >= 1, ("uncovered reset-link pair", pair, strength))
    return strength, witnesses


reset_link_edges = tuple(pair for pair in covering_pairs
                         if pair not in delayed_pairs)
core_edges = tuple(edge for edge in reset_link_edges if 14 not in edge)
core_vertices = tuple(sorted({vertex for edge in core_edges for vertex in edge}))
require((len(reset_link_edges), len(core_edges), len(core_vertices)) == (23, 22, 12),
        "core graph census")

strength_bank = {
    edge: strength_and_witnesses(edge) for edge in reset_link_edges
}
require(all(strength_and_witnesses((right, left))[0] == strength
            for (left, right), (strength, _) in strength_bank.items()),
        "orientation invariance")
core_strengths = tuple(strength_bank[edge][0] for edge in core_edges)
strength_digest = sha256(repr(tuple(
    (edge, strength_bank[edge]) for edge in core_edges
)).encode("ascii")).hexdigest()


def is_tree(indices):
    parent = {vertex: vertex for vertex in core_vertices}

    def root(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for index in indices:
        left, right = core_edges[index]
        left_root, right_root = root(left), root(right)
        if left_root == right_root:
            return False
        parent[left_root] = right_root
    return len(indices) == len(core_vertices) - 1


# Edge zero chooses one representative of each unordered complementary pair.
all_indices = set(range(len(core_edges)))
decompositions = []
for chosen in combinations(range(len(core_edges)), len(core_vertices) - 1):
    if 0 not in chosen:
        continue
    complement = tuple(sorted(all_indices - set(chosen)))
    if is_tree(chosen) and is_tree(complement):
        decompositions.append((tuple(chosen), complement))
require(len(decompositions) == 4960, "bispanning decomposition census")


def tree_weight(indices):
    answer = Fraction(1)
    for index in indices:
        answer *= core_strengths[index]
    return answer


ordered_charts = []
for first, second in decompositions:
    ordered_charts.append((tree_weight(first), first))
    ordered_charts.append((tree_weight(second), second))
ordered_charts.sort()
require(len(ordered_charts) == 9920, "ordered chart census")
best_weight, best_indices = ordered_charts[0]
best_multiplicity = sum(weight == best_weight for weight, _ in ordered_charts)
require(best_multiplicity == 1, "weighted chart not unique")
second_weight = next(weight for weight, _ in ordered_charts
                     if weight > best_weight)
second_ratio = second_weight / best_weight

EXPECTED_BEST_EDGES = (
    (2, 13), (2, 17), (2, 19), (2, 21), (3, 11), (3, 16),
    (7, 22), (10, 11), (16, 17), (18, 19), (21, 22),
)
best_edges = tuple(core_edges[index] for index in best_indices)
require(best_edges == EXPECTED_BEST_EDGES, "canonical tree drift")
complement_indices = tuple(sorted(all_indices - set(best_indices)))
complement_edges = tuple(core_edges[index] for index in complement_indices)
require(is_tree(best_indices) and is_tree(complement_indices),
        "canonical chart lost bispanning")
require(second_ratio == Fraction(111771229679330735594490655,
                                 106891942612503954628585203),
        "second-best separation drift")


# The only intrinsic core automorphism swaps rows 17 and 21.  Exact clutch
# weights break it: the image chart is the fourth ordered chart by weight.
swap = {17: 21, 21: 17}
edge_lookup = {frozenset(edge): index for index, edge in enumerate(core_edges)}
image_indices = tuple(sorted(
    edge_lookup[frozenset((swap.get(left, left), swap.get(right, right)))]
    for left, right in best_edges
))
image_weight = tree_weight(image_indices)
image_ratio = image_weight / best_weight
require(image_ratio == Fraction(719636186849558063, 684759823142197440),
        "intrinsic-C2 separation drift")
require(sum(weight < image_weight for weight, _ in ordered_charts) == 3,
        "intrinsic-C2 image rank drift")
require(sum(weight == image_weight for weight, _ in ordered_charts) == 1,
        "intrinsic-C2 image weight tie")


# The selected weak tree is itself rigid and has a unique center.  This fixes
# a canonical root, though not a cyclic ordering of its twelve vertices.
best_edge_set = {frozenset(edge) for edge in best_edges}
best_graph = {vertex: set() for vertex in core_vertices}
for left, right in best_edges:
    best_graph[left].add(right)
    best_graph[right].add(left)


def distances_from(start):
    distances = {start: 0}
    queue = [start]
    for vertex in queue:
        for neighbor in best_graph[vertex]:
            if neighbor not in distances:
                distances[neighbor] = distances[vertex] + 1
                queue.append(neighbor)
    return distances


eccentricities = {
    vertex: max(distances_from(vertex).values()) for vertex in core_vertices
}
radius = min(eccentricities.values())
diameter = max(eccentricities.values())
centers = tuple(vertex for vertex in core_vertices
                if eccentricities[vertex] == radius)
leaves = tuple(vertex for vertex in core_vertices
               if len(best_graph[vertex]) == 1)
require((radius, diameter, centers, leaves)
        == (4, 8, (17,), (7, 10, 13, 18)),
        "canonical tree geometry")

degree_classes = {}
for vertex in core_vertices:
    degree_classes.setdefault(len(best_graph[vertex]), []).append(vertex)
classes = tuple(tuple(sorted(group)) for _, group
                in sorted(degree_classes.items()))
automorphism_count = 0
for images in product(*(permutations(group) for group in classes)):
    mapping = {}
    for group, image_group in zip(classes, images):
        mapping.update(zip(group, image_group))
    image_edges = {
        frozenset((mapping[left], mapping[right])) for left, right in best_edges
    }
    automorphism_count += image_edges == best_edge_set
require(automorphism_count == 1, "canonical tree not rigid")


def determinant_bareiss(matrix):
    work = [row[:] for row in matrix]
    previous = 1
    sign = 1
    for column in range(len(work) - 1):
        pivot_row = next((row for row in range(column, len(work))
                          if work[row][column]), None)
        require(pivot_row is not None, ("singular incidence", column))
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, len(work)):
            for index in range(column + 1, len(work)):
                numerator = (work[row][index] * pivot
                             - work[row][column] * work[column][index])
                require(numerator % previous == 0, "Bareiss divisibility")
                work[row][index] = numerator // previous
        previous = pivot
        for row in range(column + 1, len(work)):
            work[row][column] = 0
    return sign * work[-1][-1]


root_vertex = 17
row_vertices = tuple(vertex for vertex in core_vertices if vertex != root_vertex)


def incidence(indices):
    matrix = [[0] * len(indices) for _ in row_vertices]
    row_lookup = {vertex: row for row, vertex in enumerate(row_vertices)}
    for column, index in enumerate(indices):
        left, right = core_edges[index]
        if left in row_lookup:
            matrix[row_lookup[left]][column] = -1
        if right in row_lookup:
            matrix[row_lookup[right]][column] = 1
    return matrix


def solve(left, right):
    size = len(left)
    work = [[Fraction(value) for value in left[row] + right[row]]
            for row in range(size)]
    for column in range(size):
        pivot_row = next(row for row in range(column, size)
                         if work[row][column])
        work[column], work[pivot_row] = work[pivot_row], work[column]
        pivot = work[column][column]
        work[column] = [value / pivot for value in work[column]]
        for row in range(size):
            if row == column:
                continue
            multiplier = work[row][column]
            work[row] = [left_value - multiplier * right_value
                         for left_value, right_value
                         in zip(work[row], work[column])]
    answer = [row[size:] for row in work]
    require(all(value.denominator == 1 for row in answer for value in row),
            "nonintegral canonical transition")
    return [[value.numerator for value in row] for row in answer]


incidence_best = incidence(best_indices)
incidence_complement = incidence(complement_indices)
incidence_determinants = (
    determinant_bareiss(incidence_best),
    determinant_bareiss(incidence_complement),
)
require(all(abs(value) == 1 for value in incidence_determinants),
        "nonunimodular canonical chart")
transition = solve(incidence_best, incidence_complement)
transition_determinant = determinant_bareiss(transition)
require(abs(transition_determinant) == 1, "canonical transition not GL11")
transition_digest = sha256("\n".join(
    ",".join(map(str, row)) for row in transition
).encode("ascii")).hexdigest()


print("THM-3269 SCALE-INVARIANT WEIGHTED BISPANNING EXACT AUDIT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("reset_link=(11 states,23 blocked edges),core=(V12,E22)")
print("clutch_strength=max_U_over_L,orientation_and_positive_row_scale_invariant")
print("strength_range=(%s,%s)" % (min(core_strengths), max(core_strengths)))
print("strength_bank_sha256=%s" % strength_digest)
print("bispanning_charts=(4960 unordered,9920 ordered)")
print("unique_minimum_tree=%s" % (best_edges,))
print("canonical_complement=%s" % (complement_edges,))
print("second_over_best=%s" % second_ratio)
print("intrinsic_C2_image_over_best=%s,strict_ordered_rank=4" % image_ratio)
print("canonical_tree=(Aut1,center17,radius4,diameter8,leaves=(7,10,13,18))")
print("incidence_determinants=%s,transition_det=%d" %
      (incidence_determinants, transition_determinant))
print("canonical_transition_sha256=%s" % transition_digest)
print("conditional_bridge=external_C12_vertex_label_only;tree_pair_sidecar_removed")
print("scope=first_link_analytic_selector,no_C12_label,no_C7_carrier,no_LRC_or_GMC_decrement")
print("all_exact_checks=PASS")
