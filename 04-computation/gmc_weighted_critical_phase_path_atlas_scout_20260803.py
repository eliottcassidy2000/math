#!/usr/bin/env python3
"""Exact hostile audit and weighted critical-phase path scout.

This companion independently rebuilds the THM-3269 clutch weights from the
upstream response formulas, enumerates every spanning tree rather than using
the theorem companion's half-atlas enumeration, and then combines the unique
weighted chart with THM-3273's intrinsic order-twelve critical quotient.
"""

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path

from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
THM3238_SCRIPT = (
    ROOT / "04-computation/gmc_complete_physical_bank_unique_reset_thm3238.py"
)
THM3254_SCRIPT = ROOT / "04-computation/gmc_first_shell_pair_clutch_thm3254.py"
DEPENDENCIES = {
    THM3238_SCRIPT:
        "201e7348cc4f1e7fe4cfd51cfda42db85b8943d8d33f2d9080f20df562ecccaa",
    ROOT / "05-knowledge/results/gmc_complete_physical_bank_unique_reset_thm3238.out":
        "77b6a45b1715e9412732e3e89103809071eab4e3225f95510b7b59b022ddc93b",
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
    ROOT / "04-computation/gmc_scale_invariant_weighted_bispanning_thm3269.py":
        "41ae9aeb01fea1384f59f3a2687b1a0482954bf202e5be2b6fc928ef579b116a",
    ROOT / "05-knowledge/results/gmc_scale_invariant_weighted_bispanning_thm3269.out":
        "65e6e1bb04f3d42d64c2a8e5322c0f3b37c9d05ed6b053670a2cd46293742e3c",
    ROOT / "01-canon/theorems/THM-3273-critical-group-c12-quotient-and-relative-c7-equivariance-boundary.md":
        "00bd27bbc69f58f4fae7d76bc0147a69660098f5bf8718b22913099336fb17b2",
    ROOT / "04-computation/gmc_critical_group_phase_carrier_thm3273.py":
        "ab50ab5f2c5ea0fd9f7c25504f38bc3982fbfa9cc8fa6892b9ee3dfdf13b738b",
    ROOT / "05-knowledge/results/gmc_critical_group_phase_carrier_thm3273.out":
        "da1d102a73ae19b3d28b54dbe6a8243b1ef587c0782bda7172fb17e6f019d7d0",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected in DEPENDENCIES.items():
    require(hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected,
            ("dependency hash drift", dependency.name))

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0, "optimization-sensitive assert")
require(float_literals == 0, "floating literal")


def literal_assignment(module, name):
    for node in module.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name) and target.id == name
                        for target in node.targets)):
            return ast.literal_eval(node.value)
    raise RuntimeError(("missing literal assignment", name))


# Load only the pre-main exact definitions of THM-3238.  This rebuilds the
# response rows from coefficient formulas; no THM-3254 or THM-3269 interval,
# cover, graph, weight, minimizer, or chart assertion is imported.
upstream_tree = ast.parse(THM3238_SCRIPT.read_text(encoding="utf-8"))
upstream_prefix = []
for upstream_node in upstream_tree.body:
    if isinstance(upstream_node, ast.FunctionDef) and upstream_node.name == "main":
        break
    upstream_prefix.append(upstream_node)
upstream_namespace = {
    "__file__": str(THM3238_SCRIPT),
    "__name__": "hostile_thm3238_exact_prefix",
}
upstream_module = ast.fix_missing_locations(
    ast.Module(body=upstream_prefix, type_ignores=[])
)
exec(compile(upstream_module, str(THM3238_SCRIPT), "exec"), upstream_namespace)

reset = upstream_namespace["RESET"]
values = upstream_namespace["VALUES"]
multiplicity = upstream_namespace["MULTIPLICITY"]
certificate = upstream_namespace["CERTIFICATE"]
upsets = upstream_namespace["UPSETS"]
bank = upstream_namespace["BANK"]
coefficient_vectors = upstream_namespace["coefficient_vectors"]
states = upstream_namespace["STATES"]
reset_counts = Counter(reset)
require(reset == (1, 3, 3, 4, 5, 6, 7, 8), "reset drift")


def reset_distance(state):
    counts = Counter(state)
    return sum(abs(counts[value] - reset_counts[value]) for value in values)


link_states = tuple(sorted(state for state in states if reset_distance(state) == 1))
require(len(link_states) == 11, "independent link census")
require(all(bool(state) and all(Counter(state)[value] <= multiplicity[value]
                                for value in values)
            for state in link_states), "nonphysical link state")


def response_values(state):
    vector = coefficient_vectors(1, bank, 1, 3, state)
    exact_values = tuple(
        sum((vector[degree][shape] for shape in upset), Fraction(0))
        for (degree, _, _, _), upset in zip(certificate, upsets)
    )
    require(all(value.denominator == 1 for value in exact_values),
            ("nonintegral reconstructed response", state))
    return tuple(value.numerator for value in exact_values)


responses = {state: response_values(state) for state in link_states}
require(response_values(reset) == (0,) * 22, "reset response not zero")

source_tree = ast.parse(THM3254_SCRIPT.read_text(encoding="utf-8"))
covering_pairs = tuple(literal_assignment(source_tree, "ROW_COVERING_PAIRS"))
delayed_pairs = set(literal_assignment(source_tree, "DELAYED_WITNESSES"))
reset_link_edges = tuple(edge for edge in covering_pairs if edge not in delayed_pairs)
core_edges = tuple(edge for edge in reset_link_edges if 14 not in edge)
core_vertices = tuple(sorted({vertex for edge in core_edges for vertex in edge}))
require((len(covering_pairs), len(delayed_pairs), len(reset_link_edges),
         len(core_edges), len(core_vertices)) == (31, 8, 23, 22, 12),
        "independent graph bank census")


def trap_interval_direct(state, left_row, right_row,
                         left_scale=Fraction(1), right_scale=Fraction(1)):
    """Solve the one reset-edge inequality directly."""
    source = responses[state]
    delta_left = -left_scale * source[left_row - 1]
    delta_right = -right_scale * source[right_row - 1]
    if delta_left > 0:
        bound = -delta_right / delta_left
        return None if bound < 0 else (Fraction(0), bound)
    if delta_left < 0:
        bound = max(Fraction(0), -delta_right / delta_left)
        return bound, None
    return None if delta_right > 0 else (Fraction(0), None)


def strength_and_witnesses(pair, left_scale=Fraction(1),
                           right_scale=Fraction(1)):
    intervals = tuple(
        (state, interval)
        for state in link_states
        if (interval := trap_interval_direct(
            state, *pair, left_scale, right_scale)) is not None
    )
    require(not any(lower == 0 and upper is None
                    for _, (lower, upper) in intervals),
            ("one-state trap", pair))
    left_endpoints = tuple(
        (state, upper) for state, (lower, upper) in intervals
        if lower == 0 and upper is not None
    )
    right_endpoints = tuple(
        (state, lower) for state, (lower, upper) in intervals
        if upper is None and lower > 0
    )
    require(left_endpoints and right_endpoints, ("missing half-line", pair))
    candidates = tuple(
        (upper / lower, left_state, right_state)
        for left_state, upper in left_endpoints
        for right_state, lower in right_endpoints
    )
    strength = max(item[0] for item in candidates)
    witnesses = tuple(sorted(
        (left_state, right_state)
        for value, left_state, right_state in candidates if value == strength
    ))
    require(strength >= 1, ("uncovered edge", pair))
    return strength, witnesses


strength_bank = {edge: strength_and_witnesses(edge) for edge in reset_link_edges}
core_strengths = tuple(strength_bank[edge][0] for edge in core_edges)
require(all(value > 1 for value in core_strengths),
        "core strength is not strictly greater than one")
require(all(strength_and_witnesses((right, left))[0] == strength
            for (left, right), (strength, _) in strength_bank.items()),
        "independent orientation check")

# A nontrivial exact rescaling checks the endpoint-covariance implementation,
# separately from the algebraic proof that both endpoints acquire beta/alpha.
alpha = Fraction(7, 11)
beta = Fraction(13, 17)
require(all(strength_and_witnesses(edge, alpha, beta)[0] == strength
            for edge, (strength, _) in strength_bank.items()),
        "positive projective row-scale check")

strength_digest = hashlib.sha256(repr(tuple(
    (edge, strength_bank[edge]) for edge in core_edges
)).encode("ascii")).hexdigest()
require(strength_digest
        == "c7b8441f54bd893f33f1767db913da0d40f392f20d74f6d177782a0490dc911a",
        "independent strength bank digest")
require((min(core_strengths), max(core_strengths)) == (
    Fraction(12654135326210, 11849988963879),
    Fraction(66479967673052544, 306432210196733),
), "independent strength range")


def is_tree(indices):
    if len(indices) != len(core_vertices) - 1:
        return False
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
    return True


def tree_weight(indices):
    answer = Fraction(1)
    for index in indices:
        answer *= core_strengths[index]
    return answer


# Hostile enumeration: inspect all C(22,11) subsets, not one representative
# from each complementary pair.  This independently recovers all 74,748
# spanning trees and then filters the ordered bispanning charts.
spanning_trees = []
for chosen in combinations(range(len(core_edges)), len(core_vertices) - 1):
    if is_tree(chosen):
        spanning_trees.append((tree_weight(chosen), chosen))
require(len(spanning_trees) == 74748, "all spanning tree census")
spanning_lookup = {indices for _, indices in spanning_trees}
all_indices = frozenset(range(len(core_edges)))
ordered_charts = sorted(
    (weight, indices)
    for weight, indices in spanning_trees
    if tuple(sorted(all_indices - set(indices))) in spanning_lookup
)
require(len(ordered_charts) == 9920, "ordered bispanning census")

best_weight, best_indices = ordered_charts[0]
require(sum(weight == best_weight for weight, _ in ordered_charts) == 1,
        "independent minimizer tie")
second_weight = next(weight for weight, _ in ordered_charts
                     if weight > best_weight)
second_ratio = second_weight / best_weight
best_edges = tuple(core_edges[index] for index in best_indices)
complement_indices = tuple(sorted(all_indices - set(best_indices)))
complement_edges = tuple(core_edges[index] for index in complement_indices)
expected_best_edges = (
    (2, 13), (2, 17), (2, 19), (2, 21), (3, 11), (3, 16),
    (7, 22), (10, 11), (16, 17), (18, 19), (21, 22),
)
expected_complement_edges = (
    (2, 10), (3, 22), (7, 18), (10, 16), (10, 22), (11, 13),
    (13, 18), (13, 22), (16, 21), (17, 22), (19, 22),
)
require(best_edges == expected_best_edges, "independent best tree")
require(complement_edges == expected_complement_edges,
        "independent best cotree")
require(second_ratio == Fraction(
    111771229679330735594490655,
    106891942612503954628585203,
), "independent second gap")

swap = {17: 21, 21: 17}
edge_index = {frozenset(edge): index for index, edge in enumerate(core_edges)}
image_indices = tuple(sorted(
    edge_index[frozenset((swap.get(left, left), swap.get(right, right)))]
    for left, right in best_edges
))
image_weight = tree_weight(image_indices)
require(image_weight / best_weight == Fraction(
    719636186849558063, 684759823142197440),
    "independent involution separation")
require(sum(weight < image_weight for weight, _ in ordered_charts) == 3,
        "involution chart rank")

best_graph = {vertex: set() for vertex in core_vertices}
for left, right in best_edges:
    best_graph[left].add(right)
    best_graph[right].add(left)


def distances(graph, start):
    answer = {start: 0}
    queue = [start]
    for vertex in queue:
        for neighbor in graph[vertex]:
            if neighbor not in answer:
                answer[neighbor] = answer[vertex] + 1
                queue.append(neighbor)
    return answer


eccentricity = {
    vertex: max(distances(best_graph, vertex).values())
    for vertex in core_vertices
}
radius = min(eccentricity.values())
diameter = max(eccentricity.values())
centers = tuple(vertex for vertex in core_vertices
                if eccentricity[vertex] == radius)
leaves = tuple(vertex for vertex in core_vertices
               if len(best_graph[vertex]) == 1)
require((radius, diameter, centers, leaves) == (4, 8, (17,), (7, 10, 13, 18)),
        "independent weak-tree geometry")

degree_classes = {}
for vertex in core_vertices:
    degree_classes.setdefault(len(best_graph[vertex]), []).append(vertex)
classes = tuple(tuple(sorted(group)) for _, group in sorted(degree_classes.items()))
best_edge_set = {frozenset(edge) for edge in best_edges}
automorphism_count = 0
for images in product(*(permutations(group) for group in classes)):
    mapping = {}
    for group, image_group in zip(classes, images):
        mapping.update(zip(group, image_group))
    image_edges = {
        frozenset((mapping[left], mapping[right])) for left, right in best_edges
    }
    automorphism_count += image_edges == best_edge_set
require(automorphism_count == 1, "independent weak-tree rigidity")


def reduced_incidence(indices, root_vertex):
    row_vertices = tuple(vertex for vertex in core_vertices
                         if vertex != root_vertex)
    row_index = {vertex: index for index, vertex in enumerate(row_vertices)}
    matrix = Matrix.zeros(len(row_vertices), len(indices))
    for column, edge_index_value in enumerate(indices):
        left, right = core_edges[edge_index_value]
        if left != root_vertex:
            matrix[row_index[left], column] = -1
        if right != root_vertex:
            matrix[row_index[right], column] = 1
    return matrix


root_vertex = centers[0]
incidence_best = reduced_incidence(best_indices, root_vertex)
incidence_complement = reduced_incidence(complement_indices, root_vertex)
incidence_determinants = (
    int(incidence_best.det()), int(incidence_complement.det())
)
require(incidence_determinants == (1, -1), "independent incidence signs")
transition_sympy = incidence_best.inv() * incidence_complement
require(all(getattr(value, "q", 1) == 1 for value in transition_sympy),
        "nonintegral transition")
transition = [[int(transition_sympy[row, column])
               for column in range(transition_sympy.cols)]
              for row in range(transition_sympy.rows)]
transition_determinant = int(Matrix(transition).det())
require(transition_determinant == -1, "independent transition determinant")
transition_digest = hashlib.sha256("\n".join(
    ",".join(map(str, row)) for row in transition
).encode("ascii")).hexdigest()
require(transition_digest
        == "93f9ae1643ad9606f15826d8f1d3809fd0b04ce1d74c11d92d354edbdd737d4f",
        "independent transition digest")


# Rebuild the full critical group rooted at the newly canonical center 17.
# A primitive adjugate row supplies an arbitrary coordinate.  The unique
# incident primitive direction then removes that unit choice, independently
# of the post-audit THM-3269 companion.
nonroot_vertices = tuple(vertex for vertex in core_vertices
                         if vertex != root_vertex)
nonroot_index = {vertex: index for index, vertex in enumerate(nonroot_vertices)}
laplacian = Matrix.zeros(len(nonroot_vertices))
for left, right in core_edges:
    if left != root_vertex:
        laplacian[nonroot_index[left], nonroot_index[left]] += 1
    if right != root_vertex:
        laplacian[nonroot_index[right], nonroot_index[right]] += 1
    if left != root_vertex and right != root_vertex:
        laplacian[nonroot_index[left], nonroot_index[right]] -= 1
        laplacian[nonroot_index[right], nonroot_index[left]] -= 1
critical_smith = smith_normal_form(laplacian, domain=ZZ)
critical_diagonal = tuple(abs(int(critical_smith[index, index]))
                          for index in range(critical_smith.rows))
require(critical_diagonal == (1,) * 10 + (74748,),
        "independent cyclic critical group")
critical_order = int(laplacian.det())
adjugate = laplacian.adjugate()

normalized_full_charts = []
primitive_vertices = []
for row, vertex in enumerate(nonroot_vertices):
    adjugate_row = tuple(int(adjugate[row, column])
                         for column in range(adjugate.cols))
    if reduce(gcd, (critical_order,) + adjugate_row) != 1:
        continue
    primitive_vertices.append(vertex)
    chart = {root_vertex: 0}
    chart.update({candidate: adjugate_row[column] % critical_order
                  for column, candidate in enumerate(nonroot_vertices)})
    require(best_graph[root_vertex] == {2, 16}, "canonical root neighbors")
    require(gcd(chart[2], critical_order) == 2,
            ("nonprimitive incident class", vertex, chart[2]))
    require(gcd(chart[16], critical_order) == 1,
            ("primitive incident class", vertex, chart[16]))
    generator_inverse = pow(chart[16], -1, critical_order)
    normalized_full_charts.append({
        candidate: (generator_inverse * value) % critical_order
        for candidate, value in chart.items()
    })

require(tuple(primitive_vertices) == (7, 16, 18, 22),
        "primitive adjugate-row census")
require(all(chart == normalized_full_charts[0]
            for chart in normalized_full_charts),
        "normalized full critical chart depends on primitive row")
full_exponent = normalized_full_charts[0]
expected_full_exponent = {
    2: 53266, 3: 60022, 7: 48259, 10: 39646,
    11: 9088, 13: 2344, 16: 1, 17: 0,
    18: 289, 19: 25012, 21: 49832, 22: 21481,
}
require(full_exponent == expected_full_exponent,
        "canonical full critical exponents")
require(len(set(full_exponent.values())) == len(core_vertices),
        "full critical vertex collision")
critical_vertex_order = tuple(sorted(core_vertices,
                                     key=full_exponent.__getitem__))
require(critical_vertex_order
        == (17, 16, 18, 13, 11, 22, 19, 10, 7, 21, 2, 3),
        "canonical full critical order")

# Use the full-Jacobian generator's image in J12 as the primary phase gauge.
# The quotient's singleton generator fibre independently selects [7-17],
# which is exactly seven times the full-Jacobian generator.
phase = {vertex: exponent % 12
         for vertex, exponent in full_exponent.items()}
expected_phase = {
    2: 10, 3: 10, 7: 7, 10: 10, 11: 4, 13: 4,
    16: 1, 17: 0, 18: 1, 19: 4, 21: 8, 22: 1,
}
require(phase == expected_phase, "full-generator J12 phase chart")
phase_counts = Counter(phase.values())
require(dict(sorted(phase_counts.items()))
        == {0: 1, 1: 3, 4: 3, 7: 1, 8: 1, 10: 3},
        "canonical rooted phase fibres")
require({value: phase_counts[value] for value in phase_counts
         if gcd(value, 12) == 1} == {1: 3, 7: 1},
        "generator-valued fibre asymmetry")
require(tuple(vertex for vertex in core_vertices if phase[vertex] == 7) == (7,),
        "singleton quotient generator vertex")
singleton_normalized_phase = {
    vertex: (7 * value) % 12 for vertex, value in phase.items()
}
require(singleton_normalized_phase == {
    2: 10, 3: 10, 7: 1, 10: 10, 11: 4, 13: 4,
    16: 7, 17: 0, 18: 7, 19: 4, 21: 8, 22: 7,
}, "singleton quotient generator normalization")


def oriented_target_set(edges):
    answer = set()
    for left, right in edges:
        answer.add((phase[right] - phase[left]) % 12)
        answer.add((phase[left] - phase[right]) % 12)
    return answer


edge_targets = oriented_target_set(core_edges)
tree_targets = oriented_target_set(best_edges)
cotree_targets = oriented_target_set(complement_edges)
require(edge_targets == set(range(12)) - {4, 8}, "oriented edge phase gap")
require(tree_targets == {1, 2, 3, 5, 6, 7, 9, 10, 11},
        "weak-tree edge target set")
require(cotree_targets == {0, 1, 3, 5, 6, 7, 9, 11},
        "strong-cotree edge target set")

weight_lookup = {frozenset(edge): strength_bank[edge][0] for edge in core_edges}


def enumerate_simple_paths(edges):
    graph = {vertex: [] for vertex in core_vertices}
    for left, right in edges:
        graph[left].append(right)
        graph[right].append(left)
    phase_records = {target: [] for target in range(12)}
    pair_best = {}
    path_count = 0
    for start in core_vertices:
        stack = [(start, (start,), Fraction(1))]
        while stack:
            vertex, path, cost = stack.pop()
            if len(path) > 1:
                path_count += 1
                target = (phase[vertex] - phase[start]) % 12
                phase_records[target].append((cost, path))
                pair = (start, vertex)
                if pair not in pair_best or cost < pair_best[pair][0]:
                    pair_best[pair] = (cost, path)
            for neighbor in graph[vertex]:
                if neighbor not in path:
                    edge = frozenset((vertex, neighbor))
                    stack.append((neighbor, path + (neighbor,),
                                  cost * weight_lookup[edge]))
    return phase_records, pair_best, path_count


def phase_minima(records):
    answer = {}
    for target, target_records in records.items():
        minimum = min(cost for cost, _ in target_records)
        paths = tuple(sorted(path for cost, path in target_records
                             if cost == minimum))
        answer[target] = minimum, paths
    return answer


full_records, full_pair_best, full_path_count = enumerate_simple_paths(core_edges)
tree_records, tree_pair_best, tree_path_count = enumerate_simple_paths(best_edges)
cotree_records, _, cotree_path_count = enumerate_simple_paths(complement_edges)
require((full_path_count, tree_path_count, cotree_path_count) == (21226, 132, 132),
        "simple path census")
full_minima = phase_minima(full_records)
tree_minima = phase_minima(tree_records)
cotree_minima = phase_minima(cotree_records)

expected_paths = {
    0: ((3, 11, 10), (10, 11, 3)),
    1: ((7, 22, 21),),
    2: ((21, 2),),
    3: ((18, 19),),
    4: ((19, 2, 21),),
    5: ((21, 22),),
    6: ((7, 22), (22, 7)),
    7: ((22, 21),),
    8: ((21, 2, 19),),
    9: ((19, 18),),
    10: ((2, 21),),
    11: ((21, 22, 7),),
}
require(all(full_minima[target][1] == expected_paths[target]
            for target in range(12)), "weighted target path atlas")
require(all(tree_minima[target] == full_minima[target]
            for target in range(12)),
        "weak tree fails to carry a global target geodesic")
require(sum(cotree_minima[target] == full_minima[target]
            for target in range(12)) == 0,
        "strong cotree carries a global target geodesic")

hop_distances = tuple(
    min(len(path) - 1 for _, path in full_records[target])
    for target in range(12)
)
require(hop_distances == (1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1),
        "unweighted target distances")

path_cost_bank = tuple(
    (target, full_minima[target][0], full_minima[target][1])
    for target in range(12)
)
path_cost_digest = hashlib.sha256(repr(path_cost_bank).encode("ascii")).hexdigest()
require(path_cost_digest
        == "74ee0d0c8a965882b3a3f8bc6254eeba83109eebbbc57f69590f903dce6219f1",
        "weighted phase path bank digest")

# The target quotient is essential: the weak tree does not preserve most
# endpoint-specific weighted geodesics, even though it preserves one global
# minimizer for every phase difference.
pairwise_geodesics_preserved = sum(
    full_pair_best[pair][0] == tree_pair_best[pair][0]
    for pair in full_pair_best
)
require((len(full_pair_best), pairwise_geodesics_preserved) == (132, 48),
        "endpoint-geodesic loss census")

atlas_edge_set = set()
for _, paths in full_minima.values():
    for path in paths:
        for index in range(len(path) - 1):
            atlas_edge_set.add(frozenset((path[index], path[index + 1])))
atlas_edges = tuple(edge for edge in core_edges
                    if frozenset(edge) in atlas_edge_set)
expected_atlas_edges = (
    (2, 19), (2, 21), (3, 11), (7, 22),
    (10, 11), (18, 19), (21, 22),
)
require(atlas_edges == expected_atlas_edges, "phase-geodesic backbone")

atlas_graph = {vertex: set() for vertex in core_vertices}
for left, right in atlas_edges:
    atlas_graph[left].add(right)
    atlas_graph[right].add(left)
seen = set()
component_sizes = []
for start in core_vertices:
    if start in seen:
        continue
    seen.add(start)
    queue = [start]
    for vertex in queue:
        for neighbor in atlas_graph[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    component_sizes.append(len(queue))
require(tuple(sorted(component_sizes)) == (1, 1, 1, 3, 6),
        "phase backbone forest components")

atlas_indices = frozenset(edge_index[frozenset(edge)] for edge in atlas_edges)
ordered_tree_sets = tuple(frozenset(indices) for _, indices in ordered_charts)
phase_preserving_charts = tuple(
    tree for tree in ordered_tree_sets if atlas_indices <= tree
)
require(len(phase_preserving_charts) == 36,
        "phase-geodesic bispanning subatlas census")
require(sum(atlas_indices <= frozenset(indices)
            for _, indices in spanning_trees) == 205,
        "all spanning phase-geodesic carrier census")

subatlas_graph = {tree: set() for tree in phase_preserving_charts}
for left_position, left_tree in enumerate(phase_preserving_charts):
    for right_tree in phase_preserving_charts[left_position + 1:]:
        if len(left_tree - right_tree) == len(right_tree - left_tree) == 1:
            subatlas_graph[left_tree].add(right_tree)
            subatlas_graph[right_tree].add(left_tree)
subatlas_edge_count = sum(len(neighbors)
                          for neighbors in subatlas_graph.values()) // 2
subatlas_degrees = tuple(len(subatlas_graph[tree])
                         for tree in phase_preserving_charts)
require((subatlas_edge_count, min(subatlas_degrees), max(subatlas_degrees))
        == (120, 4, 10), "phase-preserving exchange graph census")

subatlas_seen = {phase_preserving_charts[0]}
subatlas_queue = [phase_preserving_charts[0]]
for tree in subatlas_queue:
    for neighbor in subatlas_graph[tree]:
        if neighbor not in subatlas_seen:
            subatlas_seen.add(neighbor)
            subatlas_queue.append(neighbor)
require(len(subatlas_seen) == 36, "phase-preserving subatlas disconnected")
subatlas_diameter = 0
for start in phase_preserving_charts:
    distance = {start: 0}
    queue = [start]
    for tree in queue:
        for neighbor in subatlas_graph[tree]:
            if neighbor not in distance:
                distance[neighbor] = distance[tree] + 1
                queue.append(neighbor)
    subatlas_diameter = max(subatlas_diameter, max(distance.values()))
require(subatlas_diameter == 5, "phase-preserving subatlas diameter")
best_tree_set = frozenset(best_indices)
require(len(subatlas_graph[best_tree_set]) == 8,
        "weighted chart phase-preserving exchange degree")

require(sum(weight < best_weight for weight, _ in spanning_trees) == 59,
        "weighted tree all-spanning rank")
require(sum(weight == best_weight for weight, _ in spanning_trees) == 1,
        "weighted tree all-spanning tie")
require(sum(tree_weight(tuple(sorted(tree))) == best_weight
            for tree in phase_preserving_charts) == 1,
        "phase subatlas weighted minimum tie")


def path_bank_text():
    entries = []
    for target, cost, paths in path_cost_bank:
        entries.append("%d:%s@%s" % (target, paths, cost))
    return ";".join(entries)


print("WEIGHTED CRITICAL-PHASE RESPONSE-PATH EXACT SCOUT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("independent_THM3269_response_reconstruction=(11 states,22 rows,all_integral)")
print("independent_clutch=(23 reset edges,22 core edges,orientation_and_scale_PASS)")
print("strength_range=(%s,%s)" % (min(core_strengths), max(core_strengths)))
print("strength_bank_sha256=%s" % strength_digest)
print("independent_spanning_trees=74748,bispanning=(4960 unordered,9920 ordered)")
print("independent_unique_weighted_chart=%s" % (best_edges,))
print("independent_complement=%s" % (complement_edges,))
print("independent_second_over_best=%s" % second_ratio)
print("independent_tree=(Aut1,center17,radius4,diameter8,leaves=(7,10,13,18))")
print("independent_incidence=(det%s,transition_det=%d,digest=%s)" %
      (incidence_determinants, transition_determinant, transition_digest))
print("canonical_full_Jac_root=17,generator=[16-17],vertex_order=%s" %
      (critical_vertex_order,))
print("secondary_J12_singleton_generator=[7-17]=7*[16-17]")
print("canonical_vertex_phase=%s" % (tuple(
    (vertex, phase[vertex]) for vertex in core_vertices
),))
print("oriented_edge_targets=all_except_4_8")
print("tree_edge_targets=%s,cotree_edge_targets=%s" %
      (tuple(sorted(tree_targets)), tuple(sorted(cotree_targets))))
print("unweighted_phase_hop_distances=%s" % (hop_distances,))
print("weighted_phase_path_bank=%s" % path_bank_text())
print("weighted_phase_path_bank_sha256=%s" % path_cost_digest)
print("phase_geodesic_backbone=%s,forest_components=(6,3,1,1,1)" %
      (atlas_edges,))
print("weak_tree_global_phase_geodesics=12/12,endpoint_geodesics=48/132")
print("strong_cotree_global_phase_geodesics=0/12")
print("phase_backbone_carriers=(205 all spanning,36 ordered_bispanning)")
print("phase_exchange_subatlas=(V36,E120,degree4..10,connected,diameter5,Tstar_degree8)")
print("Tstar_all_spanning_weight_rank=60,phase_subatlas_weight_minimum=unique")
print("loss=phase_target_forgets_endpoints_and_6229_part;product_forgets_witness_bank_and_is_not_response_composition")
print("scope=canonical_internal_response_phase_only,no_vertex_torsor,no_owner_phase,no_positive_current,no_GMC_or_LRC_decrement")
print("all_exact_checks=PASS")
