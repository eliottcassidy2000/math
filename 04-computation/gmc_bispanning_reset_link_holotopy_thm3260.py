#!/usr/bin/env python3
"""Exact companion for THM-3260's reset-link holotopy atlas."""

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from math import lcm
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
    ROOT / "01-canon/theorems/THM-3255-twelve-balance-multiplicative-singer-rank-defect-and-phase-marker-boundary.md":
        "8d827a3c1c93904db6715240373e8626323802240a4831d14394f3637389609d",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


for dependency, expected in DEPENDENCIES.items():
    require(sha256(lf_bytes(dependency)).hexdigest() == expected,
            ("dependency hash drift", dependency.name))

tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(tree)
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


# Read only the two literal graph-defining banks from the pinned THM-3254
# companion.  No theorem assertion or computed result is imported.
source_tree = ast.parse(THM3254_SCRIPT.read_text(encoding="utf-8"))
covering_pairs = tuple(literal_assignment(source_tree, "ROW_COVERING_PAIRS"))
delayed_witnesses = literal_assignment(source_tree, "DELAYED_WITNESSES")
delayed_pairs = set(delayed_witnesses)
reset_link_edges = tuple(pair for pair in covering_pairs
                         if pair not in delayed_pairs)
require((len(covering_pairs), len(delayed_pairs), len(reset_link_edges))
        == (31, 8, 23), "THM-3254 pair split")


def vertices(edges):
    return tuple(sorted({vertex for edge in edges for vertex in edge}))


def adjacency(edges, supplied_vertices=None):
    answer = {vertex: set() for vertex in
              (vertices(edges) if supplied_vertices is None
               else supplied_vertices)}
    for left, right in edges:
        answer[left].add(right)
        answer[right].add(left)
    return answer


def connected(edges, supplied_vertices=None):
    nodes = vertices(edges) if supplied_vertices is None else supplied_vertices
    if not nodes:
        return True
    graph = adjacency(edges, nodes)
    seen = {nodes[0]}
    queue = [nodes[0]]
    for vertex in queue:
        for neighbor in graph[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    return len(seen) == len(nodes)


all_vertices = vertices(reset_link_edges)
require(len(all_vertices) == 13 and connected(reset_link_edges),
        "connected reset-link graph")

bridges = []
for edge in reset_link_edges:
    remainder = list(reset_link_edges)
    remainder.remove(edge)
    if not connected(tuple(remainder), all_vertices):
        bridges.append(edge)
require(tuple(bridges) == ((13, 14),), "unique reset-link bridge")
require(adjacency(reset_link_edges)[14] == {13}, "row-14 leaf")

core_vertices = tuple(vertex for vertex in all_vertices if vertex != 14)
core_edges = tuple(edge for edge in reset_link_edges if 14 not in edge)
require((len(core_vertices), len(core_edges), connected(core_edges, core_vertices))
        == (12, 22, True), "leafless core census")
cycle_rank = len(core_edges) - len(core_vertices) + 1
require(cycle_rank == 11, "core Betti number")


# The core is bipartite, but not planar.  Verify an explicit subdivision of
# K_3,3; the nine paths are internally vertex-disjoint.
kuratowski_paths = (
    (10, 11), (10, 16), (10, 22),
    (3, 11), (3, 16), (3, 22),
    (18, 13, 11), (18, 19, 2, 21, 16), (18, 7, 22),
)
edge_set = {frozenset(edge) for edge in core_edges}
branch_left = {10, 3, 18}
branch_right = {11, 16, 22}
internal_sets = []
for path in kuratowski_paths:
    require(path[0] in branch_left and path[-1] in branch_right,
            ("K33 branch endpoint", path))
    require(all(frozenset((path[index], path[index + 1])) in edge_set
                for index in range(len(path) - 1)), ("K33 path edge", path))
    internal_sets.append(set(path[1:-1]))
require(all(not (internal_sets[left] & internal_sets[right])
            for left in range(len(internal_sets))
            for right in range(left + 1, len(internal_sets))),
        "K33 internally disjoint paths")
require(all(not (inside & (branch_left | branch_right))
            for inside in internal_sets), "K33 internal branch exclusion")

colors = {}
graph = adjacency(core_edges, core_vertices)
for start in core_vertices:
    if start in colors:
        continue
    colors[start] = 0
    queue = [start]
    for vertex in queue:
        for neighbor in graph[vertex]:
            if neighbor not in colors:
                colors[neighbor] = 1 - colors[vertex]
                queue.append(neighbor)
            require(colors[neighbor] != colors[vertex], "core bipartite")
color_parts = tuple(tuple(vertex for vertex in core_vertices
                          if colors[vertex] == color)
                    for color in (0, 1))
require(tuple(map(len, color_parts)) == (5, 7), "core bipartition census")


# Determine the intrinsic automorphism group exactly by permuting only within
# degree classes.  There are just 2*7! candidates after singleton classes.
degree_classes = {}
for vertex in core_vertices:
    degree_classes.setdefault(len(graph[vertex]), []).append(vertex)
classes = tuple(tuple(sorted(group)) for _, group in sorted(degree_classes.items()))
automorphisms = []
for images in product(*(permutations(group) for group in classes)):
    mapping = {}
    for group, image_group in zip(classes, images):
        mapping.update(zip(group, image_group))
    image_edges = {frozenset((mapping[left], mapping[right]))
                   for left, right in core_edges}
    if image_edges == edge_set:
        automorphisms.append(tuple(mapping[vertex] for vertex in core_vertices))
identity = core_vertices
swap_17_21 = tuple(21 if vertex == 17 else 17 if vertex == 21 else vertex
                   for vertex in core_vertices)
require(tuple(sorted(automorphisms)) == tuple(sorted((identity, swap_17_21))),
        ("core automorphism group", automorphisms))


# Adjoining the eight delayed edges gives a second, relative-homology layer.
# It adds precisely the two vertices 9 and 12.  The relative cellular boundary
# has rank two, so the six-dimensional kernel is H_1(H,G;Z).
full_vertices = vertices(covering_pairs)
require((len(full_vertices), len(covering_pairs), connected(covering_pairs))
        == (15, 31, True), "full covering graph census")
full_cycle_rank = len(covering_pairs) - len(full_vertices) + 1
require(full_cycle_rank == 17, "full covering graph Betti number")

full_bridges = []
for edge in covering_pairs:
    remainder = list(covering_pairs)
    remainder.remove(edge)
    if not connected(tuple(remainder), full_vertices):
        full_bridges.append(edge)
require(tuple(full_bridges) == ((3, 9),), "full covering graph bridge")

new_vertices = tuple(vertex for vertex in full_vertices
                     if vertex not in all_vertices)
delayed_edges = tuple(edge for edge in covering_pairs if edge in delayed_pairs)
require(new_vertices == (9, 12), "relative vertices")


def rational_rank(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    pivot_row = 0
    for column in range(columns):
        pivot = next((row for row in range(pivot_row, rows)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [value / scale for value in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right
                         for left, right in zip(work[row], work[pivot_row])]
        pivot_row += 1
    return pivot_row


relative_boundary = []
for vertex in new_vertices:
    relative_boundary.append(tuple(
        (-1 if left == vertex else 0) + (1 if right == vertex else 0)
        for left, right in delayed_edges
    ))
relative_boundary_rank = rational_rank(relative_boundary)
relative_cycle_rank = len(delayed_edges) - relative_boundary_rank
require((relative_boundary_rank, relative_cycle_rank,
         full_cycle_rank - cycle_rank) == (2, 6, 6),
        "relative homology rank")


# The uncoloured full graph has exactly S_3 x C_2 symmetry: S_3 permutes
# {10,17,21}, C_2 swaps {14,18}, and all other vertices are fixed.  Retaining
# the reset/delayed edge colours cuts this back to the core C_2.
full_graph = adjacency(covering_pairs, full_vertices)
full_edge_set = {frozenset(edge) for edge in covering_pairs}
full_degree_classes = {}
for vertex in full_vertices:
    full_degree_classes.setdefault(len(full_graph[vertex]), []).append(vertex)
full_classes = tuple(tuple(sorted(group)) for _, group
                     in sorted(full_degree_classes.items()))
full_automorphisms = []
for images in product(*(permutations(group) for group in full_classes)):
    mapping = {}
    for group, image_group in zip(full_classes, images):
        mapping.update(zip(group, image_group))
    image_edges = {frozenset((mapping[left], mapping[right]))
                   for left, right in covering_pairs}
    if image_edges == full_edge_set:
        full_automorphisms.append(tuple(mapping[vertex]
                                        for vertex in full_vertices))

expected_full_automorphisms = []
for triple_image in permutations((10, 17, 21)):
    for swap_pair in ((14, 18), (18, 14)):
        mapping = dict(zip((10, 17, 21), triple_image))
        mapping.update(zip((14, 18), swap_pair))
        expected_full_automorphisms.append(tuple(mapping.get(vertex, vertex)
                                                 for vertex in full_vertices))
require(tuple(sorted(full_automorphisms))
        == tuple(sorted(expected_full_automorphisms)),
        "full automorphism group S3xC2")


def permutation_order(image, ordered_vertices):
    mapping = dict(zip(ordered_vertices, image))
    seen = set()
    answer = 1
    for start in ordered_vertices:
        if start in seen:
            continue
        length = 0
        vertex = start
        while vertex not in seen:
            seen.add(vertex)
            length += 1
            vertex = mapping[vertex]
        answer = lcm(answer, length)
    return answer


full_automorphism_orders = tuple(sorted({
    permutation_order(image, full_vertices) for image in full_automorphisms
}))
require(full_automorphism_orders == (1, 2, 3, 6),
        "no order-7 or order-12 full automorphism")

reset_edge_set = {frozenset(edge) for edge in reset_link_edges}
coloured_automorphisms = []
for image in full_automorphisms:
    mapping = dict(zip(full_vertices, image))
    image_reset_edges = {frozenset((mapping[left], mapping[right]))
                         for left, right in reset_link_edges}
    if image_reset_edges == reset_edge_set:
        coloured_automorphisms.append(image)
expected_coloured = []
for swap in (False, True):
    expected_coloured.append(tuple(
        (21 if vertex == 17 else 17 if vertex == 21 else vertex)
        if swap else vertex for vertex in full_vertices
    ))
require(tuple(sorted(coloured_automorphisms))
        == tuple(sorted(expected_coloured)),
        "coloured full automorphism group C2")


# Enumerate all unordered decompositions of the 22 core edges into two
# spanning trees.  Edge zero chooses one member of each complementary pair.
vertex_index = {vertex: index for index, vertex in enumerate(core_vertices)}
indexed_edges = tuple((vertex_index[left], vertex_index[right])
                      for left, right in core_edges)
edge_indices = set(range(len(core_edges)))


def is_tree(indices):
    parent = list(range(len(core_vertices)))

    def root(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for index in indices:
        left, right = indexed_edges[index]
        left, right = root(left), root(right)
        if left == right:
            return False
        parent[left] = right
    return True  # exactly V-1 acyclic edges


decomposition_tuples = []
for chosen in combinations(range(len(core_edges)), len(core_vertices) - 1):
    if 0 not in chosen or not is_tree(chosen):
        continue
    chosen_set = set(chosen)
    if is_tree(edge_indices - chosen_set):
        decomposition_tuples.append(chosen)
require(len(decomposition_tuples) == 4960, "unordered bispanning census")
decomposition_digest = sha256("\n".join(
    ",".join(map(str, chosen)) for chosen in decomposition_tuples
).encode("ascii")).hexdigest()
require(decomposition_digest ==
        "28285b3925c94703da8a0267642ffe8d20c30f5d76b443c52c146aa525c2b822",
        "bispanning decomposition digest")

full_mask = (1 << len(core_edges)) - 1
decomposition_masks = {
    sum(1 << index for index in chosen) for chosen in decomposition_tuples
}
exchange_adjacency = {mask: set() for mask in decomposition_masks}
for mask in decomposition_masks:
    left = [index for index in range(len(core_edges)) if mask >> index & 1]
    right = [index for index in range(len(core_edges)) if not (mask >> index & 1)]
    for remove in left:
        for insert in right:
            neighbor = mask ^ (1 << remove) ^ (1 << insert)
            if not (neighbor & 1):
                neighbor = full_mask ^ neighbor
            if neighbor in decomposition_masks:
                exchange_adjacency[mask].add(neighbor)
exchange_degrees = tuple(len(exchange_adjacency[mask])
                         for mask in decomposition_masks)
require((sum(exchange_degrees) // 2, min(exchange_degrees), max(exchange_degrees))
        == (43408, 13, 27), "symmetric-exchange graph census")
start_mask = next(iter(decomposition_masks))
seen = {start_mask}
queue = [start_mask]
for mask in queue:
    for neighbor in exchange_adjacency[mask]:
        if neighbor not in seen:
            seen.add(neighbor)
            queue.append(neighbor)
require(len(seen) == 4960, "connected symmetric-exchange atlas")


# The nontrivial graph automorphism acts freely on unordered tree pairs.
edge_lookup = {frozenset(edge): index for index, edge in enumerate(core_edges)}
swap_map = {17: 21, 21: 17}
edge_permutation = []
for left, right in core_edges:
    image = frozenset((swap_map.get(left, left), swap_map.get(right, right)))
    edge_permutation.append(edge_lookup[image])


def permute_mask(mask):
    return sum(1 << edge_permutation[index]
               for index in range(len(core_edges)) if mask >> index & 1)


fixed_decompositions = 0
for mask in decomposition_masks:
    image = permute_mask(mask)
    if not (image & 1):
        image = full_mask ^ image
    fixed_decompositions += image == mask
require(fixed_decompositions == 0, "no intrinsic invariant tree pair")


# Exhibit one integral chart.  Reduced incidence matrices of complementary
# trees are unimodular; their comparison is an integral GL_11 transition.
first_tree = decomposition_tuples[0]
second_tree = tuple(sorted(edge_indices - set(first_tree)))
root_vertex = 22
row_vertices = tuple(vertex for vertex in core_vertices if vertex != root_vertex)


def incidence_matrix(tree_indices):
    matrix = [[0] * len(tree_indices) for _ in row_vertices]
    row_lookup = {vertex: row for row, vertex in enumerate(row_vertices)}
    for column, edge_index in enumerate(tree_indices):
        left, right = core_edges[edge_index]
        if left in row_lookup:
            matrix[row_lookup[left]][column] = -1
        if right in row_lookup:
            matrix[row_lookup[right]][column] = 1
    return matrix


def determinant_bareiss(matrix):
    matrix = [row[:] for row in matrix]
    previous = 1
    sign = 1
    for column in range(len(matrix) - 1):
        pivot_row = next((row for row in range(column, len(matrix))
                          if matrix[row][column]), None)
        require(pivot_row is not None, ("singular incidence", column))
        if pivot_row != column:
            matrix[column], matrix[pivot_row] = matrix[pivot_row], matrix[column]
            sign = -sign
        pivot = matrix[column][column]
        for row in range(column + 1, len(matrix)):
            for index in range(column + 1, len(matrix)):
                numerator = (matrix[row][index] * pivot
                             - matrix[row][column] * matrix[column][index])
                require(numerator % previous == 0, "Bareiss divisibility")
                matrix[row][index] = numerator // previous
        previous = pivot
        for row in range(column + 1, len(matrix)):
            matrix[row][column] = 0
    return sign * matrix[-1][-1]


def solve(left, right):
    size = len(left)
    augmented = [[Fraction(value) for value in left[row] + right[row]]
                 for row in range(size)]
    for column in range(size):
        pivot_row = next(row for row in range(column, size)
                         if augmented[row][column])
        augmented[column], augmented[pivot_row] = (
            augmented[pivot_row], augmented[column])
        pivot = augmented[column][column]
        augmented[column] = [value / pivot for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            multiple = augmented[row][column]
            augmented[row] = [augmented[row][index]
                              - multiple * augmented[column][index]
                              for index in range(2 * size)]
    answer = [row[size:] for row in augmented]
    require(all(value.denominator == 1 for row in answer for value in row),
            "integral tree-pair transition")
    return [[value.numerator for value in row] for row in answer]


incidence_first = incidence_matrix(first_tree)
incidence_second = incidence_matrix(second_tree)
determinants = (determinant_bareiss(incidence_first),
                determinant_bareiss(incidence_second))
require(all(abs(value) == 1 for value in determinants),
        "unimodular tree incidences")
transition = solve(incidence_first, incidence_second)
transition_determinant = determinant_bareiss(transition)
require(abs(transition_determinant) == 1, "GL11 transition")
transition_digest = sha256("\n".join(
    ",".join(map(str, row)) for row in transition
).encode("ascii")).hexdigest()


print("THM-3260 BISPANNING RESET-LINK HOLOTOPY ATLAS EXACT AUDIT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("pair_bank=(31 total,8 delayed,23 reset-link)")
print("full_graph=(V13,E23,beta11),unique_bridge=(13,14),leaf=14")
print("core_graph=(V12,E22,beta11),bipartition=%s" % (color_parts,))
print("kuratowski=K33_subdivision,left=(10,3,18),right=(11,16,22),PASS")
print("automorphism_group=C2,nontrivial=(17 21)")
print("covering_graph=(V15,E31,beta17),unique_bridge=(3,9)")
print("relative_layer=(8 delayed edges,2 new vertices,boundary_rank2,H1_rank6)")
print("covering_automorphism_group=S3(10,17,21)xC2(14,18),orders=(1,2,3,6)")
print("edge_coloured_automorphism_group=C2,nontrivial=(17 21)")
print("unordered_two_tree_decompositions=4960,digest=%s" % decomposition_digest)
print("symmetric_exchange_graph=(V4960,E43408,min_degree13,max_degree27,connected)")
print("nontrivial_automorphism_fixed_tree_pairs=0")
print("first_tree=%s" % (tuple(core_edges[index] for index in first_tree),))
print("second_tree=%s" % (tuple(core_edges[index] for index in second_tree),))
print("incidence_determinants=%s,transition_det=%d" %
      (determinants, transition_determinant))
print("transition_digest=%s" % transition_digest)
print("conditional_bridge=cyclic_difference+C12_vertex_label+tree_pair_gives_integral_11D_polarization")
print("scope=no_faithful_C12_torsor_action,no_order7_action,no_planar_dual,no_canonical_owner_phase_map")
print("all_exact_checks=PASS")
