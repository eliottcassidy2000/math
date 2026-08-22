#!/usr/bin/env python3
"""Exact companion for THM-3273's graph-Jacobian carrier boundary."""

import ast
from collections import Counter
from functools import reduce
from hashlib import sha256
from math import gcd
from pathlib import Path

from sympy import Matrix, zeros
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


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


source_tree = ast.parse(THM3254_SCRIPT.read_text(encoding="utf-8"))
covering_edges = tuple(literal_assignment(source_tree, "ROW_COVERING_PAIRS"))
delayed_edge_set = set(literal_assignment(source_tree, "DELAYED_WITNESSES"))
delayed_edges = tuple(edge for edge in covering_edges if edge in delayed_edge_set)
reset_edges = tuple(edge for edge in covering_edges if edge not in delayed_edge_set)
core_edges = tuple(edge for edge in reset_edges if 14 not in edge)
require((len(covering_edges), len(delayed_edges), len(reset_edges), len(core_edges))
        == (31, 8, 23, 22), "edge-bank split")


def vertices(edges):
    return tuple(sorted({vertex for edge in edges for vertex in edge}))


def reduced_laplacian(edges, root):
    ordered = vertices(edges)
    require(root in ordered, ("missing root", root))
    nonroot = tuple(vertex for vertex in ordered if vertex != root)
    index = {vertex: position for position, vertex in enumerate(nonroot)}
    laplacian = zeros(len(nonroot))
    for left, right in edges:
        if left != root:
            laplacian[index[left], index[left]] += 1
        if right != root:
            laplacian[index[right], index[right]] += 1
        if left != root and right != root:
            laplacian[index[left], index[right]] -= 1
            laplacian[index[right], index[left]] -= 1
    return nonroot, laplacian


def smith_diagonal(matrix):
    smith = smith_normal_form(matrix, domain=ZZ)
    return tuple(abs(int(smith[index, index]))
                 for index in range(min(smith.rows, smith.cols)))


core_vertices, core_laplacian = reduced_laplacian(core_edges, 22)
reset_vertices, reset_laplacian = reduced_laplacian(reset_edges, 22)
full_vertices, full_laplacian = reduced_laplacian(covering_edges, 22)

core_smith = smith_diagonal(core_laplacian)
reset_smith = smith_diagonal(reset_laplacian)
full_smith = smith_diagonal(full_laplacian)
require(core_smith == (1,) * 10 + (74748,), "core critical group")
require(reset_smith == (1,) * 11 + (74748,), "reset critical group")
require(full_smith == (1,) * 11 + (4, 12, 113184),
        "full critical group")
require(74748 == 12 * 6229 and gcd(12, 6229) == 1,
        "core critical factorization")
require(5432832 == 4 * 12 * 113184,
        "full critical order")
require(5432832 == (2 ** 9) * (3 ** 4) * 131,
        "full critical factorization")
require(74748 % 7 and 74748 % 13 and 5432832 % 7 and 5432832 % 13,
        "no 7-primary or 13-primary critical torsion")


# Absence of 13-torsion does not exclude an order-13 automorphism: the
# characteristic C_131 primary factor has Aut(C_131)=C_130.  Compute the
# full-graph Abel--Jacobi chart mod 131 and show that no nontrivial element of
# that order-13 subgroup, even followed by a translation, preserves its image.
def null_vector_mod_prime(matrix, prime):
    work = [[int(matrix[row, column]) % prime
             for column in range(matrix.cols)]
            for row in range(matrix.rows)]
    pivot_columns = []
    pivot_row = 0
    for column in range(matrix.cols):
        pivot = next((row for row in range(pivot_row, matrix.rows)
                      if work[row][column]), None)
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        inverse = pow(work[pivot_row][column], -1, prime)
        work[pivot_row] = [(inverse * value) % prime
                           for value in work[pivot_row]]
        for row in range(matrix.rows):
            if row == pivot_row or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [(left - scale * right) % prime
                         for left, right in zip(work[row], work[pivot_row])]
        pivot_columns.append(column)
        pivot_row += 1
    free_columns = tuple(column for column in range(matrix.cols)
                         if column not in pivot_columns)
    require(len(free_columns) == 1, ("modular nullity", prime, free_columns))
    vector = [0] * matrix.cols
    vector[free_columns[0]] = 1
    for row, column in reversed(tuple(enumerate(pivot_columns))):
        vector[column] = (-sum(work[row][free] * vector[free]
                               for free in free_columns)) % prime
    require(all(sum(int(matrix[row, column]) * vector[column]
                    for column in range(matrix.cols)) % prime == 0
                for row in range(matrix.rows)), "modular kernel check")
    return tuple(vector)


phase131_vector = null_vector_mod_prime(full_laplacian, 131)
require(phase131_vector == (46, 117, 57, 117, 1, 59, 108,
                            44, 91, 30, 1, 91, 41, 1),
        "full Abel-Jacobi mod-131 chart")
phase131 = {22: 0}
phase131.update({vertex: phase131_vector[index]
                 for index, vertex in enumerate(full_vertices)})
phase131_image = set(phase131.values())
require(len(phase131_image) == 11, "full mod-131 vertex image")
order13_multipliers = tuple(pow(39, exponent, 131)
                            for exponent in range(13))
require(len(set(order13_multipliers)) == 13
        and order13_multipliers[0] == order13_multipliers[-1] * 39 % 131,
        "order-13 multiplier subgroup")
require(not any({(unit * value + shift) % 131
                 for value in phase131_image} == phase131_image
                for unit in order13_multipliers[1:]
                for shift in range(131)),
        "no affine order-13 preservation of vertex image")


# The delayed relative chain complex C_1(H,G)->C_0(H,G) has eight edge
# generators and the two new vertex generators 9 and 12.
new_vertices = tuple(vertex for vertex in vertices(covering_edges)
                     if vertex not in vertices(reset_edges))
require(new_vertices == (9, 12), "relative vertex bank")
relative_boundary = Matrix([
    [(-1 if left == vertex else 0) + (1 if right == vertex else 0)
     for left, right in delayed_edges]
    for vertex in new_vertices
])
require(tuple(relative_boundary) == (0, 1, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, -1, -1, 0),
        "relative boundary entries")
require(smith_diagonal(relative_boundary) == (1, 1),
        "saturated relative boundary")

relative_basis = (
    (1, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 1, 0, 0, 0, 0, 0),
    (0, 0, 0, 1, 0, 0, 0, 0),
    (0, 0, 0, 0, 1, 0, 0, 0),
    (0, 0, 0, 0, 0, 0, 0, 1),
    (0, 0, 0, 0, 0, 1, -1, 0),
)
require(all(relative_boundary * Matrix(vector) == zeros(2, 1)
            for vector in relative_basis), "relative cycle basis")
require(Matrix.hstack(*(Matrix(vector) for vector in relative_basis)).rank() == 6,
        "relative cycle independence")

# The edge-coloured graph involution swaps delayed edges (11,17),(11,21)
# and fixes the other relative basis elements.
relative_c2_trace = 4
relative_signature = ((6 + relative_c2_trace) // 2,
                      (6 - relative_c2_trace) // 2)
c7_reflection_trace = 0  # one fixed point on Q[C7], then remove constants
c7_reflection_signature = ((6 + c7_reflection_trace) // 2,
                           (6 - c7_reflection_trace) // 2)
require(relative_signature == (5, 1), "relative C2 signature")
require(c7_reflection_signature == (3, 3), "C7 reflection signature")
require(relative_signature not in ((6, 0), c7_reflection_signature),
        "relative layer is neither trivial nor affine-reflection C7 augmentation")


# A primitive adjugate row identifies the cyclic core Jacobian with Z/74748.
# Reduction mod 12 is the canonical order-12 quotient, up to a unit choice.
core_order = int(core_laplacian.det())
adjugate = core_laplacian.adjugate()
row_position = core_vertices.index(17)
primitive_row = tuple(int(adjugate[row_position, column])
                      for column in range(adjugate.cols))
require(primitive_row == (12951, 7527, 2466, 8883, 7593, 6369,
                          14988, 34229, 4932, 5961, 9313),
        "primitive adjugate row")
require(reduce(gcd, (core_order,) + primitive_row) == 1,
        "primitive cyclic coordinate")
row_times_laplacian = tuple(
    int(value) for value in (Matrix(1, len(primitive_row), primitive_row)
                             * core_laplacian)
)
expected_row_product = tuple(core_order if index == row_position else 0
                             for index in range(len(core_vertices)))
require(row_times_laplacian == expected_row_product,
        "adjugate quotient certificate")

phase = {22: 0}
phase.update({vertex: primitive_row[index] % 12
              for index, vertex in enumerate(core_vertices)})
expected_phase = {
    2: 3, 3: 3, 7: 6, 10: 3, 11: 9, 13: 9,
    16: 0, 17: 5, 18: 0, 19: 9, 21: 1, 22: 0,
}
require(phase == expected_phase, "vertex phase chart")
phase_counts = Counter(phase.values())
require(dict(sorted(phase_counts.items()))
        == {0: 3, 1: 1, 3: 3, 5: 1, 6: 1, 9: 3},
        "six-class vertex image")

pair_differences = {(phase[left] - phase[right]) % 12
                    for left in phase for right in phase}
edge_differences = set()
for left, right in core_edges:
    edge_differences.add((phase[left] - phase[right]) % 12)
    edge_differences.add((phase[right] - phase[left]) % 12)
require(pair_differences == set(range(12)), "complete pair-difference carrier")
require(edge_differences == set(range(12)) - {4, 8},
        "local-edge phase gap")

swap = {vertex: (21 if vertex == 17 else 17 if vertex == 21 else vertex)
        for vertex in phase}
require(all(phase[swap[vertex]] == (5 * phase[vertex]) % 12
            for vertex in phase), "C2 acts by multiplication by 5")

# The two delayed edges swapped by the graph C2 supply exactly the two missing
# phase increments.  Consequently they kill the one-dimensional kernel of the
# augmentation-function sampler on directed core edges.
repair_edges = (delayed_edges[3], delayed_edges[4])
repair_increments = tuple((phase[right] - phase[left]) % 12
                          for left, right in repair_edges)
require(repair_edges == ((11, 17), (11, 21)), "repair edge labels")
require(repair_increments == (8, 4), "repair edge phases")


def directed_evaluation_matrix(edges):
    rows = []
    for left, right in edges:
        for tail, head in ((left, right), (right, left)):
            row = [0] * 12
            row[(phase[head] - phase[tail]) % 12] = 1
            rows.append(row)
    return Matrix(rows)


augmentation_basis = zeros(12, 11)
for column in range(11):
    augmentation_basis[column, column] = 1
    augmentation_basis[11, column] = -1
core_sampler = directed_evaluation_matrix(core_edges) * augmentation_basis
repaired_sampler = directed_evaluation_matrix(core_edges + repair_edges) \
    * augmentation_basis
missing_vector = zeros(12, 1)
missing_vector[4, 0] = 1
missing_vector[8, 0] = -1
require(core_sampler.rank() == 10, "core augmentation sampler rank")
require(directed_evaluation_matrix(core_edges) * missing_vector == zeros(44, 1),
        "core sampler missing line")
require(repaired_sampler.rank() == 11,
        "delayed pair completes augmentation sampler")


print("THM-3273 CRITICAL-GROUP PHASE CARRIER EXACT AUDIT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("core_critical_group=Z/74748=Z/12xZ/6229")
print("reset_critical_group=Z/74748")
print("full_critical_group=Z/4xZ/12xZ/113184,order=5432832")
print("critical_primary_boundary=no_7_torsion,no_13_torsion")
print("C131_automorphism_boundary=order13_exists,no_nontrivial_affine_element_preserves_11class_vertex_image")
print("relative_boundary_snf=(1,1),relative_H1=Z^6,saturated")
print("relative_C2_signature=(5,1),C7_reflection_augmentation=(3,3)")
print("primitive_adjugate_row=%s" % (primitive_row,))
print("C12_vertex_image=(6 classes,multiplicities 3,3,3,1,1,1)")
print("C12_pair_differences=all12,core_edge_differences=all_except_4_8")
print("delayed_phase_repair=((11,17)->8,(11,21)->4),Aug_sampler_rank=10_to_11")
print("graph_C2_on_C12=multiplication_by_5")
print("scope=relative_phase_carrier_only,no_vertex_torsor,no_absolute_marker")
print("all_exact_checks=PASS")
