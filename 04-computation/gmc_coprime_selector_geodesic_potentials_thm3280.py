#!/usr/bin/env python3
"""Exact companion for THM-3280's two-potential integral completion."""

import ast
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / "01-canon/theorems/THM-3260-bispanning-reset-link-holotopy-atlas-and-nonplanar-c12-boundary.md":
        "754250ba866b2190a51540e9b4f75234079c6cd0a20bc23083d45d2a8dffe0c9",
    ROOT / "01-canon/theorems/THM-3269-scale-invariant-clutch-strength-and-canonical-weighted-bispanning-polarization.md":
        "55ca134eece22299ebc0c0e997f67a3247c6ee70281c34dd4f9f87cf631647d0",
    ROOT / "01-canon/theorems/THM-3277-weighted-critical-phase-geodesic-backbone-and-exchange-subatlas.md":
        "a10483aefd514292a5277ac6e5e426e961da55aa6f2ff2c1249bf82e679c42ec",
    ROOT / "04-computation/gmc_weighted_critical_phase_path_atlas_scout_20260803.py":
        "0beca08f9214bd6befeafdc0ccc7648be33dd8c1c2d2f12560d2e56708f2cfcf",
    ROOT / "05-knowledge/results/gmc_weighted_critical_phase_path_atlas_scout_20260803.out":
        "a55cf2bb5c130d06aaae1c84dc70d2e78e570d82253c8fdb624c4f9d86ede2e0",
    ROOT / "01-canon/theorems/THM-3278-selector-origin-bit-weighted-tree-bipartition-and-critical-character-boundary.md":
        "d6c0fe74d61b11951c9efa734fadbba0bec7f86a76621812b584ec1590eba89e",
    ROOT / "04-computation/fc3_selector_origin_bipartition_phase_bridge_thm3278.py":
        "07cf5cc1056bdd978a3a1d1146fee8139bef9a84e3b83aad08a0522b4df63a0f",
    ROOT / "05-knowledge/results/fc3_selector_origin_bipartition_phase_bridge_thm3278.out":
        "5e67eec14698f5258d17c5ba1ed295de3b9ff355aa279adb6d523ce3df683c7a",
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

selector_output = (ROOT / "05-knowledge/results/"
                   "fc3_selector_origin_bipartition_phase_bridge_thm3278.out").read_text()
geodesic_output = (ROOT / "05-knowledge/results/"
                   "gmc_weighted_critical_phase_path_atlas_scout_20260803.out").read_text()
require("rank_cut_polynomial=x+x^2+x^4+x^5+x^10" in selector_output
        and "difference_orbit=(rank11,index35_in_integral_augmentation"
        in selector_output, "selector semantic control drift")
require("weighted_phase_path_bank=" in geodesic_output
        and "phase_geodesic_backbone=" in geodesic_output,
        "geodesic semantic control drift")


ORDER = 12
SELECTOR_PHASES = (1, 2, 4, 5, 10)
GEODESIC_LENGTHS = (2, 2, 1, 1, 2, 1, 1, 1, 2, 1, 1, 2)
EVEN_GEODESIC_PHASES = tuple(
    phase for phase, length in enumerate(GEODESIC_LENGTHS) if length % 2 == 0
)
require(EVEN_GEODESIC_PHASES == (0, 1, 4, 8, 11),
        "even geodesic phase word")

f_selector = tuple(int(index in SELECTOR_PHASES) for index in range(ORDER))
f_geodesic = tuple(int(index in EVEN_GEODESIC_PHASES) for index in range(ORDER))
require(sum(f_selector) == sum(f_geodesic) == 5, "marker masses")


def cyclic_shift(vector, amount):
    return tuple(vector[(index - amount) % ORDER] for index in range(ORDER))


def cyclic_difference(vector):
    shifted = cyclic_shift(vector, 1)
    return tuple(shifted[index] - vector[index] for index in range(ORDER))


def matrix_from_columns(columns):
    rows = len(columns[0])
    return Matrix(rows, len(columns),
                  lambda row, column: columns[column][row])


def circulant(vector):
    return matrix_from_columns(tuple(cyclic_shift(vector, shift)
                                     for shift in range(ORDER)))


def bareiss_determinant(matrix):
    require(matrix.rows == matrix.cols, "Bareiss square matrix")
    work = [[int(matrix[row, column]) for column in range(matrix.cols)]
            for row in range(matrix.rows)]
    sign = 1
    previous = 1
    for pivot_index in range(matrix.rows - 1):
        pivot_row = next((row for row in range(pivot_index, matrix.rows)
                          if work[row][pivot_index] != 0), None)
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = work[pivot_row], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, matrix.rows):
            for column in range(pivot_index + 1, matrix.cols):
                numerator = (work[row][column] * pivot
                             - work[row][pivot_index] * work[pivot_index][column])
                require(numerator % previous == 0, "Bareiss nonintegral division")
                work[row][column] = numerator // previous
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def smith_invariants(matrix):
    smith = smith_normal_form(matrix, domain=ZZ)
    return tuple(abs(int(smith[index, index]))
                 for index in range(min(smith.rows, smith.cols))
                 if smith[index, index] != 0)


selector_circulant = circulant(f_selector)
geodesic_circulant = circulant(f_geodesic)
selector_resultant = bareiss_determinant(selector_circulant)
geodesic_resultant = bareiss_determinant(geodesic_circulant)
require(selector_resultant == 175, "selector resultant")
require(geodesic_resultant == 405, "geodesic resultant")
require(smith_invariants(selector_circulant) == (1,) * 10 + (5, 35),
        "selector circulant Smith form")
require(smith_invariants(geodesic_circulant) == (1,) * 10 + (9, 45),
        "geodesic circulant Smith form")

uncentered_combined = selector_circulant.row_join(geodesic_circulant)
require(smith_invariants(uncentered_combined) == (1,) * 11 + (5,),
        "uncentered combined lattice index")


# Use the saturated augmentation basis b_j=x^j-1, j=1,...,11.  A zero-sum
# coefficient vector w has b-coordinates w_1,...,w_11; deleting an arbitrary
# cyclic-difference coordinate instead would introduce a spurious factor 12.
def augmentation_coordinates(vector):
    require(sum(vector) == 0, "augmentation coordinate sum")
    return tuple(vector[index] for index in range(1, ORDER))


def augmentation_multiplication(vector):
    columns = []
    for power in range(1, ORDER):
        shifted = cyclic_shift(vector, power)
        columns.append(augmentation_coordinates(tuple(
            shifted[index] - vector[index] for index in range(ORDER)
        )))
    return matrix_from_columns(tuple(columns))


def difference_orbit(vector):
    difference = cyclic_difference(vector)
    return matrix_from_columns(tuple(
        augmentation_coordinates(cyclic_shift(difference, power))
        for power in range(ORDER)
    ))


selector_aug = augmentation_multiplication(f_selector)
geodesic_aug = augmentation_multiplication(f_geodesic)
require(abs(bareiss_determinant(selector_aug)) == 35
        and smith_invariants(selector_aug) == (1,) * 10 + (35,),
        "selector augmentation lattice")
require(abs(bareiss_determinant(geodesic_aug)) == 81
        and smith_invariants(geodesic_aug) == (1,) * 9 + (9, 9),
        "geodesic augmentation lattice")

combined_aug = selector_aug.row_join(geodesic_aug)
require(smith_invariants(combined_aug) == (1,) * 11,
        "coprime augmentation completion")

selector_orbit = difference_orbit(f_selector)
geodesic_orbit = difference_orbit(f_geodesic)
require(smith_invariants(selector_orbit) == (1,) * 10 + (35,),
        "selector difference orbit")
require(smith_invariants(geodesic_orbit) == (1,) * 9 + (9, 9),
        "geodesic difference orbit")
require(smith_invariants(selector_orbit.row_join(geodesic_orbit)) == (1,) * 11,
        "combined difference orbits")


def full_rank_index(matrix):
    invariants = smith_invariants(matrix)
    require(len(invariants) == 11, "full augmentation rank")
    index = 1
    for invariant in invariants:
        index *= invariant
    return index


# One relative geodesic potential already kills the cyclic selector defect in
# six of the eleven nonzero relative positions.  The reverse direction cannot
# work with one selector vector because the geodesic cokernel needs two
# generators at the prime 3.
single_geodesic_extension_indices = tuple(
    full_rank_index(selector_aug.row_join(geodesic_aug[:, power - 1]))
    for power in range(1, ORDER)
)
require(single_geodesic_extension_indices
        == (1, 1, 7, 5, 1, 7, 1, 5, 7, 1, 1),
        "single geodesic extension census")
single_selector_extension_indices = tuple(
    full_rank_index(geodesic_aug.row_join(selector_aug[:, power - 1]))
    for power in range(1, ORDER)
)
require(all(index > 1 for index in single_selector_extension_indices),
        "one selector vector cannot kill rank-two 3-primary defect")


# The completion is much sparser than two full orbits.  Fix one geodesic
# shift and choose ten of the twelve selector shifts.  Exactly two choices
# give an integral basis, for every fixed geodesic shift.
sparse_basis_bank = []
for geodesic_shift in range(ORDER):
    bases = []
    for selector_shifts in combinations(range(ORDER), 10):
        columns = tuple(selector_orbit[:, shift]
                        for shift in selector_shifts) + (
                            geodesic_orbit[:, geodesic_shift],
                        )
        determinant = bareiss_determinant(matrix_from_columns(tuple(
            tuple(int(value) for value in column) for column in columns
        )))
        if abs(determinant) == 1:
            bases.append((selector_shifts, determinant))
    require(len(bases) == 2, ("sparse completion count", geodesic_shift))
    sparse_basis_bank.append(tuple(bases))

require(sparse_basis_bank[0] == (
    ((0, 1, 2, 3, 5, 6, 7, 8, 10, 11), 1),
    ((0, 2, 4, 5, 6, 7, 8, 9, 10, 11), -1),
), "canonical sparse basis")
sparse_basis_digest = sha256(repr(tuple(sparse_basis_bank)).encode("ascii")).hexdigest()


# Reverse sparsification.  The geodesic cokernel needs two generators at the
# prime 3, so one selector relative vector cannot suffice.  The canonical
# pair f_A(x-1), f_A(x^2-1) plus nine geodesic orbit columns does suffice.
reverse_bases = []
selector_relative = (selector_aug[:, 0], selector_aug[:, 1])
for geodesic_shifts in combinations(range(ORDER), 9):
    columns = selector_relative + tuple(geodesic_orbit[:, shift]
                                        for shift in geodesic_shifts)
    determinant = bareiss_determinant(matrix_from_columns(tuple(
        tuple(int(value) for value in column) for column in columns
    )))
    if abs(determinant) == 1:
        reverse_bases.append((geodesic_shifts, determinant))
require(len(reverse_bases) == 7
        and reverse_bases[0]
        == ((0, 1, 2, 4, 5, 7, 8, 10, 11), -1),
        "reverse sparse basis census")
reverse_basis_digest = sha256(repr(tuple(reverse_bases)).encode("ascii")).hexdigest()


print("THM-3280 COPRIME SELECTOR-GEODESIC POTENTIAL EXACT AUDIT")
print("dependency_hash_checks=%d" % len(DEPENDENCIES))
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("coordinate_contract=(rank_C12_and_genuine_J12_both_root17_to0_primitive16_to1;identity_abstract_only)")
print("selector_word=x+x^2+x^4+x^5+x^10,mass5")
print("geodesic_even_word=1+x+x^4+x^8+x^11,mass5,lengths=%s" %
      (GEODESIC_LENGTHS,))
print("uncentered_selector=(det175,SNF_1^10_5_35)")
print("uncentered_geodesic=(det405,SNF_1^10_9_45)")
print("uncentered_combined=(index5,SNF_1^11_5,mass_mod5_obstruction)")
print("augmentation_selector=(index35,SNF_1^10_35)")
print("augmentation_geodesic=(index81,SNF_1^9_9_9)")
print("coprime_completion=(gcd_35_81=1,combined_SNF_1^11,Aug_ZC12_exact)")
print("single_geodesic_extension_indices_j1_to11=%s" %
      (single_geodesic_extension_indices,))
print("single_selector_extension_indices_j1_to11=%s" %
      (single_selector_extension_indices,))
print("sparse_basis_at_geodesic_shift0=%s" % (sparse_basis_bank[0],))
print("sparse_basis=(2_per_geodesic_shift,24_total,digest=%s)" %
      sparse_basis_digest)
print("reverse_sparse_basis=(canonical_A_relative_1_2,7_of_220_B_nines,first=%s,digest=%s)" %
      (reverse_bases[0], reverse_basis_digest))
print("scope=canonical_coordinate_identification_not_vertex_equivariant_or_physical,no_FC3_LRC_GMC_decrement")
print("all_exact_checks=PASS")
