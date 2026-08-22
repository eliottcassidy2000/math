#!/usr/bin/env python3
"""Exact selector-origin cut and abstract phase bridges for THM-3278."""

import ast
import hashlib
from itertools import combinations
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / (
        "01-canon/theorems/THM-3266-all-common-two-face-row-pairs-require-"
        "one-origin-bit.md"
    ): "e4c906cbfa68161c435285826f4e337ed0f17fbc9c9d8e2cc90b61f21001be80",
    ROOT / (
        "04-computation/fc3_common_covering_pair_face_blind_conflict_atlas_"
        "20260803.py"
    ): "c03a837f1ed2fbbadc4c9aaef8609a79b1411a1898a56a57546e5460f1fdca56",
    ROOT / (
        "05-knowledge/results/fc3_common_covering_pair_face_blind_conflict_"
        "atlas_20260803.out"
    ): "40d074ea26d838bd43edea52b2cdaffeaf27a6a2cdb0dc9abedcd6156eed0e82",
    ROOT / (
        "01-canon/theorems/THM-3269-scale-invariant-clutch-strength-and-"
        "canonical-weighted-bispanning-polarization.md"
    ): "55ca134eece22299ebc0c0e997f67a3247c6ee70281c34dd4f9f87cf631647d0",
    ROOT / "04-computation/gmc_scale_invariant_weighted_bispanning_thm3269.py":
        "41ae9aeb01fea1384f59f3a2687b1a0482954bf202e5be2b6fc928ef579b116a",
    ROOT / "05-knowledge/results/gmc_scale_invariant_weighted_bispanning_thm3269.out":
        "65e6e1bb04f3d42d64c2a8e5322c0f3b37c9d05ed6b053670a2cd46293742e3c",
    ROOT / (
        "01-canon/theorems/THM-3274-norm-fibre-constrained-phase-transfer-"
        "and-refinement-invoice.md"
    ): "8570c0bd8217a1c5d2e2682381553e38332dac54aa72c3fa629b5c995ddfb21e",
    ROOT / "04-computation/lrc_constrained_norm_phase_transfer_atlas_20260803.py":
        "24a77d063bc36b45c43d10427124a33b68e8972e63fa97ddfdeba305d0cdb523",
    ROOT / "05-knowledge/results/lrc_constrained_norm_phase_transfer_atlas_20260803.out":
        "bdd318c9cd54a7ac26aaacf6e7bd199288192aaf8d7c5413b7641deeb5915b99",
    ROOT / (
        "01-canon/theorems/THM-3275-unrestricted-twenty-two-row-face-blind-"
        "selector-obstruction.md"
    ): "3b74419b801f380e37dfd3d8abfb933a167dcae1fcc5d81b6fba6dba63debe9d",
    ROOT / (
        "04-computation/fc3_unrestricted_twenty_two_row_face_blind_"
        "selector_no_go_20260803.py"
    ): "2ad9eacf8e893f881b5616672d5fde872a80612b780372dd2318d57f75b6ea30",
    ROOT / (
        "05-knowledge/results/fc3_unrestricted_twenty_two_row_face_blind_"
        "selector_no_go_20260803.out"
    ): "664b7677e89873dda26cd004878c134e621086095d4155d307b52213471baf63",
    ROOT / (
        "01-canon/theorems/THM-3277-weighted-critical-phase-geodesic-"
        "backbone-and-exchange-subatlas.md"
    ): "a10483aefd514292a5277ac6e5e426e961da55aa6f2ff2c1249bf82e679c42ec",
    ROOT / "04-computation/gmc_weighted_critical_phase_path_atlas_scout_20260803.py":
        "0beca08f9214bd6befeafdc0ccc7648be33dd8c1c2d2f12560d2e56708f2cfcf",
    ROOT / "05-knowledge/results/gmc_weighted_critical_phase_path_atlas_scout_20260803.out":
        "a55cf2bb5c130d06aaae1c84dc70d2e78e570d82253c8fdb624c4f9d86ede2e0",
}


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


for dependency, expected in DEPENDENCIES.items():
    require(lf_sha256(dependency) == expected,
            ("dependency drift", dependency.name))

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assert_nodes == 0 and float_literals == 0,
        "optimization-sensitive or floating evidence")

THM3269_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_scale_invariant_weighted_bispanning_thm3269.out"
)
THM3266_OUTPUT = ROOT / (
    "05-knowledge/results/fc3_common_covering_pair_face_blind_conflict_atlas_"
    "20260803.out"
)
THM3274_OUTPUT = ROOT / (
    "05-knowledge/results/lrc_constrained_norm_phase_transfer_atlas_20260803.out"
)
THM3275_OUTPUT = ROOT / (
    "05-knowledge/results/fc3_unrestricted_twenty_two_row_face_blind_"
    "selector_no_go_20260803.out"
)
THM3277_OUTPUT = ROOT / (
    "05-knowledge/results/gmc_weighted_critical_phase_path_atlas_scout_"
    "20260803.out"
)


def literal_line(path, prefix):
    matches = [
        line[len(prefix):]
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(prefix)
    ]
    require(len(matches) == 1, ("transcript key drift", path.name, prefix))
    return ast.literal_eval(matches[0])


tree = tuple(literal_line(THM3269_OUTPUT, "unique_minimum_tree="))
cotree = tuple(literal_line(THM3269_OUTPUT, "canonical_complement="))
rank_label = dict(literal_line(THM3269_OUTPUT, "canonical_C12_rank_label="))
j12_phase = dict(literal_line(THM3269_OUTPUT, "canonical_J12_phase="))
phase_edges = dict(literal_line(THM3269_OUTPUT, "canonical_six_edge_sampler="))
core_edges = tuple(sorted(set(tree) | set(cotree)))
core_vertices = tuple(sorted({vertex for edge in core_edges for vertex in edge}))
require(len(core_edges) == 22 and len(core_vertices) == 12,
        "core graph census drift")

# THM-3266 formerly required the identity of a common covering pair as a
# supplied sidecar.  The unique primitive edge selected by THM-3269 is one of
# those pairs and has a single conflict orientation on all nineteen of its
# pair-specific hostile states.
common_covering_pairs = tuple(literal_line(
    THM3266_OUTPUT, "common_covering_pairs="
))
pair_conflict_records = tuple(
    ast.literal_eval(line.split("=", 1)[1])
    for line in THM3266_OUTPUT.read_text(encoding="utf-8").splitlines()
    if line.startswith("pair_conflict=")
)
primitive_pair_record = next(
    record for record in pair_conflict_records if record[0] == (16, 17)
)
require(
    len(common_covering_pairs) == 24
    and (16, 17) in common_covering_pairs
    and primitive_pair_record[3] == (19, 0)
    and "static_face_bit=NECESSARY_AND_SUFFICIENT_for_each_of_24_pairs"
    in THM3266_OUTPUT.read_text(encoding="utf-8"),
    "canonical primitive covering-pair sidecar drift",
)

# The phase-geodesic backbone intentionally omits the origin carrier.  The
# weighted tree supplies it as one of the four edges completing five forest
# components; the delayed edge is a third, direct type-four layer.
backbone_line = next(
    line for line in THM3277_OUTPUT.read_text(encoding="utf-8").splitlines()
    if line.startswith("phase_geodesic_backbone=")
)
backbone = tuple(ast.literal_eval(
    backbone_line.split("=", 1)[1].split(",forest_components=", 1)[0]
))
tree_completion = tuple(sorted(set(tree) - set(backbone)))
require(
    tree_completion == ((2, 13), (2, 17), (3, 16), (16, 17))
    and "forest_components=(6,3,1,1,1)" in backbone_line,
    "phase-backbone completion drift",
)

conflicts = tuple(
    ast.literal_eval(line.split("=", 1)[1])
    for line in THM3275_OUTPUT.read_text(encoding="utf-8").splitlines()
    if line.startswith("unrestricted_conflict=")
)
require(tuple(record[0] for record in conflicts)
        == ((3, 4, 5), (1, 3, 4, 5)), "conflict locus drift")
small_core_banks = tuple(
    frozenset(record[1]) & frozenset(core_vertices) for record in conflicts
)
full_core_banks = tuple(
    frozenset(record[2]) & frozenset(core_vertices) for record in conflicts
)
require(len(set(small_core_banks)) == len(set(full_core_banks)) == 1,
        "conflict core availability depends on state")
small_core = small_core_banks[0]
full_core = full_core_banks[0]
require(
    small_core == frozenset({2, 11, 16, 18, 22})
    and full_core == frozenset({3, 7, 10, 13, 17, 19, 21})
    and small_core.isdisjoint(full_core)
    and small_core | full_core == frozenset(core_vertices),
    "availability bipartition drift",
)
require(all((left in small_core) != (right in small_core)
            for left, right in core_edges), "core edge misses availability cut")

# Reconstruct the unique root-normalized F2 vertex cochain with derivative one
# on every core edge.  Root 17 is canonical by THM-3269.
root = 17
other_vertices = tuple(vertex for vertex in core_vertices if vertex != root)
rooted_solutions = []
for mask in range(1 << len(other_vertices)):
    colour = {root: 0}
    colour.update({
        vertex: (mask >> index) & 1
        for index, vertex in enumerate(other_vertices)
    })
    if all((colour[left] - colour[right]) % 2 == 1
           for left, right in core_edges):
        rooted_solutions.append(colour)
require(len(rooted_solutions) == 1, "rooted bipartition is not unique")
origin_colour = rooted_solutions[0]
require(frozenset(vertex for vertex, value in origin_colour.items() if value)
        == small_core, "rooted colour is not the small-face class")

# The primitive THM-3269 edge is a canonical paired selector on both conflict
# states.  Orient every canonical phase edge from the small to the full class.
require(all(16 in record[1] and 17 in record[2] for record in conflicts),
        "primitive endpoint selector lost legality")
small_to_full_sampler = []
for orbit in range(1, 7):
    left, right = phase_edges[orbit]
    require((left in small_core) != (right in small_core),
            ("phase edge misses origin cut", orbit))
    tail, head = (left, right) if left in small_core else (right, left)
    residue = (j12_phase[head] - j12_phase[tail]) % 12
    small_to_full_sampler.append((orbit, tail, head, residue))
small_to_full_sampler = tuple(small_to_full_sampler)
require(small_to_full_sampler == (
    (1, 16, 17, 11), (2, 2, 21, 10), (3, 18, 19, 3),
    (4, 11, 17, 8), (5, 22, 21, 7), (6, 22, 7, 6),
), "origin-oriented phase sampler drift")
all_sampler_residues = {
    residue
    for _, _, _, forward in small_to_full_sampler
    for residue in (forward, (-forward) % 12)
}
require(all_sampler_residues == set(range(1, 12)),
        "oriented origin sampler misses a charged residue")

# The cut is not a critical-group or C12 character.  The normalized full
# Jacobian exponents are the exact THM-3269 values in its equation (18).
full_jac_exponent = {
    17: 0, 16: 1, 18: 289, 13: 2344, 11: 9088, 22: 21481,
    19: 25012, 10: 39646, 7: 48259, 21: 49832, 2: 53266, 3: 60022,
}
require(full_jac_exponent[17] == 0 and origin_colour[17] == 0,
        "root gauge drift")
jacobian_parity_mismatches = tuple(
    vertex for vertex in core_vertices
    if full_jac_exponent[vertex] % 2 != origin_colour[vertex]
)
j12_parity_mismatches = tuple(
    vertex for vertex in core_vertices
    if j12_phase[vertex] % 2 != origin_colour[vertex]
)
rank_parity_mismatches = tuple(
    vertex for vertex in core_vertices
    if rank_label[vertex] % 2 != origin_colour[vertex]
)
require(
    jacobian_parity_mismatches == j12_parity_mismatches == (2, 7, 11)
    and rank_parity_mismatches == (2, 3, 10, 11, 13, 18, 21)
    and full_jac_exponent[11] == 9088
    and j12_phase[11] == rank_label[11] == 4,
    "critical-character hostile drift",
)

# In the nonlinear rank label, the cut indicator has every cyclic Fourier
# mode.  Its cyclic difference therefore spans a finite-index sublattice of
# the integral augmentation ideal.
x = sp.symbols("x")
rank_support = tuple(sorted(rank_label[vertex] for vertex in small_core))
require(rank_support == (1, 2, 4, 5, 10), "rank-label cut support drift")
f_coefficients = tuple(int(index in rank_support) for index in range(12))
f_polynomial = sum(f_coefficients[index] * x ** index for index in range(12))
cyclic_polynomial = x ** 12 - 1
require(sp.gcd(f_polynomial, cyclic_polynomial) == 1,
        "rank-label cut lost a cyclic Fourier mode")
resultant = int(sp.resultant(f_polynomial, cyclic_polynomial, x))
require(resultant == 175, "rank-label cut resultant drift")


def circulant(coefficients):
    matrix = sp.zeros(12)
    for column in range(12):
        for exponent, coefficient in enumerate(coefficients):
            matrix[(exponent + column) % 12, column] = coefficient
    return matrix


f_circulant = circulant(f_coefficients)
f_smith = tuple(abs(int(value)) for value in
                smith_normal_form(f_circulant, domain=ZZ).diagonal())
require(int(f_circulant.det()) == 175
        and f_smith == (1,) * 10 + (5, 35), "cut circulant Smith drift")
g_coefficients = tuple(
    f_coefficients[(index - 1) % 12] - f_coefficients[index]
    for index in range(12)
)
require(g_coefficients == (0, -1, 0, 1, -1, 0, 1, 0, 0, 0, -1, 1),
        "cyclic cut difference drift")
g_polynomial = sum(g_coefficients[index] * x ** index for index in range(12))
require(sp.gcd(g_polynomial, cyclic_polynomial) == x - 1,
        "cyclic difference charged spectrum drift")
g_circulant = circulant(g_coefficients)
g_smith = tuple(abs(int(value)) for value in
                smith_normal_form(g_circulant, domain=ZZ).diagonal())
require(g_circulant.rank() == 11
        and g_smith == (1,) * 10 + (35, 0), "cut-difference lattice drift")

# The genuine phase signed primitive edge is cyclic for THM-3274's
# irreducible norm-fibre transfer.  This is a varying-increment coefficient
# statement, not a physical response walk.
norm_quotient = sp.Matrix(literal_line(THM3274_OUTPUT, "quotient_rows="))
norm_char_high = tuple(literal_line(
    THM3274_OUTPUT, "characteristic_polynomial_high_to_low="
))
require(tuple(int(value) for value in norm_quotient.charpoly().all_coeffs())
        == norm_char_high and norm_quotient == norm_quotient.T,
        "norm-fibre quotient drift")
primitive_phase_vector = sp.zeros(12, 1)
primitive_phase_vector[j12_phase[16], 0] = 1
primitive_phase_vector[j12_phase[17], 0] = -1
krylov = sp.Matrix.hstack(*(
    (norm_quotient ** power) * primitive_phase_vector
    for power in range(12)
))
krylov_determinant = int(krylov.det())
require(krylov_determinant
        == -30396564468585450830924251521449984,
        "primitive origin-phase vector is not cyclic")

# Pull the nonlinear rank-label bit back along an abstract C12 phase label on
# the THM-3274 seam decoder.  Even this one Boolean output needs seven raw
# target-count coordinates; total degree does not lower the threshold.
P = 13


def field_add(left, right):
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def field_multiply(left, right):
    return (
        (left[0] * right[0] + 2 * left[1] * right[1]) % P,
        (left[0] * right[1] + left[1] * right[0]) % P,
    )


alpha = (1, 2)
field_points = []
point = (1, 0)
for exponent in range(168):
    field_points.append(point)
    point = field_multiply(point, alpha)
require(point == (1, 0) and len(set(field_points)) == 168,
        "seam field model drift")
field_phase = {point: exponent % 12
               for exponent, point in enumerate(field_points)}
seam_exponents = tuple(range(6)) + tuple(range(162, 168))
seam_increments = tuple(field_points[exponent] for exponent in seam_exponents)
raw_profiles = []
for source in field_points:
    profile = [0] * 12
    for increment in seam_increments:
        target = field_add(source, increment)
        if target != (0, 0):
            profile[field_phase[target]] += 1
    raw_profiles.append(tuple(profile))
raw_profiles = tuple(raw_profiles)
rank_cut_bits = tuple(
    int(exponent % 12 in rank_support) for exponent in range(168)
)


def coordinates_decode_cut(coordinates, include_degree=False):
    seen = {}
    for profile, bit in zip(raw_profiles, rank_cut_bits):
        signature = tuple(profile[index] for index in coordinates)
        if include_degree:
            signature += (sum(profile),)
        if signature in seen and seen[signature] != bit:
            return False
        seen[signature] = bit
    return True


def first_decoder_width(include_degree):
    for width in range(8):
        good = tuple(
            coordinates for coordinates in combinations(range(12), width)
            if coordinates_decode_cut(coordinates, include_degree)
        )
        if good:
            return width, good
    raise RuntimeError("rank-cut seam decoder exceeds width seven")


raw_cut_width, raw_cut_seven_sets = first_decoder_width(False)
degree_cut_width, degree_cut_seven_sets = first_decoder_width(True)
require(
    raw_cut_width == degree_cut_width == 7
    and len(raw_cut_seven_sets) == 95
    and len(degree_cut_seven_sets) == 122
    and raw_cut_seven_sets[0] == degree_cut_seven_sets[0]
    == (0, 1, 2, 3, 4, 5, 11),
    "rank-cut raw decoder width drift",
)

print("THM-3278 SELECTOR-ORIGIN BIPARTITION AND PHASE BRIDGE EXACT AUDIT")
print(f"dependency_hash_checks={len(DEPENDENCIES)}")
print(f"assert_nodes={assert_nodes},float_literals={float_literals}")
print("conflict_states=((3,4,5),(1,3,4,5))")
print("core=(V12,E22);availability_bipartition=(small5,full7)")
print("small_core=(2,11,16,18,22);full_core=(3,7,10,13,17,19,21)")
print("rooted_F2_cochain=(root17=0,dc=1_on_all22_edges,unique)")
print("primitive_endpoint_selector=(small:16,full:17,both_conflicts)")
print("global_controller=(canonical_pair_(16,17),pair_identity_sidecar_0,face_origin_bits_1)")
print("canonical_pair_conflicts=(19,one_orientation_small16_full17)")
print("phase_backbone_completion=((2,13),(2,17),(3,16),(16,17));origin_edge_in_completion")
print("small_to_full_six_edge_sampler=" + repr(small_to_full_sampler))
print("forward_reverse_phase_residues=all_nonzero_C12")
print("critical_character_mismatches=(Jac/J12:(2,7,11),rank:(2,3,10,11,13,18,21))")
print("common_character_hostile=row11:(colour1,Jac9088,J12_4,rank_4)")
print("rank_cut_polynomial=x+x^2+x^4+x^5+x^10")
print("rank_cut_resultant=175;circulant_SNF=(1^10,5,35)")
print("cyclic_difference=-x+x^3-x^4+x^6-x^10+x^11")
print("difference_orbit=(rank11,index35_in_integral_augmentation,SNF_1^10_35_0)")
print("primitive_phase_vector=e1-e0;norm_fibre_Krylov_det="
      + repr(krylov_determinant))
print("seam_rank_cut_decoder=(raw_min7_sets95,degree_min7_sets122,first_(0,1,2,3,4,5,11))")
print("scope=localized_origin_bit_internalized_not_erased;abstract_transfers_only")
print("no_memoryless_selector_no_physical_walk_no_FC3_or_LRC_or_GMC_decrement")
print("all_exact_checks=PASS")
