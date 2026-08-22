#!/usr/bin/env python3
"""Exact LRC norm-phase factorization and physical-map boundary scout.

This is a bounded discovery/type-audit artifact, not a promoted theorem
companion.  It starts with the exact THM-3246/3252/3253 owner data and the
deterministic THM-3234 Singer model, then asks how much of the twelve-valued
norm phase survives the nearest endpoint-current quotients.
"""

import ast
import contextlib
import hashlib
import io
import runpy
from collections import Counter
from math import gcd
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]


def rooted(path):
    return ROOT / path


DEPENDENCIES = {
    rooted(
        "01-canon/theorems/"
        "THM-2791-full-arm-orbit-transfer-and-lower-central-chord.md"
    ): "046b27dda461f33870ba178514801b97a2ff45dd3b56e5b9044ced4e2c2d9c6d",
    rooted(
        "01-canon/theorems/"
        "THM-2803-endpoint-current-determinant-fibre-projective-nonflatness.md"
    ): "54c0b645b509f6b66cf311f5de197a179b04d431da526ada390087602f8b7cf1",
    rooted(
        "01-canon/theorems/"
        "THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-"
        "endpoint-translation-no-go.md"
    ): "160fb4b1274c9dad3039d4b90f32cff76e2a17231d5a5ebcf93f53f1ffa7f26f",
    rooted(
        "01-canon/theorems/"
        "THM-3234-singer-owner-compactification-and-pointed-heisenberg-"
        "carrier-gate.md"
    ): "ef77a1f8fce16eb851eb38d5110a61ab73aa693f2d0ee9e11a912aa4fc302c87",
    rooted(
        "01-canon/theorems/"
        "THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word.md"
    ): "6badc0c9aba09b56d3d055a96cb8ef8b619d8492508bf21476eba5f624b13055",
    rooted(
        "01-canon/theorems/"
        "THM-3247-heisenberg-central-fourier-decomposition-and-canonical-"
        "current-cyclicity.md"
    ): "88b4ccc6d15b87c0aeffcb51f5a3d10cc9037861913e65788e2d74ea82b0252b",
    rooted(
        "01-canon/theorems/"
        "THM-3252-singer-compactified-owner-hodge-word-universal-charged-"
        "cyclicity.md"
    ): "1f8797de2d5fac74814fb78ca4f4d500de8c42eb14a6e1721e5f3e2a2810a873",
    rooted(
        "01-canon/theorems/"
        "THM-3253-positive-owner-mass-newton-cyclicity-and-maximal-common-"
        "heisenberg-module.md"
    ): "b94aea11abe97a6cc1a3826a91fab59d4c04e15f0e6acd9c924c5463b7bd63e8",
    rooted("04-computation/lrc_second_owner_all_dilation_seam_thm3246.py"):
        "e23b098b38aa2199a348f48f8ab4ac0ce5913c870ead972bd31296494fc25a4b",
    rooted("05-knowledge/results/lrc_second_owner_all_dilation_seam_thm3246.out"):
        "d7f7dd96b01c597113e78f903cad36246cb47b10e9a1758cb831aa0e83e8cebc",
    rooted("04-computation/lrc_owner_hodge_charged_cyclicity_thm3252.py"):
        "4471da35fa0fd63dc5c920ef7b695936be722ccbfbb12aeb8cbd55e2607d15c9",
    rooted("05-knowledge/results/lrc_owner_hodge_charged_cyclicity_thm3252.out"):
        "6576214230219b9759646d50f88636ec7a35eab459b19ea722813302744e9d99",
    rooted("04-computation/lrc_positive_owner_mass_newton_cyclicity_thm3253.py"):
        "89aa2a399848ae52e8dd18de9967c7ea2940c04521434ad99407f7be96bdd700",
    rooted("05-knowledge/results/lrc_positive_owner_mass_newton_cyclicity_thm3253.out"):
        "a96010c22126d391bf490f8535dcb3b93809f63e8b705fd23c6223a962bdae58",
    rooted("04-computation/lrc_twelve_balance_singer_rank_defect_thm3255.py"):
        "e1b42874879ae1057418bb9aa0f95bf8d5af2140af415e49ea8a0c7f72cfd35f",
    rooted("05-knowledge/results/lrc_twelve_balance_singer_rank_defect_thm3255.out"):
        "7a229b92d63577a3d79eb78b34418a39b42d083ab0a2f065a6bc1106978d6e45",
    rooted(
        "04-computation/"
        "lrc_multiplicative_singer_twelve_balance_discovery_20260803.py"
    ): "1066ed22d8f6e2ac8a424f62a2d8ea2649133642ebf9c133853893d7b580e143",
    rooted(
        "05-knowledge/results/"
        "lrc_multiplicative_singer_twelve_balance_discovery_20260803.out"
    ): "cb4a89fe66d7ae38dd11266b66d5a2c6509999cabcef9a7e02e252df71462d4f",
}


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n")


for dependency, expected_hash in DEPENDENCIES.items():
    require(
        hashlib.sha256(lf_bytes(dependency)).hexdigest() == expected_hash,
        ("dependency hash drift", dependency.name),
    )

syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")

# Gate the exact type boundaries used below.  This is a selected-contract
# audit, not a claim that string search proves global mathematical absence.
CLAUSES = {
    rooted(
        "01-canon/theorems/"
        "THM-2791-full-arm-orbit-transfer-and-lower-central-chord.md"
    ): (
        "same-ancestry partial translation germ on the THM-2471 rail sheet",
        "map is now from a fixed literal rail sheet to the endpoint-origin current",
    ),
    rooted(
        "01-canon/theorems/"
        "THM-3234-singer-owner-compactification-and-pointed-heisenberg-"
        "carrier-gate.md"
    ): (
        "destroyed:   interval geometry, phase endpoints, q-values, physical ancestry;",
        "sidecar:     a lawful owner-to-endpoint-origin map carrying the LRC predicate;",
    ),
    rooted(
        "01-canon/theorems/"
        "THM-3247-heisenberg-central-fourier-decomposition-and-canonical-"
        "current-cyclicity.md"
    ): (
        "canonical coefficient current J_q on the endpoint plane;",
        "a lawful physical current plus an H_13-equivariant descent/intertwiner;",
    ),
    rooted(
        "01-canon/theorems/"
        "THM-3252-singer-compactified-owner-hodge-word-universal-charged-"
        "cyclicity.md"
    ): (
        "externally relabels the coefficient plane by a primitive Singer",
        "No\nowner-to-endpoint map",
    ),
    rooted(
        "01-canon/theorems/"
        "THM-3253-positive-owner-mass-newton-cyclicity-and-maximal-common-"
        "heisenberg-module.md"
    ): (
        "placing\nthem at Singer-plane points and on the central slice `delta=0` is an abstract\nrelocation",
        "No physical LRC owner-to-plane map, canonical endpoint packet,",
    ),
}

clause_checks = 0
for document, clauses in CLAUSES.items():
    document_text = lf_bytes(document).decode("utf-8")
    for clause in clauses:
        require(clause in document_text, ("typed contract drift", document.name, clause))
        clause_checks += 1

for output_path in (
    rooted("05-knowledge/results/lrc_second_owner_all_dilation_seam_thm3246.out"),
    rooted("05-knowledge/results/lrc_owner_hodge_charged_cyclicity_thm3252.out"),
    rooted("05-knowledge/results/lrc_positive_owner_mass_newton_cyclicity_thm3253.out"),
    rooted("05-knowledge/results/lrc_twelve_balance_singer_rank_defect_thm3255.out"),
    rooted(
        "05-knowledge/results/"
        "lrc_multiplicative_singer_twelve_balance_discovery_20260803.out"
    ),
):
    require(b"all_exact_checks=PASS" in lf_bytes(output_path),
            ("stored exact companion is not PASS", output_path.name))

# Replay THM-3246 to inherit the exact all-dilation sign word.  The finite
# field, quotient, endpoint-fibre and translation computations below are new
# and independent of the inherited script.
dependency_script = rooted(
    "04-computation/lrc_second_owner_all_dilation_seam_thm3246.py"
)
dependency_output = rooted(
    "05-knowledge/results/lrc_second_owner_all_dilation_seam_thm3246.out"
)
dependency_stdout = io.StringIO()
with contextlib.redirect_stdout(dependency_stdout):
    inherited = runpy.run_path(str(dependency_script))
require(
    dependency_stdout.getvalue().encode("utf-8") == lf_bytes(dependency_output),
    "THM-3246 transcript drift",
)

P = 13
N = P * P - 1
ALPHA = (1, 2)
ZERO = (0, 0)
ONE = (1, 0)


def field_add(left, right):
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def field_multiply(left, right):
    a, b = left
    c, d = right
    return ((a * c + 2 * b * d) % P, (a * d + b * c) % P)


def field_power(base, exponent):
    result = ONE
    while exponent:
        if exponent % 2:
            result = field_multiply(result, base)
        base = field_multiply(base, base)
        exponent //= 2
    return result


def field_norm(point):
    return (point[0] * point[0] - 2 * point[1] * point[1]) % P


def determinant(left, right):
    return (left[0] * right[1] - left[1] * right[0]) % P


orbit = tuple(field_power(ALPHA, exponent) for exponent in range(N))
require(len(set(orbit)) == N and ZERO not in orbit, "primitive Singer orbit")
require(field_power(ALPHA, N) == ONE, "Singer order divides 168")
require(field_power(ALPHA, N // 2) != ONE, "Singer order is not 84")
log_alpha = {point: exponent for exponent, point in enumerate(orbit)}

norm_generator = field_norm(ALPHA)
norm_powers = tuple(pow(norm_generator, exponent, P) for exponent in range(12))
require(norm_generator == 6 and len(set(norm_powers)) == 12,
        "primitive norm generator")
require(
    all(field_norm(point) == norm_powers[exponent % 12]
        for exponent, point in enumerate(orbit)),
    "norm phase equals exponent modulo twelve",
)

# Verify the exact 14-direction x 12-scalar factorization geometrically.
for left_exponent, left in enumerate(orbit):
    for right_exponent, right in enumerate(orbit):
        same_projective_direction = determinant(left, right) == 0
        require(
            same_projective_direction
            == ((left_exponent - right_exponent) % 14 == 0),
            "projective direction/exponent mismatch",
        )

canonical_direction_supports = Counter(
    tuple(sorted({exponent % 12 for exponent in range(direction, N, 14)}))
    for direction in range(14)
)
even_phases = (0, 2, 4, 6, 8, 10)
odd_phases = (1, 3, 5, 7, 9, 11)
require(
    canonical_direction_supports
    == Counter({even_phases: 7, odd_phases: 7}),
    "projective directions retain exactly phase parity",
)

# Exhaust all primitive Singer gauges.  For each of their 14 projective
# fibres, audit both the source owner phase j mod 12 and the norm phase of the
# placed target alpha^(b+a*j).  Each fibre has six phases, twice each.
units = tuple(multiplier for multiplier in range(N) if gcd(multiplier, N) == 1)
require(len(units) == 48, "unit multiplier count")
projective_fibre_checks = 0
source_fibre_patterns = Counter()
target_fibre_patterns = Counter()
source_support_sizes = Counter()
target_support_sizes = Counter()
seam_phase_patterns = Counter()
seam_line_patterns = Counter()

q_word = tuple(inherited["q_word"])
negative_seam = tuple(owner for owner, value in enumerate(q_word) if value < 0)
expected_seam = (0, 1, 2, 3, 4, 5, 162, 163, 164, 165, 166, 167)
require(negative_seam == expected_seam, "exact negative seam")
require(tuple(sorted(owner % 12 for owner in negative_seam)) == tuple(range(12)),
        "source seam phase transversal")

for multiplier in units:
    for shift in range(N):
        fibres = [[] for _ in range(14)]
        for owner in range(N):
            placed_exponent = (shift + multiplier * owner) % N
            fibres[placed_exponent % 14].append((owner, placed_exponent))
        for fibre in fibres:
            require(len(fibre) == 12, "projective fibre size")
            source_counts = Counter(owner % 12 for owner, _ in fibre)
            target_counts = Counter(exponent % 12 for _, exponent in fibre)
            source_fibre_patterns[tuple(sorted(source_counts.values()))] += 1
            target_fibre_patterns[tuple(sorted(target_counts.values()))] += 1
            source_support_sizes[len(source_counts)] += 1
            target_support_sizes[len(target_counts)] += 1
            projective_fibre_checks += 1

        seam_target_phases = Counter(
            (shift + multiplier * owner) % 12 for owner in negative_seam
        )
        seam_lines = Counter(
            (shift + multiplier * owner) % 14 for owner in negative_seam
        )
        seam_phase_patterns[tuple(sorted(seam_target_phases.values()))] += 1
        seam_line_patterns[tuple(sorted(seam_lines.values()))] += 1

six_twice = (2, 2, 2, 2, 2, 2)
twelve_once = (1,) * 12
require(projective_fibre_checks == 48 * 168 * 14,
        "projective gauge fibre count")
require(source_fibre_patterns == Counter({six_twice: projective_fibre_checks}),
        "source phase projective mixing")
require(target_fibre_patterns == Counter({six_twice: projective_fibre_checks}),
        "target norm phase projective mixing")
require(source_support_sizes == Counter({6: projective_fibre_checks}),
        "source projective support size")
require(target_support_sizes == Counter({6: projective_fibre_checks}),
        "target projective support size")
require(seam_phase_patterns == Counter({twelve_once: 48 * 168}),
        "seam target phase transversality in every gauge")
require(seam_line_patterns == Counter({twelve_once: 48 * 168}),
        "seam projective transversality in every gauge")

canonical_seam_phase_to_line = tuple(
    (owner % 12, owner % 14) for owner in negative_seam
)
require(
    canonical_seam_phase_to_line
    == tuple((phase, phase) for phase in range(6))
    + tuple((phase, phase + 2) for phase in range(6, 12)),
    "canonical seam phase/line matching",
)
canonical_seam_missing_lines = tuple(
    sorted(set(range(14)) - {line for _, line in canonical_seam_phase_to_line})
)
require(canonical_seam_missing_lines == (6, 7), "canonical seam line gap")

# Audit the nearest canonical endpoint-current coordinate v=det(q,R).
# Across all nonzero increments q, v alone mixes all 12 phases uniformly.
# Even retaining the projective direction [q] leaves exactly the six phases
# of one parity.  Only the full q, together with the chosen F_169 norm, fixes
# the twelve-valued phase.
plane = tuple((row, column) for row in range(P) for column in range(P))
determinant_phase_counts = [[0] * 12 for _ in range(P)]
direction_determinant_phase_counts = {
    (direction, value): [0] * 12
    for direction in range(14)
    for value in range(P)
}
for exponent, increment in enumerate(orbit):
    phase = exponent % 12
    direction = exponent % 14
    for origin in plane:
        value = determinant(increment, origin)
        determinant_phase_counts[value][phase] += 1
        direction_determinant_phase_counts[(direction, value)][phase] += 1

require(
    all(tuple(counts) == (182,) * 12 for counts in determinant_phase_counts),
    "determinant fibre must be completely phase blind",
)
require(
    all(
        tuple(sorted(count for count in counts if count)) == (26,) * 6
        and sum(count != 0 for count in counts) == 6
        for counts in direction_determinant_phase_counts.values()
    ),
    "projective direction plus determinant retains parity only",
)

# The standard H_13 point action contains all affine translations.  The
# translation (r,w)->(r+1,w) gives both a full-phase and a parity-breaking
# hostile.  We also enumerate its entire nonzero-to-nonzero transition law
# and the aggregate over all 168 nonzero translations.
translation = (1, 0)
translation_differences = Counter()
for point, source_exponent in log_alpha.items():
    target = field_add(point, translation)
    if target == ZERO:
        continue
    target_exponent = log_alpha[target]
    translation_differences[(target_exponent - source_exponent) % 12] += 1
require(
    translation_differences
    == Counter({0: 13, **{difference: 14 for difference in range(1, 12)}}),
    "unit translation phase-transition census",
)
require(field_add(orbit[0], translation) == orbit[70],
        "same-parity hostile 0 to 10")
require(field_add(orbit[1], translation) == orbit[40],
        "parity-breaking hostile 1 to 4")

all_translation_differences = Counter()
for translation_point in orbit:
    for point, source_exponent in log_alpha.items():
        target = field_add(point, translation_point)
        if target == ZERO:
            continue
        target_exponent = log_alpha[target]
        all_translation_differences[
            (target_exponent - source_exponent) % 12
        ] += 1
require(
    all_translation_differences
    == Counter({0: 2184, **{difference: 2352 for difference in range(1, 12)}}),
    "all nonzero translation phase-transition census",
)

print("LRC norm-phase physical factorization boundary exact scout")
print(
    f"dependency_hash_checks={len(DEPENDENCIES)},"
    "dependency_replay=THM3246_PASS,stored_companions_3246_3252_3253_3255=PASS"
)
print(
    f"assert_nodes={assert_nodes},float_literals={float_literals},"
    f"typed_contract_clause_checks={clause_checks}"
)
print(
    "singer_field=F13[u]/(u^2-2),alpha=(1,2),orbit=168,"
    "Norm(alpha)=6,phase=exponent_mod12"
)
print(
    "factorization_ladder=full_q_plus_chosen_norm:C12_exact;"
    "projective_[q]:C2_parity_exact;determinant_v:trivial"
)
print(
    f"primitive_gauges={len(units) * N},"
    f"projective_fibre_checks={projective_fibre_checks},"
    "each_fibre=6_phases_each_twice"
)
print(
    "canonical_projective_phase_supports="
    "{0,2,4,6,8,10}*7+{1,3,5,7,9,11}*7"
)
print(
    "endpoint_determinant_fibres=13,states_each=2184,"
    "phase_counts_each=12*182"
)
print(
    "projective_direction_plus_det_fibres=182,states_each=156,"
    "phase_counts_each=6*26"
)
print(
    "negative_seam=one_per_C12_phase+12_distinct_projective_directions_"
    "in_every_8064_gauges"
)
print(
    "canonical_seam_phase_to_line="
    "0:0,1:1,2:2,3:3,4:4,5:5,6:8,7:9,8:10,9:11,10:12,11:13;"
    "missing_lines=6,7"
)
print(
    "translation_(1,0)_nonzero_transitions=167;"
    "phase_preserving=13;each_nonzero_phase_difference=14"
)
print(
    "translation_hostiles=alpha^0=(1,0)->(2,0)=alpha^70:phase_0_to_10;"
    "alpha^1=(1,2)->(2,2)=alpha^40:phase_1_to_4_parity_break"
)
print(
    "nonzero_translations=168;nonzero_to_nonzero_transitions=28056;"
    "phase_preserving=2184;"
    "each_nonzero_phase_difference=2352"
)
print(
    "nearest_physical_source=THM2791_same_ancestry_partial_rail_sheet_germ;"
    "nearest_current=THM3247_J_q(R)=P_(R+q)Q_R"
)
print(
    "first_missing_field=allocated_full_increment_q_and_endpoint_origin_R_"
    "on_the_same_literal_ancestry_sheet"
)
print(
    "additional_phase_sidecar=chosen_F169_Singer_structure_or_equivalent_"
    "absolute_norm_gauge;projective_direction_is_insufficient"
)
print(
    "map_status=pinned_nearest_contracts_leave_owner_to_endpoint_"
    "intertwiner_OPEN;no_commutative_physical_square_constructed"
)
print(
    "scope=finite_exact_abstract_factorization_and_selected_type_audit;"
    "no_physical_owner_current,row_exclusion,or_LRC14_decrement"
)
print("all_exact_checks=PASS")
