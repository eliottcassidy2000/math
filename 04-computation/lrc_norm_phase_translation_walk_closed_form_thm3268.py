#!/usr/bin/env python3
"""Exact norm-phase walk recurrence and fixed-translation hostile for THM-3268.

This companion independently reconstructs THM-3267's chosen F_169 model.
At each step it allows every nonzero increment that keeps the next point in
F_169^*, equivalently every edge of the loopless complete digraph on
F_169^*.  It does not claim that one fixed translation, a physical ancestry
path, or an LRC current has this transition law.
"""

import ast
import hashlib
from collections import Counter
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / (
        "01-canon/theorems/"
        "THM-3255-twelve-balance-multiplicative-singer-rank-defect-and-"
        "phase-marker-boundary.md"
    ): "8d827a3c1c93904db6715240373e8626323802240a4831d14394f3637389609d",
    ROOT / (
        "01-canon/theorems/"
        "THM-3267-norm-phase-factorization-ladder-and-projective-"
        "determinant-blindness.md"
    ): "49231a793e1fbac900c49d9bfcdb7f3373ea7dc29ca2da576b30468be8ecdcc4",
    ROOT / (
        "04-computation/"
        "lrc_norm_phase_physical_factorization_no_go_scout_20260803.py"
    ): "a16658e5b2ad41d9b32ebd8b7b5bd7ffedc26a913b7a0bff75d9fad2f2a617e3",
    ROOT / (
        "05-knowledge/results/"
        "lrc_norm_phase_physical_factorization_no_go_scout_20260803.out"
    ): "392d903d8451cb9df895050bab54ead3f515580950b6a9dd243061a3409a857c",
}


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


for dependency, expected_hash in DEPENDENCIES.items():
    require(
        lf_sha256(dependency) == expected_hash,
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

inherited_output = (
    ROOT
    / "05-knowledge/results/"
    / "lrc_norm_phase_physical_factorization_no_go_scout_20260803.out"
).read_text(encoding="utf-8")
require(
    "nonzero_translations=168;nonzero_to_nonzero_transitions=28056;"
    "phase_preserving=2184;each_nonzero_phase_difference=2352"
    in inherited_output,
    "THM-3267 one-step translation census drift",
)
require("all_exact_checks=PASS" in inherited_output, "THM-3267 output not PASS")

# Reconstruct F_169=F_13[u]/(u^2-2) without importing the inherited scout.
P = 13
ZERO = (0, 0)
ONE = (1, 0)
ALPHA = (1, 2)
PHASE_COUNT = 12
POINT_COUNT = P * P - 1


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


points = tuple(field_power(ALPHA, exponent) for exponent in range(POINT_COUNT))
require(
    len(set(points)) == POINT_COUNT
    and ZERO not in points
    and field_multiply(points[-1], ALPHA) == ONE,
    "primitive Singer orbit drift",
)
phase = {point: exponent % PHASE_COUNT for exponent, point in enumerate(points)}
require(
    all(
        field_norm(point) == pow(6, exponent % PHASE_COUNT, P)
        for exponent, point in enumerate(points)
    ),
    "norm phase is not exponent modulo twelve",
)
phase_fibres = Counter(phase.values())
require(
    phase_fibres == Counter({residue: 14 for residue in range(PHASE_COUNT)}),
    "norm fibres are not twelve classes of fourteen",
)

# Every nonzero increment supplies 167 nonzero-to-nonzero edges.  Conversely,
# every ordered pair of distinct nonzero points has the unique label y-x.
translation_edge_counts = Counter()
for increment in points:
    translation_edge_counts[sum(
        field_add(source, increment) != ZERO for source in points
    )] += 1
require(
    translation_edge_counts == Counter({167: 168}),
    "nonzero translation edge census drift",
)

# Derive the equitable quotient from every one of the 168 source points.
# A source sees 13 other points in its phase and all 14 points in each other
# phase, so the relative quotient row is (13,14,...,14).
relative_quotient_rows = Counter()
for source in points:
    source_phase = phase[source]
    counts = Counter(
        (phase[target] - source_phase) % PHASE_COUNT
        for target in points
        if target != source
    )
    relative_quotient_rows[tuple(
        counts[difference] for difference in range(PHASE_COUNT)
    )] += 1
expected_relative_row = (13,) + (14,) * 11
require(
    relative_quotient_rows == Counter({expected_relative_row: 168}),
    "phase partition is not equitable",
)

quotient = tuple(tuple(
    13 if source_phase == target_phase else 14
    for target_phase in range(PHASE_COUNT)
) for source_phase in range(PHASE_COUNT))
constant = (1,) * PHASE_COUNT
require(
    tuple(sum(row[column] * constant[column] for column in range(PHASE_COUNT))
          for row in quotient) == (167,) * PHASE_COUNT,
    "constant quotient eigenvalue drift",
)
centered_markers = tuple(tuple(
    11 if coordinate == marked_phase else -1
    for coordinate in range(PHASE_COUNT)
) for marked_phase in range(PHASE_COUNT))
require(
    all(
        tuple(sum(row[column] * marker[column]
                  for column in range(PHASE_COUNT)) for row in quotient)
        == tuple(-value for value in marker)
        for marker in centered_markers
    ),
    "rank-eleven centered-marker eigenspace drift",
)


def preserving_formula(length):
    return 14 * (167 ** length + 11 * ((-1) ** length))


def shifted_formula(length):
    return 14 * (167 ** length - ((-1) ** length))


# Point-level dynamic programming deliberately does not use the quotient
# matrix.  It retains the initial phase and traverses every point edge.
walk_dp = [
    [int(phase[point] == initial_phase) for point in points]
    for initial_phase in range(PHASE_COUNT)
]
preserving_terms = []
shifted_terms = []
point_dp_profiles = []
for length in range(6):
    difference_counts = [0] * PHASE_COUNT
    for initial_phase, row in enumerate(walk_dp):
        for point, multiplicity in zip(points, row):
            difference_counts[
                (phase[point] - initial_phase) % PHASE_COUNT
            ] += multiplicity
    require(
        difference_counts
        == [preserving_formula(length)]
        + [shifted_formula(length)] * (PHASE_COUNT - 1),
        ("point-level walk formula failed", length),
    )
    require(
        sum(difference_counts) == POINT_COUNT * 167 ** length,
        ("point-level total drift", length),
    )
    preserving_terms.append(difference_counts[0])
    shifted_terms.append(difference_counts[1])
    point_dp_profiles.append(tuple(difference_counts))

    next_dp = [[0] * POINT_COUNT for _ in range(PHASE_COUNT)]
    for initial_phase, row in enumerate(walk_dp):
        for source_index, multiplicity in enumerate(row):
            if not multiplicity:
                continue
            for target_index in range(POINT_COUNT):
                if target_index != source_index:
                    next_dp[initial_phase][target_index] += multiplicity
    walk_dp = next_dp

require(
    all(
        sequence[index + 2]
        == 166 * sequence[index + 1] + 167 * sequence[index]
        for sequence in (preserving_terms, shifted_terms)
        for index in range(4)
    ),
    "order-two recurrence drift",
)
require(
    preserving_terms[0] == 168
    and preserving_terms[1] - 166 * preserving_terms[0] == -168 * 153
    and shifted_terms[0] == 0
    and shifted_terms[1] - 166 * shifted_terms[0] == 2352,
    "ordinary generating function numerator drift",
)
require(
    all(
        preserving_terms[length] - shifted_terms[length]
        == 168 * ((-1) ** length)
        for length in range(6)
    ),
    "alternating phase-preserving excess drift",
)

# General finite-field formula.  For F_(q^2)/F_q the norm map is onto with
# fibres of size q+1.  Therefore the same loopless complete-graph quotient is
# (q+1)J_(q-1)-I, with eigenvalues q^2-2 and -1.  These exact algebraic checks
# exercise its initial values, total, and recurrence for many integer q; the
# universal proof is the fibre-size argument, not finite interpolation.
def general_preserving(q, length):
    return (q + 1) * (
        (q * q - 2) ** length + (q - 2) * ((-1) ** length)
    )


def general_shifted(q, length):
    return (q + 1) * (
        (q * q - 2) ** length - ((-1) ** length)
    )


general_q_checks = 0
for q in range(2, 31):
    preserve = tuple(general_preserving(q, length) for length in range(9))
    shift = tuple(general_shifted(q, length) for length in range(9))
    require(preserve[0] == q * q - 1 and shift[0] == 0,
            ("general initial condition drift", q))
    require(
        all(
            preserve[length] + (q - 2) * shift[length]
            == (q * q - 1) * (q * q - 2) ** length
            for length in range(9)
        ),
        ("general total drift", q),
    )
    require(
        all(
            sequence[index + 2]
            == (q * q - 3) * sequence[index + 1]
            + (q * q - 2) * sequence[index]
            for sequence in (preserve, shift)
            for index in range(7)
        ),
        ("general recurrence drift", q),
    )
    general_q_checks += 1

# Hostile: matching the one-step census does not authorize iterating one fixed
# translation.  Reuse t=(1,0), demand every intermediate point stay nonzero,
# and compare lengths two and thirteen with the free-translation walk.
def fixed_translation_profile(length):
    counts = Counter()
    for source in points:
        target = source
        valid = True
        for _ in range(length):
            target = field_add(target, ONE)
            if target == ZERO:
                valid = False
                break
        if valid:
            counts[(phase[target] - phase[source]) % PHASE_COUNT] += 1
    return tuple(counts[difference] for difference in range(PHASE_COUNT))


fixed_profiles = tuple(fixed_translation_profile(length) for length in range(1, 14))
require(fixed_profiles[0] == expected_relative_row,
        "fixed-translation one-step control drift")
require(fixed_profiles[1] == (12,) + (14,) * 11,
        "fixed-translation length-two hostile drift")
require(fixed_profiles[12] == (156,) + (0,) * 11,
        "fixed-translation characteristic-thirteen hostile drift")
require(
    tuple(map(sum, fixed_profiles)) == tuple(range(167, 155, -1)) + (156,),
    "fixed-translation survival totals drift",
)

profile_digest = hashlib.sha256(
    repr(tuple(point_dp_profiles)).encode("ascii")
).hexdigest()

print("LRC nonzero-translation norm-phase walk exact closed form")
print(f"dependency_hash_checks={len(DEPENDENCIES)}")
print(f"assert_nodes={assert_nodes},float_literals={float_literals}")
print("field=F13[u]/(u^2-2);alpha=(1,2);nonzero_points=168;norm_phases=12*14")
print("allowed_step=(x,t,x+t)_all_nonzero;translation_graph=loopless_complete_digraph_K168;outdegree=167")
print("equitable_quotient=14*J_12-I_12;J_12=12x12_all_ones;relative_row=(13,14,14,14,14,14,14,14,14,14,14,14)")
print("quotient_spectrum=(167:multiplicity_1,-1:multiplicity_11)")
print("centered_phase_markers=rank_11_augmentation_eigenspace_with_eigenvalue_-1")
print("definition=C_n(d)=point_walks_(x0,...,xn)_with_consecutive_unequal_and_phase(xn)-phase(x0)=d")
print("closed_form_C_n(0)=14*(167^n+11*(-1)^n)")
print("closed_form_C_n(d_nonzero)=14*(167^n-(-1)^n)")
print("recurrence=C_(n+2)=166*C_(n+1)+167*C_n")
print("ogf_C0=168*(1-153*z)/((1-167*z)*(1+z))")
print("ogf_Cnonzero=2352*z/((1-167*z)*(1+z))")
print("preserving_terms_n0_to_n5=" + repr(tuple(preserving_terms)))
print("each_nonzero_shift_terms_n0_to_n5=" + repr(tuple(shifted_terms)))
print("point_level_dp_profiles_n0_to_n5_sha256=" + profile_digest)
print("normalized_augmentation_multiplier_per_step=-1/167")
print("general_prime_power_q_quotient=(q+1)*J_(q-1)-I_(q-1);spectrum=(q^2-2,-1^(q-2))")
print("general_C_n(identity)=(q+1)*((q^2-2)^n+(q-2)*(-1)^n)")
print("general_C_n(nonidentity)=(q+1)*((q^2-2)^n-(-1)^n)")
print("general_recurrence=C_(n+2)+(3-q^2)*C_(n+1)+(2-q^2)*C_n=0")
print(f"general_symbolic_consequence_checks_q2_to_q30={general_q_checks}")
print("fixed_translation_hostile_t=(1,0):n1=(13,14*11);n2=(12,14*11);n13=(156,0*11)")
print("scope=translations_may_vary_each_step;not_repeated_fixed_translation_or_physical_ancestry;no_LRC14_decrement")
print("all_exact_checks=PASS")
