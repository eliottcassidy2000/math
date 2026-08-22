#!/usr/bin/env python3
"""Exact controls for the Rule-30/LRC14 observer-typing no-go.

This is an unnumbered, dependency-free, deterministic standard-library
companion.  It checks a deliberately small proposed bridge: sampling the
Rule 30 centre trace at the three valuation depths of the 165 first-depth-one
LRC profiles.  It also checks the elementary additive phase obstruction and
its explicit nonlinear boundary.  It proves no LRC row exclusion.
"""

from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from pathlib import Path


HORIZON = 256
P = 13
HOSTILE_SCRIPT = (
    "04-computation/lrc14_correlation_inverse_hostile_thm2344.py"
)
HOSTILE_OUTPUT = (
    "05-knowledge/results/lrc14_correlation_inverse_hostile_thm2344.out"
)
EXPECTED_HOSTILE_SCRIPT_SHA256 = (
    "ce2513ecd6e8290677d69bb5a00b8cb64216035331ebd10887353a34463271fd"
)
EXPECTED_HOSTILE_OUTPUT_SHA256 = (
    "c5e17223240b8a383c8a17a6869bb500e363e4010c2a286f62fb4310e11a4bed"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path):
    return Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def lf_sha256(path):
    return sha256(lf_bytes(path)).hexdigest()


def rule30_set(horizon):
    """Centre trace from a sparse set of black cells."""
    black = {0}
    trace = []
    for time in range(horizon + 1):
        trace.append(int(0 in black))
        black = {
            position
            for position in range(-time - 1, time + 2)
            if (
                (position - 1 in black)
                ^ ((position in black) or (position + 1 in black))
            )
        }
    return tuple(trace)


def rule30_array(horizon):
    """Independent centre trace from a fixed zero-padded dense array."""
    width = 2 * horizon + 3
    origin = horizon + 1
    row = [0] * width
    row[origin] = 1
    trace = []
    for _ in range(horizon + 1):
        trace.append(row[origin])
        next_row = [0] * width
        for position in range(1, width - 1):
            next_row[position] = row[position - 1] ^ (
                row[position] | row[position + 1]
            )
        row = next_row
    return tuple(trace)


def xor_multiple(value, multiplier):
    """Scalar multiplication in an elementary abelian two-group."""
    return value if multiplier % 2 else 0


def permutation_compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def permutation_power(permutation, exponent):
    out = tuple(range(len(permutation)))
    base = permutation
    while exponent:
        if exponent & 1:
            out = permutation_compose(out, base)
        base = permutation_compose(base, base)
        exponent >>= 1
    return out


repository_root = Path(__file__).resolve().parents[1]
script_hash = lf_sha256(Path(__file__))
hostile_script_hash = lf_sha256(repository_root / HOSTILE_SCRIPT)
hostile_output_hash = lf_sha256(repository_root / HOSTILE_OUTPUT)
require(
    hostile_script_hash == EXPECTED_HOSTILE_SCRIPT_SHA256,
    ("THM-2344 script hash changed", hostile_script_hash),
)
require(
    hostile_output_hash == EXPECTED_HOSTILE_OUTPUT_SHA256,
    ("THM-2344 output hash changed", hostile_output_hash),
)

# Two disjoint implementations of the literal single-seed Rule 30 trace.
trace_sparse = rule30_set(HORIZON)
trace_dense = rule30_array(HORIZON)
require(trace_sparse == trace_dense, "Rule 30 implementations disagree")
prefix20 = "".join(str(bit) for bit in trace_sparse[:20])
require(prefix20 == "11011100110001011001", ("prefix changed", prefix20))

# THM-2349's exact first-depth-one valuation-profile universe.
profiles = tuple(
    (1, middle, deepest)
    for deepest in range(5, 20)
    for middle in range(1, deepest)
)
require(len(profiles) == 165, "first-depth-one profile universe changed")
require(sum(middle == 1 for _, middle, _ in profiles) == 15, "repeat count changed")

signature_fibres = defaultdict(list)
for profile in profiles:
    signature = tuple(trace_sparse[depth] for depth in profile)
    signature_fibres[signature].append(profile)

signature_summary = tuple(
    (
        signature,
        len(rows),
        sum(middle == 1 for _, middle, _ in rows),
        sum(middle > 1 for _, middle, _ in rows),
    )
    for signature, rows in sorted(signature_fibres.items())
)
expected_signature_summary = (
    ((1, 0, 0), 36, 0, 36),
    ((1, 0, 1), 36, 0, 36),
    ((1, 1, 0), 51, 8, 43),
    ((1, 1, 1), 42, 7, 35),
)
require(
    signature_summary == expected_signature_summary,
    ("Rule 30 observer fibre profile changed", signature_summary),
)
mixed_owner_type_fibres = sum(
    repeated > 0 and strict > 0
    for _, _, repeated, strict in signature_summary
)
require(mixed_owner_type_fibres == 2, "owner-type mixing changed")

# Appending a deterministic signature to the full profile is graph-isomorphic
# to the original profile set: projection to the first coordinate is inverse.
fine_observer = tuple(
    (profile, tuple(trace_sparse[depth] for depth in profile))
    for profile in profiles
)
require(len(set(fine_observer)) == 165, "fine observer lost a profile")
require(
    tuple(item[0] for item in fine_observer) == profiles,
    "fine observer projection changed",
)

# Additive/XOR phase no-go.  A homomorphism F_13^2 -> F_2^m is determined by
# two basis images.  Each image v must satisfy 13v=0, but 13v=v in exponent 2.
xor_width_rows = []
candidate_images_tested = 0
for width in range(1, 17):
    torsion_images = []
    for value in range(1 << width):
        candidate_images_tested += 1
        if xor_multiple(value, P) == 0:
            torsion_images.append(value)
    require(torsion_images == [0], ("nonzero 13-torsion in XOR group", width))
    xor_width_rows.append((width, len(torsion_images), len(torsion_images) ** 2))

# Sharp scope boundary: four bits can nonlinearly label thirteen states, and
# a permutation on those labels can have order thirteen.  This does not give
# an additive phase map and it explicitly carries the missing phase labels.
nonlinear_bits = 4
state_count = 1 << nonlinear_bits
cycle13 = list(range(state_count))
for state in range(P):
    cycle13[state] = (state + 1) % P
cycle13 = tuple(cycle13)
identity16 = tuple(range(state_count))
require(permutation_power(cycle13, P) == identity16, "13-cycle power failed")
require(
    all(permutation_power(cycle13, exponent) != identity16 for exponent in range(1, P)),
    "nonlinear phase encoding has smaller order",
)
require((1 << 3) < P <= state_count, "nonlinear bit boundary changed")

# THM-2344's same-axis phase cancellation in pure residue arithmetic.  For
# b=a+m and R=0 mod 13, left, bare-conjugate, and deep phases sum to zero.
same_axis_phase_checks = 0
for inert, left_index, deep_step, target_shift in product(range(P), repeat=4):
    del inert  # The first target coordinate is inert in this aligned model.
    right_index = (left_index + deep_step) % P
    exponent = (
        left_index * target_shift
        - right_index * target_shift
        + deep_step * target_shift
    ) % P
    require(exponent == 0, "same-axis phase failed to cancel")
    same_axis_phase_checks += 1
require(same_axis_phase_checks == P**4, "same-axis phase universe changed")

semantic_object = {
    "horizon": HORIZON,
    "prefix20": prefix20,
    "profile_count": len(profiles),
    "signature_summary": signature_summary,
    "mixed_owner_type_fibres": mixed_owner_type_fibres,
    "fine_observer_count": len(fine_observer),
    "xor_width_rows": xor_width_rows,
    "candidate_images_tested": candidate_images_tested,
    "nonlinear_boundary": (nonlinear_bits, state_count, P),
    "same_axis_phase_checks": same_axis_phase_checks,
    "hostile_hashes": (hostile_script_hash, hostile_output_hash),
}
semantic_hash = sha256(
    json.dumps(semantic_object, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("UNNUMBERED Rule30/LRC14 depth-observer and additive-phase no-go")
print(f"script_sha256_lf={script_hash}")
print(f"rule30_prefix20={prefix20}")
print(f"independent_rule30_implementations_match_through={HORIZON}")
print(f"first_depth_one_profile_count={len(profiles)}")
print(f"repeated_profile_count={sum(middle == 1 for _, middle, _ in profiles)}")
print(f"depth_bit_signature_summary={signature_summary}")
print(f"mixed_repeated_strict_signature_fibres={mixed_owner_type_fibres}")
print("coarse_rule30_depth_signature_preserves_owner_multiplicity=NO")
print(f"full_profile_plus_signature_count={len(fine_observer)}")
print("full_profile_plus_signature_adds_information=NO-DETERMINISTIC-GRAPH")
print(f"xor_phase_width_rows={tuple(xor_width_rows)}")
print(f"xor_candidate_images_tested={candidate_images_tested}")
print("additive_hom_F13sq_to_F2m=ZERO_FOR_EVERY_FINITE_m")
print(
    "nonlinear_boundary="
    f"{nonlinear_bits}_bits_encode_{P}_labels_and_explicit_order_{P}_permutation"
)
print("nonlinear_encoding_excluded_by_additive_no_go=NO")
print("nonlinear_encoding_requires_explicit_phase_labels=YES")
print(f"same_axis_phase_identity_checks={same_axis_phase_checks}")
print("same_axis_rule30_word_supplies_cross_axis_phase=NO")
print(f"thm2344_script_sha256_lf={hostile_script_hash}")
print(f"thm2344_output_sha256_lf={hostile_output_hash}")
print(f"semantic_sha256={semantic_hash}")
print("status=FINITE-EXACT-OBSERVER-PROBE+ELEMENTARY-TYPED-NO-GO;UNNUMBERED")
print("scalar_rows_excluded=0")
print("lrc14_status=OPEN")
print("ALL CHECKS PASSED")
