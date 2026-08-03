#!/usr/bin/env python3
"""Exact constrained-increment transfer atlas in the chosen F_169 model.

The script separates two semantics throughout:

* varying-set semantics: at every step a fresh increment may be chosen from
  the declared set, provided the next point is nonzero;
* fixed-increment semantics: one declared increment is repeated.

It studies a multiplicative norm fibre, an additive projective scalar line,
and the twelve THM-3246 seam points.  The seam points are used here as an
abstract increment alphabet; no physical ancestry or LRC transition is
claimed.
"""

import ast
import hashlib
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT / (
        "01-canon/theorems/"
        "THM-3246-all-dilation-second-owner-seam-stabilization-and-sign-word.md"
    ): "6badc0c9aba09b56d3d055a96cb8ef8b619d8492508bf21476eba5f624b13055",
    ROOT / (
        "01-canon/theorems/"
        "THM-3267-norm-phase-factorization-ladder-and-projective-"
        "determinant-blindness.md"
    ): "49231a793e1fbac900c49d9bfcdb7f3373ea7dc29ca2da576b30468be8ecdcc4",
    ROOT / (
        "01-canon/theorems/"
        "THM-3268-nonzero-translation-norm-phase-walk-closed-form-and-"
        "rank-eleven-mixing-mode.md"
    ): "d39b6216320ed0d6c4b4a03934cf02c9eb983d253237d9a99669beaaf5c75bda",
    ROOT / "04-computation/lrc_norm_phase_translation_walk_closed_form_thm3268.py":
        "fd26af6cba1cad5019556fd2235b319f7d089273b64d7f0dc01d2877ca0a4985",
    ROOT / "05-knowledge/results/lrc_norm_phase_translation_walk_closed_form_thm3268.out":
        "a9d8f74c6f9a101b81fd1c6a652c8f9379d945d2883c502452252fdbfc1fe844",
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

inherited_walk_output = (
    ROOT / "05-knowledge/results/lrc_norm_phase_translation_walk_closed_form_thm3268.out"
).read_text(encoding="utf-8")
require(
    "equitable_quotient=14*J_12-I_12" in inherited_walk_output
    and "all_exact_checks=PASS" in inherited_walk_output,
    "THM-3268 control transcript drift",
)

# Chosen THM-3267 model.
P = 13
ZERO = (0, 0)
ONE = (1, 0)
ALPHA = (1, 2)
POINT_COUNT = 168
PHASE_COUNT = 12


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
point_index = {point: index for index, point in enumerate(points)}
phase = {point: exponent % PHASE_COUNT for exponent, point in enumerate(points)}
require(
    len(point_index) == POINT_COUNT
    and ZERO not in point_index
    and field_multiply(points[-1], ALPHA) == ONE,
    "primitive Singer model drift",
)
require(
    all(
        field_norm(point) == pow(6, exponent % PHASE_COUNT, P)
        for exponent, point in enumerate(points)
    ),
    "norm phase drift",
)
require(
    Counter(phase.values())
    == Counter({residue: 14 for residue in range(PHASE_COUNT)}),
    "phase-fibre sizes drift",
)


def point_adjacency(increments):
    matrix = sp.zeros(POINT_COUNT)
    for source_index, source in enumerate(points):
        for increment in increments:
            target = field_add(source, increment)
            if target != ZERO:
                matrix[source_index, point_index[target]] += 1
    return matrix


def outgoing_phase_profile(source, increments):
    counts = Counter()
    for increment in increments:
        target = field_add(source, increment)
        if target != ZERO:
            counts[phase[target]] += 1
    return tuple(counts[target_phase] for target_phase in range(PHASE_COUNT))


def canonical_partition(labels):
    label_to_colour = {}
    colours = []
    for label in labels:
        if label not in label_to_colour:
            label_to_colour[label] = len(label_to_colour)
        colours.append(label_to_colour[label])
    return tuple(colours)


def same_partition(left, right):
    require(len(left) == len(right), "partition length mismatch")
    return all(
        (left[i] == left[j]) == (right[i] == right[j])
        for i in range(len(left))
        for j in range(len(left))
    )


def quotient_from_labels(increments, labels):
    cells = defaultdict(list)
    for index, label in enumerate(labels):
        cells[label].append(index)
    ordered_labels = tuple(sorted(cells))
    label_index = {label: index for index, label in enumerate(ordered_labels)}
    quotient = sp.zeros(len(ordered_labels))
    for source_label in ordered_labels:
        reference = None
        for source_index in cells[source_label]:
            counts = Counter()
            source = points[source_index]
            for increment in increments:
                target = field_add(source, increment)
                if target != ZERO:
                    counts[labels[point_index[target]]] += 1
            row = tuple(counts[target_label] for target_label in ordered_labels)
            if reference is None:
                reference = row
            require(row == reference, ("partition is not equitable", source_label))
        for target_label, count in zip(ordered_labels, reference):
            quotient[label_index[source_label], label_index[target_label]] = count
    sizes = tuple(len(cells[label]) for label in ordered_labels)
    return ordered_labels, sizes, quotient


def relative_profiles_from_quotient(ordered_labels, sizes, quotient, lengths):
    label_index = {label: index for index, label in enumerate(ordered_labels)}
    state_count = len(ordered_labels)
    distributions = []
    for initial_phase in range(PHASE_COUNT):
        row = [0] * state_count
        for label, size in zip(ordered_labels, sizes):
            if label[0] == initial_phase:
                row[label_index[label]] = size
        distributions.append(sp.Matrix(1, state_count, row))
    profiles = []
    for _ in range(lengths):
        counts = [0] * PHASE_COUNT
        for initial_phase, row in enumerate(distributions):
            for target_index, target_label in enumerate(ordered_labels):
                counts[(target_label[0] - initial_phase) % PHASE_COUNT] += int(
                    row[target_index]
                )
        profiles.append(tuple(counts))
        distributions = [row * quotient for row in distributions]
    return tuple(profiles)


def direct_relative_profiles(increments, lengths):
    distributions = [
        [int(phase[point] == initial_phase) for point in points]
        for initial_phase in range(PHASE_COUNT)
    ]
    profiles = []
    for _ in range(lengths):
        counts = [0] * PHASE_COUNT
        for initial_phase, row in enumerate(distributions):
            for point, multiplicity in zip(points, row):
                counts[(phase[point] - initial_phase) % PHASE_COUNT] += multiplicity
        profiles.append(tuple(counts))
        next_distributions = [[0] * POINT_COUNT for _ in range(PHASE_COUNT)]
        for initial_phase, row in enumerate(distributions):
            for source_index, multiplicity in enumerate(row):
                if not multiplicity:
                    continue
                source = points[source_index]
                for increment in increments:
                    target = field_add(source, increment)
                    if target != ZERO:
                        next_distributions[initial_phase][point_index[target]] += multiplicity
        distributions = next_distributions
    return tuple(profiles)


# Positive control: all increments reproduce THM-3268.
all_increment_profiles = Counter(
    outgoing_phase_profile(source, points) for source in points
)
expected_all_increment_profiles = Counter()
for source_phase in range(PHASE_COUNT):
    expected_all_increment_profiles[
        tuple(
            13 if target_phase == source_phase else 14
            for target_phase in range(PHASE_COUNT)
        )
    ] = 14
require(
    all_increment_profiles == expected_all_increment_profiles,
    "all-increment positive control drift",
)

# ---------------------------------------------------------------------------
# 1. A fixed multiplicative norm fibre: phase is equitable, but order twelve
#    is irreducible.
# ---------------------------------------------------------------------------

norm_fibres = tuple(
    tuple(point for point in points if phase[point] == residue)
    for residue in range(PHASE_COUNT)
)
norm_zero = norm_fibres[0]
require(len(norm_zero) == 14, "norm-zero phase fibre size drift")
require(
    {((-point[0]) % P, (-point[1]) % P) for point in norm_zero}
    == set(norm_zero),
    "norm fibre is not inverse-closed additively",
)

norm_rows = []
for source_phase in range(PHASE_COUNT):
    profiles = {
        outgoing_phase_profile(source, norm_zero)
        for source in points
        if phase[source] == source_phase
    }
    require(len(profiles) == 1, ("norm phase is not equitable", source_phase))
    norm_rows.append(next(iter(profiles)))
norm_quotient = sp.Matrix(norm_rows)
require(norm_quotient == norm_quotient.T, "norm quotient is not symmetric")

EXPECTED_NORM_ROWS = (
    (0, 0, 2, 2, 2, 0, 2, 2, 0, 2, 1, 0),
    (0, 0, 2, 0, 2, 0, 0, 2, 2, 2, 2, 2),
    (2, 2, 1, 2, 0, 2, 1, 2, 2, 0, 0, 0),
    (2, 0, 2, 2, 2, 0, 0, 2, 0, 0, 2, 2),
    (2, 2, 0, 2, 0, 2, 2, 0, 1, 2, 1, 0),
    (0, 0, 2, 0, 2, 2, 2, 2, 2, 0, 2, 0),
    (2, 0, 1, 0, 2, 2, 2, 0, 1, 0, 2, 2),
    (2, 2, 2, 2, 0, 2, 0, 0, 0, 2, 0, 2),
    (0, 2, 2, 0, 1, 2, 1, 0, 2, 2, 0, 2),
    (2, 2, 0, 0, 2, 0, 0, 2, 2, 2, 0, 2),
    (1, 2, 0, 2, 1, 2, 2, 0, 0, 0, 2, 2),
    (0, 2, 0, 2, 0, 0, 2, 2, 2, 2, 2, 0),
)
require(tuple(norm_rows) == EXPECTED_NORM_ROWS, "norm quotient row drift")

NORM_CHAR_HIGH = (
    1,
    -13,
    -77,
    884,
    1508,
    -19811,
    -2325,
    160510,
    -96880,
    -341392,
    243616,
    27968,
    -23104,
)
require(
    tuple(int(coefficient) for coefficient in norm_quotient.charpoly().all_coeffs())
    == NORM_CHAR_HIGH,
    "norm quotient characteristic polynomial drift",
)

# Every other absolute norm fibre is a phase-shift conjugate of this quotient.
for increment_phase, increment_fibre in enumerate(norm_fibres):
    rows = []
    for source_phase in range(PHASE_COUNT):
        profiles = {
            outgoing_phase_profile(source, increment_fibre)
            for source in points
            if phase[source] == source_phase
        }
        require(len(profiles) == 1, ("rotated norm fibre not equitable", increment_phase))
        rows.append(next(iter(profiles)))
    expected = sp.zeros(PHASE_COUNT)
    for source_phase in range(PHASE_COUNT):
        for target_phase in range(PHASE_COUNT):
            expected[source_phase, target_phase] = norm_quotient[
                (source_phase - increment_phase) % PHASE_COUNT,
                (target_phase - increment_phase) % PHASE_COUNT,
            ]
    require(sp.Matrix(rows) == expected, ("norm-fibre conjugacy drift", increment_phase))


# Elementary finite-field polynomial arithmetic gives a compact exact
# irreducibility certificate modulo 331.
def trim_mod(poly, modulus):
    values = [value % modulus for value in poly]
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values)


def divmod_mod(numerator, denominator, modulus):
    numerator = list(trim_mod(numerator, modulus))
    denominator = trim_mod(denominator, modulus)
    require(denominator != (0,), "zero polynomial divisor")
    quotient = [0] * max(1, len(numerator) - len(denominator) + 1)
    inverse_lead = pow(denominator[-1], -1, modulus)
    while len(numerator) >= len(denominator) and any(numerator):
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] * inverse_lead % modulus
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            numerator[shift + index] = (
                numerator[shift + index] - coefficient * value
            ) % modulus
        while len(numerator) > 1 and numerator[-1] == 0:
            numerator.pop()
    return trim_mod(quotient, modulus), trim_mod(numerator, modulus)


def gcd_mod(left, right, modulus):
    left = trim_mod(left, modulus)
    right = trim_mod(right, modulus)
    while right != (0,):
        _, remainder = divmod_mod(left, right, modulus)
        left, right = right, remainder
    inverse_lead = pow(left[-1], -1, modulus)
    return trim_mod(tuple(inverse_lead * value for value in left), modulus)


def multiply_mod(left, right, divisor, modulus):
    product = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            product[left_index + right_index] = (
                product[left_index + right_index] + left_value * right_value
            ) % modulus
    _, remainder = divmod_mod(tuple(product), divisor, modulus)
    return remainder


def power_mod(base, exponent, divisor, modulus):
    result = (1,)
    base = trim_mod(base, modulus)
    while exponent:
        if exponent % 2:
            result = multiply_mod(result, base, divisor, modulus)
        base = multiply_mod(base, base, divisor, modulus)
        exponent //= 2
    return result


irreducibility_prime = 331
norm_char_low_mod = trim_mod(tuple(reversed(NORM_CHAR_HIGH)), irreducibility_prime)
x_polynomial = (0, 1)
require(
    power_mod(
        x_polynomial,
        irreducibility_prime ** 12,
        norm_char_low_mod,
        irreducibility_prime,
    )
    == x_polynomial,
    "norm characteristic polynomial fails Frobenius closure modulo 331",
)
for proper_degree in (6, 4):
    residue = power_mod(
        x_polynomial,
        irreducibility_prime ** proper_degree,
        norm_char_low_mod,
        irreducibility_prime,
    )
    difference = list(residue) + [0] * max(0, 2 - len(residue))
    difference[1] = (difference[1] - 1) % irreducibility_prime
    require(
        gcd_mod(tuple(difference), norm_char_low_mod, irreducibility_prime) == (1,),
        ("norm characteristic polynomial has a proper Frobenius factor", proper_degree),
    )

# Relative-phase path profiles and their exact minimal recurrence.
norm_labels = tuple((residue,) for residue in range(PHASE_COUNT))
norm_sizes = (14,) * PHASE_COUNT
norm_profiles = relative_profiles_from_quotient(
    norm_labels, norm_sizes, norm_quotient, 24
)
require(
    direct_relative_profiles(norm_zero, 6) == norm_profiles[:6],
    "norm-fibre point DP disagrees with quotient DP",
)
for difference in range(PHASE_COUNT):
    sequence = [profile[difference] for profile in norm_profiles]
    for start in range(12):
        require(
            sequence[start + 12]
            + sum(
                NORM_CHAR_HIGH[offset] * sequence[start + 12 - offset]
                for offset in range(1, 13)
            )
            == 0,
            ("norm recurrence drift", difference, start),
        )

norm_hankel_determinants = []
for difference in range(PHASE_COUNT):
    sequence = [profile[difference] for profile in norm_profiles]
    hankel = sp.Matrix(
        12,
        12,
        lambda row, column: sequence[row + column],
    )
    determinant = int(hankel.det())
    require(determinant != 0, ("norm scalar recurrence drops below order 12", difference))
    norm_hankel_determinants.append(determinant)
norm_hankel_digest = hashlib.sha256(
    repr(tuple(norm_hankel_determinants)).encode("ascii")
).hexdigest()

# Fixed-increment hostile inside the same norm fibre.
require(ONE in norm_zero, "fixed norm-fibre hostile missing")
fixed_increment_witness = (
    phase[field_add(points[0], ONE)],
    phase[field_add(points[12], ONE)],
)
require(
    phase[points[0]] == phase[points[12]] == 0
    and fixed_increment_witness == (10, 9),
    "fixed-increment non-equitability witness drift",
)


def repeated_increment_profile(increment, length):
    counts = Counter()
    for source in points:
        target = source
        valid = True
        for _ in range(length):
            target = field_add(target, increment)
            if target == ZERO:
                valid = False
                break
        if valid:
            counts[(phase[target] - phase[source]) % PHASE_COUNT] += 1
    return tuple(counts[difference] for difference in range(PHASE_COUNT))


fixed_profiles = {
    length: repeated_increment_profile(ONE, length)
    for length in (1, 2, 12, 13)
}
require(fixed_profiles[1] == (13,) + (14,) * 11, "fixed length-one drift")
require(fixed_profiles[2] == (12,) + (14,) * 11, "fixed length-two drift")
require(fixed_profiles[13] == (156,) + (0,) * 11, "fixed length-thirteen drift")

# ---------------------------------------------------------------------------
# 2. One additive projective direction: phase fails, but squared transverse
#    determinant is the exact coarsest outgoing-equitable repair.
# ---------------------------------------------------------------------------

scalar_line = tuple((scalar, 0) for scalar in range(1, P))
line_phase_profile_counts = {}
for source_phase in range(PHASE_COUNT):
    line_phase_profile_counts[source_phase] = Counter(
        outgoing_phase_profile(source, scalar_line)
        for source in points
        if phase[source] == source_phase
    )
require(
    all(len(line_phase_profile_counts[source_phase]) == 4
        for source_phase in range(PHASE_COUNT)),
    "projective line does not split every phase four ways",
)
require(
    outgoing_phase_profile(points[0], scalar_line)
    == (1, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0)
    and outgoing_phase_profile(points[12], scalar_line)
    == (1, 2, 0, 0, 0, 1, 2, 0, 2, 2, 0, 2),
    "projective-line phase witness drift",
)

line_labels = tuple((phase[point], point[1] * point[1] % P) for point in points)
line_signature_labels = tuple(
    (phase[point], outgoing_phase_profile(point, scalar_line)) for point in points
)
require(
    same_partition(
        canonical_partition(line_labels),
        canonical_partition(line_signature_labels),
    ),
    "one-step phase refinement is not (phase,b^2)",
)
line_ordered_labels, line_sizes, line_quotient = quotient_from_labels(
    scalar_line, line_labels
)
require(len(line_ordered_labels) == 48, "projective-line refinement size drift")
require(
    Counter(line_sizes) == Counter({4: 36, 2: 12}),
    "projective-line refined cell sizes drift",
)

square_values = (0, 1, 3, 4, 9, 10, 12)
line_blocks = tuple(
    tuple(
        (source_phase, size)
        for (source_phase, transverse_square), size
        in zip(line_ordered_labels, line_sizes)
        if transverse_square == square_value
    )
    for square_value in square_values
)
require(
    tuple(len(block) for block in line_blocks) == (6, 7, 7, 7, 7, 7, 7),
    "projective-line block sizes drift",
)
require(
    all(sum(size for _, size in block) == (12 if index == 0 else 26)
        for index, block in enumerate(line_blocks)),
    "projective-line merged-coset populations drift",
)

line_symbol = sp.symbols("lambda")
expected_line_char = sp.expand(
    (line_symbol - 12) ** 6
    * (line_symbol - 11)
    * (line_symbol + 1) ** 41
)
require(
    sp.expand(line_quotient.charpoly(line_symbol).as_expr()) == expected_line_char,
    "projective-line quotient spectrum drift",
)
line_identity = (
    line_quotient ** 3
    - 22 * line_quotient ** 2
    + 109 * line_quotient
    + 132 * sp.eye(48)
)
require(line_identity == sp.zeros(48), "projective-line cubic identity drift")

line_profiles = relative_profiles_from_quotient(
    line_ordered_labels, line_sizes, line_quotient, 13
)
require(
    direct_relative_profiles(scalar_line, 6) == line_profiles[:6],
    "projective-line point DP disagrees with refined quotient",
)


def projective_line_formula(length, difference):
    if difference == 0:
        numerator = (
            300 * 12 ** length
            + 26 * 11 ** length
            + 1858 * ((-1) ** length)
        )
    elif difference % 2:
        numerator = 168 * 12 ** length - 168 * ((-1) ** length)
    else:
        numerator = (
            144 * 12 ** length
            + 26 * 11 ** length
            - 170 * ((-1) ** length)
        )
    require(numerator % 13 == 0, "projective-line formula lost integrality")
    return numerator // 13


for length, profile in enumerate(line_profiles):
    require(
        profile
        == tuple(
            projective_line_formula(length, difference)
            for difference in range(PHASE_COUNT)
        ),
        ("projective-line closed form drift", length),
    )
    require(
        sum(profile) == 156 * 12 ** length + 12 * 11 ** length,
        ("projective-line total drift", length),
    )
for difference in range(PHASE_COUNT):
    sequence = [profile[difference] for profile in line_profiles]
    require(
        all(
            sequence[index + 3]
            == 22 * sequence[index + 2]
            - 109 * sequence[index + 1]
            - 132 * sequence[index]
            for index in range(10)
        ),
        ("projective-line recurrence drift", difference),
    )

# ---------------------------------------------------------------------------
# 3. THM-3246 seam alphabet: phase refinement is discrete.  Restoring zero
#    gives an additive Cayley graph; deleting zero differentiates its exact
#    characteristic polynomial.
# ---------------------------------------------------------------------------

seam_exponents = tuple(range(6)) + tuple(range(162, 168))
seam = tuple(points[exponent] for exponent in seam_exponents)
require(
    {field_multiply(points[6], increment) for increment in seam}
    == {points[exponent] for exponent in range(12)},
    "seam is not the alpha^-6..alpha^5 arc",
)
require(
    Counter(phase[increment] for increment in seam)
    == Counter({residue: 1 for residue in range(PHASE_COUNT)}),
    "seam is not one point per norm phase",
)

seam_signatures = tuple(
    (phase[source], outgoing_phase_profile(source, seam)) for source in points
)
require(len(set(seam_signatures)) == POINT_COUNT, "seam phase refinement is not discrete")
require(
    all(
        len({
            outgoing_phase_profile(source, seam)
            for source in points
            if phase[source] == source_phase
        })
        == 14
        for source_phase in range(PHASE_COUNT)
    ),
    "seam does not split each phase into fourteen singletons",
)

# The source phase is not needed: the raw twelve target-phase counts already
# locate the source among all 168 punctured points.  Record the exact nearest
# separation and affine dimension, then verify the scalar-normalized straight
# arc by its conjugacy rather than by a second discovery path.
seam_raw_profiles = tuple(
    outgoing_phase_profile(source, seam) for source in points
)
require(len(set(seam_raw_profiles)) == POINT_COUNT,
        "seam raw phase histogram is not a point decoder")
seam_decoder_l1 = min(
    sum(abs(left - right) for left, right in zip(first, second))
    for first, second in combinations(seam_raw_profiles, 2)
)
seam_decoder_hamming = min(
    sum(left != right for left, right in zip(first, second))
    for first, second in combinations(seam_raw_profiles, 2)
)
seam_decoder_affine_rank = sp.Matrix([
    [value - seam_raw_profiles[0][index]
     for index, value in enumerate(profile)]
    for profile in seam_raw_profiles[1:]
]).rank()
require((seam_decoder_l1, seam_decoder_hamming,
         seam_decoder_affine_rank) == (2, 2, 12),
        "seam decoder metric drift")


def projected_signature_count(coordinates, include_degree=False):
    signatures = set()
    for profile in seam_raw_profiles:
        projection = tuple(profile[index] for index in coordinates)
        signatures.add((projection, sum(profile)) if include_degree else projection)
    return len(signatures)


six_projection_counts = {
    coordinates: projected_signature_count(coordinates)
    for coordinates in combinations(range(PHASE_COUNT), 6)
}
best_six_count = max(six_projection_counts.values())
best_six_sets = tuple(coordinates for coordinates, count
                      in six_projection_counts.items()
                      if count == best_six_count)
require((best_six_count, best_six_sets)
        == (166, ((0, 1, 3, 5, 6, 8),)),
        "sharp six-coordinate decoder boundary")

degree_six_counts = {
    coordinates: projected_signature_count(coordinates, include_degree=True)
    for coordinates in combinations(range(PHASE_COUNT), 6)
}
require(max(degree_six_counts.values()) == 166
        and all(count < POINT_COUNT for count in degree_six_counts.values()),
        "degree sidecar unexpectedly repairs six coordinates")

injective_seven_sets = tuple(
    coordinates for coordinates in combinations(range(PHASE_COUNT), 7)
    if projected_signature_count(coordinates) == POINT_COUNT
)
require(len(injective_seven_sets) == 21
        and injective_seven_sets[0] == (0, 1, 2, 5, 8, 9, 11),
        "sharp seven-coordinate decoder census")

cyclic_window_starts = {}
for length in range(1, PHASE_COUNT + 1):
    starts = tuple(
        start for start in range(PHASE_COUNT)
        if projected_signature_count(tuple(
            (start + offset) % PHASE_COUNT for offset in range(length)
        )) == POINT_COUNT
    )
    if starts:
        cyclic_window_starts[length] = starts
require(cyclic_window_starts[min(cyclic_window_starts)]
        == (4, 7, 8, 9, 10, 11)
        and min(cyclic_window_starts) == 8,
        "cyclic window decoder threshold")

straight_arc = tuple(points[exponent] for exponent in range(12))
straight_profiles = tuple(
    outgoing_phase_profile(source, straight_arc) for source in points
)
require(len(set(straight_profiles)) == POINT_COUNT,
        "straight Singer arc is not a raw point decoder")
scale = points[6]
scale_phase = phase[scale]
for source, profile in zip(points, seam_raw_profiles):
    scaled_source = field_multiply(scale, source)
    expected = tuple(profile[(residue - scale_phase) % PHASE_COUNT]
                     for residue in range(PHASE_COUNT))
    require(outgoing_phase_profile(scaled_source, straight_arc) == expected,
            "Singer arc scalar conjugacy drift")

# The seam also has no nonidentity F_13-linear set stabilizer.
def matrix_action(matrix, point):
    return (
        (matrix[0] * point[0] + matrix[1] * point[1]) % P,
        (matrix[2] * point[0] + matrix[3] * point[1]) % P,
    )


seam_linear_stabilizer = []
seam_set = set(seam)
for a in range(P):
    for b in range(P):
        for c in range(P):
            for d in range(P):
                if (a * d - b * c) % P == 0:
                    continue
                matrix = (a, b, c, d)
                if {matrix_action(matrix, point) for point in seam} == seam_set:
                    seam_linear_stabilizer.append(matrix)
require(
    seam_linear_stabilizer == [(1, 0, 0, 1)],
    "seam linear stabilizer drift",
)

seam_adjacency = point_adjacency(seam)
require(
    Counter(sum(int(value) for value in seam_adjacency.row(row))
            for row in range(POINT_COUNT))
    == Counter({12: 156, 11: 12}),
    "seam outdegree census drift",
)
require(
    Counter(sum(int(value) for value in seam_adjacency.col(column))
            for column in range(POINT_COUNT))
    == Counter({12: 156, 11: 12}),
    "seam indegree census drift",
)

seam_char_high = tuple(
    int(coefficient) for coefficient in seam_adjacency.charpoly().all_coeffs()
)

# Cyclotomic-ring reconstruction of the full 169-vertex additive Cayley
# spectrum.  Elements of Z[zeta_13] use the basis 1,z,...,z^11 and
# z^12=-(1+...+z^11).
CYCLOTOMIC_DEGREE = 12
CYCLOTOMIC_ZERO = (0,) * CYCLOTOMIC_DEGREE
CYCLOTOMIC_ONE = (1,) + (0,) * (CYCLOTOMIC_DEGREE - 1)


def cyclotomic_add(left, right):
    return tuple(left[index] + right[index] for index in range(CYCLOTOMIC_DEGREE))


def cyclotomic_negate(value):
    return tuple(-coefficient for coefficient in value)


def cyclotomic_multiply(left, right):
    temporary = [0] * 23
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            temporary[left_index + right_index] += left_value * right_value
    for exponent in range(22, 11, -1):
        coefficient = temporary[exponent]
        if coefficient:
            temporary[exponent] = 0
            for residue in range(12):
                temporary[exponent - 12 + residue] -= coefficient
    return tuple(temporary[:12])


def cyclotomic_power_of_zeta(exponent):
    exponent %= 13
    if exponent < 12:
        return tuple(int(index == exponent) for index in range(12))
    return (-1,) * 12


def character_sum(histogram, galois_scalar):
    result = CYCLOTOMIC_ZERO
    for residue, multiplicity in histogram.items():
        for _ in range(multiplicity):
            result = cyclotomic_add(
                result,
                cyclotomic_power_of_zeta(galois_scalar * residue),
            )
    return result


def multiply_by_linear_factor(polynomial, root):
    result = [CYCLOTOMIC_ZERO] * (len(polynomial) + 1)
    for index, coefficient in enumerate(polynomial):
        result[index] = cyclotomic_add(
            result[index],
            cyclotomic_negate(cyclotomic_multiply(coefficient, root)),
        )
        result[index + 1] = cyclotomic_add(result[index + 1], coefficient)
    return result


def orbit_polynomial(histogram):
    polynomial = [CYCLOTOMIC_ONE]
    for galois_scalar in range(1, 13):
        polynomial = multiply_by_linear_factor(
            polynomial,
            character_sum(histogram, galois_scalar),
        )
    integer_coefficients = []
    for coefficient in polynomial:
        require(
            all(value == 0 for value in coefficient[1:]),
            "Galois orbit polynomial is not rational",
        )
        integer_coefficients.append(coefficient[0])
    return tuple(integer_coefficients)


character_representatives = tuple((1, slope) for slope in range(P)) + ((0, 1),)
character_histograms = []
fourier_orbit_factors_low = []
for character in character_representatives:
    histogram = Counter(
        (character[0] * increment[0] + character[1] * increment[1]) % P
        for increment in seam
    )
    character_histograms.append(
        tuple(histogram[residue] for residue in range(P))
    )
    fourier_orbit_factors_low.append(orbit_polynomial(histogram))
fourier_orbit_factors_low = tuple(fourier_orbit_factors_low)
require(
    Counter(fourier_orbit_factors_low).most_common(1)[0][1] == 2
    and len(set(fourier_orbit_factors_low)) == 13,
    "seam Fourier orbit multiplicities drift",
)
duplicate_factor = next(
    factor for factor, multiplicity in Counter(fourier_orbit_factors_low).items()
    if multiplicity == 2
)
duplicate_directions = tuple(
    representative
    for representative, factor
    in zip(character_representatives, fourier_orbit_factors_low)
    if factor == duplicate_factor
)
require(
    duplicate_directions == ((1, 7), (0, 1)),
    "seam duplicate Fourier directions drift",
)


def integer_polynomial_multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        for right_index, right_value in enumerate(right):
            result[left_index + right_index] += left_value * right_value
    return tuple(result)


full_cayley_char_low = (-12, 1)
for factor in fourier_orbit_factors_low:
    full_cayley_char_low = integer_polynomial_multiply(
        full_cayley_char_low, factor
    )
require(len(full_cayley_char_low) == 170, "full Cayley polynomial degree drift")

deleted_char_low = []
for exponent in range(1, len(full_cayley_char_low)):
    numerator = exponent * full_cayley_char_low[exponent]
    require(numerator % 169 == 0, ("deleted cofactor not integral", exponent))
    deleted_char_low.append(numerator // 169)
deleted_char_low = tuple(deleted_char_low)
require(
    tuple(reversed(seam_char_high)) == deleted_char_low,
    "direct seam characteristic polynomial disagrees with Cayley derivative",
)

lambda_symbol = sp.symbols("lambda")
seam_char_poly = sp.Poly.from_list(seam_char_high, gens=lambda_symbol)
seam_char_derivative = seam_char_poly.diff()
require(
    sp.gcd(seam_char_poly, seam_char_derivative).degree() == 0,
    "seam characteristic polynomial is not squarefree",
)
duplicate_poly = sp.Poly(
    sum(coefficient * lambda_symbol ** exponent
        for exponent, coefficient in enumerate(duplicate_factor)),
    lambda_symbol,
)
seam_cofactor, seam_remainder = sp.div(seam_char_poly, duplicate_poly)
require(seam_remainder.is_zero and seam_cofactor.degree() == 156,
        "seam duplicate-orbit factor did not survive deletion")
rational_factorization = sp.factor_list(seam_char_poly.as_expr())[1]
require(
    tuple(sorted((sp.degree(factor, lambda_symbol), multiplicity)
                 for factor, multiplicity in rational_factorization))
    == ((12, 1), (156, 1)),
    "seam rational factor degrees drift",
)

fourier_factor_digest = hashlib.sha256(
    repr(fourier_orbit_factors_low).encode("ascii")
).hexdigest()
full_cayley_char_digest = hashlib.sha256(
    repr(tuple(reversed(full_cayley_char_low))).encode("ascii")
).hexdigest()
seam_char_digest = hashlib.sha256(
    repr(seam_char_high).encode("ascii")
).hexdigest()
seam_cofactor_digest = hashlib.sha256(
    repr(tuple(int(value) for value in seam_cofactor.all_coeffs())).encode("ascii")
).hexdigest()

# Phase-matched hostile: another one-per-phase transversal has the same
# immediate discrete refinement but a different exact spectrum.
alternate_exponents = tuple(residue + 12 * (residue % 2) for residue in range(12))
alternate_transversal = tuple(points[exponent] for exponent in alternate_exponents)
require(
    Counter(phase[increment] for increment in alternate_transversal)
    == Counter({residue: 1 for residue in range(PHASE_COUNT)}),
    "alternate transversal phase census drift",
)
alternate_signatures = tuple(
    (phase[source], outgoing_phase_profile(source, alternate_transversal))
    for source in points
)
require(len(set(alternate_signatures)) == POINT_COUNT,
        "alternate transversal refinement is not discrete")
alternate_char_high = tuple(
    int(coefficient)
    for coefficient in point_adjacency(alternate_transversal).charpoly().all_coeffs()
)
alternate_char_digest = hashlib.sha256(
    repr(alternate_char_high).encode("ascii")
).hexdigest()
require(alternate_char_digest != seam_char_digest,
        "phase-matched hostile unexpectedly has seam spectrum")

# A single seam increment is the same characteristic-thirteen hostile as in
# the norm-fibre and projective-line experiments.
require(ONE in seam, "fixed seam hostile missing")
require(
    repeated_increment_profile(ONE, 13) == (156,) + (0,) * 11,
    "fixed seam increment hostile drift",
)

print("LRC constrained norm-phase transfer atlas (raw exact research artifact)")
print(f"dependency_hash_checks={len(DEPENDENCIES)}")
print(f"assert_nodes={assert_nodes},float_literals={float_literals},sympy_version={sp.__version__}")
print("field=F13[u]/(u^2-2);alpha=(1,2);phase(alpha^j)=j_mod_12;X=F169_nonzero")
print("semantics_varying=fresh_increment_from_declared_set_each_step_with_nonzero_target")
print("semantics_fixed=repeat_one_declared_increment_and_kill_on_first_zero")
print("positive_control_all_increments=THM3268_quotient_14J12-I12")
print()
print("[fixed absolute norm fibre S0={t:phase(t)=0}, size 14]")
print("phase_partition=equitable_for_varying_S0;all_12_norm_fibres_are_phase_shift_conjugate")
print("quotient_rows=" + repr(EXPECTED_NORM_ROWS))
print("quotient_symmetric=yes;row_sums=(13,14*11)")
print("characteristic_polynomial_high_to_low=" + repr(NORM_CHAR_HIGH))
print("irreducibility_certificate=mod_331_degree_12_Rabin_Frobenius_tests_PASS")
print("spectrum=12_distinct_real_algebraic_roots;minimal_matrix_polynomial=displayed_degree_12_polynomial")
print("relative_path_profiles_n0_to_n5=" + repr(norm_profiles[:6]))
print("relative_sequence_minimal_order=12_for_each_d_by_nonzero_12x12_Hankel_determinant")
print("norm_hankel_determinants_sha256=" + norm_hankel_digest)
print("recurrence=C[n+12]+sum_(i=1)^12(char_coeff[i]*C[n+12-i])=0")
print("fixed_t=(1,0)_same_phase_sources_alpha^0_and_alpha^12_target_phases=" + repr(fixed_increment_witness))
print("fixed_t_profiles_n1_n2_n12_n13=" + repr(fixed_profiles))
print()
print("[projective scalar line L*=F13^*(1,0), size 12]")
print("phase_partition=not_equitable;each_phase_has_4_distinct_outgoing_phase_profiles")
print("witness_alpha^0=" + repr(outgoing_phase_profile(points[0], scalar_line)))
print("witness_alpha^12=" + repr(outgoing_phase_profile(points[12], scalar_line)))
print("coarsest_outgoing_equitable_refinement_of_phase=(phase,b^2)=48_cells")
print("refined_cell_sizes=(36_cells_of_4,12_cells_of_2)")
print("transverse_square_blocks_s=0,1,3,4,9,10,12:" + repr(line_blocks))
print("quotient_spectrum=(12:multiplicity_6,11:multiplicity_1,-1:multiplicity_41)")
print("minimal_polynomial=(lambda-12)(lambda-11)(lambda+1)")
print("recurrence=C[n+3]=22*C[n+2]-109*C[n+1]-132*C[n]")
print("relative_path_profiles_n0_to_n5=" + repr(line_profiles[:6]))
print("closed_form_d0=(300*12^n+26*11^n+1858*(-1)^n)/13")
print("closed_form_d_odd=168*(12^n-(-1)^n)/13")
print("closed_form_d_even_nonzero=(144*12^n+26*11^n-170*(-1)^n)/13")
print("total_paths=156*12^n+12*11^n")
print("general_direction=conjugate_by_multiplication_d_inverse;normalized_sidecar=(phase(x)-phase(d),(det(d,x)/Norm(d))^2)")
print("fixed_single_scalar_increment=characteristic_13_hostile_not_the_48_state_transfer")
print()
print("[THM3246 seam N=alpha^{0..5,162..167}, size 12]")
print("interpretation=abstract_varying_increment_alphabet_only;not_a_physical_seam_transition")
print("Singer_arc_normalization=alpha^6*N={alpha^0,...,alpha^11}")
print("phase_census=one_increment_per_C12_phase")
print("phase_partition=not_equitable;first_outgoing_refinement_has_168_singletons")
print("coarsest_outgoing_equitable_refinement_of_phase=full_Singer_point_or_exponent")
print("raw_target_phase_histogram=168_distinct_without_source_phase")
print("raw_decoder_metrics=(min_L1_2,min_Hamming_2,affine_rank_12)")
print("sharp_coordinate_decoder=(min_7,injective_7sets_21,first_(0,1,2,5,8,9,11),max_6_signatures_166,best_6_(0,1,3,5,6,8))")
print("degree_sidecar=six_coordinates_still_max_166")
print("cyclic_window_decoder=(min_8,starts_(4,7,8,9,10,11))")
print("straight_arc_raw_decoder=168_distinct_by_alpha^6_scalar_conjugacy")
print("GL2_F13_set_stabilizer=identity_only")
print("punctured_transfer_degrees=out_(12*156,11*12);in_(12*156,11*12)")
print("full_additive_Cayley_character_representatives=" + repr(character_representatives))
print("full_additive_Cayley_character_histograms=" + repr(tuple(character_histograms)))
print("fourier_orbit_factors_low_to_high=" + repr(fourier_orbit_factors_low))
print("fourier_orbit_factor_count=14_degree_12_factors;distinct=13;duplicate_directions=" + repr(duplicate_directions))
print("full_Cayley_charpoly=(lambda-12)*product_14(F_direction);sha256=" + full_cayley_char_digest)
print("punctured_charpoly=full_Cayley_charpoly_derivative/169")
print("punctured_charpoly_sha256=" + seam_char_digest)
print("punctured_charpoly_first8_high_to_low=" + repr(seam_char_high[:8]))
print("punctured_charpoly_last8_high_to_low=" + repr(seam_char_high[-8:]))
print("punctured_rational_factor_degrees=(12,156);degree12_factor_high_to_low=" + repr(tuple(reversed(duplicate_factor))))
print("degree156_cofactor_sha256=" + seam_cofactor_digest)
print("punctured_charpoly_squarefree=yes;minimal_matrix_recurrence_order=168")
print("phase_matched_alternate_transversal_exponents=" + repr(alternate_exponents))
print("alternate_transversal_also_discrete_but_charpoly_sha256=" + alternate_char_digest)
print("fixed_one_seam_increment_n13=(156,0*11)_characteristic_13_hostile")
print("scope=no_physical_ancestry,current,safety,row_exclusion,or_LRC14_decrement")
print("all_exact_checks=PASS")
