#!/usr/bin/env python3
"""Exact deterministic companion for THM-3462.

Five ordered points have four positive gaps and ten labelled interval-sum
distances.  This companion derives every comparison from those interval
indicators and certifies the complete preorder atlas without using a bounded
integer grid for completeness.

Route 1 enumerates the vertices of every possible strict/zero wall stratum.
After homogeneous rescaling, a feasible stratum is the pointed polyhedron

    x_i >= 1,
    sigma_j h_j(x) >= 1  when sigma_j != 0,
    h_j(x) = 0           when sigma_j = 0.

Every nonempty pointed rational polyhedron has a vertex.  Four independent
active equations at such a vertex are among x_i=1 and h_j=-1,0,1, giving 49
affine equations and C(49,4) exact Cramer systems.

Route 2 independently inserts the 15 walls recursively and decides each
partial sign prescription by exact Fourier--Motzkin elimination.  The two
routes must return the same covector set.  Bounded positive-integer sweeps are
used only afterwards to optimize representatives, never to discover or prove
completeness of the atlas.

Only Python's standard library is used.  The companion performs no writes or
network calls.  Every mathematical gate raises an explicit exception and is
therefore active under ``python -O``.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from json import dumps
from math import comb, factorial, gcd
from pathlib import Path
import ast
import sys


THEOREM_ID = "THM-3462"
STATUS = "PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED"
DEPENDENCY_RELATIVE_PATH = (
    "01-canon/theorems/"
    "THM-3457-four-point-line-metric-turnpike-preorder-atlas.md"
)
EXPECTED_DEPENDENCY_LF_SHA256 = "1e0acb0775d6533ec29309e9835c22021b4d5c9eb41652f0be60c3565b7c01d5"
EXPECTED_COVECTOR_SHA256 = "38c1e9f857ebfe1ab5473bcab15695500fa5de7293cba436924199b4f0024eff"
EXPECTED_SEMANTIC_SHA256 = "90323a5738e47db6643fea438d3a6b7e3d744b6119568f1e35b2bd03b3297d88"

GAP_NAMES = ("a", "b", "c", "d")
EXPECTED_COMPARISON_COUNTS = {"forced": 25, "mixed": 20}
EXPECTED_CRAMER_COUNTS = {
    "affine_equations": 49,
    "systems": 211876,
    "nonsingular": 123706,
    "positive_candidates": 17143,
    "covectors": 477,
}
EXPECTED_FM_DECISIONS = 7053
EXPECTED_TIE_CENSUS = {
    "strict": 114,
    "2": 162,
    "2+2": 124,
    "2+2+2": 36,
    "2+2+2+2": 5,
    "3+2": 14,
    "3+2+2": 16,
    "3+2+2+2": 3,
    "3+3+2": 2,
    "4+3+2": 1,
}
EXPECTED_F_VECTOR = (20, 125, 218, 114)

# Human order for the theorem statement.  The exact algorithms use the
# lexicographically sorted coefficient rows derived from the comparisons.
EXPECTED_WALL_ROWS_HUMAN = (
    (1, -1, 0, 0),
    (1, -1, -1, 0),
    (1, -1, -1, -1),
    (1, 0, -1, 0),
    (1, 0, -1, -1),
    (1, 0, 0, -1),
    (1, 1, -1, 0),
    (1, 1, -1, -1),
    (1, 1, 0, -1),
    (0, 1, -1, 0),
    (0, 1, -1, -1),
    (0, 1, 0, -1),
    (1, 1, 1, -1),
    (0, 1, 1, -1),
    (0, 0, 1, -1),
)

TOTAL_HOSTILE = (6, 5, 4, 8)
HEIGHT_HOSTILE = (2, 3, 4, 10)
REVERSAL = (4, 3, 2, 1, 0)
IDENTITY = (0, 1, 2, 3, 4)


def require(condition, message):
    """Raise an optimization-safe gate failure."""

    if not condition:
        raise RuntimeError(message)


def sign(value):
    """Return -1, 0, or 1 exactly."""

    return (value > 0) - (value < 0)


def dot(left, right):
    return sum((x * y for x, y in zip(left, right)), 0)


def vector_neg(vector):
    return tuple(-entry for entry in vector)


def vector_sub(left, right):
    return tuple(x - y for x, y in zip(left, right))


def primitive_tuple(values):
    common = 0
    for value in values:
        common = gcd(common, abs(value))
    require(common > 0, f"cannot primitive-normalize {values}")
    return tuple(value // common for value in values)


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def lf_sha256(path):
    return sha256(lf_bytes(path)).hexdigest()


def derive_intervals():
    """Derive the ten nonempty contiguous interval indicators."""

    rows = []
    for span in range(1, 5):
        for left in range(0, 5 - span):
            right = left + span
            label = f"{left}{right}"
            coefficients = tuple(1 if left <= gap < right else 0 for gap in range(4))
            rows.append((label, coefficients))
    require(len(rows) == comb(5, 2) == 10, f"interval count drift: {len(rows)}")
    require(len({row for _, row in rows}) == 10, "duplicate interval indicator")
    return tuple(rows)


INTERVALS = derive_intervals()
EDGE_LABELS = tuple(label for label, _ in INTERVALS)
EDGE_COEFFICIENTS = dict(INTERVALS)
EDGE_INDEX = {label: index for index, label in enumerate(EDGE_LABELS)}


def normalize_mixed_wall(row):
    """Orient an interval-difference wall by its first nonzero coordinate."""

    common = 0
    for value in row:
        common = gcd(common, abs(value))
    require(common > 0, f"zero mixed wall: {row}")
    reduced = tuple(value // common for value in row)
    first = next(value for value in reduced if value)
    if first < 0:
        reduced = vector_neg(reduced)
    positive = tuple(index for index, value in enumerate(reduced) if value > 0)
    negative = tuple(index for index, value in enumerate(reduced) if value < 0)
    require(positive and negative, f"wall is not mixed: {reduced}")
    require(max(positive) < min(negative), f"wall blocks are not ordered/disjoint: {reduced}")
    require(
        positive == tuple(range(min(positive), max(positive) + 1)),
        f"positive wall block is not contiguous: {reduced}",
    )
    require(
        negative == tuple(range(min(negative), max(negative) + 1)),
        f"negative wall block is not contiguous: {reduced}",
    )
    require(set(reduced).issubset({-1, 0, 1}), f"non-unit wall row: {reduced}")
    return reduced


def derive_comparisons():
    forced = []
    mixed = []
    walls = set()
    for (left_label, left_row), (right_label, right_row) in combinations(INTERVALS, 2):
        difference = vector_sub(left_row, right_row)
        nonnegative = all(value >= 0 for value in difference)
        nonpositive = all(value <= 0 for value in difference)
        if nonnegative or nonpositive:
            require(nonnegative != nonpositive, f"duplicate distance rows {left_label},{right_label}")
            forced.append(
                (left_label, right_label, 1 if nonnegative else -1, difference)
            )
            continue
        wall = normalize_mixed_wall(difference)
        orientation = 1 if difference == wall else -1
        mixed.append((left_label, right_label, orientation, wall))
        walls.add(wall)

    require(
        {"forced": len(forced), "mixed": len(mixed)} == EXPECTED_COMPARISON_COUNTS,
        f"comparison partition drift: forced={len(forced)}, mixed={len(mixed)}",
    )
    require(walls == set(EXPECTED_WALL_ROWS_HUMAN), f"15-wall set drift: {sorted(walls)}")
    return tuple(forced), tuple(mixed), tuple(sorted(walls))


FORCED_COMPARISONS, MIXED_COMPARISONS, WALLS = derive_comparisons()
WALL_INDEX = {wall: index for index, wall in enumerate(WALLS)}


def block_name(indices):
    return "+".join(GAP_NAMES[index] for index in indices)


def wall_name(row):
    positive = tuple(index for index, value in enumerate(row) if value > 0)
    negative = tuple(index for index, value in enumerate(row) if value < 0)
    return f"{block_name(positive)}={block_name(negative)}"


def determinant3(matrix):
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def determinant4(matrix):
    total = 0
    for column in range(4):
        minor = tuple(
            tuple(matrix[row][other] for other in range(4) if other != column)
            for row in range(1, 4)
        )
        total += (-1 if column % 2 else 1) * matrix[0][column] * determinant3(minor)
    return total


def cramer_numerators(equations):
    matrix = tuple(row for row, _ in equations)
    right = tuple(value for _, value in equations)
    determinant = determinant4(matrix)
    if determinant == 0:
        return None
    numerators = []
    for column in range(4):
        replaced = tuple(
            tuple(right[row] if index == column else matrix[row][index] for index in range(4))
            for row in range(4)
        )
        numerators.append(determinant4(replaced))
    if determinant < 0:
        determinant = -determinant
        numerators = [-value for value in numerators]
    return tuple(numerators), determinant


def covector_from_numerators(numerators):
    return tuple(sign(dot(wall, numerators)) for wall in WALLS)


def covector_from_point(point):
    return tuple(sign(dot(wall, point)) for wall in WALLS)


def affine_equations():
    equations = []
    for index in range(4):
        row = tuple(1 if index == other else 0 for other in range(4))
        equations.append((row, 1))
    for wall in WALLS:
        for value in (-1, 0, 1):
            equations.append((wall, value))
    require(len(equations) == 49, f"affine equation count drift: {len(equations)}")
    return tuple(equations)


AFFINE_EQUATIONS = affine_equations()


def cramer_atlas():
    """Enumerate all global active-equation vertices exactly."""

    systems = 0
    nonsingular = 0
    positive_candidates = 0
    representatives = {}
    for indices in combinations(range(len(AFFINE_EQUATIONS)), 4):
        systems += 1
        solution = cramer_numerators(tuple(AFFINE_EQUATIONS[index] for index in indices))
        if solution is None:
            continue
        nonsingular += 1
        numerators, _ = solution
        if not all(value > 0 for value in numerators):
            continue
        positive_candidates += 1
        covector = covector_from_numerators(numerators)
        candidate = primitive_tuple(numerators)
        old = representatives.get(covector)
        if old is None or (sum(candidate), max(candidate), candidate) < (
            sum(old),
            max(old),
            old,
        ):
            representatives[covector] = candidate

    counts = {
        "affine_equations": len(AFFINE_EQUATIONS),
        "systems": systems,
        "nonsingular": nonsingular,
        "positive_candidates": positive_candidates,
        "covectors": len(representatives),
    }
    require(counts == EXPECTED_CRAMER_COUNTS, f"Cramer census drift: {counts}")
    return frozenset(representatives), representatives, counts


def canonical_inequalities(rows):
    """Normalize positive scalar multiples and retain the strongest duplicate."""

    strongest = {}
    for coefficients, right in rows:
        coefficients = tuple(Fraction(value) for value in coefficients)
        right = Fraction(right)
        first = next((value for value in coefficients if value), None)
        if first is None:
            if right > 0:
                return None
            continue
        scale = abs(first)
        normalized = tuple(value / scale for value in coefficients)
        bound = right / scale
        if normalized not in strongest or bound > strongest[normalized]:
            strongest[normalized] = bound
    return tuple(sorted(strongest.items()))


def eliminate_variable(rows, variable):
    positive = []
    negative = []
    zero = []
    for coefficients, right in rows:
        coefficient = coefficients[variable]
        reduced = coefficients[:variable] + coefficients[variable + 1 :]
        item = (coefficients, reduced, right, coefficient)
        if coefficient > 0:
            positive.append(item)
        elif coefficient < 0:
            negative.append(item)
        else:
            zero.append((reduced, right))

    combined = list(zero)
    for _, positive_reduced, positive_right, positive_coefficient in positive:
        for _, negative_reduced, negative_right, negative_coefficient in negative:
            coefficients = tuple(
                (-negative_coefficient) * left + positive_coefficient * right
                for left, right in zip(positive_reduced, negative_reduced)
            )
            bound = (
                (-negative_coefficient) * positive_right
                + positive_coefficient * negative_right
            )
            combined.append((coefficients, bound))
    return canonical_inequalities(combined)


def fourier_motzkin_feasible(covector_prefix):
    """Decide a partial wall-sign prescription by exact elimination."""

    rows = []
    for index in range(4):
        row = tuple(1 if index == other else 0 for other in range(4))
        rows.append((row, 1))
    for wall, wall_sign in zip(WALLS, covector_prefix):
        if wall_sign > 0:
            rows.append((wall, 1))
        elif wall_sign < 0:
            rows.append((vector_neg(wall), 1))
        else:
            rows.append((wall, 0))
            rows.append((vector_neg(wall), 0))

    canonical = canonical_inequalities(rows)
    if canonical is None:
        return False
    variable_count = 4
    while variable_count:
        # Choosing the smallest positive/negative pairing count changes only
        # efficiency, not the Fourier--Motzkin certificate.
        costs = []
        for variable in range(variable_count):
            positive = sum(coefficients[variable] > 0 for coefficients, _ in canonical)
            negative = sum(coefficients[variable] < 0 for coefficients, _ in canonical)
            costs.append((positive * negative, positive + negative, variable))
        variable = min(costs)[2]
        canonical = eliminate_variable(canonical, variable)
        if canonical is None:
            return False
        variable_count -= 1
    return True


def fourier_motzkin_atlas():
    """Insert walls recursively; this path does not use Cramer candidates."""

    prefixes = [()]
    level_counts = [1]
    decisions = 0
    for _ in WALLS:
        next_prefixes = []
        for prefix in prefixes:
            for wall_sign in (-1, 0, 1):
                decisions += 1
                extended = prefix + (wall_sign,)
                if fourier_motzkin_feasible(extended):
                    next_prefixes.append(extended)
        prefixes = next_prefixes
        level_counts.append(len(prefixes))
    require(decisions == EXPECTED_FM_DECISIONS, f"FM decision count drift: {decisions}")
    require(len(prefixes) == 477, f"FM terminal covector count drift: {len(prefixes)}")
    return frozenset(prefixes), tuple(level_counts), decisions


def distance_signature(point):
    values = [(dot(coefficients, point), label) for label, coefficients in INTERVALS]
    values.sort(key=lambda item: (item[0], EDGE_INDEX[item[1]]))
    blocks = []
    current_value = None
    current_labels = []
    for value, label in values:
        if current_value is None or value == current_value:
            current_value = value
            current_labels.append(label)
            continue
        blocks.append(tuple(current_labels))
        current_value = value
        current_labels = [label]
    blocks.append(tuple(current_labels))
    return tuple(blocks)


def signature_text(signature):
    return "<".join("=".join(block) for block in signature)


def signature_json(signature):
    return [list(block) for block in signature]


def covector_text(covector):
    symbols = {-1: "-", 0: "0", 1: "+"}
    return "".join(symbols[value] for value in covector)


def tie_type(signature):
    sizes = tuple(sorted((len(block) for block in signature if len(block) > 1), reverse=True))
    return "strict" if not sizes else "+".join(str(size) for size in sizes)


def rational_rank(rows):
    matrix = [list(map(Fraction, row)) for row in rows if any(row)]
    if not matrix:
        return 0
    rank = 0
    column = 0
    while rank < len(matrix) and column < len(matrix[0]):
        pivot = next((row for row in range(rank, len(matrix)) if matrix[row][column]), None)
        if pivot is None:
            column += 1
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for row in range(len(matrix)):
            if row == rank or not matrix[row][column]:
                continue
            multiple = matrix[row][column]
            matrix[row] = [
                value - multiple * pivot_entry
                for value, pivot_entry in zip(matrix[row], matrix[rank])
            ]
        rank += 1
        column += 1
    return rank


def edge_under_permutation(label, permutation):
    left = permutation[int(label[0])]
    right = permutation[int(label[1])]
    if left > right:
        left, right = right, left
    return f"{left}{right}"


def permute_signature(signature, permutation):
    return tuple(
        tuple(
            sorted(
                (edge_under_permutation(label, permutation) for label in block),
                key=EDGE_INDEX.__getitem__,
            )
        )
        for block in signature
    )


def block_rank(signature, label):
    for index, block in enumerate(signature):
        if label in block:
            return index
    raise RuntimeError(f"missing edge label {label} from signature")


def atlas_records(covectors, representatives):
    records = []
    signatures = {}
    for index, covector in enumerate(sorted(covectors)):
        representative = representatives[covector]
        require(
            covector_from_point(representative) == covector,
            f"Cramer representative sign drift at {covector}",
        )
        signature = distance_signature(representative)
        require(signature not in signatures, f"two covectors share {signature_text(signature)}")
        signatures[signature] = covector
        zero_rows = [wall for wall, wall_sign in zip(WALLS, covector) if wall_sign == 0]
        dimension = 3 - rational_rank(zero_rows)
        require(0 <= dimension <= 3, f"invalid projective dimension {dimension}")
        records.append(
            {
                "id": f"C{index:03d}",
                "covector": covector,
                "cramer_representative": representative,
                "signature": signature,
                "tie_type": tie_type(signature),
                "dimension": dimension,
            }
        )
    require(len(records) == len(signatures) == 477, "covector/signature bijection failed")
    return records, signatures


def reversal_and_s5_audit(records, signature_to_covector):
    record_by_covector = {record["covector"]: record for record in records}
    covector_to_signature = {record["covector"]: record["signature"] for record in records}
    reversal_map = {}
    for record in records:
        reversed_signature = permute_signature(record["signature"], REVERSAL)
        require(reversed_signature in signature_to_covector, f"reversal left atlas at {record['id']}")
        reversal_map[record["covector"]] = signature_to_covector[reversed_signature]

    fixed = tuple(sorted(covector for covector, image in reversal_map.items() if covector == image))
    require(len(fixed) == 5, f"reversal fixed count drift: {len(fixed)}")
    a_equals_d = WALL_INDEX[(1, 0, 0, -1)]
    b_equals_c = WALL_INDEX[(0, 1, -1, 0)]
    symmetry_slice = tuple(
        sorted(
            covector
            for covector in record_by_covector
            if covector[a_equals_d] == 0 and covector[b_equals_c] == 0
        )
    )
    require(symmetry_slice == fixed, "reversal-fixed signatures are not exactly a=d,b=c")

    ratio_controls = (
        (1, 2, 2, 1),
        (1, 1, 1, 1),
        (3, 2, 2, 3),
        (2, 1, 1, 2),
        (3, 1, 1, 3),
    )
    ratio_covectors = tuple(sorted(covector_from_point(point) for point in ratio_controls))
    require(ratio_covectors == fixed, "five p/q reversal strata drift")

    reversal_orbits = {
        frozenset((covector, reversal_map[covector])) for covector in record_by_covector
    }
    require(len(reversal_orbits) == 241, f"reversal orbit count drift: {len(reversal_orbits)}")
    representatives = tuple(sorted(min(orbit) for orbit in reversal_orbits))
    require(len(representatives) == 241, "reversal representative count drift")

    vertex_permutations = tuple(permutations(range(5)))
    require(len(vertex_permutations) == factorial(5) == 120, "S5 enumeration failed")
    positional_signatures = set(signature_to_covector)

    # The unique diameter and the strict incident chains recover line order
    # from either endpoint.  Exact slice intersections audit that only the
    # identity and full reversal return an S5 relabelling to the positional
    # atlas.
    for record in records:
        signature = record["signature"]
        require(signature[-1] == ("04",), f"diameter is not unique at {record['id']}")
        require(
            [block_rank(signature, label) for label in ("01", "02", "03", "04")]
            == sorted(block_rank(signature, label) for label in ("01", "02", "03", "04")),
            f"left endpoint chain failed at {record['id']}",
        )
        require(
            [block_rank(signature, label) for label in ("34", "24", "14", "04")]
            == sorted(block_rank(signature, label) for label in ("34", "24", "14", "04")),
            f"right endpoint chain failed at {record['id']}",
        )
        slice_returns = {
            permute_signature(signature, permutation)
            for permutation in vertex_permutations
            if permute_signature(signature, permutation) in positional_signatures
        }
        expected = {signature, covector_to_signature[reversal_map[record["covector"]]]}
        require(slice_returns == expected, f"S5 returns beyond reversal at {record['id']}")

    all_labelled = set()
    orbit_sizes = Counter()
    for covector in representatives:
        signature = covector_to_signature[covector]
        orbit = {permute_signature(signature, permutation) for permutation in vertex_permutations}
        stabilizer = tuple(
            permutation
            for permutation in vertex_permutations
            if permute_signature(signature, permutation) == signature
        )
        require(len(orbit) * len(stabilizer) == 120, f"orbit-stabilizer failed at {covector}")
        expected_size = 60 if covector in fixed else 120
        require(len(orbit) == expected_size, f"S5 orbit size drift at {covector}: {len(orbit)}")
        require(all_labelled.isdisjoint(orbit), f"two reversal classes collapse under S5 at {covector}")
        all_labelled.update(orbit)
        orbit_sizes[len(orbit)] += 1

    require(orbit_sizes == Counter({120: 236, 60: 5}), f"S5 orbit census drift: {orbit_sizes}")
    require(len(all_labelled) == 28620, f"labelled signature total drift: {len(all_labelled)}")
    return reversal_map, fixed, representatives, dict(orbit_sizes), len(all_labelled)


def positive_compositions(total):
    for a in range(1, total - 2):
        for b in range(1, total - a - 1):
            for c in range(1, total - a - b):
                d = total - a - b - c
                if d > 0:
                    yield (a, b, c, d)


def representative_audit(covectors):
    """Optimize representatives only after exact completeness is known."""

    minimum_total = {}
    total_coverage = {}
    for total in range(4, 24):
        for point in positive_compositions(total):
            covector = covector_from_point(point)
            require(covector in covectors, f"integer point left certified atlas: {point}")
            if covector not in minimum_total:
                minimum_total[covector] = point
        total_coverage[total] = len(minimum_total)
    require(total_coverage[21] == 475, f"denominator-21 coverage drift: {total_coverage[21]}")
    require(total_coverage[22] == 475, f"denominator-22 coverage drift: {total_coverage[22]}")
    require(total_coverage[23] == 477, f"denominator-23 coverage drift: {total_coverage[23]}")

    minimum_height = {}
    height_coverage = {}
    for height in range(1, 11):
        for point in product(range(1, height + 1), repeat=4):
            if max(point) != height:
                continue
            covector = covector_from_point(point)
            require(covector in covectors, f"height point left certified atlas: {point}")
            if covector not in minimum_height:
                minimum_height[covector] = point
        height_coverage[height] = len(minimum_height)
    require(height_coverage[10] == 477, f"height-10 coverage drift: {height_coverage[10]}")
    require(height_coverage[9] < 477, "height 9 unexpectedly covers the atlas")

    for name, representatives in (("total", minimum_total), ("height", minimum_height)):
        require(set(representatives) == set(covectors), f"{name} representatives incomplete")
        for covector, point in representatives.items():
            require(gcd(gcd(point[0], point[1]), gcd(point[2], point[3])) == 1, f"nonprimitive {name} representative at {covector}: {point}")

    total_covector = covector_from_point(TOTAL_HOSTILE)
    reversed_total_covector = covector_from_point(tuple(reversed(TOTAL_HOSTILE)))
    final_total = {
        covector for covector, point in minimum_total.items() if sum(point) == 23
    }
    require(
        final_total == {total_covector, reversed_total_covector},
        f"final total-23 signatures drift: {final_total}",
    )
    require(minimum_total[total_covector] == TOTAL_HOSTILE, f"total hostile representative drift: {minimum_total[total_covector]}")

    height_covector = covector_from_point(HEIGHT_HOSTILE)
    require(max(minimum_height[height_covector]) == 10, f"height hostile is not sharp: {minimum_height[height_covector]}")

    total_signature = distance_signature(TOTAL_HOSTILE)
    require(
        block_rank(total_signature, "23") < block_rank(total_signature, "12")
        < block_rank(total_signature, "01") < block_rank(total_signature, "34"),
        "total hostile does not force c<b<a<d",
    )
    require(block_rank(total_signature, "34") < block_rank(total_signature, "13"), "total hostile does not force d<b+c")
    require(block_rank(total_signature, "02") < block_rank(total_signature, "24"), "total hostile does not force a+b<c+d")

    height_signature = distance_signature(HEIGHT_HOSTILE)
    require(
        block_rank(height_signature, "01") < block_rank(height_signature, "12")
        < block_rank(height_signature, "23") < block_rank(height_signature, "02"),
        "height hostile does not force a<b<c<a+b",
    )
    require(block_rank(height_signature, "03") < block_rank(height_signature, "34"), "height hostile does not force a+b+c<d")

    return (
        minimum_total,
        minimum_height,
        total_coverage,
        height_coverage,
        total_covector,
        height_covector,
    )


def covector_digest(covectors):
    payload = dumps(
        sorted(covectors),
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    return sha256(payload).hexdigest()


def signature_digest(records):
    payload = (
        "\n".join(
            f"{covector_text(record['covector'])}|{signature_text(record['signature'])}"
            for record in records
        )
        + "\n"
    ).encode("ascii")
    return sha256(payload).hexdigest()


def source_self_audit(source_path):
    source = lf_bytes(source_path).decode("utf-8")
    tree = ast.parse(source, filename=str(source_path))
    allowed_imports = {
        "ast",
        "collections",
        "fractions",
        "hashlib",
        "itertools",
        "json",
        "math",
        "pathlib",
        "sys",
    }
    imports = set()
    banned_names = {"eval", "exec", "compile", "open", "input", "breakpoint", "__import__"}
    banned_attributes = {
        "write",
        "writelines",
        "write_text",
        "write_bytes",
        "touch",
        "mkdir",
        "unlink",
        "rmdir",
        "rename",
    }
    assert_count = 0
    banned_calls = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            require(node.module is not None, "relative import in companion")
            imports.add(node.module.split(".")[0])
        elif isinstance(node, ast.Assert):
            assert_count += 1
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id in banned_names:
                banned_calls.append(node.func.id)
            elif isinstance(node.func, ast.Attribute) and node.func.attr in banned_attributes:
                banned_calls.append(node.func.attr)
    require(imports <= allowed_imports, f"non-stdlib/unapproved import roots: {sorted(imports - allowed_imports)}")
    require(assert_count == 0, f"source contains {assert_count} assert statements")
    require(not banned_calls, f"source contains banned write/dynamic calls: {sorted(banned_calls)}")
    return {
        "imports": tuple(sorted(imports)),
        "asserts": assert_count,
        "banned_calls": tuple(banned_calls),
    }


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    source_path = Path(__file__).resolve()
    repo_root = source_path.parents[1]
    dependency_path = repo_root / DEPENDENCY_RELATIVE_PATH
    require(dependency_path.is_file(), f"missing dependency {DEPENDENCY_RELATIVE_PATH}")
    source_hash = lf_sha256(source_path)
    dependency_hash = lf_sha256(dependency_path)
    if EXPECTED_DEPENDENCY_LF_SHA256 != "UNFROZEN":
        require(
            dependency_hash == EXPECTED_DEPENDENCY_LF_SHA256,
            f"dependency LF hash drift: {dependency_hash}",
        )
    self_audit = source_self_audit(source_path)

    wall_multiplicities = Counter(wall for _, _, _, wall in MIXED_COMPARISONS)
    require(sum(wall_multiplicities.values()) == 20, "mixed comparison multiplicity drift")

    cramer_covectors, cramer_representatives, cramer_counts = cramer_atlas()
    fm_covectors, fm_level_counts, fm_decisions = fourier_motzkin_atlas()
    require(cramer_covectors == fm_covectors, "Cramer and Fourier--Motzkin atlases disagree")

    covector_hash = covector_digest(cramer_covectors)
    if EXPECTED_COVECTOR_SHA256 != "UNFROZEN":
        require(covector_hash == EXPECTED_COVECTOR_SHA256, f"covector digest drift: {covector_hash}")

    records, signature_to_covector = atlas_records(cramer_covectors, cramer_representatives)
    tie_census = dict(Counter(record["tie_type"] for record in records))
    require(tie_census == EXPECTED_TIE_CENSUS, f"tie census drift: {tie_census}")
    face_vector = tuple(Counter(record["dimension"] for record in records)[dimension] for dimension in range(4))
    require(face_vector == EXPECTED_F_VECTOR, f"face-vector drift: {face_vector}")
    require(face_vector[0] - face_vector[1] + face_vector[2] - face_vector[3] == -1, "Euler control failed")

    maximal = tuple(record for record in records if record["tie_type"] == "4+3+2")
    require(len(maximal) == 1, f"maximal tie stratum drift: {len(maximal)}")
    require(covector_from_point((1, 1, 1, 1)) == maximal[0]["covector"], "equal-gap control missed maximal stratum")

    reversal_map, reversal_fixed, reversal_representatives, s5_orbit_sizes, labelled_total = reversal_and_s5_audit(
        records, signature_to_covector
    )
    (
        minimum_total,
        minimum_height,
        total_coverage,
        height_coverage,
        total_covector,
        height_covector,
    ) = representative_audit(cramer_covectors)

    record_by_covector = {record["covector"]: record for record in records}
    for record in records:
        covector = record["covector"]
        record["total_representative"] = minimum_total[covector]
        record["height_representative"] = minimum_height[covector]
        record["reversal_id"] = record_by_covector[reversal_map[covector]]["id"]

    sig_hash = signature_digest(records)
    semantic_payload = {
        "theorem_id": THEOREM_ID,
        "status": STATUS,
        "scope": "five ordered real line points; ten labelled interval-sum distances",
        "completeness_route": "49 affine equations, exact four-by-four Cramer vertices; no bounded discovery scan",
        "independent_route": "recursive exact Fourier-Motzkin insertion of 15 walls",
        "intervals": [
            {"label": label, "coefficients": list(coefficients)}
            for label, coefficients in INTERVALS
        ],
        "comparisons": {
            "forced": len(FORCED_COMPARISONS),
            "mixed": len(MIXED_COMPARISONS),
            "walls_lexicographic": [
                {"name": wall_name(wall), "coefficients": list(wall), "multiplicity": wall_multiplicities[wall]}
                for wall in WALLS
            ],
        },
        "cramer": cramer_counts,
        "fourier_motzkin": {
            "decisions": fm_decisions,
            "level_counts": list(fm_level_counts),
            "terminal_covectors": len(fm_covectors),
        },
        "covector_sha256": covector_hash,
        "signature_sha256": sig_hash,
        "tie_census": tie_census,
        "f_vector": list(face_vector),
        "reversal": {
            "fixed_ids": [record_by_covector[covector]["id"] for covector in reversal_fixed],
            "classes": len(reversal_representatives),
        },
        "s5": {
            "orbit_sizes": {str(size): count for size, count in sorted(s5_orbit_sizes.items())},
            "labelled_total": labelled_total,
        },
        "representatives": {
            "total_coverage_21_22_23": [total_coverage[21], total_coverage[22], total_coverage[23]],
            "height_coverage_9_10": [height_coverage[9], height_coverage[10]],
            "total_hostile": list(TOTAL_HOSTILE),
            "total_hostile_id": record_by_covector[total_covector]["id"],
            "height_hostile": list(HEIGHT_HOSTILE),
            "height_hostile_id": record_by_covector[height_covector]["id"],
        },
        "atlas": [
            {
                "id": record["id"],
                "covector": list(record["covector"]),
                "signature": signature_json(record["signature"]),
                "tie_type": record["tie_type"],
                "dimension": record["dimension"],
                "total_representative": list(record["total_representative"]),
                "height_representative": list(record["height_representative"]),
                "reversal_id": record["reversal_id"],
            }
            for record in records
        ],
    }
    semantic_bytes = dumps(
        semantic_payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    semantic_hash = sha256(semantic_bytes).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, f"semantic digest drift: {semantic_hash}")

    print("THM-3462 EXACT DETERMINISTIC COMPANION")
    print(f"status={STATUS}")
    print("arithmetic=stdlib integer/Fraction only; explicit exceptions remain active under -O")
    print("side_effects=stdout and read-only source/dependency hashing; no writes; no network")
    print(
        "self_audit="
        f"imports:{','.join(self_audit['imports'])};asserts:{self_audit['asserts']};"
        f"banned_calls:{len(self_audit['banned_calls'])}"
    )
    print(f"script_lf_sha256={source_hash}")
    print(f"dependency={DEPENDENCY_RELATIVE_PATH}")
    print(f"dependency_lf_sha256={dependency_hash}")
    print(
        "interval_indicators="
        + ";".join(
            f"{label}:{''.join(str(value) for value in coefficients)}"
            for label, coefficients in INTERVALS
        )
    )
    print("comparison_partition=total:45,forced:25,mixed:20")
    print("WALLS_LEXICOGRAPHIC_BEGIN")
    for index, wall in enumerate(WALLS):
        print(
            f"W{index:02d}|{wall_name(wall)}|coefficients={wall}|"
            f"mixed_pair_multiplicity={wall_multiplicities[wall]}"
        )
    print("WALLS_LEXICOGRAPHIC_END")
    print(
        "cramer_certificate="
        f"affine_equations:{cramer_counts['affine_equations']},systems:{cramer_counts['systems']},"
        f"nonsingular:{cramer_counts['nonsingular']},positive_candidates:{cramer_counts['positive_candidates']},"
        f"covectors:{cramer_counts['covectors']}"
    )
    print(
        "fourier_motzkin_certificate="
        f"decisions:{fm_decisions},terminal_covectors:{len(fm_covectors)},"
        f"level_counts:{','.join(str(value) for value in fm_level_counts)}"
    )
    print(f"covector_sha256={covector_hash}")
    print(f"signature_sha256={sig_hash}")
    print("tie_census=" + ",".join(f"{key}:{tie_census[key]}" for key in EXPECTED_TIE_CENSUS))
    print(f"projective_f_vector={face_vector};euler_compact_support=-1")
    print(
        "reversal="
        f"fixed:5,classes:241,fixed_ids:{','.join(record_by_covector[covector]['id'] for covector in reversal_fixed)}"
    )
    print("s5_no_extra_collapse=236_orbits*120+5_orbits*60=28620")
    print(
        "primitive_total_bound="
        f"coverage21:{total_coverage[21]},coverage22:{total_coverage[22]},coverage23:{total_coverage[23]},"
        f"sharp_hostile:{TOTAL_HOSTILE},id:{record_by_covector[total_covector]['id']}"
    )
    print(
        "primitive_height_bound="
        f"coverage9:{height_coverage[9]},coverage10:{height_coverage[10]},"
        f"sharp_hostile:{HEIGHT_HOSTILE},id:{record_by_covector[height_covector]['id']}"
    )
    print("ATLAS_BEGIN")
    for record in records:
        print(
            f"{record['id']}|covector={covector_text(record['covector'])}|dim={record['dimension']}|"
            f"tie={record['tie_type']}|total_rep={record['total_representative']}|"
            f"height_rep={record['height_representative']}|reverse={record['reversal_id']}|"
            f"signature={signature_text(record['signature'])}"
        )
    print("ATLAS_END")
    print(f"semantic_sha256={semantic_hash}")
    print(
        "VERDICT=PASS: exact 477-preorder atlas, f-vector (20,125,218,114), "
        "241 reversal classes, and 28620 labelled signatures certified"
    )


if __name__ == "__main__":
    main()
