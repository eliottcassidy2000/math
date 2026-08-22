#!/usr/bin/env python3
"""Exact deterministic companion for THM-3457.

The normalized positive gap simplex is cut by five and only five possible
distance-equality forms.  Completeness is certified without a bounded grid:
for each of the 3^5 sign prescriptions, this program forms the closed rational
polytope in a+b+c=1, enumerates every possible vertex as an intersection of
two of the eight boundary lines (three coordinate lines plus five equality
lines), and averages all closure candidates.  Because every closure vertex is
present and every open constraint has nonnegative slack on the closure, that
average realizes the open stratum exactly iff the stratum is nonempty.

Only Python's standard library is used.  Every mathematical gate raises an
explicit exception and remains active under ``python -O``.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations, product
from json import dumps
from math import factorial
from pathlib import Path
import sys


THEOREM_ID = "THM-3457"
DEPENDENCY_RELATIVE_PATH = (
    "01-canon/theorems/"
    "THM-3454-fibonacci-selected-u-spine-farey-lorentz-isometry-and-one-tie-edge-order.md"
)
EXPECTED_THM3454_LF_SHA256 = (
    "2f55187055b8158f5a99fd154df14beec5c2b6e7a0a4f65bc5995d495be7d058"
)
EXPECTED_SEMANTIC_SHA256 = "a04174188df0edc25bef47c3f67dd5309f28e835519e4efbfa26e3df03c2cb0a"

EDGE_LABELS = ("01", "12", "23", "02", "13", "03")
EDGE_INDEX = {label: index for index, label in enumerate(EDGE_LABELS)}
EDGE_COEFFICIENTS = {
    "01": (1, 0, 0),
    "12": (0, 1, 0),
    "23": (0, 0, 1),
    "02": (1, 1, 0),
    "13": (0, 1, 1),
    "03": (1, 1, 1),
}

EQUALITY_FORMS = (
    ("a-b", (1, -1, 0)),
    ("b-c", (0, 1, -1)),
    ("a-c", (1, 0, -1)),
    ("a-b-c", (1, -1, -1)),
    ("c-a-b", (-1, -1, 1)),
)
COORDINATE_FORMS = (
    ("a", (1, 0, 0)),
    ("b", (0, 1, 0)),
    ("c", (0, 0, 1)),
)
BOUNDARY_FORMS = COORDINATE_FORMS + EQUALITY_FORMS
SUM_FORM = (1, 1, 1)

IDENTITY = (0, 1, 2, 3)
REVERSAL = (3, 2, 1, 0)
ENDPOINT_SWAP = (3, 1, 2, 0)
INTERIOR_SWAP = (0, 2, 1, 3)


def blocks(*items):
    """Make a literal signature compactly while retaining frozen edge order."""

    return tuple(tuple(item.split("=")) for item in items)


FROZEN_ATLAS = (
    ("S1", "strict", (1, 2, 4), blocks("01", "12", "02", "23", "13", "03")),
    ("S1r", "strict", (4, 2, 1), blocks("23", "12", "13", "01", "02", "03")),
    ("S2", "strict", (2, 3, 4), blocks("01", "12", "23", "02", "13", "03")),
    ("S2r", "strict", (4, 3, 2), blocks("23", "12", "01", "13", "02", "03")),
    ("S3", "strict", (1, 3, 2), blocks("01", "23", "12", "02", "13", "03")),
    ("S3r", "strict", (2, 3, 1), blocks("23", "01", "12", "13", "02", "03")),
    ("S4", "strict", (2, 1, 4), blocks("12", "01", "02", "23", "13", "03")),
    ("S4r", "strict", (4, 1, 2), blocks("12", "23", "13", "01", "02", "03")),
    ("S5", "strict", (3, 2, 4), blocks("12", "01", "23", "02", "13", "03")),
    ("S5r", "strict", (4, 2, 3), blocks("12", "23", "01", "13", "02", "03")),
    ("T1", "one_tie", (1, 2, 3), blocks("01", "12", "23=02", "13", "03")),
    ("T1r", "one_tie", (3, 2, 1), blocks("23", "12", "01=13", "02", "03")),
    ("T2", "one_tie", (2, 1, 3), blocks("12", "01", "23=02", "13", "03")),
    ("T2r", "one_tie", (3, 1, 2), blocks("12", "23", "01=13", "02", "03")),
    ("T3", "one_tie", (2, 2, 1), blocks("23", "01=12", "13", "02", "03")),
    ("T3r", "one_tie", (1, 2, 2), blocks("01", "12=23", "02", "13", "03")),
    ("T4", "one_tie", (2, 2, 3), blocks("01=12", "23", "02", "13", "03")),
    ("T4r", "one_tie", (3, 2, 2), blocks("12=23", "01", "13", "02", "03")),
    ("T5", "one_tie", (1, 1, 3), blocks("01=12", "02", "23", "13", "03")),
    ("T5r", "one_tie", (3, 1, 1), blocks("12=23", "13", "01", "02", "03")),
    ("D1", "two_tie", (1, 1, 2), blocks("01=12", "23=02", "13", "03")),
    ("D1r", "two_tie", (2, 1, 1), blocks("12=23", "01=13", "02", "03")),
    ("D2", "two_tie", (1, 2, 1), blocks("01=23", "12", "02=13", "03")),
    ("D3", "two_tie", (2, 1, 2), blocks("12", "01=23", "02=13", "03")),
    ("M", "maximal_tie", (1, 1, 1), blocks("01=12=23", "02=13", "03")),
)

EXPECTED_CATEGORY_CENSUS = {
    "strict": 10,
    "one_tie": 10,
    "two_tie": 4,
    "maximal_tie": 1,
}
EXPECTED_REVERSAL_MAP = {
    "S1": "S1r", "S1r": "S1",
    "S2": "S2r", "S2r": "S2",
    "S3": "S3r", "S3r": "S3",
    "S4": "S4r", "S4r": "S4",
    "S5": "S5r", "S5r": "S5",
    "T1": "T1r", "T1r": "T1",
    "T2": "T2r", "T2r": "T2",
    "T3": "T3r", "T3r": "T3",
    "T4": "T4r", "T4r": "T4",
    "T5": "T5r", "T5r": "T5",
    "D1": "D1r", "D1r": "D1",
    "D2": "D2",
    "D3": "D3",
    "M": "M",
}
REVERSAL_REPRESENTATIVE_IDS = (
    "S1", "S2", "S3", "S4", "S5",
    "T1", "T2", "T3", "T4", "T5",
    "D1", "D2", "D3", "M",
)
EXPECTED_FIXED_IDS = ("D2", "D3", "M")


def require(condition, detail):
    if not condition:
        raise RuntimeError(f"{THEOREM_ID} exact companion failure: {detail}")


def lf_sha256(path):
    data = Path(path).read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


def vector_neg(vector):
    return tuple(-entry for entry in vector)


def vector_sub(left, right):
    return tuple(x - y for x, y in zip(left, right))


def dot(coefficients, point):
    return sum(
        (Fraction(coefficient) * coordinate for coefficient, coordinate in zip(coefficients, point)),
        Fraction(0),
    )


def signum(value):
    if value < 0:
        return -1
    if value > 0:
        return 1
    return 0


def normalize_gaps(gaps):
    require(len(gaps) == 3, f"gap triple has length {len(gaps)}")
    require(all(value > 0 for value in gaps), f"nonpositive gap in {gaps}")
    total = sum(gaps)
    return tuple(Fraction(value, total) for value in gaps)


def sign_pattern(point):
    return tuple(signum(dot(coefficients, point)) for _, coefficients in EQUALITY_FORMS)


def distance_values(point):
    a, b, c = point
    return {
        "01": a,
        "12": b,
        "23": c,
        "02": a + b,
        "13": b + c,
        "03": a + b + c,
    }


def distance_signature(point):
    values = distance_values(point)
    levels = {}
    for label in EDGE_LABELS:
        levels.setdefault(values[label], []).append(label)
    return tuple(
        tuple(sorted(levels[value], key=EDGE_INDEX.__getitem__))
        for value in sorted(levels)
    )


def signature_text(signature):
    return " < ".join("=".join(block) for block in signature)


def fraction_text(value):
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def point_text(point):
    return "(" + ",".join(fraction_text(value) for value in point) + ")"


def signs_text(signs):
    symbols = {-1: "-", 0: "0", 1: "+"}
    return "(" + ",".join(symbols[value] for value in signs) + ")"


def determinant3(matrix):
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def solve_intersection(first, second):
    matrix = (SUM_FORM, first, second)
    denominator = determinant3(matrix)
    if denominator == 0:
        return None
    right = (Fraction(1), Fraction(0), Fraction(0))
    coordinates = []
    for column in range(3):
        replaced = tuple(
            tuple(right[row] if index == column else matrix[row][index] for index in range(3))
            for row in range(3)
        )
        coordinates.append(Fraction(determinant3(replaced), denominator))
    point = tuple(coordinates)
    require(sum(point, Fraction(0)) == 1, f"intersection off normalized plane: {point}")
    return point


def arrangement_candidates():
    candidates = set()
    for (_, first), (_, second) in combinations(BOUNDARY_FORMS, 2):
        point = solve_intersection(first, second)
        if point is not None:
            candidates.add(point)
    result = tuple(sorted(candidates))
    for simplex_vertex in (
        (Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1)),
    ):
        require(simplex_vertex in result, f"missing simplex vertex {simplex_vertex}")
    return result


def closure_accepts(point, signs):
    if any(coordinate < 0 for coordinate in point):
        return False
    for sign, (_, coefficients) in zip(signs, EQUALITY_FORMS):
        value = dot(coefficients, point)
        if sign == 0 and value != 0:
            return False
        if sign < 0 and value > 0:
            return False
        if sign > 0 and value < 0:
            return False
    return True


def exact_stratum_witness(signs, candidates):
    require(len(signs) == len(EQUALITY_FORMS), f"bad sign-vector length: {signs}")
    require(all(sign in (-1, 0, 1) for sign in signs), f"bad sign value: {signs}")
    closure_points = tuple(point for point in candidates if closure_accepts(point, signs))
    if not closure_points:
        return None

    count = len(closure_points)
    average = tuple(
        sum((point[index] for point in closure_points), Fraction(0)) / count
        for index in range(3)
    )
    require(closure_accepts(average, signs), f"closure average escaped: {signs}, {average}")

    open_slacks = list(average)
    slack_families = [[point[index] for point in closure_points] for index in range(3)]
    for sign, (_, coefficients) in zip(signs, EQUALITY_FORMS):
        if sign == 0:
            require(dot(coefficients, average) == 0, f"zero wall lost at {signs}")
            continue
        multiplier = Fraction(sign)
        open_slacks.append(multiplier * dot(coefficients, average))
        slack_families.append(
            [multiplier * dot(coefficients, point) for point in closure_points]
        )

    require(all(slack >= 0 for slack in open_slacks), f"negative averaged slack at {signs}")
    if all(slack > 0 for slack in open_slacks):
        require(sign_pattern(average) == signs, f"witness has wrong signs at {signs}: {average}")
        return average

    for slack, family in zip(open_slacks, slack_families):
        if slack == 0:
            require(
                all(value == 0 for value in family),
                f"zero average from nonzero nonnegative closure slacks at {signs}",
            )
    return None


def tie_statistics(signature):
    missing_edges = sum(len(block) * (len(block) - 1) // 2 for block in signature)
    strict_arcs = 15 - missing_edges
    strict_pairs = {
        (left, right)
        for left_index, left_block in enumerate(signature)
        for right_block in signature[left_index + 1 :]
        for left in left_block
        for right in right_block
    }
    extension_count = 0
    for order in permutations(EDGE_LABELS):
        positions = {label: index for index, label in enumerate(order)}
        if all(positions[left] < positions[right] for left, right in strict_pairs):
            extension_count += 1
    factorial_count = 1
    for block in signature:
        factorial_count *= factorial(len(block))
    require(
        extension_count == factorial_count,
        f"linear-extension mismatch for {signature}: {extension_count} != {factorial_count}",
    )
    return missing_edges, strict_arcs, extension_count


def category_from_signature(signature):
    missing_edges, _, _ = tie_statistics(signature)
    categories = {0: "strict", 1: "one_tie", 2: "two_tie", 4: "maximal_tie"}
    require(missing_edges in categories, f"unexpected tie multiplicity {missing_edges}: {signature}")
    return categories[missing_edges]


def transformed_edge(label, vertex_permutation):
    first = vertex_permutation[int(label[0])]
    second = vertex_permutation[int(label[1])]
    low, high = sorted((first, second))
    return f"{low}{high}"


def permute_signature(signature, vertex_permutation):
    require(
        tuple(sorted(vertex_permutation)) == IDENTITY,
        f"not a vertex permutation: {vertex_permutation}",
    )
    return tuple(
        tuple(
            sorted(
                (transformed_edge(label, vertex_permutation) for label in block),
                key=EDGE_INDEX.__getitem__,
            )
        )
        for block in signature
    )


def block_rank(signature, label):
    for rank, block in enumerate(signature):
        if label in block:
            return rank
    raise RuntimeError(f"{THEOREM_ID} edge label {label} missing from {signature}")


def rational_json(point):
    return [[value.numerator, value.denominator] for value in point]


def signature_json(signature):
    return [list(block) for block in signature]


def matrix_multiply(left, right):
    require(bool(left) and bool(right), "empty matrix multiplication")
    inner = len(left[0])
    require(all(len(row) == inner for row in left), "ragged left matrix")
    require(all(len(row) == len(right[0]) for row in right), "ragged right matrix")
    require(inner == len(right), f"matrix shape mismatch: {len(left)}x{inner}, {len(right)}x{len(right[0])}")
    return tuple(
        tuple(
            sum((left[i][k] * right[k][j] for k in range(inner)), Fraction(0))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def matrix_vector_multiply(matrix, vector):
    require(all(len(row) == len(vector) for row in matrix), "matrix/vector shape mismatch")
    return tuple(
        sum((entry * coordinate for entry, coordinate in zip(row, vector)), Fraction(0))
        for row in matrix
    )


def matrix_json(matrix):
    return [rational_json(row) for row in matrix]


def centred_gram_audit(row_id, point):
    a, b, c = point
    coordinates = (Fraction(0), a, a + b, a + b + c)
    require(coordinates[-1] == 1, f"{row_id} coordinate normalization failed: {coordinates}")

    centering = tuple(
        tuple(Fraction(1 if i == j else 0) - Fraction(1, 4) for j in range(4))
        for i in range(4)
    )
    centred = matrix_vector_multiply(centering, coordinates)
    direct_centred = tuple(
        coordinate - sum(coordinates, Fraction(0)) / 4 for coordinate in coordinates
    )
    require(centred == direct_centred, f"{row_id} Hx disagrees with direct centering")
    require(sum(centred, Fraction(0)) == 0, f"{row_id} centred coordinates do not sum to zero")

    squared_distances = tuple(
        tuple((coordinates[i] - coordinates[j]) ** 2 for j in range(4))
        for i in range(4)
    )
    left_product = matrix_multiply(centering, squared_distances)
    double_centred = matrix_multiply(left_product, centering)
    gram_from_distances = tuple(
        tuple(-Fraction(1, 2) * value for value in row)
        for row in double_centred
    )
    gram_outer = tuple(
        tuple(centred[i] * centred[j] for j in range(4))
        for i in range(4)
    )
    require(
        gram_from_distances == gram_outer,
        f"{row_id} exact identity -H D^2 H/2 = (Hx)(Hx)^T failed",
    )

    for i, j in combinations(range(4), 2):
        reconstructed = gram_from_distances[i][i] + gram_from_distances[j][j] - 2 * gram_from_distances[i][j]
        require(
            reconstructed == squared_distances[i][j],
            f"{row_id} squared-distance reconstruction failed at ({i},{j})",
        )

    for first_row, second_row in combinations(range(4), 2):
        for first_column, second_column in combinations(range(4), 2):
            minor = (
                gram_from_distances[first_row][first_column]
                * gram_from_distances[second_row][second_column]
                - gram_from_distances[first_row][second_column]
                * gram_from_distances[second_row][first_column]
            )
            require(minor == 0, f"{row_id} Gram matrix has a nonzero 2x2 minor")
    require(
        sum((gram_from_distances[i][i] for i in range(4)), Fraction(0)) > 0,
        f"{row_id} Gram matrix is zero",
    )
    return {
        "id": row_id,
        "coordinates": coordinates,
        "centred": centred,
        "gram": gram_from_distances,
        "rank": 1,
        "reconstructed_pairs": 6,
    }


def audit_equality_forms():
    form_lookup = {coefficients: name for name, coefficients in EQUALITY_FORMS}
    dynamic_pairs = []
    universal_pairs = []
    for left, right in combinations(EDGE_LABELS, 2):
        difference = vector_sub(EDGE_COEFFICIENTS[left], EDGE_COEFFICIENTS[right])
        if difference in form_lookup:
            dynamic_pairs.append((left, right, form_lookup[difference], 1))
        elif vector_neg(difference) in form_lookup:
            dynamic_pairs.append((left, right, form_lookup[vector_neg(difference)], -1))
        else:
            uniformly_nonnegative = all(value >= 0 for value in difference) and any(
                value > 0 for value in difference
            )
            uniformly_nonpositive = all(value <= 0 for value in difference) and any(
                value < 0 for value in difference
            )
            require(
                uniformly_nonnegative or uniformly_nonpositive,
                f"unclassified distance comparison {left}-{right}: {difference}",
            )
            universal_pairs.append((left, right, 1 if uniformly_nonnegative else -1))

    expected_dynamic = (
        ("01", "12", "a-b", 1),
        ("01", "23", "a-c", 1),
        ("01", "13", "a-b-c", 1),
        ("12", "23", "b-c", 1),
        ("23", "02", "c-a-b", 1),
        ("02", "13", "a-c", 1),
    )
    require(tuple(dynamic_pairs) == expected_dynamic, f"equality-form incidence drift: {dynamic_pairs}")
    require(len(universal_pairs) == 9, f"expected 9 universal comparisons, got {universal_pairs}")
    return tuple(dynamic_pairs), tuple(universal_pairs)


def fibonacci_audit(signature_to_id):
    fibonacci = [0, 1]
    for _ in range(2, 65):
        fibonacci.append(fibonacci[-1] + fibonacci[-2])

    boundary = (fibonacci[1], fibonacci[2], fibonacci[3])
    require(boundary == (1, 1, 2), f"Fibonacci boundary drift: {boundary}")
    boundary_id = signature_to_id[distance_signature(tuple(Fraction(x) for x in boundary))]
    require(boundary_id == "D1", f"k=3 boundary is {boundary_id}, not D1")

    audited = []
    for k in range(4, 63):
        gaps = (fibonacci[k - 2], fibonacci[k - 1], fibonacci[k])
        a, b, c = gaps
        require(c == a + b and 0 < a < b, f"Fibonacci wall failed at k={k}: {gaps}")
        row_id = signature_to_id[distance_signature(tuple(Fraction(x) for x in gaps))]
        require(row_id == "T1", f"Fibonacci k={k} landed in {row_id}")
        pell = b * b - a * b - a * a
        require(pell in (-1, 1), f"Cassini-Pell value failed at k={k}: {pell}")
        audited.append((k, gaps, pell))

    hostile = (1, 3, 4)
    hostile_id = signature_to_id[distance_signature(tuple(Fraction(x) for x in hostile))]
    hostile_pell = hostile[1] ** 2 - hostile[0] * hostile[1] - hostile[0] ** 2
    require(hostile[2] == hostile[0] + hostile[1], f"hostile left wall: {hostile}")
    require(hostile[0] < hostile[1], f"hostile left a<b chamber: {hostile}")
    require(hostile_id == "T1", f"hostile landed in {hostile_id}")
    require(hostile_pell == 5, f"hostile Pell value drift: {hostile_pell}")
    require(hostile_pell not in (-1, 1), "hostile unexpectedly passes Cassini-Pell")
    return boundary, boundary_id, tuple(audited), hostile, hostile_id, hostile_pell


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    repo_root = Path(__file__).resolve().parents[1]
    dependency_path = repo_root / DEPENDENCY_RELATIVE_PATH
    require(dependency_path.is_file(), f"missing dependency {dependency_path}")
    dependency_hash = lf_sha256(dependency_path)
    require(
        dependency_hash == EXPECTED_THM3454_LF_SHA256,
        f"THM-3454 LF hash drift: {dependency_hash}",
    )

    dynamic_pairs, universal_pairs = audit_equality_forms()
    candidates = arrangement_candidates()

    rows = []
    signature_to_id = {}
    sign_to_id = {}
    for row_id, declared_category, gaps, frozen_signature in FROZEN_ATLAS:
        point = normalize_gaps(gaps)
        computed_signature = distance_signature(point)
        require(
            computed_signature == frozen_signature,
            f"frozen signature mismatch for {row_id}: {signature_text(computed_signature)}",
        )
        computed_category = category_from_signature(computed_signature)
        require(
            computed_category == declared_category,
            f"category mismatch for {row_id}: {computed_category} != {declared_category}",
        )
        signs = sign_pattern(point)
        require(computed_signature not in signature_to_id, f"duplicate signature at {row_id}")
        require(signs not in sign_to_id, f"duplicate sign stratum at {row_id}")
        signature_to_id[computed_signature] = row_id
        sign_to_id[signs] = row_id
        missing_edges, strict_arcs, extension_count = tie_statistics(computed_signature)
        rows.append(
            {
                "id": row_id,
                "category": declared_category,
                "gaps": gaps,
                "point": point,
                "signs": signs,
                "signature": computed_signature,
                "missing_edges": missing_edges,
                "strict_arcs": strict_arcs,
                "linear_extensions": extension_count,
            }
        )

    require(len(rows) == 25, f"frozen atlas has {len(rows)} rows")
    category_census = dict(Counter(row["category"] for row in rows))
    require(category_census == EXPECTED_CATEGORY_CENSUS, f"category census drift: {category_census}")
    gram_audits = tuple(centred_gram_audit(row["id"], row["point"]) for row in rows)

    stratum_decisions = []
    feasible_by_sign = {}
    for signs in product((-1, 0, 1), repeat=len(EQUALITY_FORMS)):
        witness = exact_stratum_witness(signs, candidates)
        if witness is None:
            stratum_decisions.append({"signs": signs, "feasible": False})
            continue
        signature = distance_signature(witness)
        require(signature in signature_to_id, f"unfrozen feasible signature {signature_text(signature)}")
        row_id = signature_to_id[signature]
        require(signs == rows[[row["id"] for row in rows].index(row_id)]["signs"], f"sign mismatch at {row_id}")
        feasible_by_sign[signs] = (witness, row_id, signature)
        stratum_decisions.append(
            {"signs": signs, "feasible": True, "witness": witness, "row_id": row_id}
        )

    require(len(stratum_decisions) == 243, f"stratum decision count {len(stratum_decisions)}")
    require(len(feasible_by_sign) == 25, f"feasible stratum count {len(feasible_by_sign)}")
    require(set(feasible_by_sign) == set(sign_to_id), "feasible and frozen sign strata differ")
    for signs, (_, row_id, _) in feasible_by_sign.items():
        require(sign_to_id[signs] == row_id, f"stratum-row mismatch: {signs}, {row_id}")

    id_to_row = {row["id"]: row for row in rows}
    reversal_map = {}
    for row in rows:
        reversed_signature = permute_signature(row["signature"], REVERSAL)
        require(reversed_signature in signature_to_id, f"reversal left atlas at {row['id']}")
        reversed_id = signature_to_id[reversed_signature]
        reversal_map[row["id"]] = reversed_id
        gap_reversal_signature = distance_signature(tuple(reversed(row["point"])))
        require(gap_reversal_signature == reversed_signature, f"gap/label reversal mismatch at {row['id']}")
        expected_reversed_gaps = tuple(reversed(row["gaps"]))
        require(
            id_to_row[reversed_id]["gaps"] == expected_reversed_gaps,
            f"frozen reversed representative mismatch at {row['id']}",
        )
    require(reversal_map == EXPECTED_REVERSAL_MAP, f"reversal map drift: {reversal_map}")
    fixed_ids = tuple(row["id"] for row in rows if reversal_map[row["id"]] == row["id"])
    require(fixed_ids == EXPECTED_FIXED_IDS, f"reversal fixed rows drift: {fixed_ids}")
    reversal_orbits = {
        frozenset((row["id"], reversal_map[row["id"]]))
        for row in rows
    }
    require(len(reversal_orbits) == 14, f"reversal orbit count {len(reversal_orbits)}")
    require(
        all(rep in id_to_row for rep in REVERSAL_REPRESENTATIVE_IDS),
        "missing reversal representative",
    )
    require(
        {frozenset((rep, reversal_map[rep])) for rep in REVERSAL_REPRESENTATIVE_IDS}
        == reversal_orbits,
        "14 representative list does not cover reversal orbits",
    )

    all_vertex_permutations = tuple(permutations(range(4)))
    require(len(all_vertex_permutations) == 24, "S4 enumeration did not have 24 permutations")
    atlas_signatures = set(signature_to_id)
    s4_orbits = {}
    stabilizers = {}
    for representative_id in REVERSAL_REPRESENTATIVE_IDS:
        signature = id_to_row[representative_id]["signature"]
        orbit = {permute_signature(signature, permutation) for permutation in all_vertex_permutations}
        stabilizer = tuple(
            permutation
            for permutation in all_vertex_permutations
            if permute_signature(signature, permutation) == signature
        )
        require(len(orbit) * len(stabilizer) == 24, f"orbit-stabilizer failed at {representative_id}")
        s4_orbits[representative_id] = orbit
        stabilizers[representative_id] = stabilizer

    expected_orbit_sizes = {
        representative_id: (12 if representative_id in EXPECTED_FIXED_IDS else 24)
        for representative_id in REVERSAL_REPRESENTATIVE_IDS
    }
    actual_orbit_sizes = {key: len(value) for key, value in s4_orbits.items()}
    require(actual_orbit_sizes == expected_orbit_sizes, f"S4 orbit-size drift: {actual_orbit_sizes}")

    arbitrary_labelled_signatures = set()
    for representative_id in REVERSAL_REPRESENTATIVE_IDS:
        orbit = s4_orbits[representative_id]
        require(
            arbitrary_labelled_signatures.isdisjoint(orbit),
            f"S4 representative {representative_id} collapses an earlier orbit",
        )
        arbitrary_labelled_signatures.update(orbit)
    require(
        len(arbitrary_labelled_signatures) == 11 * 24 + 3 * 12 == 300,
        f"arbitrary-labelled signature count {len(arbitrary_labelled_signatures)}",
    )

    for row in rows:
        orbit_intersection = {
            candidate for candidate in arbitrary_labelled_signatures
            if candidate in {
                permute_signature(row["signature"], permutation)
                for permutation in all_vertex_permutations
            }
            and candidate in atlas_signatures
        }
        expected_intersection = {
            row["signature"],
            id_to_row[reversal_map[row["id"]]]["signature"],
        }
        require(
            orbit_intersection == expected_intersection,
            f"S4 collapses beyond reversal at {row['id']}: {orbit_intersection}",
        )

    diameter_stabilizer = {
        permutation
        for permutation in all_vertex_permutations
        if {permutation[0], permutation[3]} == {0, 3}
    }
    require(
        diameter_stabilizer == {IDENTITY, REVERSAL, ENDPOINT_SWAP, INTERIOR_SWAP},
        f"diameter stabilizer drift: {diameter_stabilizer}",
    )
    for row in rows:
        signature = row["signature"]
        require(signature[-1] == ("03",), f"non-singleton/non-diameter maximum at {row['id']}")
        require(block_rank(signature, "01") < block_rank(signature, "02"), f"01 !< 02 at {row['id']}")
        require(block_rank(signature, "23") < block_rank(signature, "13"), f"23 !< 13 at {row['id']}")
        for permutation in all_vertex_permutations:
            transformed = permute_signature(signature, permutation)
            if transformed in atlas_signatures:
                require(
                    permutation in (IDENTITY, REVERSAL),
                    f"non-geometric S4 slice return at {row['id']}: {permutation}",
                )
        require(
            permute_signature(signature, ENDPOINT_SWAP) not in atlas_signatures,
            f"endpoint swap returned to slice at {row['id']}",
        )
        require(
            permute_signature(signature, INTERIOR_SWAP) not in atlas_signatures,
            f"interior swap returned to slice at {row['id']}",
        )

    boundary, boundary_id, fib_audited, hostile, hostile_id, hostile_pell = fibonacci_audit(signature_to_id)

    semantic_payload = {
        "theorem_id": THEOREM_ID,
        "scope": "FINITE-EXACT companion; theorem is PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED",
        "dependency": {
            "path": DEPENDENCY_RELATIVE_PATH,
            "lf_sha256": dependency_hash,
        },
        "edge_labels": list(EDGE_LABELS),
        "equality_forms": [
            {"name": name, "coefficients": list(coefficients)}
            for name, coefficients in EQUALITY_FORMS
        ],
        "dynamic_pairs": [list(item) for item in dynamic_pairs],
        "universal_pairs": [list(item) for item in universal_pairs],
        "arrangement_candidate_count": len(candidates),
        "atlas": [
            {
                "id": row["id"],
                "category": row["category"],
                "gaps": list(row["gaps"]),
                "normalized": rational_json(row["point"]),
                "signs": list(row["signs"]),
                "signature": signature_json(row["signature"]),
                "reversal": reversal_map[row["id"]],
                "missing_edges": row["missing_edges"],
                "strict_arcs": row["strict_arcs"],
                "linear_extensions": row["linear_extensions"],
            }
            for row in rows
        ],
        "category_census": EXPECTED_CATEGORY_CENSUS,
        "centred_rank_one_gram": {
            "identity": "G=-H*D_squared*H/2=(Hx)(Hx)^T",
            "distance_reconstruction": "d_ij^2=G_ii+G_jj-2G_ij",
            "arxiv_2505_09814_role": "typing boundary for XX^T only; not an atlas dependency or proof",
            "rows": [
                {
                    "id": audit["id"],
                    "coordinates": rational_json(audit["coordinates"]),
                    "centred": rational_json(audit["centred"]),
                    "gram": matrix_json(audit["gram"]),
                    "rank": audit["rank"],
                    "reconstructed_pairs": audit["reconstructed_pairs"],
                }
                for audit in gram_audits
            ],
        },
        "sign_strata": [
            {
                "signs": list(decision["signs"]),
                "feasible": decision["feasible"],
                **(
                    {
                        "witness": rational_json(decision["witness"]),
                        "row_id": decision["row_id"],
                    }
                    if decision["feasible"]
                    else {}
                ),
            }
            for decision in stratum_decisions
        ],
        "reversal": {
            "representatives": list(REVERSAL_REPRESENTATIVE_IDS),
            "fixed": list(fixed_ids),
            "orbit_count": len(reversal_orbits),
        },
        "s4": {
            "representatives": [
                {
                    "id": representative_id,
                    "orbit_size": len(s4_orbits[representative_id]),
                    "stabilizer_size": len(stabilizers[representative_id]),
                }
                for representative_id in REVERSAL_REPRESENTATIVE_IDS
            ],
            "arbitrary_labelled_signature_count": len(arbitrary_labelled_signatures),
        },
        "fibonacci": {
            "symbolic_wall": "(F[k-2],F[k-1],F[k]); c=a+b and a<b for k>=4",
            "boundary_k3": {"gaps": list(boundary), "row_id": boundary_id},
            "exact_control_range": [fib_audited[0][0], fib_audited[-1][0]],
            "hostile": {
                "gaps": list(hostile),
                "row_id": hostile_id,
                "cassini_pell": hostile_pell,
            },
        },
    }
    semantic_bytes = dumps(
        semantic_payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("ascii")
    semantic_hash = sha256(semantic_bytes).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(
            semantic_hash == EXPECTED_SEMANTIC_SHA256,
            f"semantic hash drift: {semantic_hash}",
        )

    print("THM-3457 EXACT DETERMINISTIC COMPANION")
    print("status=PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED")
    print("arithmetic=stdlib Fraction/integer only; every gate uses explicit exceptions")
    print(f"dependency={DEPENDENCY_RELATIVE_PATH}")
    print(f"dependency_lf_sha256={dependency_hash}")
    print("edge_labels=" + ",".join(EDGE_LABELS))
    print("equality_forms=" + ",".join(name for name, _ in EQUALITY_FORMS))
    print(f"comparison_partition=dynamic:{len(dynamic_pairs)},universal:{len(universal_pairs)}")
    print(
        "strata_certificate="
        f"sign_vectors:243,candidate_intersections:{len(candidates)},feasible:25,infeasible:218"
    )
    print("atlas_census=strict:10,one_tie:10,two_tie:4,maximal_tie:1")
    print("ATLAS_BEGIN")
    for row in rows:
        print(
            f"{row['id']}|{row['category']}|gaps={row['gaps']}|"
            f"normalized={point_text(row['point'])}|signs={signs_text(row['signs'])}|"
            f"signature={signature_text(row['signature'])}|reverse={reversal_map[row['id']]}|"
            f"missing={row['missing_edges']}|arcs={row['strict_arcs']}|"
            f"extensions={row['linear_extensions']}"
        )
    print("ATLAS_END")
    print("reversal_orbits=14;fixed=3;category_orbits=strict:5,one_tie:5,two_tie:3,maximal_tie:1")
    print("s4_no_further_collapse=proved_by_exact_slice_intersections_and_diameter_stabilizer")
    print("s4_arbitrary_labelled_signatures=300=11*24+3*12")
    print("S4_REPRESENTATIVES_BEGIN")
    for representative_id in REVERSAL_REPRESENTATIVE_IDS:
        row = id_to_row[representative_id]
        print(
            f"{representative_id}|orbit={len(s4_orbits[representative_id])}|"
            f"stabilizer={len(stabilizers[representative_id])}|"
            f"signature={signature_text(row['signature'])}"
        )
    print("S4_REPRESENTATIVES_END")
    print("tie_metrics=strict:(missing0,arcs15,extensions1);one_tie:(1,14,2);two_tie:(2,13,4);maximal_tie:(4,11,12)")
    print("centred_gram_gate=25/25:rank1;G=-H*D_squared*H/2=(Hx)(Hx)^T;distance_pairs=150/150")
    print("arxiv_2505_09814_boundary=XX^T_typing_only;not_an_atlas_dependency_or_proof")
    print("fibonacci_wall=k>=4:gaps=(F[k-2],F[k-1],F[k]),c=a+b,a<b,row=T1")
    print("fibonacci_boundary=k=3:gaps=(1,1,2),row=D1")
    print("non_fibonacci_hostile=gaps=(1,3,4),row=T1,cassini_pell=5")
    print(f"semantic_sha256={semantic_hash}")
    print("VERDICT=PASS: exact 25-row atlas, 14 reversal/S4 classes, and 300 labelled signatures certified")


if __name__ == "__main__":
    main()
