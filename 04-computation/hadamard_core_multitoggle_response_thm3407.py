#!/usr/bin/env python3
"""Exact finite controls for THM-3407.

Standard library only; integer/Fraction arithmetic throughout.  The program
constructs the Paley-I controls independently and compares three determinant
paths: direct toggling, the event matrix, and both support compressions.  It
also freezes the first genuinely oriented three-event hostiles and pins the
THM-3394 order-668 certificate without dynamically importing its renderer.

The universal determinant, trade-floor, and rank-one rectangle statements are
proved in the theorem.  These computations are finite controls, not substitutes
for those proofs.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_SEMANTIC_SHA256 = "6333b6b8292916e6917b9d04c82bb7b4bf86d6fe2c9dfda45acfce5e3eb98a09"

DEPENDENCY_PINS = (
    (
        "01-canon/theorems/THM-3403-hadamard-core-maxdet-smith-and-circuit-descent.md",
        "645abd68beb3c4e461f51de99b0690cd5cbb0f18164b44ea6f416bf22a56fa44",
    ),
    (
        "04-computation/hadamard_core_descent_thm3403.py",
        "4f75f8be55840747790d29f54bbf9d202a445dcbaf09c09bca38f385d01c4894",
    ),
    (
        "05-knowledge/results/hadamard_core_descent_thm3403.out",
        "0a428ad558077b95f76724cd51f68d0a50a06b14a66400d50a4b37b1240e6322",
    ),
    (
        "01-canon/theorems/THM-3394-twelve-formerly-missing-hadamard-orders-through-2000.md",
        "63e3e841b609f2e9a38bbaa13dbbd665d6c05d7594fdb92351653207c3cadf86",
    ),
    (
        "04-computation/hadamard_twelve_order_bank_thm3394.py",
        "7ae931b3cf268550287bd0621b9b85b8ea167126fadfb90d57b5106d0f82fb2d",
    ),
    (
        "04-computation/hadamard_twelve_order_signword_thm3394.b85",
        "68f7ceebb67005bf1b968171f7e6897cc33bde68adbd63f14bd45edfeb7b3f06",
    ),
    (
        "05-knowledge/results/hadamard_twelve_order_bank_thm3394.out",
        "d8efee90947015a7e6fc28a1685cc3d378357a85e1d4814953b32b17c5cd76a9",
    ),
)

ORDER668_NORMALIZED_SHA256 = (
    "73f1de1539849e1dc7e6085cc69c563fd2965c44970263e8203384bd1a46aa63"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def matmul(left, right):
    right_t = transpose(right)
    return [
        [sum(x * y for x, y in zip(row, column)) for column in right_t]
        for row in left
    ]


def scalar_eye(size, scalar):
    return [
        [scalar if row == column else 0 for column in range(size)]
        for row in range(size)
    ]


def determinant_bareiss(matrix):
    """Fraction-free exact determinant, including the empty determinant."""
    work = [row[:] for row in matrix]
    size = len(work)
    require(all(len(row) == size for row in work), "determinant input is not square")
    if size == 0:
        return 1
    sign = 1
    previous = 1
    for column in range(size - 1):
        pivot_row = next(
            (row for row in range(column, size) if work[row][column]), None
        )
        if pivot_row is None:
            return 0
        if pivot_row != column:
            work[column], work[pivot_row] = work[pivot_row], work[column]
            sign = -sign
        pivot = work[column][column]
        for row in range(column + 1, size):
            for index in range(column + 1, size):
                numerator = (
                    work[row][index] * pivot
                    - work[row][column] * work[column][index]
                )
                require(
                    numerator % previous == 0,
                    "Bareiss division lost exactness",
                )
                work[row][index] = numerator // previous
        previous = pivot
        for row in range(column + 1, size):
            work[row][column] = 0
    return sign * work[-1][-1]


def is_prime(number):
    if number < 2:
        return False
    if number % 2 == 0:
        return number == 2
    divisor = 3
    while divisor * divisor <= number:
        if number % divisor == 0:
            return False
        divisor += 2
    return True


def paley_type_i(prime):
    require(
        is_prime(prime) and prime % 4 == 3,
        "Paley-I control requires a prime congruent to three modulo four",
    )

    def character(value):
        value %= prime
        if value == 0:
            return 0
        return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1

    jacobsthal = [
        [character(row - column) for column in range(prime)]
        for row in range(prime)
    ]
    core = [
        [
            jacobsthal[row][column] - (1 if row == column else 0)
            for column in range(prime)
        ]
        for row in range(prime)
    ]
    return [[1] * (prime + 1)] + [[1] + row for row in core]


def core_data(hadamard):
    order = len(hadamard)
    require(all(len(row) == order for row in hadamard), "Hadamard control is ragged")
    require(all(value == 1 for value in hadamard[0]), "top row is not normalized")
    require(all(row[0] == 1 for row in hadamard), "left column is not normalized")
    sign_core = [row[1:] for row in hadamard[1:]]
    binary_core = [[(1 - value) // 2 for value in row] for row in sign_core]
    return sign_core, binary_core


def check_hadamard_and_core(hadamard):
    order = len(hadamard)
    m = order // 4
    sign_core, binary_core = core_data(hadamard)
    size = order - 1
    require(
        matmul(hadamard, transpose(hadamard)) == scalar_eye(order, order),
        "Hadamard Gram identity failed",
    )
    target = [
        [m * (1 + (1 if row == column else 0)) for column in range(size)]
        for row in range(size)
    ]
    require(
        matmul(binary_core, transpose(binary_core)) == target,
        "binary core row Gram identity failed",
    )
    require(
        matmul(transpose(binary_core), binary_core) == target,
        "binary core column Gram identity failed",
    )
    determinant = determinant_bareiss(binary_core)
    require(abs(determinant) == 2 * m ** (2 * m), "binary maxdet control failed")
    return sign_core, binary_core, determinant


def toggled(binary_core, events):
    changed = [row[:] for row in binary_core]
    for row, column in events:
        changed[row][column] ^= 1
    return changed


def direct_response(binary_core, determinant, events):
    return Fraction(determinant_bareiss(toggled(binary_core, events)), determinant)


def event_matrix(sign_core, events):
    return [
        [
            sign_core[row_a][column_a] * sign_core[row_b][column_a]
            for row_b, _ in events
        ]
        for row_a, column_a in events
    ]


def event_response(sign_core, m, events):
    size = len(events)
    q_matrix = event_matrix(sign_core, events)
    numerator = [
        [
            2 * m * (1 if row == column else 0) - q_matrix[row][column]
            for column in range(size)
        ]
        for row in range(size)
    ]
    return Fraction(determinant_bareiss(numerator), (2 * m) ** size)


def support_responses(sign_core, m, events):
    if not events:
        return Fraction(1), Fraction(1)
    rows = sorted({row for row, _ in events})
    columns = sorted({column for _, column in events})
    row_at = {value: index for index, value in enumerate(rows)}
    column_at = {value: index for index, value in enumerate(columns)}
    l_matrix = [
        [sign_core[row][column] for column in columns]
        for row in rows
    ]
    d_matrix = [[0 for _ in columns] for _ in rows]
    for row, column in events:
        d_matrix[row_at[row]][column_at[column]] = sign_core[row][column]
    dlt = matmul(d_matrix, transpose(l_matrix))
    ltd = matmul(transpose(l_matrix), d_matrix)
    row_numerator = [
        [
            2 * m * (1 if row == column else 0) - dlt[row][column]
            for column in range(len(rows))
        ]
        for row in range(len(rows))
    ]
    column_numerator = [
        [
            2 * m * (1 if row == column else 0) - ltd[row][column]
            for column in range(len(columns))
        ]
        for row in range(len(columns))
    ]
    return (
        Fraction(determinant_bareiss(row_numerator), (2 * m) ** len(rows)),
        Fraction(
            determinant_bareiss(column_numerator),
            (2 * m) ** len(columns),
        ),
    )


def signed_minor(sign_core, events):
    if not events:
        return 1
    rows = [row for row, _ in events]
    columns = [column for _, column in events]
    sign = 1
    for row, column in events:
        sign *= sign_core[row][column]
    minor = [
        [sign_core[row][column] for column in columns]
        for row in rows
    ]
    return sign * determinant_bareiss(minor)


def pair_signature(sign_core, events):
    q_matrix = event_matrix(sign_core, events)
    return tuple(
        q_matrix[left][right] * q_matrix[right][left]
        for left, right in combinations(range(len(events)), 2)
    )


def submask_iter(mask):
    current = mask
    while True:
        yield current
        if current == 0:
            break
        current = (current - 1) & mask


def mask_events(candidates, mask):
    return tuple(
        event for index, event in enumerate(candidates) if mask & (1 << index)
    )


def window_audit(order, sign_core, binary_core, determinant, m):
    candidates = tuple((row, column) for row in range(3) for column in range(3))
    values = {}
    for mask in range(1 << len(candidates)):
        events = mask_events(candidates, mask)
        direct = direct_response(binary_core, determinant, events)
        event = event_response(sign_core, m, events)
        support_row, support_column = support_responses(sign_core, m, events)
        require(
            direct == event == support_row == support_column,
            ("three-path response mismatch", order, mask),
        )
        values[mask] = direct

    nonzero_by_degree = [0] * (len(candidates) + 1)
    matching_nonzero_by_degree = [0] * (len(candidates) + 1)
    for mask in range(1 << len(candidates)):
        degree = mask.bit_count()
        coefficient = sum(
            (-1 if (degree - submask.bit_count()) % 2 else 1) * values[submask]
            for submask in submask_iter(mask)
        )
        events = mask_events(candidates, mask)
        predicted = Fraction(((-1) ** degree) * signed_minor(sign_core, events), (2 * m) ** degree)
        require(coefficient == predicted, ("Boolean Mobius mismatch", order, mask))
        rows = [row for row, _ in events]
        columns = [column for _, column in events]
        is_matching = len(rows) == len(set(rows)) and len(columns) == len(set(columns))
        if coefficient:
            nonzero_by_degree[degree] += 1
            require(is_matching, ("nonmatching interaction survived", order, mask))
            matching_nonzero_by_degree[degree] += 1

    equality_sizes = {}
    for mask, response in values.items():
        if mask and abs(response) == 1:
            equality_sizes[mask.bit_count()] = equality_sizes.get(mask.bit_count(), 0) + 1
    if order == 4:
        require(set(equality_sizes) == {4, 6}, "order-four equality boundary changed")
        require(min(equality_sizes) == order, "order-four trade floor changed")
    else:
        require(not equality_sizes, ("fixed window acquired an equality trade", order))
    return {
        "nonzero_mobius_by_degree": nonzero_by_degree,
        "equality_sizes": sorted(equality_sizes.items()),
    }


def pair_kind(sign_core, first, second):
    row, column = first
    other_row, other_column = second
    if row == other_row:
        return "shared_row", 1
    if column == other_column:
        return "shared_column", 1
    chi = (
        sign_core[row][column]
        * sign_core[other_row][other_column]
        * sign_core[other_row][column]
        * sign_core[row][other_column]
    )
    return ("positive_plaquette" if chi == 1 else "negative_plaquette"), chi


def two_toggle_audit(order, sign_core, binary_core, determinant, m):
    size = order - 1
    positions = tuple((row, column) for row in range(size) for column in range(size))
    counts = {
        "shared_row": 0,
        "shared_column": 0,
        "positive_plaquette": 0,
        "negative_plaquette": 0,
    }
    representatives = {}
    direct_checks = 0
    low = Fraction(m - 1, m)
    high = Fraction(2 * m * m - 2 * m + 1, 2 * m * m)
    for first, second in combinations(positions, 2):
        kind, chi = pair_kind(sign_core, first, second)
        counts[kind] += 1
        representatives.setdefault(kind, (first, second))
        expected = low if chi == 1 else high
        events = (first, second)
        event = event_response(sign_core, m, events)
        support_row, support_column = support_responses(sign_core, m, events)
        require(
            event == support_row == support_column == expected,
            ("two-toggle compiler mismatch", order, first, second),
        )
        if order <= 12:
            require(
                direct_response(binary_core, determinant, events) == expected,
                ("two-toggle direct mismatch", order, first, second),
            )
            direct_checks += 1

    if order == 20:
        for kind, events in representatives.items():
            expected = low if kind != "negative_plaquette" else high
            require(
                direct_response(binary_core, determinant, events) == expected,
                ("order-twenty representative mismatch", kind),
            )
            direct_checks += 1

    row_or_column = size * comb(size, 2)
    choose_rows = comb(size, 2)
    expected_counts = {
        "shared_row": row_or_column,
        "shared_column": row_or_column,
        "positive_plaquette": 2 * choose_rows * (2 * m - 1) ** 2,
        "negative_plaquette": 2 * choose_rows * 2 * m * (2 * m - 1),
    }
    require(counts == expected_counts, ("two-toggle multiplicity mismatch", order))
    require(sum(counts.values()) == comb(size * size, 2), "two-toggle universe mismatch")

    maximum = abs(determinant)
    low_determinant = 2 * (m - 1) * m ** (2 * m - 1)
    high_determinant = (2 * m * m - 2 * m + 1) * m ** (2 * m - 2)
    require(low * maximum == low_determinant, "low shell magnitude mismatch")
    require(high * maximum == high_determinant, "high shell magnitude mismatch")
    require(high_determinant - low_determinant == m ** (2 * m - 2), "shell gap mismatch")
    if m % 2:
        require(low_determinant % 4 == 0, "odd-quarter low shell lost mod-four divisibility")
        require(high_determinant % 2 == 1, "odd-quarter high shell lost oddness")
    return {
        "counts": counts,
        "direct_checks": direct_checks,
        "low": str(low),
        "high": str(high),
    }


def sign_rank_one(rectangle):
    if not rectangle or not rectangle[0]:
        return True
    return all(
        rectangle[0][0] * rectangle[row][column]
        == rectangle[row][0] * rectangle[0][column]
        for row in range(len(rectangle))
        for column in range(len(rectangle[0]))
    )


def rectangle_audit(sign_core, binary_core, determinant, m, rows, columns, label):
    rectangle = [[sign_core[row][column] for column in columns] for row in rows]
    require(sign_rank_one(rectangle), ("rectangle is not rank one", label))
    events = tuple((row, column) for row in rows for column in columns)
    expected = Fraction(1) - Fraction(len(rows) * len(columns), 2 * m)
    direct = direct_response(binary_core, determinant, events)
    event = event_response(sign_core, m, events)
    support_row, support_column = support_responses(sign_core, m, events)
    require(
        direct == event == support_row == support_column == expected,
        ("rank-one rectangle law failed", label),
    )
    return str(expected), len(events)


def trade_and_rectangle_audit(order, hadamard, sign_core, binary_core, determinant, m):
    size = order - 1
    for count in range(size + 1):
        events = tuple((0, column) for column in range(count))
        expected = Fraction(1) - Fraction(count, 2 * m)
        require(
            direct_response(binary_core, determinant, events) == expected,
            ("same-row tariff failed", order, count),
        )
        column_events = tuple((row, 0) for row in range(count))
        require(
            direct_response(binary_core, determinant, column_events) == expected,
            ("same-column tariff failed", order, count),
        )

    row_pair = (0, 1)
    difference_columns = tuple(
        column
        for column in range(size)
        if sign_core[row_pair[0]][column] != sign_core[row_pair[1]][column]
    )
    require(len(difference_columns) == 2 * m, "row difference set has wrong size")
    row_response, row_trade_size = rectangle_audit(
        sign_core,
        binary_core,
        determinant,
        m,
        row_pair,
        difference_columns,
        (order, "row_swap"),
    )
    require((row_response, row_trade_size) == ("-1", order), "row-swap trade changed")

    column_pair = (0, 1)
    difference_rows = tuple(
        row
        for row in range(size)
        if sign_core[row][column_pair[0]] != sign_core[row][column_pair[1]]
    )
    require(len(difference_rows) == 2 * m, "column difference set has wrong size")
    column_response, column_trade_size = rectangle_audit(
        sign_core,
        binary_core,
        determinant,
        m,
        difference_rows,
        column_pair,
        (order, "column_swap"),
    )
    require(
        (column_response, column_trade_size) == ("-1", order),
        "column-swap trade changed",
    )

    if order == 8:
        response, size_eight = rectangle_audit(
            sign_core,
            binary_core,
            determinant,
            m,
            (1, 2, 3, 6),
            (1, 4),
            (8, "closed_quadruple_field"),
        )
        require((response, size_eight) == ("-1", 8), "closed-quadruple trade changed")

    pair_signatures = {
        tuple(left * right for left, right in zip(hadamard[first], hadamard[second]))
        for first, second in combinations(range(order), 2)
    }
    column_vectors = transpose(hadamard)
    column_pair_signatures = {
        tuple(left * right for left, right in zip(column_vectors[first], column_vectors[second]))
        for first, second in combinations(range(order), 2)
    }
    if m % 2 and m > 1:
        require(
            len(pair_signatures) == comb(order, 2),
            ("row pair-product collision in odd-quarter control", order),
        )
        require(
            len(column_pair_signatures) == comb(order, 2),
            ("column pair-product collision in odd-quarter control", order),
        )
    return {
        "same_line_counts": size + 1,
        "row_swap_size": row_trade_size,
        "column_swap_size": column_trade_size,
        "row_pair_signature_count": len(pair_signatures),
        "column_pair_signature_count": len(column_pair_signatures),
    }


def abstract_triple_palette_audit():
    """Exhaust the 64 oriented sign matrices with diagonal one."""
    counts = {}
    allowed = {
        0: {0, -4},
        1: {0},
        2: {0, 4},
        3: {4},
    }
    for values in product((-1, 1), repeat=6):
        q_matrix = [
            [1, values[0], values[1]],
            [values[2], 1, values[3]],
            [values[4], values[5], 1],
        ]
        chis = (
            q_matrix[0][1] * q_matrix[1][0],
            q_matrix[0][2] * q_matrix[2][0],
            q_matrix[1][2] * q_matrix[2][1],
        )
        gamma = q_matrix[0][1] * q_matrix[1][2] * q_matrix[2][0]
        gamma_reverse = q_matrix[0][2] * q_matrix[2][1] * q_matrix[1][0]
        require(
            gamma * gamma_reverse == chis[0] * chis[1] * chis[2],
            "triple cycle-product identity failed",
        )
        negative_pairs = sum(value == -1 for value in chis)
        minor = determinant_bareiss(q_matrix)
        require(minor in allowed[negative_pairs], "triple palette escaped")
        key = "%d:%+d:%+d" % (negative_pairs, gamma, minor)
        counts[key] = counts.get(key, 0) + 1
    expected = {
        "0:-1:-4": 4,
        "0:+1:+0": 4,
        "1:-1:+0": 12,
        "1:+1:+0": 12,
        "2:-1:+0": 12,
        "2:+1:+4": 12,
        "3:-1:+4": 4,
        "3:+1:+4": 4,
    }
    require(counts == expected, "abstract triple palette multiplicities changed")
    print(
        "TRIPLE_PALETTE offdiagonal_states=64 "
        "minor_by_negative_pairs=0:{-4,0};1:{0};2:{0,4};3:{4} "
        "orientation_needed_exactly_for_even_pair_parity"
    )
    return counts


TRIPLE_HOSTILES = {
    8: (
        (((0, 0), (1, 2), (2, 1)), ((0, 0), (1, 1), (2, 3))),
        (-1, 1, -1),
        (4, 0),
        (Fraction(7, 16), Fraction(1, 2)),
        (14, 16),
    ),
    12: (
        (((0, 2), (1, 1), (2, 0)), ((0, 1), (1, 0), (2, 3))),
        (-1, -1, 1),
        (4, 0),
        (Fraction(16, 27), Fraction(11, 18)),
        (864, 891),
    ),
    20: (
        (((0, 2), (1, 1), (2, 0)), ((0, 0), (1, 1), (2, 3))),
        (-1, -1, 1),
        (4, 0),
        (Fraction(92, 125), Fraction(37, 50)),
        (14375000, 14453125),
    ),
}


def triple_hostile_audit(order, sign_core, binary_core, determinant, m):
    event_pair, expected_pairs, expected_minors, expected_responses, expected_determinants = (
        TRIPLE_HOSTILES[order]
    )
    actual_pairs = tuple(pair_signature(sign_core, events) for events in event_pair)
    actual_minors = tuple(signed_minor(sign_core, events) for events in event_pair)
    actual_responses = tuple(event_response(sign_core, m, events) for events in event_pair)
    direct_responses = tuple(
        direct_response(binary_core, determinant, events) for events in event_pair
    )
    actual_determinants = tuple(
        abs(determinant_bareiss(toggled(binary_core, events))) for events in event_pair
    )
    require(actual_pairs == (expected_pairs, expected_pairs), "triple pair data changed")
    require(actual_minors == expected_minors, "triple signed minors changed")
    require(actual_responses == expected_responses, "triple event responses changed")
    require(direct_responses == expected_responses, "triple direct responses changed")
    require(actual_determinants == expected_determinants, "triple determinants changed")
    return {
        "events": event_pair,
        "pair_data": expected_pairs,
        "signed_minors": expected_minors,
        "responses": tuple(map(str, expected_responses)),
        "determinants": expected_determinants,
    }


def dependency_and_order668_audit():
    for relative, expected in DEPENDENCY_PINS:
        path = ROOT / relative
        require(path.is_file(), "missing dependency " + relative)
        require(lf_hash(path) == expected, "dependency hash drift " + relative)
        print("PIN " + relative + " " + expected)

    output = (ROOT / "05-knowledge/results/hadamard_twelve_order_bank_thm3394.out").read_text()
    required_fragment = (
        "shape=668x668\n"
        "   contract period=166 rowsums=2,0,0,0 offpeak_PAF_sum=-4\n"
        "   raw_sha256=bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0 "
        "normalized_sha256=" + ORDER668_NORMALIZED_SHA256
        + " pairs=222778 hamming_min=334 hamming_max=334"
    )
    require(required_fragment in output, "order-668 output pin is absent")

    m = 167
    size = 667
    low_count = 2 * comb(size, 2) * (size + (2 * m - 1) ** 2)
    high_count = 2 * comb(size, 2) * 2 * m * (2 * m - 1)
    require(
        (low_count, high_count)
        == (49_555_629_432, 49_407_259_284),
        "order-668 shell counts changed",
    )
    require(low_count + high_count == comb(size * size, 2), "order-668 pair total changed")
    data = {
        "normalized_sha256": ORDER668_NORMALIZED_SHA256,
        "pair_product_plus_minus": (333, 334),
        "one_response": str(Fraction(333, 334)),
        "low_response": str(Fraction(166, 167)),
        "high_response": str(Fraction(55_445, 55_778)),
        "shell_counts": (low_count, high_count),
        "row_swap_shape": (2, 334),
        "row_swap_response": "-1",
    }
    print(
        "ORDER668 normalized_sha256=%s product_counts=333,334 "
        "responses(one,low,high)=(333/334,166/167,55445/55778) "
        "shell_counts=%d,%d row_swap=2x334:rho=-1"
        % (ORDER668_NORMALIZED_SHA256, low_count, high_count)
    )
    return data


def main():
    print("THM-3407 HADAMARD CORE MULTI-TOGGLE RESPONSE: EXACT CONTROLS")
    order668 = dependency_and_order668_audit()
    triple_palette = abstract_triple_palette_audit()
    semantic = {
        "order668": order668,
        "triple_palette": triple_palette,
        "orders": {},
    }

    for prime in (3, 7, 11, 19):
        hadamard = paley_type_i(prime)
        order = prime + 1
        m = order // 4
        sign_core, binary_core, determinant = check_hadamard_and_core(hadamard)
        pair_data = two_toggle_audit(
            order, sign_core, binary_core, determinant, m
        )
        window_data = window_audit(
            order, sign_core, binary_core, determinant, m
        )
        trade_data = trade_and_rectangle_audit(
            order, hadamard, sign_core, binary_core, determinant, m
        )
        triple_data = None
        if order in TRIPLE_HOSTILES:
            triple_data = triple_hostile_audit(
                order, sign_core, binary_core, determinant, m
            )
        semantic["orders"][str(order)] = {
            "pairs": pair_data,
            "window": window_data,
            "trades": trade_data,
            "triple": triple_data,
        }
        counts = pair_data["counts"]
        print(
            "PALEY N=%d m=%d two_shells(low,high)=(%s,%s) "
            "counts(low,high)=(%d,%d) direct_pair_checks=%d"
            % (
                order,
                m,
                pair_data["low"],
                pair_data["high"],
                counts["shared_row"]
                + counts["shared_column"]
                + counts["positive_plaquette"],
                counts["negative_plaquette"],
                pair_data["direct_checks"],
            )
        )
        print(
            "   window masks=512 mobius_nonzero=%s equality_sizes=%s "
            "same_line_tests=%d row_column_swap=%d,%d"
            % (
                ",".join(map(str, window_data["nonzero_mobius_by_degree"])),
                window_data["equality_sizes"],
                trade_data["same_line_counts"],
                trade_data["row_swap_size"],
                trade_data["column_swap_size"],
            )
        )
        if triple_data is not None:
            print(
                "   triple_hostile pair_data=%s signed_minors=%s responses=%s determinants=%s"
                % (
                    triple_data["pair_data"],
                    triple_data["signed_minors"],
                    triple_data["responses"],
                    triple_data["determinants"],
                )
            )

    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(digest == EXPECTED_SEMANTIC_SHA256, "semantic digest drift")
    print("semantic_sha256=" + digest)
    print("FINITE CONTROLS ONLY; UNIVERSAL FORMULAS AND TRADE FLOOR USE THE PROOF")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
