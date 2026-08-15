#!/usr/bin/env python3
"""Exact controls for THM-3411's pairwise-independent toggle filter.

The proof is symbolic Boolean-to-Walsh algebra.  This standard-library
companion supplies independent direct-determinant controls, reconstructs all
twenty-six extremal OA laws, exhausts the nonattacking four-event universe of
the Paley order-eight core, and checks a fixed-row Paley order-twelve window.
All gates survive ``python -O``.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
from collections import Counter
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINS = (
    (
        "THM-3396",
        ROOT / "01-canon/theorems/THM-3396-four-bit-pairwise-independent-fourier-cone.md",
        "db3fb5fe0abe5a9f549c0f3c95a6fe17b0815030dce68e26f9ff9df922ca005f",
    ),
    (
        "THM-3396-script",
        ROOT / "04-computation/four_bit_pairwise_independent_fourier_cone_thm3396.py",
        "43bc4e9bdc99bf0723ace6d124128a00fd8f89bd43e8738cf7b4bb8ce2139bcc",
    ),
    (
        "THM-3396-output",
        ROOT / "05-knowledge/results/four_bit_pairwise_independent_fourier_cone_thm3396.out",
        "47bc7bda869b6b5ce799fee770e8ca56dc1fbbcb1c9a33145678d411f3edfa0b",
    ),
    (
        "THM-3407",
        ROOT / "01-canon/theorems/THM-3407-hadamard-core-multitoggle-response-plaquette-shells-and-trade-distance.md",
        "d0d277966d34782f9272da5ea8b6200378e1b80abf3c54de1e2d796870ea0479",
    ),
    (
        "THM-3407-script",
        ROOT / "04-computation/hadamard_core_multitoggle_response_thm3407.py",
        "83f38703b126cdc9bf9358bce20803ed890a1f19d685c7644e909c1a8df45e90",
    ),
    (
        "THM-3407-output",
        ROOT / "05-knowledge/results/hadamard_core_multitoggle_response_thm3407.out",
        "46b6b6176a23ebb0554c0359315f5584227576c753b5e3e1d1f6cc8ce1f1d31f",
    ),
)

HOSTILES = (
    (((0, 1), (1, 0), (2, 3), (3, 2)), (0, 0, 0, 0), 8),
    (((0, 0), (1, 1), (2, 2), (3, 4)), (0, 4, 0, 4), 0),
    (((0, 1), (1, 2), (2, 5), (3, 0)), (0, 4, 0, 4), 8),
    (((0, 0), (1, 1), (2, 3), (5, 6)), (0, 0, 0, 0), 0),
)

EXPECTED_ORDER8_HISTOGRAM = (
    (Fraction(0), 504),
    (Fraction(5, 1536), 672),
    (Fraction(3, 512), 8064),
    (Fraction(5, 768), 2016),
    (Fraction(11, 1536), 5376),
    (Fraction(1, 128), 8064),
    (Fraction(13, 1536), 2016),
    (Fraction(5, 512), 2688),
)

EXPECTED_ORDER12_WINDOW = (
    7_920,
    (
        (Fraction(0), 480),
        (Fraction(5, 7776), 96),
        (Fraction(5, 2592), 1512),
        (Fraction(17, 7776), 1008),
        (Fraction(1, 432), 3792),
        (Fraction(7, 2592), 928),
        (Fraction(23, 7776), 104),
    ),
    77,
    480,
    "960bf2469f49da3fa61b78b82f33c243f81ab71ac96dbec95a6677d519805548",
    158,
)
EXPECTED_SEMANTIC_SHA256 = "e702bae75e08cb5bba1e103a5eeb20f21ba540b3300edea9010f9f9c3bad9e02"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def sha256_lf(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant, including the empty matrix."""
    size = len(matrix)
    require(all(len(row) == size for row in matrix), ("nonsquare", matrix))
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
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
                require(numerator % previous == 0, ("Bareiss", matrix))
                work[row][index] = numerator // previous
        previous = pivot
        for row in range(column + 1, size):
            work[row][column] = 0
    return sign * work[-1][-1]


def paley_core(prime: int) -> tuple[list[list[int]], list[list[int]]]:
    require(prime % 4 == 3, prime)

    def character(value: int) -> int:
        value %= prime
        if value == 0:
            return 0
        return 1 if pow(value, (prime - 1) // 2, prime) == 1 else -1

    sign_core = [
        [character(row - column) - int(row == column) for column in range(prime)]
        for row in range(prime)
    ]
    binary_core = [[(1 - value) // 2 for value in row] for row in sign_core]
    hadamard = [[1] * (prime + 1)] + [[1] + row for row in sign_core]
    gram = [
        [sum(left * right for left, right in zip(row, other)) for other in hadamard]
        for row in hadamard
    ]
    require(
        gram
        == [
            [(prime + 1) * int(row == column) for column in range(prime + 1)]
            for row in range(prime + 1)
        ],
        ("Paley Gram", prime),
    )
    return sign_core, binary_core


def signed_minor(
    sign_core: list[list[int]], events: tuple[tuple[int, int], ...], mask: int
) -> int:
    indices = tuple(index for index in range(4) if mask & (1 << index))
    rows = tuple(events[index][0] for index in indices)
    columns = tuple(events[index][1] for index in indices)
    if len(set(rows)) < len(rows) or len(set(columns)) < len(columns):
        return 0
    sign = 1
    for index in indices:
        row, column = events[index]
        sign *= sign_core[row][column]
    matrix = [
        [sign_core[events[left][0]][events[right][1]] for right in indices]
        for left in indices
    ]
    return sign * determinant(matrix)


def interactions(sign_core, events):
    return tuple(signed_minor(sign_core, events, mask) for mask in range(16))


def response_values(m: int, packet: tuple[int, ...]) -> tuple[Fraction, ...]:
    values = []
    for toggle_mask in range(16):
        value = Fraction(0)
        for mask, minor in enumerate(packet):
            if mask & ~toggle_mask:
                continue
            value += Fraction((-1) ** mask.bit_count() * minor, (2 * m) ** mask.bit_count())
        values.append(value)
    return tuple(values)


def walsh(values: tuple[Fraction, ...], mask: int) -> Fraction:
    return sum(
        (value if (toggle_mask & mask).bit_count() % 2 == 0 else -value)
        for toggle_mask, value in enumerate(values)
    ) / 16


def coefficient_packet(m: int, packet: tuple[int, ...]) -> tuple[Fraction, ...]:
    four_minor = packet[15]
    return tuple(
        Fraction(4 * m * packet[15 ^ (1 << index)] - four_minor, 256 * m**4)
        for index in range(4)
    ) + (Fraction(four_minor, 256 * m**4),)


def sharp_norm(coefficients: tuple[Fraction, ...]) -> Fraction:
    return max(
        max(abs(value) for value in coefficients),
        sum(abs(value) for value in coefficients) / 3,
    )


def oa_vertices():
    vertices = []
    for index in range(5):
        for sign in (-1, 1):
            point = tuple(Fraction(sign * int(position == index)) for position in range(5))
            vertices.append(("H8", point))
    for signs in itertools.product((-1, 1), repeat=5):
        product = 1
        for sign in signs:
            product *= sign
        if product == 1:
            vertices.append(("H12", tuple(Fraction(sign, 3) for sign in signs)))
    require(len(vertices) == 26, len(vertices))
    records = []
    for kind, point in vertices:
        weights = []
        for toggle_mask in range(16):
            x = tuple(-1 if toggle_mask & (1 << index) else 1 for index in range(4))
            parity = x[0] * x[1] * x[2] * x[3]
            weight = Fraction(
                1 + parity * (point[4] + sum(point[index] * x[index] for index in range(4))),
                16,
            )
            require(weight >= 0, (kind, point, toggle_mask, weight))
            weights.append(weight)
        require(sum(weights) == 1, (kind, point, weights))
        run_count = 8 if kind == "H8" else 12
        cells = tuple(weight * run_count for weight in weights)
        require(all(cell.denominator == 1 for cell in cells), (kind, point, cells))
        records.append((kind, point, tuple(cell.numerator for cell in cells), tuple(weights)))
    return tuple(records)


def filter_audit(m, packet, vertices):
    values = response_values(m, packet)
    predicted = coefficient_packet(m, packet)
    masks = tuple(15 ^ (1 << index) for index in range(4)) + (15,)
    actual = tuple(walsh(values, mask) for mask in masks)
    require(actual == predicted, ("Walsh packet", m, packet, actual, predicted))
    uniform = sum(values) / 16
    deviations = []
    for _, point, _, weights in vertices:
        expectation = sum(weight * value for weight, value in zip(weights, values))
        contracted = sum(left * right for left, right in zip(point, predicted))
        require(expectation - uniform == contracted, ("OA contraction", point))
        deviations.append(abs(contracted))
    require(max(deviations) == sharp_norm(predicted), ("sharp vertex", predicted))
    return predicted, sharp_norm(predicted), values


def direct_values(binary_core, events):
    determinant_base = determinant(binary_core)
    values = []
    for toggle_mask in range(16):
        changed = [row[:] for row in binary_core]
        for index, (row, column) in enumerate(events):
            if toggle_mask & (1 << index):
                changed[row][column] ^= 1
        values.append(Fraction(determinant(changed), determinant_base))
    return tuple(values)


def hostile_audit(sign_core, binary_core, vertices):
    records = []
    for events, expected_three, expected_four in HOSTILES:
        packet = interactions(sign_core, events)
        three = tuple(packet[15 ^ (1 << index)] for index in range(4))
        require((three, packet[15]) == (expected_three, expected_four), (
            "hostile interaction", events, three, packet[15],
        ))
        coefficients, norm, values = filter_audit(2, packet, vertices)
        require(values == direct_values(binary_core, events), ("direct determinant", events))
        numerator = tuple(value * (256 * 2**4) for value in coefficients)
        require(all(value.denominator == 1 for value in numerator), numerator)
        records.append(
            (
                events,
                three,
                packet[15],
                tuple(value.numerator for value in numerator),
                norm,
                digest(values),
            )
        )
    require(records[1][1] == records[2][1] and records[1][2] != records[2][2], records)
    return tuple(records)


def nonattacking_census(sign_core, m, row_bank, column_size, vertices):
    histogram = Counter()
    interaction_profiles = Counter()
    total = 0
    full_walsh_checks = 0
    for rows in row_bank:
        for columns_unordered in itertools.combinations(range(column_size), 4):
            for columns in itertools.permutations(columns_unordered):
                events = tuple(zip(rows, columns))
                three = tuple(
                    signed_minor(sign_core, events, 15 ^ (1 << index))
                    for index in range(4)
                )
                four = signed_minor(sign_core, events, 15)
                high_packet = [0] * 16
                for index, value in enumerate(three):
                    high_packet[15 ^ (1 << index)] = value
                high_packet[15] = four
                coefficients = coefficient_packet(m, tuple(high_packet))
                norm = sharp_norm(coefficients)
                profile = (three, four)
                if profile not in interaction_profiles or total % 97 == 0:
                    packet = interactions(sign_core, events)
                    audited_coefficients, audited_norm, _ = filter_audit(
                        m, packet, vertices
                    )
                    require(
                        (audited_coefficients, audited_norm) == (coefficients, norm),
                        ("sampled full response", m, events),
                    )
                    full_walsh_checks += 1
                histogram[norm] += 1
                interaction_profiles[profile] += 1
                total += 1
                require(
                    norm <= Fraction(4 * m + 5, 48 * m**4),
                    ("universal bound", m, events, coefficients, norm),
                )
    return (
        total,
        tuple(sorted(histogram.items())),
        len(interaction_profiles),
        sum(count for (three, four), count in interaction_profiles.items() if three == (0, 0, 0, 0) and four == 0),
        digest(tuple(sorted(interaction_profiles.items()))),
        full_walsh_checks,
    )


def abstract_bound_audit():
    checks = 0
    worst = Fraction(0)
    for m in range(1, 13):
        bound = Fraction(4 * m + 5, 48 * m**4)
        for three in itertools.product((-4, 0, 4), repeat=4):
            for four in (-16, -8, 0, 8, 16):
                coefficients = tuple(
                    Fraction(4 * m * value - four, 256 * m**4) for value in three
                ) + (Fraction(four, 256 * m**4),)
                norm = sharp_norm(coefficients)
                require(norm <= bound, (m, three, four, norm, bound))
                worst = max(worst, norm / bound)
                checks += 1
    return checks, worst


def main() -> None:
    for name, path, expected in PINS:
        require(sha256_lf(path) == expected, ("dependency changed", name, path))
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    vertices = oa_vertices()
    vertex_histogram = Counter(kind for kind, _, _, _ in vertices)
    require(vertex_histogram == Counter({"H8": 10, "H12": 16}), vertex_histogram)

    sign8, binary8 = paley_core(7)
    hostile_records = hostile_audit(sign8, binary8, vertices)
    census8 = nonattacking_census(
        sign8,
        2,
        itertools.combinations(range(7), 4),
        7,
        vertices,
    )
    require(census8[0] == 29_400, census8)
    require(census8[1] == EXPECTED_ORDER8_HISTOGRAM, census8[1])
    require(census8[3] == 504, census8)

    sign12, _ = paley_core(11)
    census12 = nonattacking_census(sign12, 3, ((0, 1, 2, 3),), 11, vertices)
    require(census12[0] == 7_920, census12)
    if EXPECTED_ORDER12_WINDOW is not None:
        require(census12 == EXPECTED_ORDER12_WINDOW, census12)

    bound_audit = abstract_bound_audit()
    require(bound_audit == (4_860, Fraction(1)), bound_audit)
    vertex_record = tuple((kind, point, cells) for kind, point, cells, _ in vertices)
    semantic_record = (
        tuple((name, expected) for name, _, expected in PINS),
        vertex_record,
        hostile_records,
        census8,
        census12,
        bound_audit,
    )
    semantic = digest(semantic_record)
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, (semantic, EXPECTED_SEMANTIC_SHA256))

    print("THM-3411 PAIRWISE-INDEPENDENT HADAMARD TOGGLE FILTER")
    print(f"source_sha256_lf={sha256_lf(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINS)}")
    print("status=PROVED_CANDIDATE_plus_VERIFIED_EXACT;independent_audit_requested")
    print(f"OA_vertices=(H8,H12,total)={(vertex_histogram['H8'],vertex_histogram['H12'],len(vertices))};vertex_sha256={digest(vertex_record)}")
    print(f"hostiles=(events,three_minors,four_minor,Walsh_numerators,sharp_bias,response_sha256)={hostile_records}")
    print(f"order8_nonattacking=(total,bias_histogram,interaction_profiles,blind,response_checks)={(census8[0],census8[1],census8[2],census8[3],census8[5])};profile_sha256={census8[4]}")
    print(f"order12_fixed_rows_0123=(total,bias_histogram,interaction_profiles,blind,response_checks)={(census12[0],census12[1],census12[2],census12[3],census12[5])};profile_sha256={census12[4]}")
    print(f"universal_bound_abstract_controls=(checks,worst_ratio)={bound_audit}")
    print("scope=four_distinct_candidate_positions;rows_or_columns_may_repeat;unbiased_pairwise_independent_filter;normalized_determinant_average;no_Hadamard_completion_or_existence_claim")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
