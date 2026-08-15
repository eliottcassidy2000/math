#!/usr/bin/env python3
"""Exact controls for THM-3413's strength-k toggle filtration.

The theorem is a formal Boolean--Walsh transform.  This standard-library
companion checks the transform and its inverse on abstract packets, exhausts
all distinct toggle-position sets in the Paley order-four core, replays direct
determinants in larger Paley cores, and verifies every parity-halfcube detector.
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
    (
        "THM-3411",
        ROOT / "01-canon/theorems/THM-3411-pairwise-independent-toggle-filter-and-sharp-high-minor-norm.md",
        "624443ef3614fdc8993ef05e84e85d9e57f1f6fd97535c36d1e6c26767853c13",
    ),
    (
        "THM-3411-script",
        ROOT / "04-computation/hadamard_pairwise_independent_toggle_filter_thm3411.py",
        "376c8d5bf5a9f68672bd680ea5df367b5eb752fb30228b2154abdfc37e4d3415",
    ),
    (
        "THM-3411-output",
        ROOT / "05-knowledge/results/hadamard_pairwise_independent_toggle_filter_thm3411.out",
        "10fdb3d5cce9fd98e64b8ea78087ec4e66d06097a1c52af0bf518de0388ee93b",
    ),
)

HOSTILES_4 = (
    ((0, 1), (1, 0), (2, 3), (3, 2)),
    ((0, 0), (1, 1), (2, 2), (3, 4)),
    ((0, 1), (1, 2), (2, 5), (3, 0)),
    ((0, 0), (1, 1), (2, 3), (5, 6)),
)

LARGER_CONTROLS = (
    (7, ((0, 1), (1, 2), (2, 0))),
    (7, HOSTILES_4[0]),
    (7, ((0, 0), (0, 1), (1, 0), (2, 3), (3, 2))),
    (7, ((0, 0), (0, 1), (1, 0), (1, 2), (2, 1), (3, 3))),
    (11, ((0, 0), (1, 2), (2, 4), (3, 6), (4, 8), (5, 10), (6, 1))),
)

EXPECTED_SEMANTIC_SHA256 = "0db6229ea8922fb2f9c570103dd8aff7e30da424bcc196f29de19b54e2cab329"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def sha256_lf(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def determinant(matrix: list[list[int]]) -> int:
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
    indices = tuple(index for index in range(len(events)) if mask & (1 << index))
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
    return tuple(
        signed_minor(sign_core, events, mask) for mask in range(1 << len(events))
    )


def response_values(m: int, packet: tuple[int, ...]) -> tuple[Fraction, ...]:
    size = packet.__len__().bit_length() - 1
    require(1 << size == len(packet), len(packet))
    values = []
    for toggle_mask in range(1 << size):
        value = Fraction(0)
        for mask, minor in enumerate(packet):
            if mask & ~toggle_mask:
                continue
            value += Fraction(
                (-1) ** mask.bit_count() * minor,
                (2 * m) ** mask.bit_count(),
            )
        values.append(value)
    return tuple(values)


def direct_values(binary_core, events):
    determinant_base = determinant(binary_core)
    values = []
    for toggle_mask in range(1 << len(events)):
        changed = [row[:] for row in binary_core]
        for index, (row, column) in enumerate(events):
            if toggle_mask & (1 << index):
                changed[row][column] ^= 1
        values.append(Fraction(determinant(changed), determinant_base))
    return tuple(values)


def walsh(values: tuple[Fraction, ...], mask: int) -> Fraction:
    return sum(
        value if (toggle_mask & mask).bit_count() % 2 == 0 else -value
        for toggle_mask, value in enumerate(values)
    ) / len(values)


def high_coefficients(m: int, packet: tuple[int, ...]) -> tuple[Fraction, ...]:
    coefficients = []
    for target in range(len(packet)):
        value = Fraction(0)
        for source, minor in enumerate(packet):
            if target & ~source:
                continue
            value += Fraction(
                (-1) ** (target.bit_count() + source.bit_count()) * minor,
                (4 * m) ** source.bit_count(),
            )
        coefficients.append(value)
    return tuple(coefficients)


def inverse_packet(m: int, coefficients: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    recovered = []
    for source in range(len(coefficients)):
        value = sum(
            coefficient
            for target, coefficient in enumerate(coefficients)
            if source & ~target == 0
        )
        recovered.append((4 * m) ** source.bit_count() * value)
    return tuple(recovered)


def parity_expectation(values, target, sign):
    selected = tuple(
        value
        for toggle_mask, value in enumerate(values)
        if (-1) ** (toggle_mask & target).bit_count() == sign
    )
    require(len(selected) * 2 == len(values), (len(values), target, sign))
    return sum(selected) / len(selected)


def packet_audit(m: int, packet: tuple[int, ...], direct=None):
    values = response_values(m, packet)
    if direct is not None:
        require(values == direct, ("direct response", m, packet, values, direct))
    coefficients = high_coefficients(m, packet)
    actual = tuple(walsh(values, mask) for mask in range(len(packet)))
    require(coefficients == actual, ("Walsh transform", m, packet))
    require(inverse_packet(m, coefficients) == packet, ("inverse", m, packet))
    uniform = coefficients[0]
    parity_checks = 0
    for target in range(1, len(packet)):
        for sign in (-1, 1):
            expectation = parity_expectation(values, target, sign)
            require(
                expectation - uniform == sign * coefficients[target],
                ("parity detector", target, sign),
            )
            parity_checks += 1
    for strength in range(len(packet).bit_length()):
        high_masks = tuple(
            mask for mask in range(1, len(packet)) if mask.bit_count() > strength
        )
        coefficient_blind = all(coefficients[mask] == 0 for mask in high_masks)
        minor_blind = all(packet[mask] == 0 for mask in high_masks)
        require(coefficient_blind == minor_blind, ("blind converse", strength))
    return coefficients, values, parity_checks


def walsh_transform(values: list[int]) -> tuple[int, ...]:
    transformed = values[:]
    stride = 1
    while stride < len(transformed):
        for start in range(0, len(transformed), 2 * stride):
            for offset in range(stride):
                left = transformed[start + offset]
                right = transformed[start + stride + offset]
                transformed[start + offset] = left + right
                transformed[start + stride + offset] = left - right
        stride *= 2
    return tuple(transformed)


def parity_array_audit(max_size: int):
    checks = 0
    records = []
    for size in range(1, max_size + 1):
        length = 1 << size
        count = length // 2
        for target in range(1, length):
            for sign in (-1, 1):
                indicator = [
                    int((-1) ** (toggle_mask & target).bit_count() == sign)
                    for toggle_mask in range(length)
                ]
                moments = walsh_transform(indicator)
                expected = tuple(
                    count if mask == 0 else sign * count if mask == target else 0
                    for mask in range(length)
                )
                require(moments == expected, ("parity OA", size, target, sign))
                checks += length
        records.append((size, digest(tuple(
            tuple(
                toggle_mask
                for toggle_mask in range(length)
                if (-1) ** (toggle_mask & target).bit_count() == 1
            )
            for target in range(1, length)
        ))))
    return checks, tuple(records)


def abstract_transform_audit():
    records = []
    total_checks = 0
    for size in range(1, 9):
        for seed in range(5):
            packet = tuple(
                1
                if mask == 0
                else ((-1) ** ((mask + seed).bit_count()))
                * (((mask * mask + 3 * seed + 5 * mask.bit_count()) % 19) - 9)
                for mask in range(1 << size)
            )
            coefficients, values, parity_checks = packet_audit(seed + 1, packet)
            total_checks += len(packet) + parity_checks
            records.append((size, seed, digest((packet, coefficients, values))))
    return total_checks, digest(tuple(records))


def paley_four_exhaustion():
    sign_core, binary_core = paley_core(3)
    positions = tuple(itertools.product(range(3), repeat=2))
    counts = Counter()
    response_digest = []
    total_parity = 0
    for size in range(1, len(positions) + 1):
        for events in itertools.combinations(positions, size):
            packet = interactions(sign_core, events)
            coefficients, values, parity_checks = packet_audit(
                1, packet, direct_values(binary_core, events)
            )
            total_parity += parity_checks
            for strength in range(size + 1):
                blind = all(
                    packet[mask] == 0
                    for mask in range(1, 1 << size)
                    if mask.bit_count() > strength
                )
                counts[(size, strength, blind)] += 1
            response_digest.append((events, packet, coefficients, values))
    require(sum(count for (size, strength, blind), count in counts.items() if strength == size) == (1 << 9) - 1, counts)
    return (
        (1 << 9) - 1,
        total_parity,
        tuple(sorted(counts.items())),
        digest(tuple(response_digest)),
    )


def larger_paley_controls():
    cores = {}
    records = []
    total_parity = 0
    for prime, events in LARGER_CONTROLS:
        if prime not in cores:
            cores[prime] = paley_core(prime)
        sign_core, binary_core = cores[prime]
        packet = interactions(sign_core, events)
        coefficients, values, parity_checks = packet_audit(
            (prime + 1) // 4,
            packet,
            direct_values(binary_core, events),
        )
        total_parity += parity_checks
        records.append(
            (
                prime + 1,
                events,
                tuple(packet[mask] for mask in range(len(packet)) if mask.bit_count() >= max(1, len(events) - 1)),
                digest(coefficients),
                digest(values),
            )
        )
    return len(records), total_parity, tuple(records)


def thm3411_regression():
    sign_core, binary_core = paley_core(7)
    records = []
    for events in HOSTILES_4:
        packet = interactions(sign_core, events)
        coefficients, values, _ = packet_audit(
            2, packet, direct_values(binary_core, events)
        )
        four_minor = packet[15]
        expected = tuple(
            Fraction(8 * packet[15 ^ (1 << index)] - four_minor, 4096)
            for index in range(4)
        ) + (Fraction(four_minor, 4096),)
        observed = tuple(coefficients[15 ^ (1 << index)] for index in range(4)) + (
            coefficients[15],
        )
        require(observed == expected, ("THM-3411 regression", events))
        records.append((events, observed, digest(values)))
    return tuple(records)


def convolution_audit():
    sign_core, binary_core = paley_core(7)
    events = LARGER_CONTROLS[2][1]
    packet = interactions(sign_core, events)
    coefficients, values, _ = packet_audit(
        2, packet, direct_values(binary_core, events)
    )
    size = len(events)
    length = 1 << size
    moments = {
        0b00111: Fraction(1, 5),
        0b11001: Fraction(-1, 7),
        0b11110: Fraction(1, 11),
        0b11111: Fraction(-1, 13),
    }
    require(all(mask.bit_count() > 2 for mask in moments), moments)
    require(sum(abs(value) for value in moments.values()) < 1, moments)

    deltas = []
    law_digests = []
    for power in range(1, 13):
        powered = {mask: value**power for mask, value in moments.items()}
        weights = tuple(
            Fraction(1, length)
            * (
                1
                + sum(
                    moment
                    * (-1) ** (toggle_mask & mask).bit_count()
                    for mask, moment in powered.items()
                )
            )
            for toggle_mask in range(length)
        )
        require(all(weight >= 0 for weight in weights), (power, weights))
        require(sum(weights) == 1, (power, weights))
        recovered_moments = tuple(
            sum(
                weight * (-1) ** (toggle_mask & mask).bit_count()
                for toggle_mask, weight in enumerate(weights)
            )
            for mask in range(length)
        )
        expected_moments = tuple(
            1 if mask == 0 else powered.get(mask, 0) for mask in range(length)
        )
        require(recovered_moments == expected_moments, (power, recovered_moments))
        expectation = sum(weight * value for weight, value in zip(weights, values))
        delta = expectation - coefficients[0]
        predicted = sum(
            coefficients[mask] * moment for mask, moment in powered.items()
        )
        require(delta == predicted, ("convolution response", power, delta, predicted))
        deltas.append(delta)
        law_digests.append(digest(weights))

    roots = tuple(moments.values())
    characteristic = [Fraction(1)]
    for root in roots:
        updated = [Fraction(0)] * (len(characteristic) + 1)
        for index, coefficient in enumerate(characteristic):
            updated[index] -= root * coefficient
            updated[index + 1] += coefficient
        characteristic = updated
    for start in range(len(deltas) - len(roots)):
        recurrence = sum(
            characteristic[index] * deltas[start + index]
            for index in range(len(characteristic))
        )
        require(recurrence == 0, ("C-finite recurrence", start, recurrence))

    parity_mask = 0b00111
    parity_coefficient = coefficients[parity_mask]
    parity_words = tuple(
        sign**power * parity_coefficient
        for sign in (-1, 1)
        for power in range(1, 7)
    )
    return (
        events,
        tuple(sorted(moments.items())),
        tuple(characteristic),
        tuple(deltas),
        tuple(law_digests),
        parity_words,
    )


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

    abstract = abstract_transform_audit()
    parity_arrays = parity_array_audit(9)
    paley4 = paley_four_exhaustion()
    larger = larger_paley_controls()
    regression = thm3411_regression()
    convolution = convolution_audit()
    semantic_record = (
        tuple((name, expected) for name, _, expected in PINS),
        abstract,
        parity_arrays,
        paley4,
        larger,
        regression,
        convolution,
    )
    semantic = digest(semantic_record)
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, (semantic, EXPECTED_SEMANTIC_SHA256))

    print("THM-3413 STRENGTH-K ORTHOGONAL-ARRAY TOGGLE FILTRATION")
    print(f"source_sha256_lf={sha256_lf(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINS)}")
    print("status=PROVED_CANDIDATE_plus_VERIFIED_EXACT;independent_audit_requested")
    print(f"abstract_transform=(checks,digest)={abstract}")
    print(f"parity_halfcubes=(moment_checks,size_digest)={(parity_arrays[0], digest(parity_arrays[1]))}")
    print(f"paley_order4_exhaustion=(position_sets,parity_tests,blind_table_digest,response_digest)={(paley4[0], paley4[1], digest(paley4[2]), paley4[3])}")
    print(f"larger_paley_controls=(packets,parity_tests,records_digest)={(larger[0], larger[1], digest(larger[2]))}")
    print(f"thm3411_regression=(packets,digest)={(len(regression), digest(regression))}")
    print(f"convolution_spectrum=(powers,recurrence_degree,digest)={(len(convolution[3]), len(convolution[2]) - 1, digest(convolution))}")
    print("scope=distinct_toggle_positions;rows_or_columns_may_repeat;strength_k_Rademacher_laws;normalized_determinant_average;no_completion_or_existence_claim")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
