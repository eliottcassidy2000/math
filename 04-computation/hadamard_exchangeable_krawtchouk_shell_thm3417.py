#!/usr/bin/env python3
"""Exact controls for THM-3417's exchangeable Krawtchouk shell theorem.

The companion independently checks both radial inversions, the symmetrized
parity detector bank, direct Paley determinants, strict information-loss
hostiles, the fixed-shell convolution boundary, and the low-degree Hermite
normalization.  It uses only the Python standard library and exact arithmetic.
"""

from __future__ import annotations

import hashlib
import itertools
from fractions import Fraction
from math import comb
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINS = (
    (
        "THM-3413",
        ROOT
        / "01-canon/theorems/THM-3413-strength-k-orthogonal-array-toggle-filtration-and-high-minor-converse.md",
        "0b15fb0be8de43a67f8f330a9259febf0106614e6f553a668ffea4bf6d71164e",
    ),
    (
        "THM-3413-script",
        ROOT / "04-computation/hadamard_strength_k_toggle_filtration_thm3413.py",
        "dbdf97d32e4114abfd71393791ac8b956d4f1e46eef9a98ec14532107dd0a1a6",
    ),
    (
        "THM-3413-output",
        ROOT / "05-knowledge/results/hadamard_strength_k_toggle_filtration_thm3413.out",
        "55bbcec2736c72dbcadac16450fe326317971990693ab47e91676b3e51ec4172",
    ),
)

EXPECTED_SEMANTIC_SHA256 = "d2471a54b365e186d605bdc5640434219013efe1cff787a18fbf66989c4e355b"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def sha256_lf(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest(value: object) -> str:
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def choose(n: int, r: int) -> int:
    return comb(n, r) if 0 <= r <= n else 0


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


def interactions(
    sign_core: list[list[int]], events: tuple[tuple[int, int], ...]
) -> tuple[int, ...]:
    return tuple(
        signed_minor(sign_core, events, mask) for mask in range(1 << len(events))
    )


def response_values(m: int, packet: tuple[int, ...]) -> tuple[Fraction, ...]:
    size = len(packet).bit_length() - 1
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


def direct_values(
    binary_core: list[list[int]], events: tuple[tuple[int, int], ...]
) -> tuple[Fraction, ...]:
    determinant_base = determinant(binary_core)
    require(determinant_base != 0, "singular base")
    values = []
    for toggle_mask in range(1 << len(events)):
        changed = [row[:] for row in binary_core]
        for index, (row, column) in enumerate(events):
            if toggle_mask & (1 << index):
                changed[row][column] ^= 1
        values.append(Fraction(determinant(changed), determinant_base))
    return tuple(values)


def walsh_coefficients(m: int, packet: tuple[int, ...]) -> tuple[Fraction, ...]:
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


def krawtchouk(t: int, degree: int, weight: int) -> int:
    return sum(
        (-1) ** overlap
        * choose(weight, overlap)
        * choose(t - weight, degree - overlap)
        for overlap in range(degree + 1)
    )


def shell_average(
    values: tuple[Fraction, ...], t: int, weight: int
) -> Fraction:
    shell = tuple(value for mask, value in enumerate(values) if mask.bit_count() == weight)
    require(len(shell) == choose(t, weight), (t, weight, len(shell)))
    return sum(shell, Fraction()) / len(shell)


def radial_data(
    m: int, packet: tuple[int, ...], values: tuple[Fraction, ...], label: object
) -> tuple[
    tuple[Fraction, ...],
    tuple[Fraction, ...],
    tuple[int, ...],
    tuple[Fraction, ...],
]:
    t = len(packet).bit_length() - 1
    require(len(values) == 1 << t, (label, "value count"))
    coefficients = walsh_coefficients(m, packet)
    responses = tuple(shell_average(values, t, weight) for weight in range(t + 1))
    walsh_shells = tuple(
        sum(
            (
                coefficients[mask]
                for mask in range(1 << t)
                if mask.bit_count() == degree
            ),
            Fraction(),
        )
        for degree in range(t + 1)
    )
    minor_shells = tuple(
        sum(
            packet[mask]
            for mask in range(1 << t)
            if mask.bit_count() == degree
        )
        for degree in range(t + 1)
    )

    for weight, response in enumerate(responses):
        predicted = sum(
            (
                Fraction(krawtchouk(t, degree, weight), choose(t, degree))
                * walsh_shells[degree]
                for degree in range(t + 1)
            ),
            Fraction(),
        )
        require(predicted == response, (label, "Krawtchouk forward", weight))
        predicted_minor = sum(
            (
                Fraction(
                    (-1) ** degree
                    * choose(weight, degree)
                    * minor_shells[degree],
                    (2 * m) ** degree * choose(t, degree),
                )
                for degree in range(weight + 1)
            ),
            Fraction(),
        )
        require(
            predicted_minor == response,
            (label, "minor-shell forward", weight),
        )

    recovered_walsh = tuple(
        Fraction(1, 2**t)
        * sum(
            choose(t, weight)
            * krawtchouk(t, degree, weight)
            * responses[weight]
            for weight in range(t + 1)
        )
        for degree in range(t + 1)
    )
    require(recovered_walsh == walsh_shells, (label, "Krawtchouk inverse"))

    recovered_minors = tuple(
        (2 * m) ** degree
        * choose(t, degree)
        * sum(
            (-1) ** weight * choose(degree, weight) * responses[weight]
            for weight in range(degree + 1)
        )
        for degree in range(t + 1)
    )
    require(recovered_minors == minor_shells, (label, "finite-difference inverse"))

    for degree in range(t + 1):
        predicted = sum(
            Fraction(
                (-1) ** (source + degree)
                * choose(source, degree)
                * minor_shells[source],
                (4 * m) ** source,
            )
            for source in range(degree, t + 1)
        )
        require(predicted == walsh_shells[degree], (label, "radial Mobius", degree))
        recovered = (4 * m) ** degree * sum(
            choose(source, degree) * walsh_shells[source]
            for source in range(degree, t + 1)
        )
        require(recovered == minor_shells[degree], (label, "radial inverse", degree))

    return responses, walsh_shells, minor_shells, coefficients


def walsh_transform(values: tuple[int, ...]) -> tuple[int, ...]:
    data = list(values)
    step = 1
    while step < len(data):
        for start in range(0, len(data), 2 * step):
            for offset in range(step):
                left = data[start + offset]
                right = data[start + step + offset]
                data[start + offset] = left + right
                data[start + step + offset] = left - right
        step *= 2
    return tuple(data)


def orthogonality_audit(max_t: int) -> int:
    checks = 0
    for t in range(max_t + 1):
        for left in range(t + 1):
            for right in range(t + 1):
                value = sum(
                    choose(t, weight)
                    * krawtchouk(t, left, weight)
                    * krawtchouk(t, right, weight)
                    for weight in range(t + 1)
                )
                expected = 2**t * choose(t, left) * int(left == right)
                require(value == expected, ("orthogonality", t, left, right))
                checks += 1
    return checks


def abstract_audit(max_t: int) -> tuple[int, str]:
    records = []
    for t in range(1, max_t + 1):
        packet = tuple(
            1
            if mask == 0
            else (-1) ** (mask.bit_count() + mask % 3)
            * (1 + mask.bit_count() + 2 * mask)
            for mask in range(1 << t)
        )
        values = response_values(3, packet)
        records.append(radial_data(3, packet, values, ("abstract", t))[:3])
    return max_t, digest(records)


def symmetrized_parity_audit(max_t: int) -> tuple[int, int, str]:
    cases = 0
    moment_checks = 0
    records = []
    for t in range(1, max_t + 1):
        for degree in range(1, t + 1):
            level_size = choose(t, degree)
            for sign in (-1, 1):
                row_counts = []
                for mask in range(1 << t):
                    radial_sum = krawtchouk(t, degree, mask.bit_count())
                    numerator = level_size + sign * radial_sum
                    require(numerator % 2 == 0, ("parity count", t, degree, mask))
                    count = numerator // 2
                    require(0 <= count <= level_size, ("parity positivity", t, degree))
                    row_counts.append(count)
                moments = walsh_transform(tuple(row_counts))
                runs = level_size * 2 ** (t - 1)
                require(moments[0] == runs, ("parity runs", t, degree, sign))
                for target, moment in enumerate(moments[1:], start=1):
                    expected = (
                        sign * 2 ** (t - 1)
                        if target.bit_count() == degree
                        else 0
                    )
                    require(
                        moment == expected,
                        ("parity moment", t, degree, sign, target),
                    )
                    moment_checks += 1
                records.append((t, degree, sign, min(row_counts), max(row_counts), runs))
                cases += 1
    return cases, moment_checks, digest(records)


def paley_four_exhaustion() -> tuple[int, int, int, str]:
    sign_core, binary_core = paley_core(3)
    positions = tuple(itertools.product(range(3), repeat=2))
    records = []
    repeated_row_or_column = 0
    shell_checks = 0
    for selection in range(1, 1 << len(positions)):
        events = tuple(
            position for index, position in enumerate(positions) if selection & (1 << index)
        )
        packet = interactions(sign_core, events)
        event_values = response_values(1, packet)
        direct = direct_values(binary_core, events)
        require(event_values == direct, ("Paley four direct", events))
        radial = radial_data(1, packet, direct, ("Paley four", events))
        records.append((events, radial[:3]))
        shell_checks += len(events) + 1
        if len({row for row, _ in events}) < len(events) or len(
            {column for _, column in events}
        ) < len(events):
            repeated_row_or_column += 1
    return len(records), repeated_row_or_column, shell_checks, digest(records)


def level_profile(coefficients: tuple[Fraction, ...], degree: int) -> tuple[Fraction, ...]:
    return tuple(
        sorted(
            coefficient
            for mask, coefficient in enumerate(coefficients)
            if mask.bit_count() == degree
        )
    )


def h4_collision_audit() -> tuple[object, ...]:
    sign_core, binary_core = paley_core(3)
    first = ((0, 0), (0, 1), (0, 2), (2, 0))
    second = ((0, 0), (0, 1), (1, 0), (1, 1))
    records = []
    for label, events in (("first", first), ("second", second)):
        packet = interactions(sign_core, events)
        direct = direct_values(binary_core, events)
        require(direct == response_values(1, packet), ("H4 collision direct", label))
        radial = radial_data(1, packet, direct, ("H4 collision", label))
        records.append((events, radial))
    first_radial = records[0][1]
    second_radial = records[1][1]
    expected = (Fraction(1), Fraction(1, 2), Fraction(1, 6), Fraction(0), Fraction(0))
    require(first_radial[0] == second_radial[0] == expected, "H4 equal shell means")
    first_profile = level_profile(first_radial[3], 1)
    second_profile = level_profile(second_radial[3], 1)
    require(
        first_profile
        == (Fraction(0), Fraction(1, 8), Fraction(1, 8), Fraction(1, 4)),
        "H4 first profile",
    )
    require(second_profile == (Fraction(1, 8),) * 4, "H4 second profile")
    require(first_profile != second_profile, "H4 collision must lose localization")
    return first, second, expected, first_profile, second_profile


def h8_strict_loss_audit() -> tuple[object, ...]:
    sign_core, binary_core = paley_core(7)
    events = ((0, 0), (1, 2), (2, 1), (3, 4))
    packet = interactions(sign_core, events)
    direct = direct_values(binary_core, events)
    require(direct == response_values(2, packet), "H8 strict-loss direct")
    responses, walsh_shells, minor_shells, coefficients = radial_data(
        2, packet, direct, "H8 strict loss"
    )
    expected_responses = (
        Fraction(1),
        Fraction(3, 4),
        Fraction(13, 24),
        Fraction(3, 8),
        Fraction(1, 4),
    )
    expected_walsh = (
        Fraction(9, 16),
        Fraction(3, 8),
        Fraction(1, 16),
        Fraction(0),
        Fraction(0),
    )
    expected_high = (
        (7, Fraction(1, 128)),
        (11, Fraction(0)),
        (13, Fraction(-1, 128)),
        (14, Fraction(0)),
        (15, Fraction(0)),
    )
    high = tuple(
        (mask, coefficients[mask])
        for mask in range(1 << 4)
        if mask.bit_count() > 2
    )
    require(responses == expected_responses, "H8 strict-loss responses")
    require(walsh_shells == expected_walsh, "H8 strict-loss radial Walsh")
    require(minor_shells == (1, 4, 4, 0, 0), "H8 strict-loss minors")
    require(high == expected_high, "H8 strict-loss labelled high packet")
    require(any(value for _, value in high), "H8 strict-loss nontriviality")
    return events, responses, walsh_shells, minor_shells, high


def xor_convolution(
    left: tuple[Fraction, ...], right: tuple[Fraction, ...]
) -> tuple[Fraction, ...]:
    require(len(left) == len(right), "convolution dimensions")
    size = len(left)
    return tuple(
        sum((left[index] * right[target ^ index] for index in range(size)), Fraction())
        for target in range(size)
    )


def convolution_audit() -> tuple[object, ...]:
    sign_core, binary_core = paley_core(7)
    events = ((0, 1), (1, 0), (2, 3), (3, 2))
    packet = interactions(sign_core, events)
    values = direct_values(binary_core, events)
    require(values == response_values(2, packet), "H8 convolution direct")
    _, walsh_shells, _, _ = radial_data(2, packet, values, "H8 convolution")
    require(walsh_shells[4] == Fraction(1, 512), "H8 top coefficient")
    records = []
    for weight in (1, 2):
        one_step = tuple(
            Fraction(int(mask.bit_count() == weight), choose(4, weight))
            for mask in range(1 << 4)
        )
        distribution = (Fraction(1),) + (Fraction(0),) * 15
        moments = tuple(
            Fraction(krawtchouk(4, degree, weight), choose(4, degree))
            for degree in range(5)
        )
        rho = max(abs(moment) for moment in moments[1:4])
        require(rho < 1 and moments[4] == (-1) ** weight, ("parity boundary", weight))
        responses = []
        for power in range(1, 13):
            distribution = xor_convolution(distribution, one_step)
            require(sum(distribution, Fraction()) == 1, ("convolution mass", weight, power))
            require(
                all(
                    probability == 0
                    or mask.bit_count() % 2 == (power * weight) % 2
                    for mask, probability in enumerate(distribution)
                ),
                ("convolution parity support", weight, power),
            )
            observed = sum(
                (probability * values[mask] for mask, probability in enumerate(distribution)),
                Fraction(),
            )
            predicted = sum(
                walsh_shells[degree] * moments[degree] ** power
                for degree in range(5)
            )
            require(observed == predicted, ("convolution response", weight, power))
            error = observed - walsh_shells[0] - (-1) ** (power * weight) * walsh_shells[4]
            envelope = rho**power * sum(abs(value) for value in walsh_shells[1:4])
            require(abs(error) <= envelope, ("convolution envelope", weight, power))
            responses.append(observed)
        records.append((weight, moments, tuple(responses), distribution))

    boundary_checks = 0
    for t in range(2, 13):
        for weight in range(1, t):
            moments = tuple(
                Fraction(krawtchouk(t, degree, weight), choose(t, degree))
                for degree in range(t + 1)
            )
            require(moments[t] == (-1) ** weight, ("top parity", t, weight))
            require(
                all(abs(moments[degree]) < 1 for degree in range(1, t)),
                ("strict internal mixing", t, weight),
            )
            boundary_checks += t - 1
    return events, walsh_shells, boundary_checks, digest(records)


def larger_paley_control() -> tuple[object, ...]:
    sign_core, binary_core = paley_core(11)
    events = ((0, 0), (1, 2), (2, 4), (3, 6), (4, 8), (5, 10), (6, 1))
    packet = interactions(sign_core, events)
    direct = direct_values(binary_core, events)
    require(direct == response_values(3, packet), "H12 direct")
    radial = radial_data(3, packet, direct, "H12 seven-toggle")
    return events, radial[:3]


def duplicate_position_boundary() -> tuple[object, ...]:
    sign_core, binary_core = paley_core(3)
    events = ((0, 0), (0, 0))
    packet = interactions(sign_core, events)
    additive = response_values(1, packet)
    xor_values = direct_values(binary_core, events)
    require(additive[:3] == xor_values[:3], "duplicate one-event faces")
    require(additive[3] == 0 and xor_values[3] == 1, "duplicate XOR boundary")
    return events, packet, additive, xor_values


def hermite(degree: int, value: int) -> int:
    if degree == 0:
        return 1
    previous, current = 1, value
    for index in range(1, degree):
        previous, current = current, value * current - index * previous
    return current


def hermite_audit() -> tuple[int, str]:
    checks = 0
    for t in range(3, 41):
        for weight in range(t + 1):
            y = t - 2 * weight
            exact = (
                Fraction(y, t),
                Fraction(y * y - t, t * (t - 1)),
                Fraction(y**3 - (3 * t - 2) * y, t * (t - 1) * (t - 2)),
            )
            observed = tuple(
                Fraction(krawtchouk(t, degree, weight), choose(t, degree))
                for degree in range(1, 4)
            )
            require(observed == exact, ("low Hermite normalization", t, weight))
            checks += 3

    records = []
    for root in (4, 8, 16, 32):
        t = root * root
        for x in (-2, -1, 0, 1, 2):
            weight = (t - x * root) // 2
            require(2 * weight == t - x * root, ("central lattice", root, x))
            errors = tuple(
                Fraction(root**degree * krawtchouk(t, degree, weight), choose(t, degree))
                - hermite(degree, x)
                for degree in range(7)
            )
            records.append((root, x, errors))
    return checks, digest(records)


def main() -> None:
    for label, path, expected in PINS:
        require(path.exists(), ("missing dependency", label, path))
        require(sha256_lf(path) == expected, ("dependency drift", label, path))

    orthogonality = orthogonality_audit(12)
    abstract = abstract_audit(9)
    parity = symmetrized_parity_audit(9)
    paley_four = paley_four_exhaustion()
    h4_collision = h4_collision_audit()
    h8_loss = h8_strict_loss_audit()
    convolution = convolution_audit()
    h12 = larger_paley_control()
    duplicate = duplicate_position_boundary()
    hermite_checks = hermite_audit()

    semantic_object = (
        orthogonality,
        abstract,
        parity,
        paley_four,
        h4_collision,
        h8_loss,
        convolution,
        h12,
        duplicate,
        hermite_checks,
    )
    semantic = digest(semantic_object)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic drift", semantic))

    print("THM-3417 exchangeable Krawtchouk shell exact controls")
    print(f"dependency_pins={len(PINS)}")
    print(f"krawtchouk_orthogonality_checks={orthogonality}")
    print(f"abstract_radial_audit=(max_t,digest)={abstract}")
    print(f"symmetrized_parity=(laws,moment_checks,digest)={parity}")
    print(
        "paley_order4=(position_sets,repeated_row_or_column,shell_checks,digest)="
        f"{paley_four}"
    )
    print(f"h4_radial_fibre_collision={h4_collision}")
    print(f"h8_exchangeable_strength2_strict_loss={h8_loss}")
    print(
        "fixed_shell_convolution=(events,walsh_shells,boundary_checks,digest)="
        f"{convolution}"
    )
    print(f"paley_order12_seven_toggle_digest={digest(h12)}")
    print(f"duplicate_position_boundary={duplicate}")
    print(f"hermite_low_degree_checks_and_limit_digest={hermite_checks}")
    print(f"semantic_sha256={semantic}")
    print("status=PASS")


if __name__ == "__main__":
    main()
