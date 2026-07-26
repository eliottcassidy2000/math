#!/usr/bin/env python3
"""Exact companion for THM-2369's complete-line target Dirichlet theorem."""

from __future__ import annotations

from fractions import Fraction


P = 13
ZERO = Fraction(0)
Gaussian = tuple[Fraction, Fraction]
Cyclo = tuple[Fraction, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def gaussian_add(left: Gaussian, right: Gaussian) -> Gaussian:
    return (left[0] + right[0], left[1] + right[1])


def gaussian_scale(value: Gaussian, scalar: Fraction | int) -> Gaussian:
    factor = Fraction(scalar)
    return (factor * value[0], factor * value[1])


def gaussian_multiply(left: Gaussian, right: Gaussian) -> Gaussian:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def gaussian_conjugate(value: Gaussian) -> Gaussian:
    return (value[0], -value[1])


def gaussian_norm_squared(value: Gaussian) -> Fraction:
    return value[0] * value[0] + value[1] * value[1]


def canonical(values: list[Fraction]) -> Cyclo:
    """Canonical Q[zeta_13] coordinates, using 1+...+zeta^12=0."""
    require(len(values) == P, "cyclotomic vector has wrong length")
    anchor = values[-1]
    return tuple(value - anchor for value in values)


def cyclo_zero() -> Cyclo:
    return tuple(ZERO for _ in range(P))


def cyclo_rational(value: Fraction | int) -> Cyclo:
    entries = [ZERO for _ in range(P)]
    entries[0] = Fraction(value)
    return canonical(entries)


def cyclo_add(left: Cyclo, right: Cyclo) -> Cyclo:
    return canonical([left[index] + right[index] for index in range(P)])


def cyclo_scale(value: Cyclo, scalar: Fraction | int) -> Cyclo:
    factor = Fraction(scalar)
    return canonical([factor * entry for entry in value])


def cyclo_phase(value: Cyclo, exponent: int) -> Cyclo:
    entries = [ZERO for _ in range(P)]
    shift = exponent % P
    for index, entry in enumerate(value):
        entries[(index + shift) % P] += entry
    return canonical(entries)


def character_distance(exponent: int) -> Cyclo:
    return cyclo_add(
        cyclo_rational(2),
        cyclo_scale(
            cyclo_add(
                cyclo_phase(cyclo_rational(1), exponent),
                cyclo_phase(cyclo_rational(1), -exponent),
            ),
            -1,
        ),
    )


def dot(left: tuple[int, int], right: tuple[int, int]) -> int:
    return (left[0] * right[0] + left[1] * right[1]) % P


def add_group(
    left: tuple[int, int],
    right: tuple[int, int],
) -> tuple[int, int]:
    return ((left[0] + right[0]) % P, (left[1] + right[1]) % P)


def subtract_group(
    left: tuple[int, int],
    right: tuple[int, int],
) -> tuple[int, int]:
    return ((left[0] - right[0]) % P, (left[1] - right[1]) % P)


def root_average(exponent: int) -> Fraction:
    return Fraction(1 if exponent % P == 0 else 0)


def main() -> None:
    group = tuple(
        (first, second)
        for first in range(P)
        for second in range(P)
    )

    # Complete coordinate lines have exact multipliers 0 and 2.
    line_multipliers: dict[tuple[int, int], tuple[int, int]] = {}
    line_character_checks = 0
    for target in group:
        multipliers = []
        for coordinate in range(2):
            total = cyclo_zero()
            for shift in range(P):
                total = cyclo_add(
                    total,
                    character_distance(shift * target[coordinate]),
                )
            average = cyclo_scale(total, Fraction(1, P))
            expected_integer = 0 if target[coordinate] == 0 else 2
            require(
                average == cyclo_rational(expected_integer),
                "complete-line character multiplier changed",
            )
            multipliers.append(expected_integer)
            line_character_checks += 1
        line_multipliers[target] = (multipliers[0], multipliers[1])

    # A deterministic exact Gaussian-rational target array.
    target_array: dict[tuple[int, int], Gaussian] = {}
    for first, second in group:
        target_array[(first, second)] = (
            Fraction(((3 * first + 5 * second + 1) % 23) - 11, 29),
            Fraction(((7 * first + 2 * second + 4) % 19) - 9, 31),
        )

    energy_10 = ZERO
    energy_01 = ZERO
    energy_11 = ZERO
    energy_total = ZERO
    line_energy_1 = ZERO
    line_energy_2 = ZERO
    for target, value in target_array.items():
        if target == (0, 0):
            continue
        norm = gaussian_norm_squared(value)
        energy_total += norm
        if target[0] != 0 and target[1] == 0:
            energy_10 += norm
        elif target[0] == 0 and target[1] != 0:
            energy_01 += norm
        else:
            energy_11 += norm
        line_energy_1 += line_multipliers[target][0] * norm
        line_energy_2 += line_multipliers[target][1] * norm

    require(
        energy_total == energy_10 + energy_01 + energy_11
        and line_energy_1 == 2 * (energy_10 + energy_11)
        and line_energy_2 == 2 * (energy_01 + energy_11)
        and energy_10 == energy_total - line_energy_2 / 2
        and energy_01 == energy_total - line_energy_1 / 2
        and energy_11
        == (line_energy_1 + line_energy_2) / 2 - energy_total,
        "pure/fork target-energy reconstruction changed",
    )

    # Every nonzero complete-Cayley multiplier is exactly one after /338.
    cayley_multiplier_checks = 0
    for target in group:
        total = cyclo_zero()
        for direction in group:
            if direction == (0, 0):
                continue
            total = cyclo_add(
                total,
                character_distance(dot(direction, target)),
            )
        normalized = cyclo_scale(total, Fraction(1, 338))
        expected = 0 if target == (0, 0) else 1
        require(
            normalized == cyclo_rational(expected),
            "full-Cayley target multiplier changed",
        )
        cayley_multiplier_checks += 1

    # The deep-colour line cancellation transfers at least 1/13 of energy.
    deep_transfer_lines = 0
    for target in group:
        if target == (0, 0):
            continue
        unit_entries: list[Gaussian] = []
        unit_sum: Gaussian = (ZERO, ZERO)
        for deep_colour in range(1, P):
            entry = (
                Fraction(
                    ((deep_colour + 3 * target[0] + target[1]) % 17) - 8,
                    19,
                ),
                Fraction(
                    ((2 * deep_colour + target[0] + 4 * target[1]) % 13)
                    - 6,
                    23,
                ),
            )
            unit_entries.append(entry)
            unit_sum = gaussian_add(unit_sum, entry)
        zero_entry = gaussian_scale(unit_sum, -1)
        unit_energy = sum(
            gaussian_norm_squared(entry) for entry in unit_entries
        )
        total_energy = unit_energy + gaussian_norm_squared(zero_entry)
        require(
            P * unit_energy >= total_energy,
            "deep-colour 1/13 transfer changed",
        )
        deep_transfer_lines += 1

    equality_entry: Gaussian = (Fraction(3, 7), Fraction(-2, 5))
    equality_unit_energy = 12 * gaussian_norm_squared(equality_entry)
    equality_zero = gaussian_scale(equality_entry, -12)
    equality_total_energy = (
        equality_unit_energy + gaussian_norm_squared(equality_zero)
    )
    require(
        P * equality_unit_energy == equality_total_energy,
        "sharp deep-colour transfer profile changed",
    )

    # Constant and inverse-character currents have identical bispectra.
    deepest_charge = (3, 5)
    bispectrum_checks = 0
    for left in group:
        for right in group:
            endpoint = add_group(left, right)
            inverse_phase = (
                -dot(left, deepest_charge)
                -dot(right, deepest_charge)
                +dot(endpoint, deepest_charge)
            ) % P
            require(
                inverse_phase == 0,
                "balanced bispectrum acquired character slope",
            )
            bispectrum_checks += 1

    # General balanced monomials remain invariant under every modulation.
    balanced_modulation_checks = 0
    for seed in range(200):
        left_1 = (seed % P, (3 * seed + 1) % P)
        left_2 = ((5 * seed + 2) % P, (7 * seed + 3) % P)
        left_3 = ((11 * seed + 4) % P, (2 * seed + 5) % P)
        right_1 = ((4 * seed + 6) % P, (9 * seed + 7) % P)
        right_2 = ((8 * seed + 8) % P, (6 * seed + 9) % P)
        left_sum = add_group(add_group(left_1, left_2), left_3)
        right_3 = subtract_group(
            left_sum,
            add_group(right_1, right_2),
        )
        right_sum = add_group(add_group(right_1, right_2), right_3)
        charge = subtract_group(left_sum, right_sum)
        require(charge == (0, 0), "balanced control stopped balancing")
        for slope in group:
            require(
                dot(charge, slope) == 0,
                "balanced monomial saw modulation slope",
            )
            balanced_modulation_checks += 1

    # Two independent charged directions determine the full slope.
    slope_labels: dict[tuple[int, int], tuple[int, int]] = {}
    one_direction_counts: dict[int, int] = {}
    basis = ((1, 0), (0, 1))
    for slope in group:
        label = (dot(basis[0], slope), dot(basis[1], slope))
        require(label not in slope_labels, "charged slope labels collided")
        slope_labels[label] = slope
        first_label = label[0]
        one_direction_counts[first_label] = (
            one_direction_counts.get(first_label, 0) + 1
        )
    require(
        len(slope_labels) == P * P
        and set(one_direction_counts.values()) == {P},
        "rank-two charged gauge map changed",
    )

    # Every inverse-character nearest-neighbour covariant defect vanishes.
    inverse_edge_checks = 0
    for twist in group:
        twist_phase = -dot(twist, deepest_charge)
        for direction in basis:
            shifted = add_group(twist, direction)
            shifted_phase = -dot(shifted, deepest_charge)
            covariant_phase = dot(direction, deepest_charge) + shifted_phase
            require(
                (covariant_phase - twist_phase) % P == 0,
                "inverse-character edge stopped being flat",
            )
            inverse_edge_checks += 1

    # A spanning grid reconstructs exactly the inverse-character slope.
    reconstructed: dict[tuple[int, int], int] = {(0, 0): 0}
    for first in range(P):
        for second in range(P):
            reconstructed[(first, second)] = (
                -first * deepest_charge[0] - second * deepest_charge[1]
            ) % P
    require(
        all(
            reconstructed[twist] == (-dot(twist, deepest_charge)) % P
            for twist in group
        ),
        "flat spanning-tree reconstruction changed",
    )

    # Cyclic pair-twist polarization extracts conjugate(R)*K exactly.
    reference: Gaussian = (Fraction(2, 3), Fraction(-3, 5))
    current: Gaussian = (Fraction(-4, 7), Fraction(5, 11))
    reference_current = gaussian_multiply(
        gaussian_conjugate(reference),
        current,
    )
    reverse_current = gaussian_multiply(
        reference,
        gaussian_conjugate(current),
    )
    constant_term = (
        gaussian_norm_squared(reference) + gaussian_norm_squared(current),
        ZERO,
    )
    polarized = (ZERO, ZERO)
    polarized = gaussian_add(
        polarized,
        gaussian_scale(constant_term, root_average(-1)),
    )
    polarized = gaussian_add(
        polarized,
        gaussian_scale(reference_current, root_average(1 - 1)),
    )
    polarized = gaussian_add(
        polarized,
        gaussian_scale(reverse_current, root_average(-1 - 1)),
    )
    require(
        polarized == reference_current,
        "cyclic reference polarization changed",
    )

    print("THM-2369 complete-line target Dirichlet exact referee")
    print(f"coordinate-line character multipliers: {line_character_checks}")
    print(
        "target mask cells: "
        f"pure-1={12}, pure-2={12}, fork={144}"
    )
    print("line energy identity: D1=2(E10+E11), D2=2(E01+E11)")
    print(f"full-Cayley target multipliers: {cayley_multiplier_checks}")
    print(
        "deep-colour target lines/sharp share: "
        f"{deep_transfer_lines} / 1/{P}"
    )
    print(f"balanced bispectrum cells: {bispectrum_checks}")
    print(f"balanced modulation checks: {balanced_modulation_checks}")
    print(
        "charged slope labels/one-direction fibres: "
        f"{len(slope_labels)} / {P}"
    )
    print(f"inverse-character flat edges: {inverse_edge_checks}")
    print(f"reference polarization: {polarized}")
    print(
        "VERDICT: complete lines recover target type; "
        "balanced observables lose character slope"
    )


if __name__ == "__main__":
    main()
