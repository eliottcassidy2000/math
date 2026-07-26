#!/usr/bin/env python3
"""Exact companion for THM-2374's Boolean subcube Dirichlet spectrum."""

from __future__ import annotations

from fractions import Fraction


ZERO = Fraction(0)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def parity(value: int) -> int:
    return value.bit_count() & 1


def submasks(mask: int) -> tuple[int, ...]:
    values: list[int] = []
    current = mask
    while True:
        values.append(current)
        if current == 0:
            return tuple(values)
        current = (current - 1) & mask


def walsh_coefficients(
    values: tuple[Fraction, ...],
    dimension: int,
) -> dict[int, Fraction]:
    size = 1 << dimension
    require(len(values) == size, "Boolean table has wrong size")
    result: dict[int, Fraction] = {}
    for frequency in range(size):
        result[frequency] = sum(
            value * (-1 if parity(frequency & point) else 1)
            for point, value in enumerate(values)
        ) / size
    return result


def subcube_dirichlet(
    values: tuple[Fraction, ...],
    dimension: int,
    coordinate_mask: int,
) -> Fraction:
    shifts = submasks(coordinate_mask)
    size = 1 << dimension
    return sum(
        (values[point ^ shift] - values[point]) ** 2
        for point in range(size)
        for shift in shifts
    ) / (size * len(shifts))


def spectral_dirichlet(
    coefficients: dict[int, Fraction],
    coordinate_mask: int,
) -> Fraction:
    return 2 * sum(
        coefficient * coefficient
        for frequency, coefficient in coefficients.items()
        if frequency & coordinate_mask
    )


def type_support(frequency: int, type_masks: tuple[int, ...]) -> int:
    support = 0
    for index, mask in enumerate(type_masks):
        if frequency & mask:
            support |= 1 << index
    return support


def main() -> None:
    dimension = 4
    size = 1 << dimension
    all_coordinates = size - 1
    type_masks = (0b0001, 0b0010, 0b1100)
    require(
        type_masks[0] | type_masks[1] | type_masks[2] == all_coordinates,
        "type partition changed",
    )

    trials = 1054
    subcube_checks = 0
    pair_checks = 0
    inversion_checks = 0
    difference_profiles: set[tuple[Fraction, ...]] = set()
    for trial in range(trials):
        values = tuple(
            Fraction(
                (trial + 1) * (point + 1) * (point + 1)
                + (point + 1) ** 3
            )
            for point in range(size)
        )
        difference_profiles.add(
            tuple(values[point] - values[0] for point in range(1, size))
        )
        coefficients = walsh_coefficients(values, dimension)
        dirichlet = {
            mask: subcube_dirichlet(values, dimension, mask)
            for mask in range(size)
        }

        for mask in range(size):
            require(
                dirichlet[mask]
                == spectral_dirichlet(coefficients, mask),
                "complete-subcube spectral identity changed",
            )
            subcube_checks += 1

        for left_index in range(len(type_masks)):
            for right_index in range(left_index + 1, len(type_masks)):
                left = type_masks[left_index]
                right = type_masks[right_index]
                defect = (
                    dirichlet[left]
                    + dirichlet[right]
                    - dirichlet[left | right]
                )
                expected = 2 * sum(
                    coefficient * coefficient
                    for frequency, coefficient in coefficients.items()
                    if frequency & left and frequency & right
                )
                require(
                    defect == expected and defect >= 0,
                    "pairwise type defect changed",
                )
                pair_checks += 1

        mixed_energy = sum(
            coefficient * coefficient
            for frequency, coefficient in coefficients.items()
            if type_support(frequency, type_masks).bit_count() >= 2
        )
        complement_formula = (
            sum(
                dirichlet[all_coordinates ^ type_mask]
                for type_mask in type_masks
            )
            - (len(type_masks) - 1) * dirichlet[all_coordinates]
        ) / 2
        require(
            mixed_energy == complement_formula,
            "unweighted mixed-energy formula changed",
        )

        for frequency in range(1, size):
            inverted = -sum(
                (-1) ** (frequency.bit_count() - subset.bit_count())
                * dirichlet[all_coordinates ^ subset]
                for subset in submasks(frequency)
            ) / 2
            require(
                inverted == coefficients[frequency] ** 2,
                "squared Walsh-support inversion changed",
            )
            inversion_checks += 1
    require(
        len(difference_profiles) == trials,
        "exact table bank repeated a difference profile",
    )

    # THM-2348's exact two-by-two fixed-optimum hostile.
    mass = Fraction(7, 3)
    hostile = (ZERO, mass, mass, mass)
    hostile_coefficients = walsh_coefficients(hostile, 2)
    d_first = subcube_dirichlet(hostile, 2, 0b01)
    d_second = subcube_dirichlet(hostile, 2, 0b10)
    d_full = subcube_dirichlet(hostile, 2, 0b11)
    hostile_defect = d_first + d_second - d_full
    require(
        tuple(hostile_coefficients[index] for index in range(4))
        == (
            3 * mass / 4,
            -mass / 4,
            -mass / 4,
            -mass / 4,
        )
        and d_first == mass * mass / 4
        and d_second == mass * mass / 4
        and d_full == 3 * mass * mass / 8
        and hostile_defect == mass * mass / 8,
        "two-by-two hostile constants changed",
    )

    rectangle_checks = 0
    rectangle_contrasts: set[Fraction] = set()
    for seed in range(373):
        first = Fraction(seed, 11)
        second = Fraction(2 * seed + 1, 13)
        third = Fraction(3 * seed + 2, 17)
        prescribed_contrast = Fraction(seed + 1, 379)
        fourth = prescribed_contrast - first + second + third
        table = (
            first,
            second,
            third,
            fourth,
        )
        contrast = table[0] - table[1] - table[2] + table[3]
        rectangle_contrasts.add(contrast)
        defect = (
            subcube_dirichlet(table, 2, 0b01)
            + subcube_dirichlet(table, 2, 0b10)
            - subcube_dirichlet(table, 2, 0b11)
        )
        require(
            defect == contrast * contrast / 8,
            "rectangle-contrast normalization changed",
        )
        rectangle_checks += 1
    require(
        len(rectangle_contrasts) == rectangle_checks,
        "rectangle bank repeated a contrast",
    )

    # Swapping the two actual block labels translates the whole cube.
    translated = tuple(
        hostile[point ^ 0b11]
        for point in range(4)
    )
    require(
        all(
            subcube_dirichlet(hostile, 2, mask)
            == subcube_dirichlet(translated, 2, mask)
            for mask in range(4)
        ),
        "global two-block relabelling changed the Dirichlet bank",
    )

    # A selected binary face of a three-block table is not globally decisive.
    colours = range(3)
    larger_table = {
        (left, right): Fraction(1 if left == 2 and right == 2 else 0)
        for left in colours
        for right in colours
    }
    binary_face = {
        larger_table[(left, right)]
        for left in (0, 1)
        for right in (0, 1)
    }
    larger_contrast = (
        larger_table[(0, 0)]
        - larger_table[(0, 2)]
        - larger_table[(2, 0)]
        + larger_table[(2, 2)]
    )
    require(
        binary_face == {ZERO} and larger_contrast == 1,
        "larger-block face hostile changed",
    )

    print("THM-2374 binary-allocation Dirichlet exact referee")
    print(f"distinct exact rational tables: {trials}")
    print(f"complete-subcube spectral checks: {subcube_checks}")
    print(f"prime-type pair defects: {pair_checks}")
    print(f"squared-support inversion cells: {inversion_checks}")
    print(
        "two-by-two hostile Dp,Dq,Dpq,J: "
        f"{d_first}, {d_second}, {d_full}, {hostile_defect}"
    )
    print(f"distinct rectangle-contrast checks: {rectangle_checks}")
    print("global swap / larger-block face hostiles: PASS / PASS")
    print(
        "VERDICT: complete subcubes recover squared Walsh support; "
        "mixed-type defects are nonnegative"
    )


if __name__ == "__main__":
    main()
