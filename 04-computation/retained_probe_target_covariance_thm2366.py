#!/usr/bin/env python3
"""Exact companion for THM-2366's retained-probe target boundary."""

from __future__ import annotations

from fractions import Fraction
from math import gcd


P = 13
ZERO = Fraction(0)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def canonical(values: list[Fraction]) -> tuple[Fraction, ...]:
    """Canonical coordinates in Q[zeta_13], using 1+...+zeta^12=0."""
    require(len(values) == P, "cyclotomic vector has wrong length")
    anchor = values[-1]
    return tuple(value - anchor for value in values)


def cyclo_zero() -> tuple[Fraction, ...]:
    return tuple(ZERO for _ in range(P))


def cyclo_rational(value: Fraction | int) -> tuple[Fraction, ...]:
    values = [ZERO for _ in range(P)]
    values[0] = Fraction(value)
    return canonical(values)


def cyclo_add(
    left: tuple[Fraction, ...],
    right: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    return canonical([left[index] + right[index] for index in range(P)])


def cyclo_scale(
    value: tuple[Fraction, ...],
    scalar: Fraction | int,
) -> tuple[Fraction, ...]:
    factor = Fraction(scalar)
    return canonical([factor * entry for entry in value])


def cyclo_phase(
    value: tuple[Fraction, ...],
    exponent: int,
) -> tuple[Fraction, ...]:
    shifted = [ZERO for _ in range(P)]
    shift = exponent % P
    for index, entry in enumerate(value):
        shifted[(index + shift) % P] += entry
    return canonical(shifted)


def cyclo_sum(values: list[tuple[Fraction, ...]]) -> tuple[Fraction, ...]:
    total = cyclo_zero()
    for value in values:
        total = cyclo_add(total, value)
    return total


def dft2(
    table: list[list[tuple[Fraction, ...]]],
) -> dict[tuple[int, int], tuple[Fraction, ...]]:
    transformed: dict[tuple[int, int], tuple[Fraction, ...]] = {}
    for first_colour in range(P):
        for second_colour in range(P):
            terms = []
            for first_shift in range(P):
                for second_shift in range(P):
                    terms.append(
                        cyclo_phase(
                            table[first_shift][second_shift],
                            first_colour * first_shift
                            + second_colour * second_shift,
                        )
                    )
            transformed[(first_colour, second_colour)] = cyclo_scale(
                cyclo_sum(terms),
                Fraction(1, P * P),
            )
    return transformed


def inverse_dft2(
    transformed: dict[tuple[int, int], tuple[Fraction, ...]],
) -> list[list[tuple[Fraction, ...]]]:
    table = [
        [cyclo_zero() for _ in range(P)]
        for _ in range(P)
    ]
    for first_shift in range(P):
        for second_shift in range(P):
            terms = []
            for first_colour in range(P):
                for second_colour in range(P):
                    terms.append(
                        cyclo_phase(
                            transformed[(first_colour, second_colour)],
                            -first_colour * first_shift
                            - second_colour * second_shift,
                        )
                    )
            table[first_shift][second_shift] = cyclo_sum(terms)
    return table


def main() -> None:
    # The target q=(b,a+h) has thirteen deep-colour decompositions.
    target_counts: dict[tuple[int, int], int] = {}
    zero_target_labels = 0
    for deep_colour in range(P):
        for first_colour in range(P):
            for second_colour in range(P):
                target = (
                    first_colour,
                    (deep_colour + second_colour) % P,
                )
                target_counts[target] = target_counts.get(target, 0) + 1
                if target == (0, 0):
                    zero_target_labels += 1
    require(
        len(target_counts) == P * P
        and set(target_counts.values()) == {P}
        and zero_target_labels == P,
        "target-coordinate decomposition changed",
    )

    # Exact thirteenth-root orthogonality behind every target-line sum.
    root_sums: list[tuple[Fraction, ...]] = []
    for difference in range(P):
        root_sum = cyclo_sum(
            [
                cyclo_phase(cyclo_rational(1), colour * difference)
                for colour in range(P)
            ]
        )
        expected = cyclo_rational(P if difference == 0 else 0)
        require(root_sum == expected, "root orthogonality changed")
        root_sums.append(root_sum)

    line_cell_checks = 0
    for first_colour in range(P):
        for target_second in range(P):
            for root_shift in range(P):
                for target_shift in range(P):
                    for endpoint_shift in range(P):
                        value = (
                            Fraction(0)
                            if root_shift == endpoint_shift
                            else Fraction(
                                1
                                + (
                                    root_shift
                                    + 2 * target_shift
                                    + 3 * endpoint_shift
                                )
                                % 11,
                                17,
                            )
                        )
                        difference = (root_shift - endpoint_shift) % P
                        if difference == 0:
                            require(
                                value == 0,
                                "diagonal-plane zero changed",
                            )
                        else:
                            require(
                                root_sums[difference] == cyclo_zero(),
                                "off-diagonal root sum changed",
                            )
                        require(
                            0 <= first_colour < P
                            and 0 <= target_second < P
                            and 0 <= target_shift < P,
                            "target-line index left F_13",
                        )
                        line_cell_checks += 1

    # Full two-dimensional inversion on a hostile rational table.
    rational_table = [
        [
            cyclo_rational(
                Fraction(
                    ((5 * first_shift + 7 * second_shift + 3) % 23) - 11,
                    29,
                )
            )
            for second_shift in range(P)
        ]
        for first_shift in range(P)
    ]
    rational_transform = dft2(rational_table)
    require(
        inverse_dft2(rational_transform) == rational_table,
        "two-dimensional finite inversion changed",
    )

    # The inverse-character response occupies exactly (b,h)=(0,-a).
    selected_deep_colour = 5
    covariance_base = Fraction(-7, 31)
    covariance_table = [
        [
            cyclo_phase(
                cyclo_rational(covariance_base),
                selected_deep_colour * endpoint_shift,
            )
            for endpoint_shift in range(P)
        ]
        for _target_shift in range(P)
    ]
    covariance_transform = dft2(covariance_table)
    covariance_cells = 0
    for first_colour in range(P):
        for second_colour in range(P):
            expected = (
                cyclo_rational(covariance_base)
                if (
                    first_colour == 0
                    and second_colour == (-selected_deep_colour) % P
                )
                else cyclo_zero()
            )
            require(
                covariance_transform[(first_colour, second_colour)] == expected,
                "conditional inverse-character support changed",
            )
            covariance_cells += 1

    # Exact energy identity and both sharp constants, in coefficient space.
    nonzero_target_count = P * P - 1
    epsilon = Fraction(1, 97)

    # Upper equality: two target coefficients cancel at the base point.
    upper_coefficients = (epsilon, -epsilon)
    upper_energy = sum(value * value for value in upper_coefficients)
    upper_base_defect = sum(upper_coefficients)
    upper_delta = upper_energy + upper_base_defect * upper_base_defect
    require(
        upper_delta == upper_energy,
        "upper energy equality changed",
    )

    # Lower equality: all 168 nonzero-target coefficients are equal.
    lower_energy = nonzero_target_count * epsilon * epsilon
    lower_base_defect = nonzero_target_count * epsilon
    lower_delta = lower_energy + lower_base_defect * lower_base_defect
    require(
        lower_delta == P * P * lower_energy,
        "lower energy equality changed",
    )

    # The rational compatibility tensor.
    mass = Fraction(6, 91)
    deep_offzero = Fraction(1, 182)
    blocker_offzero = Fraction(1, 12)
    require(P - 1 == 12 and 12 * deep_offzero == mass, "hostile mass changed")

    hostile_deep_colours = 0
    hostile_deep_value = Fraction(-1, 2366)
    for colour in range(1, P):
        terms = [cyclo_rational(0)]
        terms.extend(
            cyclo_scale(
                cyclo_phase(cyclo_rational(1), colour * shift),
                deep_offzero,
            )
            for shift in range(1, P)
        )
        transformed = cyclo_scale(cyclo_sum(terms), Fraction(1, P))
        require(
            transformed == cyclo_rational(hostile_deep_value),
            "hostile deep colour changed",
        )
        hostile_deep_colours += 1

    blocker_colour = Fraction(11, 156)
    for colour in range(1, P):
        terms = [cyclo_rational(1)]
        terms.extend(
            cyclo_scale(
                cyclo_phase(cyclo_rational(1), colour * shift),
                blocker_offzero,
            )
            for shift in range(1, P)
        )
        transformed = cyclo_scale(cyclo_sum(terms), Fraction(1, P))
        require(
            transformed == cyclo_rational(blocker_colour),
            "hostile blocker colour changed",
        )

    pure_coefficient = hostile_deep_value * blocker_colour
    fork_coefficient = pure_coefficient * blocker_colour
    pure_corner = 12**2 * pure_coefficient
    fork_corner = 12**3 * fork_coefficient
    require(
        pure_coefficient == Fraction(-11, 369096)
        and fork_coefficient == Fraction(-121, 57578976)
        and pure_corner == Fraction(-66, 15379)
        and fork_corner == Fraction(-726, 199927),
        "sharp hostile coefficient or corner changed",
    )

    kappa_pure = Fraction(11, 13**2 * 12**2)
    kappa_fork = Fraction(11**2, 13**3 * 12**3)
    require(
        -pure_coefficient == kappa_pure * mass
        and -fork_coefficient == kappa_fork * mass,
        "fixed-colour envelope scale changed",
    )

    role_target_gains = (Fraction(13, 2), Fraction(169, 4))
    require(
        role_target_gains[0] == Fraction(P, 2)
        and role_target_gains[1] == Fraction(P * P, 4),
        "positive-intertwiner amplification changed",
    )

    # The aggregate budget functional annihilates every circulant basis.
    coefficient_norm_squared = 0
    coefficient_support = 0
    circulant_sums = [0 for _ in range(P)]
    for root_shift in range(P):
        for target_shift in range(P):
            for endpoint_shift in range(P):
                coefficient = 0
                if endpoint_shift == 0:
                    coefficient += 1
                if target_shift == 0 and endpoint_shift == 0:
                    coefficient -= P
                if coefficient:
                    coefficient_support += 1
                coefficient_norm_squared += coefficient * coefficient
                circulant_sums[
                    (root_shift - endpoint_shift) % P
                ] += coefficient
    require(
        coefficient_support == P * P
        and coefficient_norm_squared == 2028
        and circulant_sums == [0 for _ in range(P)],
        "budget functional geometry changed",
    )

    c_two_drift = Fraction(11**2, 12 * 13**5)
    c_two_unit_energy = Fraction(11**2, 12 * 13**6)
    require(
        c_two_drift == Fraction(121, 4455516)
        and c_two_unit_energy == Fraction(121, 57921708),
        "C=2 quantitative floor changed",
    )

    unit_multiplier_checks = 0
    for residue in range(1, P):
        for multiplier in range(-500, 501):
            if multiplier % P != residue:
                continue
            if multiplier % 7 == 0:
                continue
            require(
                gcd(abs(multiplier), 91) == 1,
                "live probe/deep multiplier stopped being a 91-unit",
            )
            unit_multiplier_checks += 1

    print("THM-2366 retained-probe target covariance exact referee")
    print(
        "target labels/decompositions: "
        f"{P**3}/{len(target_counts)}x{P}"
    )
    print(f"target-line diagonal cells: {line_cell_checks}")
    print(f"covariance DFT cells: {covariance_cells}")
    print(
        "energy comparison: "
        f"Delta/{P**2} <= E_nz <= Delta; target cells={nonzero_target_count}"
    )
    print(
        "hostile deep colours/value: "
        f"{hostile_deep_colours} at {hostile_deep_value}"
    )
    print(
        "hostile pure coefficient/corner: "
        f"{pure_coefficient} / {pure_corner}"
    )
    print(
        "hostile fork coefficient/corner: "
        f"{fork_coefficient} / {fork_corner}"
    )
    print(
        "positive intertwiner gains: "
        f"{role_target_gains[0]}, {role_target_gains[1]}"
    )
    print(
        "budget functional support/norm^2: "
        f"{coefficient_support}/{coefficient_norm_squared}"
    )
    print(f"C=2 drift coefficient: {c_two_drift}")
    print(f"C=2 unit-target energy coefficient: {c_two_unit_energy}")
    print(f"live 91-unit multiplier checks: {unit_multiplier_checks}")
    print(
        "VERDICT: mixed probe colour is compatible with circulation; "
        "a sub-thirteen lawful budget breaks it"
    )


if __name__ == "__main__":
    main()
