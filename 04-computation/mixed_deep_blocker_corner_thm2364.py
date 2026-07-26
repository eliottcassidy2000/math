#!/usr/bin/env python3
"""Exact companion for THM-2364's mixed thirteen-shift corner."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from typing import Dict, List, Tuple


P = 13
NONZERO = P - 1
HALF_DANGER = Fraction(1, 14)
Matrix = List[List[Fraction]]
Bits = Tuple[int, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def fractional_part(x: Fraction) -> Fraction:
    return x - x.numerator // x.denominator


def centered_distance(x: Fraction) -> Fraction:
    value = fractional_part(x)
    return min(value, 1 - value)


def danger(x: Fraction) -> int:
    return int(centered_distance(x) < HALF_DANGER)


def target_shift_danger(y: Fraction, shift: int) -> int:
    return danger(y + Fraction(shift, P))


def corner_sum(matrix: Matrix) -> Fraction:
    g_total = sum(matrix[r][0] for r in range(P))
    h_total = sum(sum(row) for row in matrix)
    return -g_total / P + h_total / (P * P)


def interaction_energy(matrix: Matrix) -> Fraction:
    row_means = [sum(row) / P for row in matrix]
    column_means = [
        sum(matrix[r][s] for r in range(P)) / P for s in range(P)
    ]
    grand_mean = sum(row_means) / P
    residual_square = sum(
        (
            matrix[r][s]
            - row_means[r]
            - column_means[s]
            + grand_mean
        )
        ** 2
        for r in range(P)
        for s in range(P)
    )
    # Normalized two-dimensional Parseval.
    return residual_square / (P * P)


def check_profile(matrix: Matrix, rho: Fraction, label: str) -> None:
    require(
        len(matrix) == P and all(len(row) == P for row in matrix),
        f"{label}: profile stopped being 13 by 13",
    )
    require(
        all(value >= 0 for row in matrix for value in row),
        f"{label}: profile lost nonnegativity",
    )
    require(
        all(matrix[0][s] == 0 for s in range(P)),
        f"{label}: deepest anchored-zero row changed",
    )

    g = [matrix[r][0] for r in range(P)]
    g_total = sum(g)
    h_total = sum(sum(row) for row in matrix)
    require(
        2 * rho >= g_total >= rho > 0,
        f"{label}: deep successor/cover mass bounds failed",
    )
    deep_successor = 2 * rho - g_total
    require(
        0 <= deep_successor <= rho
        and g_total == 2 * rho - deep_successor,
        f"{label}: exact deep successor identity failed",
    )
    for r in range(P):
        require(
            sum(matrix[r]) <= 2 * g[r],
            f"{label}: blocker translate count exceeded two",
        )

    successor_correction = 2 * g_total - h_total
    require(
        successor_correction >= 0,
        f"{label}: successor correction became negative",
    )
    mixed_sum = corner_sum(matrix)
    require(
        mixed_sum == -(11 * g_total + successor_correction) / (P * P),
        f"{label}: exact corner identity changed",
    )
    require(
        mixed_sum <= -Fraction(11, P * P) * rho,
        f"{label}: mixed-colour sum floor failed",
    )

    energy = interaction_energy(matrix)
    require(
        energy >= mixed_sum**2 / (NONZERO * NONZERO),
        f"{label}: mixed Parseval/Cauchy floor failed",
    )
    require(
        energy
        >= Fraction(121, P**4 * NONZERO**2) * rho**2,
        f"{label}: uniform mixed-energy floor failed",
    )
def deterministic_profile(seed: int) -> tuple[Matrix, Fraction]:
    matrix = [[Fraction(0) for _ in range(P)] for _ in range(P)]
    for r in range(1, P):
        g = Fraction(1 + ((7 * seed + 5 * r) % 19), 1000 + 13 * seed)
        matrix[r][0] = g
        extra = g * Fraction((seed + 3 * r) % 13, 12)
        weights = [1 + ((seed + r * s + 2 * s) % 11) for s in range(1, P)]
        weight_sum = sum(weights)
        for s, weight in enumerate(weights, start=1):
            matrix[r][s] = extra * Fraction(weight, weight_sum)
    g_total = sum(matrix[r][0] for r in range(P))
    rho = g_total * Fraction(7 + seed % 6, 12)
    return matrix, rho


def product_int(values: Tuple[int, ...]) -> int:
    result = 1
    for value in values:
        result *= value
    return result


def successor_corner(weights: Dict[Bits, Fraction], blockers: int) -> Fraction:
    """Integrate the exact pointwise fully mixed corner factor."""
    denominator = P ** (blockers + 1)
    total = Fraction(0)
    for bits, weight in weights.items():
        require(
            len(bits) == blockers + 1 and all(bit in (0, 1) for bit in bits),
            "invalid successor-bit word",
        )
        deep_factor = 2 - bits[0]
        blocker_factor = product_int(
            tuple(11 + bit for bit in bits[1:])
        )
        total -= weight * Fraction(deep_factor * blocker_factor, denominator)
    return total


def check_successor_mixture(
    weights: Dict[Bits, Fraction], blockers: int, label: str
) -> None:
    rho = sum(weights.values())
    require(rho > 0, f"{label}: successor mixture lost mass")
    corner = successor_corner(weights, blockers)
    floor = Fraction(11**blockers, P ** (blockers + 1)) * rho
    require(corner <= -floor, f"{label}: full-corner sign floor failed")
    colour_count = NONZERO ** (blockers + 1)
    energy = corner**2 / colour_count
    require(
        energy
        >= Fraction(
            11 ** (2 * blockers),
            P ** (2 * (blockers + 1)) * colour_count,
        )
        * rho**2,
        f"{label}: fully mixed energy floor failed",
    )
    if blockers == 1:
        g_total = sum(
            weight * (2 - bits[0]) for bits, weight in weights.items()
        )
        correction = sum(
            weight * (2 - bits[0]) * bits[1]
            for bits, weight in weights.items()
        )
        blocker_successor = sum(
            weight * bits[1] for bits, weight in weights.items()
        )
        require(
            correction >= blocker_successor
            and corner == -(11 * g_total + correction) / (P * P),
            f"{label}: t=1 successor sharpening failed",
        )


def main() -> None:
    # Exhaust the common open cells of the shifted danger comb and its
    # successor bit.
    boundaries = {
        fractional_part(sign * HALF_DANGER - Fraction(shift, P))
        for shift in range(P)
        for sign in (-1, 1)
    }
    require(len(boundaries) == 26, "wrong shifted-danger boundary count")
    ordered = sorted(boundaries)
    samples = []
    for index, left in enumerate(ordered):
        if index + 1 < len(ordered):
            right = ordered[index + 1]
        else:
            right = ordered[0] + 1
        samples.append(fractional_part((left + right) / 2))
    for y in samples:
        require(
            sum(target_shift_danger(y, s) for s in range(P))
            == 2 - danger(P * y),
            f"successor count changed at {y}",
        )

    for seed in range(64):
        matrix, rho = deterministic_profile(seed)
        check_profile(matrix, rho, f"control-{seed}")

    successor_mixtures = 0
    for blockers in (1, 2):
        bit_words = tuple(product((0, 1), repeat=blockers + 1))
        for seed in range(32):
            weights = {
                bits: Fraction(
                    1
                    + (
                        seed
                        + sum((index + 2) * bit for index, bit in enumerate(bits))
                    )
                    % 17,
                    1000 + 19 * seed,
                )
                for bits in bit_words
            }
            check_successor_mixture(
                weights, blockers, f"successor-{blockers}-{seed}"
            )
            successor_mixtures += 1

    # Sharpness for the abstract profile facts.  The deep shift is uniform
    # over the twelve nonzero rows, and the blocker row has one anchored
    # entry plus a uniform second unit of mass.
    g = Fraction(1, 182)
    sharp = [[Fraction(0) for _ in range(P)] for _ in range(P)]
    for r in range(1, P):
        sharp[r][0] = g
        for s in range(1, P):
            sharp[r][s] = g / NONZERO
    rho = NONZERO * g
    check_profile(sharp, rho, "sharp")
    sharp_corner = corner_sum(sharp)
    sharp_coefficient = -Fraction(11, P * P * NONZERO * NONZERO) * rho
    sharp_energy = interaction_energy(sharp)
    require(
        sharp_corner == NONZERO * NONZERO * sharp_coefficient,
        "sharp mixed coefficients stopped summing to the corner value",
    )
    require(
        sharp_energy == NONZERO * NONZERO * sharp_coefficient**2,
        "sharp profile stopped attaining the interaction-energy floor",
    )

    eta = Fraction(2593, 90090)
    all_row_coefficient = Fraction(11, P * P * NONZERO * NONZERO) * eta / 6
    fork_coefficient = (
        Fraction(11**2, P**3 * NONZERO**3) * eta / 6
    )
    require(
        all_row_coefficient == Fraction(2593, 1195871040),
        "all-row coefficient constant changed",
    )
    require(
        fork_coefficient == Fraction(28523, 186555882240),
        "all-row fork coefficient constant changed",
    )

    print("THM-2364 mixed deep/blocker corner exact companion")
    print("shifted-danger open cells checked: 26")
    print("deterministic anchored profiles checked: 64")
    print(f"one/two-blocker successor mixtures checked: {successor_mixtures}")
    print("fully mixed colours: pure 144; fork 1728")
    print("91-unit residue tuples: pure 5184; fork 373248")
    print(f"sharp profile mass: {rho}")
    print(f"sharp mixed-colour sum: {sharp_corner}")
    print(f"each sharp mixed coefficient: {sharp_coefficient}")
    print(f"sharp mixed energy: {sharp_energy}")
    print(f"all-row pure coefficient floor: {all_row_coefficient} * e_j")
    print(f"all-row fork coefficient floor: {fork_coefficient} * e_j")
    print("VERDICT: anchored zero/full word forces a fully mixed colour")


if __name__ == "__main__":
    main()
