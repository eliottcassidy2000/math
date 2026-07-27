#!/usr/bin/env python3
"""Exact companion for THM-2503.

Enumerates every strict thirteen-shift profile of the active dipole sweep,
checks the pointwise count identity and sharp service census, and verifies the
rational cyclotomic kernel and exact Fourier-energy invoices.
"""

from collections import Counter
from fractions import Fraction as F


P = 13


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mod_one(value: F) -> F:
    return value % 1


def circle_distance(value: F) -> F:
    value = mod_one(value)
    return min(value, 1 - value)


def danger(value: F, width: int) -> bool:
    return circle_distance(value) < F(width, 14)


def profile_samples(
    width: int,
    plus_shift: bool,
    support,
) -> dict[tuple[int, ...], F]:
    """One exact sample from every open profile cell."""

    direction = -1 if plus_shift else 1
    boundaries = sorted(
        {
            mod_one(direction * F(ell, P) + sign * F(width, 14))
            for ell in range(P)
            for sign in (-1, 1)
        }
    )
    profiles: dict[tuple[int, ...], F] = {}
    for index, left in enumerate(boundaries):
        right = boundaries[(index + 1) % len(boundaries)]
        if index + 1 == len(boundaries):
            right += 1
        sample = mod_one((left + right) / 2)
        if not support(sample):
            continue
        profile = tuple(
            ell
            for ell in range(P)
            if danger(
                sample + (F(ell, P) if plus_shift else -F(ell, P)),
                width,
            )
        )
        profiles.setdefault(profile, sample)
    return profiles


def matrix_rank(matrix: list[list[F]]) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0]) if rows else 0
    rank = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(rank, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][column]
        work[rank] = [entry / scale for entry in work[rank]]
        for row in range(rows):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                work[row][j] - factor * work[rank][j]
                for j in range(columns)
            ]
        rank += 1
        if rank == rows:
            break
    return rank


def root_vector(exponent: int) -> tuple[F, ...]:
    """zeta^exponent in the basis 1,...,zeta^11 of Q(zeta_13)."""

    exponent %= P
    if exponent < P - 1:
        return tuple(F(int(index == exponent)) for index in range(P - 1))
    return tuple(F(-1) for _ in range(P - 1))


def evaluation_matrix(color: int) -> list[list[F]]:
    columns = [root_vector(-color * ell) for ell in range(P)]
    return [[columns[column][row] for column in range(P)] for row in range(P - 1)]


def evaluate(coefficients: list[F], color: int) -> tuple[F, ...]:
    matrix = evaluation_matrix(color)
    return tuple(
        sum(matrix[row][column] * coefficients[column] for column in range(P))
        for row in range(P - 1)
    )


def expected_v_profiles(width: int) -> set[tuple[int, ...]]:
    if width == 1:
        return {(0,), (0, 1), (0, 12)}
    return {
        (0, 1, 2),
        (0, 1, 12),
        (0, 11, 12),
        (0, 1, 2, 3),
        (0, 1, 2, 12),
        (0, 1, 11, 12),
        (0, 10, 11, 12),
    }


def expected_a_profiles() -> set[tuple[int, ...]]:
    return {
        *((ell,) for ell in range(1, P)),
        *((ell, ell + 1) for ell in range(1, P - 1)),
    }


def active_profile(
    v_profile: tuple[int, ...],
    a_profile: tuple[int, ...],
) -> tuple[int, ...]:
    v_set = set(v_profile)
    a_set = set(a_profile)
    return tuple(int(ell not in v_set and ell not in a_set) for ell in range(P))


def weighted_mixture(
    certificate: tuple[tuple[tuple[int, ...], tuple[int, ...], F], ...],
) -> tuple[list[F], F]:
    mixture = [F(0) for _ in range(P)]
    mass = F(0)
    for v_profile, a_profile, weight in certificate:
        profile = active_profile(v_profile, a_profile)
        mass += weight
        for ell in range(P):
            mixture[ell] += weight * profile[ell]
    return mixture, mass


def audit_width(
    width: int,
    a_profiles: dict[tuple[int, ...], F],
) -> tuple[Counter, list[F], F, F, int]:
    v_profiles = profile_samples(
        width,
        plus_shift=False,
        support=lambda value: danger(value, width),
    )
    require(set(v_profiles) == expected_v_profiles(width), "V profile atlas")

    census: Counter = Counter()
    active_words: list[tuple[int, ...]] = []
    for v_profile, v_sample in v_profiles.items():
        for a_profile, a_sample in a_profiles.items():
            active = tuple(
                ell
                for ell in range(P)
                if not danger(v_sample - F(ell, P), width)
                and not danger(a_sample + F(ell, P), 1)
            )
            v_set = set(v_profile)
            a_set = set(a_profile)
            count = len(active)
            formula = (
                11
                - 2 * width
                + int(danger(P * v_sample, width))
                + int(danger(P * a_sample, 1))
                + len(v_set & a_set)
            )
            require(count == P - len(v_set | a_set), "union count")
            require(count == formula, "pointwise active-dipole identity")
            require(0 not in active, "anchored zero shift")
            census[count] += 1
            active_words.append(active)

    expected = (
        Counter({9: 20, 10: 35, 11: 14})
        if width == 1
        else Counter({7: 32, 8: 69, 9: 52, 10: 8})
    )
    require(census == expected, "sharp profile census")

    hitting_set = (2, 4, 6, 8, 10) if width == 1 else (4, 6, 8)
    hitting_floor = 4 if width == 1 else 2
    for active in active_words:
        require(
            sum(int(ell in active) for ell in hitting_set) >= hitting_floor,
            "universal sharp maximum certificate",
        )

    number_words = len(active_words)
    mixture = [
        sum(F(int(ell in active), number_words) for active in active_words)
        for ell in range(P)
    ]
    total = sum(mixture)
    lower = F(11 - 2 * width)
    upper = F(12 - width)
    require(lower <= total <= upper, "integrated service bounds")
    require(mixture[0] == 0, "mixture zero shift")
    require(sum(bool(entry) for entry in mixture) >= 11 - 2 * width, "support floor")
    require(max(mixture) >= lower / 12, "maximum-shift floor")
    require(
        max(mixture) >= F(hitting_floor, len(hitting_set)),
        "sharp maximum-shift floor",
    )

    energy = F(1, P) * sum(entry * entry for entry in mixture) - total * total / (P * P)
    energy_denominator = 2028 if width == 1 else 1690
    require(
        energy >= total * total / energy_denominator,
        "nonzero Fourier energy floor",
    )
    for color in range(1, P):
        require(any(evaluate(mixture, color)), "nonzero rational root color")

    return census, mixture, total, energy, energy_denominator


def main() -> None:
    a_profiles = profile_samples(
        1,
        plus_shift=True,
        support=lambda value: not danger(value, 1),
    )
    require(set(a_profiles) == expected_a_profiles(), "A profile atlas")
    require(Counter(map(len, a_profiles)) == Counter({1: 12, 2: 11}), "A sizes")

    for color in range(1, P):
        matrix = evaluation_matrix(color)
        require(matrix_rank(matrix) == 12, "cyclotomic evaluation rank")
        require(all(sum(row) == 0 for row in matrix), "constant kernel")

    results = {}
    for width in (1, 2):
        results[width] = audit_width(width, a_profiles)

    ordinary_max_certificate = (
        ((0,), (4, 5), F(1, 5)),
        ((0,), (6, 7), F(1, 5)),
        ((0,), (10, 11), F(1, 5)),
        ((0, 1), (2, 3), F(1, 5)),
        ((0, 12), (8, 9), F(1, 5)),
    )
    ordinary_sharp, ordinary_mass = weighted_mixture(ordinary_max_certificate)
    require(ordinary_mass == 1, "ordinary sharp mass")
    require(ordinary_sharp == [F(0), *([F(4, 5)] * 12)], "ordinary sharp ray")
    ordinary_total = sum(ordinary_sharp)
    ordinary_energy = (
        F(1, P) * sum(entry * entry for entry in ordinary_sharp)
        - ordinary_total * ordinary_total / (P * P)
    )
    require(
        ordinary_energy == ordinary_total * ordinary_total / 2028,
        "ordinary sharp energy",
    )

    guard_max_certificate = (
        ((0, 1, 2, 3), (8, 9), F(1, 3)),
        ((0, 10, 11, 12), (6, 7), F(1, 3)),
        ((0, 1, 11, 12), (4, 5), F(1, 3)),
    )
    guard_max, guard_mass = weighted_mixture(guard_max_certificate)
    require(guard_mass == 1, "guard max sharp mass")
    require(max(guard_max) == F(2, 3), "guard sharp maximum")

    guard_dual = [F(value, 130) for value in (0, 9, 9, 11, 12, 12, 12, 12, 12, 12, 11, 9, 9)]
    require(sum(entry * entry for entry in guard_dual) == F(11, 130), "guard dual norm")
    gap_census: Counter = Counter()
    for v_profile in expected_v_profiles(2):
        for a_profile in expected_a_profiles():
            profile = active_profile(v_profile, a_profile)
            gap = sum(guard_dual[ell] * profile[ell] for ell in range(P))
            gap -= F(11, 130) * sum(profile)
            gap_census[gap] += 1
    require(
        gap_census
        == Counter(
            {
                F(0): 25,
                F(1, 130): 40,
                F(2, 130): 38,
                F(3, 130): 16,
                F(4, 130): 32,
                F(6, 130): 10,
            }
        ),
        "guard dual gap census",
    )

    guard_energy_certificate = (
        ((0, 10, 11, 12), (4, 5), F(1, 130)),
        ((0, 1, 2), (4, 5), F(1, 65)),
        ((0, 11, 12), (6, 7), F(1, 65)),
        ((0, 1, 2, 3), (4, 5), F(3, 130)),
        ((0, 1, 2, 3), (6, 7), F(2, 65)),
        ((0, 10, 11, 12), (8, 9), F(3, 65)),
    )
    guard_sharp, _ = weighted_mixture(guard_energy_certificate)
    require(guard_sharp == guard_dual, "guard sharp primal decomposition")
    guard_total = sum(guard_sharp)
    guard_energy = (
        F(1, P) * sum(entry * entry for entry in guard_sharp)
        - guard_total * guard_total / (P * P)
    )
    require(guard_total == 1, "guard sharp total")
    require(guard_energy == F(1, 1690), "guard sharp energy")

    print("THM-2503 ACTIVE-DIPOLE SWEEP AUDIT")
    print("a_safe_profiles=23;sizes=1:12,2:11")
    for width in (1, 2):
        census, mixture, total, energy, energy_denominator = results[width]
        census_text = ",".join(f"{key}:{census[key]}" for key in sorted(census))
        print(
            f"L={width};v_profiles={len(expected_v_profiles(width))};"
            f"profile_census={census_text};pointwise={min(census)}..{max(census)}"
        )
        print(
            f"L={width};uniform_mixture_total={total};"
            f"nonzero_shifts={sum(bool(entry) for entry in mixture)};"
            f"max_shift={max(mixture)};energy={energy};"
            f"energy_to_floor={energy / (total * total / energy_denominator)}"
        )
    print("cyclotomic_evaluation_ranks=12x12;kernel=constants")
    print("sharp_max_floors=L1:4/5,L2:2/3")
    print("negative_real_floor=-T/156;sharp_energy=L1:T^2/2028,L2:T^2/1690")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
