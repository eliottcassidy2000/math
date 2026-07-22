#!/usr/bin/env python3
"""Exact arithmetic referee for THM-2073's LRC(14) dyadic tower."""

from __future__ import annotations

from math import floor, gcd


def require(condition: bool, message: str = "exact referee failure") -> None:
    if not condition:
        raise AssertionError(message)


def abs_residue(k: int, modulus: int) -> int:
    r = k % modulus
    return min(r, modulus - r)


def nearest_integer(numerator: int, denominator: int) -> int:
    """Unique nearest integer when the caller has distance <1/2."""
    return (2 * numerator + denominator) // (2 * denominator)


def check_capacities(limit: int = 100_000) -> tuple[list[int], list[int]]:
    alpha_equal: list[int] = []
    beta_equal: list[int] = []
    for d in range(2, limit + 1):
        alpha_num = floor(d / 7) + 1
        require(2 * alpha_num <= d, f"alpha exceeds 1/2 at d={d}")
        if 2 * alpha_num == d:
            alpha_equal.append(d)

        beta_num = floor(2 * d / 7) + 1
        require(2 * beta_num <= d, f"beta exceeds 1/2 at d={d}")
        if 2 * beta_num == d:
            beta_equal.append(d)
    require(alpha_equal == [2])
    require(beta_equal == [2, 4])
    return alpha_equal, beta_equal


def universal_odd_classes(modulus: int, multiplier: int) -> list[int]:
    """Odd z mod 2m with 7|za|_m < multiplier*m for all units a."""
    answer = []
    for z in range(1, 2 * modulus, 2):
        if all(
            7 * abs_residue(z * a, modulus) < multiplier * modulus
            for a in range(1, modulus)
            if gcd(a, modulus) == 1
        ):
            answer.append(z)
    return answer


def check_residue_shells() -> tuple[dict[int, list[int]], dict[int, list[int]], dict[int, list[int]]]:
    guard = {m: universal_odd_classes(m, 1) for m in range(2, 15)}
    expected_guard = {
        m: ([] if m % 2 == 0 else [m])
        for m in range(2, 15)
    }
    require(guard == expected_guard, f"guard shell mismatch: {guard}")

    singleton = {m: universal_odd_classes(m, 2) for m in range(3, 15, 2)}
    expected_singleton = {m: [m] for m in range(3, 15, 2)}
    require(singleton == expected_singleton, f"singleton shell mismatch: {singleton}")

    doubled = {2 * m: universal_odd_classes(2 * m, 1) for m in range(3, 15, 2)}
    expected_doubled = {2 * m: [] for m in range(3, 15, 2)}
    require(doubled == expected_doubled, f"doubled shell mismatch: {doubled}")
    return guard, singleton, doubled


def check_four_lift_formulas(max_denominator: int = 80, max_odd_speed: int = 49) -> tuple[int, int]:
    guard_checks = 0
    singleton_checks = 0
    for q in range(2, max_denominator + 1):
        modulus = 4 * q
        for p in range(q):
            for z in range(1, max_odd_speed + 1, 2):
                residue = abs_residue(z * p, q)
                nearest = nearest_integer(z * p, q)

                direct_guard = {
                    j for j in range(4)
                    if 14 * abs_residue(2 * z * (p + j * q), modulus) < modulus
                }
                if 7 * residue < q:
                    predicted_guard = {j for j in range(4) if j % 2 == nearest % 2}
                    require(direct_guard == predicted_guard, (
                        "guard ownership mismatch", q, p, z,
                        direct_guard, predicted_guard,
                    ))
                    guard_checks += 1
                else:
                    require(len(direct_guard) == 0)

                direct_singleton = {
                    j for j in range(4)
                    if 14 * abs_residue(z * (p + j * q), modulus) < modulus
                }
                if 7 * residue < 2 * q:
                    predicted_singleton = {(-pow(z, -1, 4) * nearest) % 4}
                    require(direct_singleton == predicted_singleton, (
                        "singleton ownership mismatch", q, p, z,
                        direct_singleton, predicted_singleton,
                    ))
                    singleton_checks += 1
                else:
                    require(len(direct_singleton) == 0)
    return guard_checks, singleton_checks


def check_antipodal_exclusion(max_denominator: int = 200, max_odd_speed: int = 99) -> int:
    checks = 0
    for q in range(2, max_denominator + 1):
        modulus = 2 * q
        for p in range(q):
            for z in range(1, max_odd_speed + 1, 2):
                first = abs_residue(z * p, q)
                shifted = abs_residue(z * (2 * p + q), modulus)
                # shifted/modulus is ||z(p/q+1/2)||.
                require(not (7 * first < q and 7 * shifted < modulus))
                checks += 1
    return checks


def main() -> None:
    print("LRC14 DYADIC DELETION-TOWER AUDIT -- exact arithmetic")
    alpha_equal, beta_equal = check_capacities()
    print(f"capacity equality: alpha={alpha_equal}, beta={beta_equal} PASS")

    guard, singleton, doubled = check_residue_shells()
    print("1/7 all-unit guard shell:", ", ".join(
        f"{m}:{guard[m]}" for m in guard
    ))
    print("2/7 all-unit singleton shell:", ", ".join(
        f"{m}:{singleton[m]}" for m in singleton
    ))
    print("doubled-denominator guard shell:", ", ".join(
        f"{m}:{doubled[m]}" for m in doubled
    ))

    guard_checks, singleton_checks = check_four_lift_formulas()
    print(f"four-lift formulas: {guard_checks} guard and {singleton_checks} singleton eligible cases PASS")
    antipodal_checks = check_antipodal_exclusion()
    print(f"antipodal odd-tail exclusions: {antipodal_checks} PASS")
    print("PASS")


if __name__ == "__main__":
    main()
