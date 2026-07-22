#!/usr/bin/env python3
"""Exact finite referee for THM-2095's p-adic valuation wheel."""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import gcd, prod


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primes_through(n: int) -> list[int]:
    primes: list[int] = []
    for value in range(2, n + 1):
        if all(value % p for p in primes if p * p <= value):
            primes.append(value)
    return primes


def capacity(p: int, f: int, profile: tuple[int, ...]) -> int:
    return sum(
        p**v * (p ** (f - v) // 7 + 1)
        for v in profile
        if v < f
    )


def profile_is_feasible(p: int, e: int, profile: tuple[int, ...]) -> bool:
    require(len(profile) == 6, "profile length")
    require(tuple(sorted(profile)) == profile, "profile order")
    require(profile[0] == profile[1] == 0, "two-unit condition")
    require(profile[-1] <= e, "profile cap")
    return all(capacity(p, f, profile) >= p**f for f in range(1, e + 1))


def feasible_profiles(p: int, e: int) -> list[tuple[int, ...]]:
    profiles: list[tuple[int, ...]] = []
    for tail in combinations_with_replacement(range(e + 1), 4):
        profile = (0, 0) + tail
        if profile_is_feasible(p, e, profile):
            profiles.append(profile)
    return profiles


def first_layer_allowed_primes() -> list[int]:
    allowed: list[int] = []
    for p in primes_through(41):
        if any(p <= k * (p // 7 + 1) for k in range(2, 7)):
            allowed.append(p)
    return allowed


def main() -> None:
    expected = {
        2: (8, (0, 0, 1, 2, 2, 5), [5, 13, 25, 36, 36, 29, 13, 2]),
        3: (4, (0, 0, 0, 1, 2, 2), [4, 8, 7, 2]),
        5: (2, (0, 0, 0, 0, 0, 1), [2, 1]),
        7: (2, (0, 0, 0, 0, 0, 1), [3, 2]),
        11: (1, (0, 0, 0, 0, 0, 0), [1]),
        17: (1, (0, 0, 0, 0, 0, 0), [1]),
        23: (1, (0, 0, 0, 0, 0, 0), [1]),
        29: (1, (0, 0, 0, 0, 0, 0), [1]),
    }

    allowed = first_layer_allowed_primes()
    require(allowed == list(expected), f"first-layer prime list: {allowed}")
    # The excluded small primes are load-bearing regressions, not a generic
    # p>42 consequence.
    require(13 > 6 * (13 // 7 + 1), "p=13 should fail")
    require(19 > 6 * (19 // 7 + 1), "p=19 should fail")
    require(all(p > 42 for p in primes_through(200) if p > 41), "prime cutoff regression")

    rows: list[str] = []
    for p, (max_e, witness, expected_counts) in expected.items():
        counts: list[int] = []
        last_profiles: list[tuple[int, ...]] = []
        for e in range(1, max_e + 2):
            profiles = feasible_profiles(p, e)
            counts.append(len(profiles))
            if e == max_e:
                last_profiles = profiles
        require(counts[:-1] == expected_counts, f"count ledger p={p}: {counts}")
        require(counts[-1] == 0, f"first infeasible level p={p}")
        require(witness in last_profiles, f"missing witness p={p}")

        # If a higher profile existed, capping it at max_e+1 would contradict
        # the just-checked empty level. This direct truncation audit checks the
        # capacity identity on all last-level witnesses.
        for profile in last_profiles:
            require(profile_is_feasible(p, max_e, profile), f"last profile p={p}")
        rows.append(
            f"p={p} max_e={max_e} counts={','.join(map(str, expected_counts))},0 "
            f"witness={','.join(map(str, witness))}"
        )

    general_scale = prod(p**data[0] for p, data in expected.items())
    live_scale = 3**4 * 5**2 * 11 * 17 * 23 * 29
    live_pair_bound = 57 * live_scale
    require(general_scale == 2**8 * 3**4 * 5**2 * 7**2 * 11 * 17 * 23 * 29, "general scale")
    require(live_scale == 252576225, "live scale")
    require(live_pair_bound == 14396844825, "live pair bound")

    divisor_count = (4 + 1) * (2 + 1) * 2**4
    ratio_pairs = [
        (r, s)
        for s in range(1, 58, 2)
        if s % 7 != 0
        for r in range(1, 58)
        if gcd(r, s) == 1 and (r, s) != (1, 1)
    ]
    require(divisor_count == 240, "live divisor count")
    require(len(ratio_pairs) == 1165, "live ratio count")
    require(divisor_count * len(ratio_pairs) == 279600, "marked pair ledger")

    diagonal_sum = F(1, 7) + 4 * F(1, 35) + F(1, 42) + F(2, 77)
    diagonal_margin = diagonal_sum - F(2, 7)
    require(diagonal_sum == F(709, 2310), "diagonal overlap sum")
    require(diagonal_margin == F(7, 330), "diagonal margin")

    print("THM-2095 P-ADIC GUARD-RATIO REFEREE")
    print("prime_cutoff=42")
    print("first_layer_primes=" + ",".join(map(str, allowed)))
    for row in rows:
        print(row)
    print(f"general_scale={general_scale}")
    print(f"live_scale={live_scale}")
    print(f"live_pair_bound={live_pair_bound}")
    print(f"live_scale_divisors={divisor_count}")
    print(f"live_ratio_pairs={len(ratio_pairs)}")
    print(f"marked_pair_ledger={divisor_count * len(ratio_pairs)}")
    print(f"diagonal_margin={diagonal_margin}")
    print("PASS")


if __name__ == "__main__":
    main()
