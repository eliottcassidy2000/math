#!/usr/bin/env python3
"""Exact companion for THM-2413.

Verifies the affine-drift representation of prime-index lines, the complete
slope-two/twin-center structure, reciprocal term identities, and a small
independent exact line-cover prefix.
"""

from functools import lru_cache
from fractions import Fraction as F
from math import gcd, isqrt


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primes_through(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[0:2] = b"\x00\x00"
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start::prime] = b"\x00" * (((limit - start) // prime) + 1)
    return [n for n in range(2, limit + 1) if sieve[n]]


def normalized_line(i: int, pi: int, j: int, pj: int) -> tuple[int, int, int]:
    """Return coprime (a,b,c), b>0, for by-ax=c through two points."""
    a = pj - pi
    b = j - i
    common = gcd(abs(a), abs(b))
    a //= common
    b //= common
    if b < 0:
        a, b = -a, -b
    return a, b, b * pi - a * i


def maximal_line_masks(prime_prefix: tuple[int, ...]) -> tuple[int, ...]:
    n = len(prime_prefix)
    masks: set[int] = set()
    for i in range(1, n + 1):
        # A horizontal line is a lawful singleton replacement.
        masks.add(1 << (i - 1))
        for j in range(i + 1, n + 1):
            a, b, c = normalized_line(i, prime_prefix[i - 1], j, prime_prefix[j - 1])
            mask = 0
            for k, prime in enumerate(prime_prefix, start=1):
                if b * prime - a * k == c:
                    mask |= 1 << (k - 1)
            masks.add(mask)
    # A proper subset is never preferable in an unweighted cover.
    ordered = sorted(masks, key=int.bit_count, reverse=True)
    maximal: list[int] = []
    for mask in ordered:
        if not any(mask & existing == mask for existing in maximal):
            maximal.append(mask)
    return tuple(maximal)


def minimum_line_cover(prime_prefix: tuple[int, ...]) -> int:
    n = len(prime_prefix)
    masks = maximal_line_masks(prime_prefix)
    point_masks = tuple(
        tuple(mask for mask in masks if mask & (1 << point))
        for point in range(n)
    )

    @lru_cache(maxsize=None)
    def solve(uncovered: int) -> int:
        if uncovered == 0:
            return 0

        max_gain = max((mask & uncovered).bit_count() for mask in masks)
        lower = (uncovered.bit_count() + max_gain - 1) // max_gain

        live_points = [
            point for point in range(n) if uncovered & (1 << point)
        ]
        pivot = min(
            live_points,
            key=lambda point: sum(
                1 for mask in point_masks[point] if mask & uncovered
            ),
        )
        branches = sorted(
            point_masks[pivot],
            key=lambda mask: (mask & uncovered).bit_count(),
            reverse=True,
        )

        best = uncovered.bit_count()
        for mask in branches:
            gain = (mask & uncovered).bit_count()
            if gain == 0:
                continue
            candidate = 1 + solve(uncovered & ~mask)
            if candidate < best:
                best = candidate
            if best == lower:
                break
        return best

    return solve((1 << n) - 1)


primes = primes_through(200_000)
require(len(primes) >= 10_000, "prime sieve too short")


# ---------------------------------------------------------------------------
# 1. Rational lines are exactly affine-drift level sets.
# ---------------------------------------------------------------------------

drift_incidence_checks = 0
line_keys: set[tuple[int, int, int]] = set()
for i in range(1, 61):
    for j in range(i + 1, 61):
        a, b, c = normalized_line(i, primes[i - 1], j, primes[j - 1])
        require(gcd(abs(a), b) == 1 and b > 0, "line normalization failed")
        members = tuple(
            k
            for k in range(1, 81)
            if b * primes[k - 1] - a * k == c
        )
        require(i in members and j in members, "drift level lost defining point")
        for k in range(1, 81):
            determinant = (
                (j - i) * (primes[k - 1] - primes[i - 1])
                - (k - i) * (primes[j - 1] - primes[i - 1])
            )
            require(
                (k in members) == (determinant == 0),
                "drift incidence and collinearity disagree",
            )
            drift_incidence_checks += 1
        line_keys.add((a, b, c))

# Equally spaced collinearity is exactly a vanishing second difference.
second_difference_checks = 0
for step in range(1, 13):
    for i in range(1, 81 - 2 * step):
        second = (
            primes[i + 2 * step - 1]
            - 2 * primes[i + step - 1]
            + primes[i - 1]
        )
        left = normalized_line(
            i,
            primes[i - 1],
            i + step,
            primes[i + step - 1],
        )
        right = normalized_line(
            i + step,
            primes[i + step - 1],
            i + 2 * step,
            primes[i + 2 * step - 1],
        )
        require((second == 0) == (left[:2] == right[:2]), "Delta_h^2 test")
        second_difference_checks += 1


# ---------------------------------------------------------------------------
# 2. The slope-two drift and twin-prime center weld.
# ---------------------------------------------------------------------------

drift_two = [prime - 2 * n for n, prime in enumerate(primes, start=1)]
require(drift_two[:6] == [0, -1, -1, -1, 1, 1], "initial H_2 drift")
require(drift_two[1] < drift_two[0], "the initial descent was lost")

twin_centers: list[int] = []
plateau_edges: list[int] = []
for n in range(2, 10_000):
    difference = drift_two[n] - drift_two[n - 1]
    prime_gap = primes[n] - primes[n - 1]
    require(difference == prime_gap - 2, "H_2 derivative law failed")
    require(difference >= 0, "H_2 was not nondecreasing after n=2")
    require((difference == 0) == (prime_gap == 2), "plateau/twin mismatch")
    if prime_gap == 2:
        center = primes[n - 1] + 1
        plateau_edges.append(n)
        twin_centers.append(center)
        require(primes[n - 1] + primes[n] == 2 * center, "summand weld")
        require(primes[n - 1] * primes[n] == center * center - 1,
                "multiplicand weld")
        if center != 4:
            require(center % 6 == 0, "nonexceptional twin center not 0 mod 6")

require(
    twin_centers[:12] == [4, 6, 12, 18, 30, 42, 60, 72, 102, 108, 138, 150],
    "A014574 prefix mismatch",
)

# A nondecreasing drift can repeat at three indices only across two
# consecutive twin gaps. The unique prime triple p,p+2,p+4 is 3,5,7.
level_members: dict[int, list[int]] = {}
for n, value in enumerate(drift_two[1:10_000], start=2):
    level_members.setdefault(value, []).append(n)
largest_slope_two_fibre = max(map(len, level_members.values()))
triple_fibres = {
    value: members for value, members in level_members.items() if len(members) >= 3
}
require(largest_slope_two_fibre == 3, "unexpected long slope-two line")
require(triple_fibres == {-1: [2, 3, 4]}, "unique 3,5,7 slope-two line failed")

# A014574 is the additive fixed-gap pattern cut out by multiplicative atoms.
atom_pattern_checks = 0
prime_set = set(primes)
twin_center_set = set(twin_centers)
for center in range(3, 20_000):
    is_center = center - 1 in prime_set and center + 1 in prime_set
    require(
        is_center == (center in twin_center_set),
        "adjacent multiplicative-atom characterization failed",
    )
    atom_pattern_checks += 1


# ---------------------------------------------------------------------------
# 3. Reciprocal support identity and the non-twin slope boundary.
# ---------------------------------------------------------------------------

reciprocal_checks = 0
twin_half_sum = F(0)
center_sum = F(0)
correction_sum = F(0)
for center in twin_centers[:200]:
    half_pair = F(1, 2) * (F(1, center - 1) + F(1, center + 1))
    correction = F(1, center * (center * center - 1))
    require(half_pair == F(1, center) + correction, "reciprocal weld")
    twin_half_sum += half_pair
    center_sum += F(1, center)
    correction_sum += correction
    reciprocal_checks += 1
require(twin_half_sum == center_sum + correction_sum, "partial reciprocal sum")

# The prime graph has substantially richer affine patches than slope two.
slope_four_indices = tuple(
    n for n in range(1, 10_001) if primes[n - 1] == 4 * n - 11
)
require(
    slope_four_indices[:8] == (6, 7, 10, 12, 13, 16, 18, 21),
    "slope-four eight-point control failed",
)
require(len(slope_four_indices[:8]) > largest_slope_two_fibre,
        "slope-two hostile boundary disappeared")


# ---------------------------------------------------------------------------
# 4. Independent exact small A373813 cover.
# ---------------------------------------------------------------------------

expected_cover_prefix = (
    1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4,
    4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7,
)
computed_cover_prefix = tuple(
    minimum_line_cover(tuple(primes[:n]))
    for n in range(1, len(expected_cover_prefix) + 1)
)
require(computed_cover_prefix == expected_cover_prefix, "A373813 exact prefix")
require(
    all(
        computed_cover_prefix[n] - computed_cover_prefix[n - 1] in (0, 1)
        for n in range(1, len(computed_cover_prefix))
    ),
    "line-cover increments were not binary",
)


print("theorem=THM-2413")
print("status=PROVED+FINITE-EXACT;LITERATURE-CONSEQUENCES-CITED")
print(f"prime_count={len(primes)};drift_line_keys={len(line_keys)}")
print(f"drift_incidence_checks={drift_incidence_checks}")
print(f"second_difference_checks={second_difference_checks}")
print("H2_initial=0,-1,-1,-1,1,1;unique_initial_descent=n1_to_n2")
print("H2_delta=prime_gap-2;post_n2=nonnegative;zero=exactly-twin")
print("A014574_prefix=4,6,12,18,30,42,60,72,102,108,138,150")
print(f"twin_centers_checked={len(twin_centers)}")
print(f"atom_pattern_checks={atom_pattern_checks}")
print("unique_slope2_triple=indices_2,3,4|primes_3,5,7|line_y=2x-1")
print(f"reciprocal_term_checks={reciprocal_checks}")
print("reciprocal_identity=half_twin_pair=1/center+1/(center*(center^2-1))")
print("slope4_control=6,7,10,12,13,16,18,21|line_y=4x-11")
print("A373813_exact_prefix=" + ",".join(map(str, computed_cover_prefix)))
print("hostiles=singleton-line;slope2-isolated-pair;initial-H2-descent")
print("all_checks=PASS")
