#!/usr/bin/env python3
"""Exact hostile controls for the THM-2148 finite-kernel argument.

This is a referee, not a dependency of the proof.  It exhausts every
rank-at-most-two invariant-factor group C_m x C_n with m | n and n <= 18.
For every triple of nontrivial characters it computes the *closed*
radius-1/14 centered masks exactly and checks:

* no two masks cover;
* every three-mask cover consists of the three nonzero order-two characters;
* such a cover occurs exactly when the group has a C_2 x C_2 quotient.

It separately enumerates the integral polar diamond
2|m|+|n| <= 2 and verifies that its nonzero points omit parity (1,1).

Reproduce:
    python3 04-computation/lrc14_three_target_fano_kernel_referee_thm2148.py
"""

from __future__ import annotations

from itertools import combinations
from math import gcd, lcm


MAX_INVARIANT_FACTOR = 18


def character_order(a: int, b: int, m: int, n: int) -> int:
    left = 1 if m == 1 else m // gcd(a, m)
    right = n // gcd(b, n)
    return lcm(left, right)


def centered_danger_mask(a: int, b: int, m: int, n: int) -> int:
    """Closed ||a*x/m+b*y/n|| <= 1/14, using denominator m*n."""
    mask = 0
    denominator = m * n
    for x in range(m):
        for y in range(n):
            residue = (a * x * n + b * y * m) % denominator
            distance = min(residue, denominator - residue)
            if 14 * distance <= denominator:
                mask |= 1 << (x * n + y)
    return mask


def audit_group(m: int, n: int) -> tuple[int, int]:
    size = m * n
    full = (1 << size) - 1
    characters: list[tuple[int, int, int, int]] = []
    for a in range(m):
        for b in range(n):
            if a == 0 and b == 0:
                continue
            characters.append(
                (
                    a,
                    b,
                    character_order(a, b, m, n),
                    centered_danger_mask(a, b, m, n),
                )
            )

    pair_covers = 0
    for left, right in combinations(characters, 2):
        if left[3] | right[3] == full:
            pair_covers += 1
    assert pair_covers == 0, (m, n, pair_covers)

    triple_covers = 0
    for first_index, first in enumerate(characters):
        for second_index in range(first_index + 1, len(characters)):
            second = characters[second_index]
            missing = full ^ (first[3] | second[3])
            if missing == 0:
                raise AssertionError(("unexpected pair cover", m, n, first, second))
            for third in characters[second_index + 1 :]:
                if missing & ~third[3]:
                    continue
                triple_covers += 1
                triple = (first, second, third)
                assert all(character[2] == 2 for character in triple)
                assert sum(character[0] for character in triple) % m == 0
                assert sum(character[1] for character in triple) % n == 0

    expected = int(m % 2 == 0 and n % 2 == 0)
    assert triple_covers == expected, (m, n, triple_covers, expected)
    return pair_covers, triple_covers


def main() -> None:
    groups = [
        (m, n)
        for n in range(1, MAX_INVARIANT_FACTOR + 1)
        for m in range(1, n + 1)
        if n % m == 0
    ]
    total_pair_covers = 0
    total_triple_covers = 0
    positive_groups = []
    for m, n in groups:
        pair_count, triple_count = audit_group(m, n)
        total_pair_covers += pair_count
        total_triple_covers += triple_count
        if triple_count:
            positive_groups.append((m, n))

    polar_points = sorted(
        (m, n)
        for m in range(-4, 5)
        for n in range(-4, 5)
        if (m, n) != (0, 0) and 2 * abs(m) + abs(n) <= 2
    )
    polar_parities = sorted({(m % 2, n % 2) for m, n in polar_points})
    assert polar_points == [
        (-1, 0),
        (0, -2),
        (0, -1),
        (0, 1),
        (0, 2),
        (1, 0),
    ]
    assert polar_parities == [(0, 0), (0, 1), (1, 0)]
    assert (1, 1) not in polar_parities

    print("THM-2148 THREE-TARGET FANO KERNEL REFEREE")
    print(f"invariant_factor_limit={MAX_INVARIANT_FACTOR}")
    print(f"groups_checked={len(groups)}")
    print(f"total_pair_covers={total_pair_covers}")
    print(f"total_triple_covers={total_triple_covers}")
    print(f"positive_groups={positive_groups}")
    print(f"polar_points={polar_points}")
    print(f"polar_parities={polar_parities}")
    print("missing_required_parity=(1, 1)")
    print("status=PASS")


if __name__ == "__main__":
    main()
