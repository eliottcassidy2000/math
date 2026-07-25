#!/usr/bin/env python3
"""Exact hostile controls for THM-2174.

The script verifies two target-mixed rank-two carry fibres, enumerates the
complete fixed-higher-quotient base-5 fibre, and computes exact-denominator
phase-danger clutters.  All checks remain active under ``python -O``.
"""

from collections import Counter, defaultdict
from itertools import permutations
from math import gcd


R = (1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1)
S = (1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1)

A = (6, 12, 18, 13, 23, 21, 26, 41, 51, 22, 32, 37, 48)
B = (6, 12, 18, 13, 23, 20, 27, 41, 51, 22, 32, 37, 48)

U = (6, 11, 16, 12, 21, 22, 26, 42, 51, 23, 32, 37, 49)
W = (6, 11, 16, 12, 21, 22, 26, 41, 51, 23, 31, 37, 49)

Z1 = (1, 2, 3, 2, 4, 4, 5, 8, 10, 4, 6, 7, 9)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(x * y for x, y in zip(left, right))


def circular_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def exact_margin(row: tuple[int, ...], modulus: int) -> tuple[int, tuple[int, ...]]:
    values = []
    for numerator in range(1, modulus):
        if gcd(numerator, modulus) == 1:
            values.append(
                (
                    min(circular_residue(numerator * value, modulus) for value in row),
                    numerator,
                )
            )
    maximum = max(value for value, _ in values)
    return maximum, tuple(numerator for value, numerator in values if value == maximum)


def unrestricted_margin(row: tuple[int, ...], modulus: int) -> int:
    return max(
        min(circular_residue(numerator * value, modulus) for value in row)
        for numerator in range(1, modulus)
    )


def danger_clutter(row: tuple[int, ...], modulus: int) -> tuple[tuple[int, ...], ...]:
    bad_sets = set()
    for numerator in range(1, modulus):
        if gcd(numerator, modulus) != 1:
            continue
        bad_sets.add(
            frozenset(
                i + 1
                for i, value in enumerate(row)
                if 14 * circular_residue(numerator * value, modulus) < modulus
            )
        )
    minimal = {bad for bad in bad_sets if not any(other < bad for other in bad_sets)}
    return tuple(sorted((tuple(sorted(bad)) for bad in minimal), key=lambda x: (len(x), x)))


def radix_path(
    row: tuple[int, ...], base: int = 5
) -> tuple[tuple[tuple[int, int], tuple[int, ...], tuple[int, ...]], ...]:
    path = []
    power = 1
    while True:
        quotient = tuple(value // power for value in row)
        digit = tuple((value // power) % base for value in row)
        path.append(
            (
                (-dot(R, quotient), -dot(S, quotient)),
                tuple(i + 1 for i, value in enumerate(quotient) if value),
                tuple(i + 1 for i, value in enumerate(digit) if value == 0),
            )
        )
        if not any(quotient):
            return tuple(path)
        power *= base


def rank_two_mod_five() -> bool:
    return any(
        (R[i] * S[j] - R[j] * S[i]) % 5
        for i in range(13)
        for j in range(i + 1, 13)
    )


def enumerate_fixed_quotient_fibre() -> tuple[Counter[int], int]:
    """Count primitive distinct v=5*Z1+D satisfying both relations."""
    blocks: dict[int, list[int]] = {}
    for i, quotient in enumerate(Z1):
        blocks.setdefault(quotient, []).append(i)

    block_options = []
    for block in blocks.values():
        options = []
        digit_words = (
            permutations(range(5), len(block))
            if len(block) > 1
            else ((digit,) for digit in range(5))
        )
        for digits in digit_words:
            delta_r = sum(R[i] * digit for i, digit in zip(block, digits))
            delta_s = sum(S[i] * digit for i, digit in zip(block, digits))
            mask = sum(1 << i for i, digit in zip(block, digits) if digit == 0)
            common_gcd = 0
            for i, digit in zip(block, digits):
                common_gcd = gcd(common_gcd, 5 * Z1[i] + digit)
            options.append((delta_r, delta_s, mask, common_gcd))
        block_options.append(options)

    dynamic = {(0, 0, 0, 0): 1}
    for options in block_options:
        next_dynamic: defaultdict[tuple[int, int, int, int], int] = defaultdict(int)
        for (sum_r, sum_s, mask, common_gcd), count in dynamic.items():
            for delta_r, delta_s, new_mask, block_gcd in options:
                next_dynamic[
                    (
                        sum_r + delta_r,
                        sum_s + delta_s,
                        mask | new_mask,
                        gcd(common_gcd, block_gcd),
                    )
                ] += count
        dynamic = next_dynamic

    masks: Counter[int] = Counter()
    for (sum_r, sum_s, mask, common_gcd), count in dynamic.items():
        if (sum_r, sum_s) == (-5, 5) and common_gcd == 1:
            masks[mask] += count
    return masks, len(dynamic)


def main() -> None:
    require(rank_two_mod_five(), "relations lost rank two modulo five")

    for name, row in (("A", A), ("B", B), ("U", U), ("W", W)):
        require(dot(R, row) == dot(S, row) == 0, f"{name} lost a relation")
        require(gcd(*row) == 1, f"{name} is not primitive")
        require(len(set(row)) == 13, f"{name} is not distinct")

    path_a = radix_path(A)
    path_b = radix_path(B)
    require(
        tuple((carry, owner) for carry, owner, _ in path_a)
        == tuple((carry, owner) for carry, owner, _ in path_b),
        "A/B carry-owner paths differ",
    )
    require(unrestricted_margin(A, 5) == 1, "wrong A margin modulo five")
    require(unrestricted_margin(B, 5) == 0, "wrong B margin modulo five")

    masks, dynamic_states = enumerate_fixed_quotient_fibre()
    safe = masks[0]
    singleton_counts = tuple(masks[1 << i] for i in range(13))
    singleton = sum(singleton_counts)
    plural = sum(count for mask, count in masks.items() if mask.bit_count() >= 2)
    total = sum(masks.values())
    require((safe, singleton, plural, total, len(masks)) == (76787, 346271, 1783379, 2206437, 1836), "fixed-fibre census changed")
    require(all(singleton_counts), "one-deletion sidecar missed a coordinate")

    path_u = radix_path(U)
    path_w = radix_path(W)
    require(path_u == path_w, "U/W full carry-owner-zero-mask paths differ")
    require(exact_margin(U, 25) == (1, (1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19, 21, 22, 23, 24)), "wrong exact U margin")
    require(exact_margin(W, 25) == (2, (3, 7, 18, 22)), "wrong exact W margin")
    require(unrestricted_margin(U, 25) == unrestricted_margin(W, 25) == 5, "wrong unrestricted denominator-25 control")

    clutter_u = danger_clutter(U, 25)
    clutter_w = danger_clutter(W, 25)
    expected_u = (
        (1,), (2,), (3,), (5,), (6,), (8,), (10,), (11,), (4, 12), (7, 9, 13)
    )
    require(clutter_u == expected_u, "wrong U danger clutter")
    require(clutter_w == ((),), "wrong W danger clutter")

    print("THM-2174 PHASE-DANGER CLUTTER REFEREE")
    print("base5_relation_rank=2")
    print("AB_carries=" + str(tuple(carry for carry, _, _ in path_a)))
    print("AB_owners=" + str(tuple(owner for _, owner, _ in path_a)))
    print("AB_m5=(1,0)")
    print(f"fixed_quotient_dp_states={dynamic_states}")
    print(f"fixed_quotient_counts={(safe, singleton, plural, total)}")
    print(f"fixed_quotient_distinct_masks={len(masks)}")
    print(f"fixed_quotient_singletons={singleton_counts}")
    print("UW_zero_masks=" + str(tuple(mask for _, _, mask in path_u)))
    print(f"UW_exact_m25={(exact_margin(U, 25)[0], exact_margin(W, 25)[0])}")
    print(f"UW_unrestricted_m25={(unrestricted_margin(U, 25), unrestricted_margin(W, 25))}")
    print(f"U_clutter25={clutter_u}")
    print(f"W_clutter25={clutter_w}")
    print("first_metric_depths=((2,4),(3,3),(5,2),(7,2),(11,2),(13,2))")
    print("all assertions passed")


if __name__ == "__main__":
    main()
