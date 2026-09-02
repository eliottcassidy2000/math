#!/usr/bin/env python3
"""Exact audit of one joint clock/capacity hostile in the LRC(14) minority lane.

The row is safe; this script only verifies failure of three sufficient methods:
the complete integer half-turn grid, the p=8..14 unit banks, and MC7 adaptive
divisor capacity.  All arithmetic is integer or fractions.Fraction arithmetic.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt


H = 420
ANCHOR = 2 * H
M = 28 * H
D = tuple(11 + 1680 * k for k in range(7))
C = (525, 945, 1365, 1575)
PAD = 1287
W = D + C + (PAD,)
S = (ANCHOR,) + W


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def circular_numerator(v: int, numerator: int, denominator: int) -> int:
    residue = (v * numerator) % denominator
    return min(residue, denominator - residue)


def trial_divisors(n: int) -> set[int]:
    result = {1, n}
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            result.add(d)
            result.add(n // d)
    return result


def factorization(n: int) -> list[tuple[int, int]]:
    factors = []
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            factors.append((p, exponent))
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors.append((n, 1))
    return factors


def product_divisors(n: int) -> set[int]:
    result = {1}
    for p, exponent in factorization(n):
        result = {d * p**e for d in result for e in range(exponent + 1)}
    return result


def mc7_capacity(a: int) -> Fraction:
    total = Fraction()
    for v in S:
        if v % a:
            q = a // gcd(a, v)
            total += Fraction((q + 6) // 7, q)
    return total


def sign_class(v: int, p: int) -> int | None:
    modulus = 2 * p
    if gcd(v, modulus) != 1:
        return None
    residue = v % modulus
    return min(residue, modulus - residue)


def odd_c_mask(u: int, odd_s: tuple[int, ...]) -> int:
    # For j=7s and speed 7u, danger is 103 <= su (mod 240) <= 137.
    return sum(
        1 << index
        for index, s in enumerate(odd_s)
        if 103 <= (s * u) % 240 <= 137
    )


def main() -> None:
    require(len(W) == len(set(W)) == 12, "W is not twelve distinct speeds")
    require(all(w > 0 and w % 2 == 1 for w in W), "W is not positive odd")
    require(gcd(*S) == 1, "row is not primitive")

    missing_gates = [m for m in range(2, 15) if not any(v % m == 0 for v in S)]
    require(not missing_gates, "THM-366 gate is missing")

    # Symbolic partition checks behind the full grid proof.
    nonseven = tuple(j for j in range(M) if j % 7)
    d_multiplicity = Counter(
        sum(5040 < (j * d) % M < 6720 for d in D) for j in nonseven
    )
    require(d_multiplicity == Counter({1: len(nonseven)}), "D is not a partition")

    odd_s_full = tuple(range(1, 1680, 2))
    c_multiplicity = Counter(
        sum(5040 < ((7 * s) * c) % M < 6720 for c in C)
        for s in odd_s_full
    )
    require(c_multiplicity == Counter({1: len(odd_s_full)}), "C is not a partition")

    # Full direct replay, including anchor-unsafe j.
    grid_safe = []
    grid_odd_blocked = 0
    grid_anchor_blocked = 0
    for j in range(M):
        numerator = M // 2 + j
        distances = [circular_numerator(v, numerator, M) for v in S]
        if 14 * min(distances) >= M:
            grid_safe.append(j)
        if j % 14 == 0:
            require(distances[0] == 0, "anchor did not block a multiple of 14")
            grid_anchor_blocked += 1
        else:
            require(min(distances[1:]) * 14 < M, "odd block is missing")
            grid_odd_blocked += 1
    require(not grid_safe, "half-turn grid has a witness")

    # Unit sign-class theorem and direct bank replay.
    bank_counts = []
    for p in range(8, 15):
        modulus = 2 * p
        units = tuple(a for a in range(1, modulus) if gcd(a, modulus) == 1)
        survivors = []
        for a in units:
            distances = [circular_numerator(v, a, modulus) for v in S]
            if 14 * min(distances) >= modulus:
                survivors.append(a)
        require(not survivors, f"p={p} unit bank has a witness")
        represented = {sign_class(w, p) for w in W} - {None}
        target = {min(a, modulus - a) for a in units}
        if H % p:
            require(represented == target, f"p={p} sign classes are incomplete")
        bank_counts.append((p, len(units), len(represented), H % p == 0))

    # Two independent divisor constructors feed the exact MC7 audit.
    divisors_trial = {d for v in S for d in trial_divisors(v) if d >= 2}
    divisors_product = {d for v in S for d in product_divisors(v) if d >= 2}
    require(divisors_trial == divisors_product, "divisor constructors disagree")
    capacities = {a: mc7_capacity(a) for a in sorted(divisors_trial)}
    minimum = min(capacities.values())
    minimizers = [a for a, value in capacities.items() if value == minimum]
    failures = [a for a, value in capacities.items() if value < 1]
    require(not failures, "MC7 closes the row")
    require(minimum == Fraction(8, 7), "MC7 minimum changed")
    require(minimizers == [7, 21, 35, 105], "MC7 minimizers changed")

    # Exhaust all C candidates modulo 240.  Four is minimal, and the cover is
    # unique up to independently changing the signs of four residue classes.
    odd_s = tuple(range(1, 240, 2))
    odd_u = tuple(range(1, 240, 2))
    masks_by_u = {u: odd_c_mask(u, odd_s) for u in odd_u}
    users_by_mask: dict[int, list[int]] = defaultdict(list)
    for u, mask in masks_by_u.items():
        users_by_mask[mask].append(u)
    masks = tuple(users_by_mask)
    full_mask = (1 << len(odd_s)) - 1
    cover_counts = {}
    cover_examples = {}
    for size in range(1, 5):
        expanded_count = 0
        first_example = None
        for chosen in combinations(masks, size):
            union = 0
            multiplicity = 1
            for mask in chosen:
                union |= mask
                multiplicity *= len(users_by_mask[mask])
            if union == full_mask:
                expanded_count += multiplicity
                if first_example is None:
                    first_example = tuple(tuple(users_by_mask[mask]) for mask in chosen)
        cover_counts[size] = expanded_count
        cover_examples[size] = first_example
    require(cover_counts == {1: 0, 2: 0, 3: 0, 4: 16}, "C cover census changed")
    expected_classes = {
        frozenset((15, 225)),
        frozenset((45, 195)),
        frozenset((75, 165)),
        frozenset((105, 135)),
    }
    actual_classes = {frozenset(group) for group in cover_examples[4]}
    require(actual_classes == expected_classes, "minimal C classes changed")

    # Positive control: this is a method hostile, not an unsafe row.
    witness_distances = [circular_numerator(v, 6, 17) for v in S]
    require(min(witness_distances) == 2, "positive-control clearance changed")
    require(14 * min(witness_distances) > 17, "positive control is not strict")
    tight = [v for v, distance in zip(S, witness_distances) if distance == 2]
    require(tight == [11, 1575], "positive-control tight set changed")

    print("LRC14 MINORITY ANCHOR JOINT CLOCK/CAPACITY HOSTILE -- EXACT AUDIT")
    print(f"h={H} anchor={ANCHOR} modulus={M}")
    print("W=" + ",".join(map(str, W)))
    print(f"primitive={gcd(*S)} distinct_odds={len(set(W))} missing_thm366={missing_gates}")
    print(f"D_partition_universe={len(nonseven)} multiplicity={dict(d_multiplicity)}")
    print(f"C_partition_universe={len(odd_s_full)} multiplicity={dict(c_multiplicity)}")
    print(
        f"grid_total={M} anchor_blocked={grid_anchor_blocked} "
        f"odd_blocked={grid_odd_blocked} safe={len(grid_safe)}"
    )
    for p, unit_count, class_count, anchor_killed in bank_counts:
        print(
            f"unit_bank_p={p} clocks={unit_count} represented_sign_classes={class_count} "
            f"anchor_killed={anchor_killed} survivors=0"
        )
    print(
        f"mc7_represented_divisors={len(capacities)} failures={len(failures)} "
        f"minimum={minimum} minimizers={','.join(map(str, minimizers))}"
    )
    print(
        "C_cover_unique_masks="
        f"{len(masks)} mask_size_hist={dict(sorted(Counter(m.bit_count() for m in masks).items()))}"
    )
    print(
        "C_cover_counts="
        + ",".join(f"size{size}:{cover_counts[size]}" for size in sorted(cover_counts))
    )
    print(
        "C_minimal_classes="
        + ";".join(
            ",".join(map(str, sorted(group))) for group in sorted(expected_classes, key=min)
        )
    )
    print(f"positive_control=6/17 minimum=2/17 tight={','.join(map(str, tight))}")
    print("VERDICT PASS_FINITE_EXACT_METHOD_HOSTILE_NOT_LRC_COUNTEREXAMPLE")


if __name__ == "__main__":
    main()
