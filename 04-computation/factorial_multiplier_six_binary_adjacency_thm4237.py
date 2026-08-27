#!/usr/bin/env python3
"""Exact arithmetic controls for THM-4237.

The theorem is a symbolic corollary of THM-3474.  This companion checks the
new multiplier-six bit lemma on a declared exhaustive universe, checks the
prime-power specialization, and freezes both sides of the sharp boundary.
It never constructs the large factorial-moment coefficients.
"""

from hashlib import sha256


ODD_H_BOUND = 1 << 16
PRIME_BOUND = 500
EXPONENT_BOUND = 25


def primes_below(bound: int):
    primes = []
    for candidate in range(2, bound):
        prime = True
        divisor = 2
        while divisor * divisor <= candidate:
            if candidate % divisor == 0:
                prime = False
                break
            divisor += 1
        if prime:
            primes.append(candidate)
    return primes


def is_submask(left: int, right: int) -> bool:
    return left & ~right == 0


def has_adjacent_ones(integer: int) -> bool:
    return integer & (integer >> 1) != 0


def reset_limit(prime: int) -> int:
    return min(5, (prime - 1) // 2)


def surviving_multipliers(prime: int, exponent: int):
    height = prime**exponent
    return tuple(
        multiplier
        for multiplier in range(1, reset_limit(prime) + 1)
        if is_submask(multiplier * height, 6 * height)
    )


def direct_survivors(height: int, limit: int):
    return tuple(
        multiplier
        for multiplier in range(1, limit + 1)
        if is_submask(multiplier * height, 6 * height)
    )


def exhaustive_odd_height_check():
    cells = 0
    passed = True
    histograms = {}
    for limit in (3, 5):
        histogram = {}
        for height in range(1, ODD_H_BOUND, 2):
            survivors = direct_survivors(height, limit)
            expected = () if has_adjacent_ones(height) else (
                (2,) if limit == 3 else (2, 4)
            )
            passed &= survivors == expected
            histogram[survivors] = histogram.get(survivors, 0) + 1
            cells += 1
        histograms[limit] = tuple(sorted(histogram.items()))
    return cells, passed, histograms


def prime_power_check():
    cells = 0
    classification = True
    congruence_family = True
    rows = []
    for prime in primes_below(PRIME_BOUND):
        if prime < 7:
            continue
        for exponent in range(1, EXPONENT_BOUND + 1):
            height = prime**exponent
            survivors = surviving_multipliers(prime, exponent)
            classification &= (not survivors) == has_adjacent_ones(height)
            if prime % 4 == 3 and exponent % 2 == 1:
                congruence_family &= height % 4 == 3 and not survivors
            rows.append(
                (
                    prime,
                    exponent,
                    height & 3,
                    has_adjacent_ones(height),
                    survivors,
                )
            )
            cells += 1
    digest = sha256(repr(rows).encode()).hexdigest()
    return cells, classification, congruence_family, digest


def main():
    cells, passed, histograms = exhaustive_odd_height_check()
    print(
        f"odd_height_universe=1<=H<{ODD_H_BOUND},H_odd "
        f"limits=3,5 cells={cells} exact_adjacency_classification={passed}"
    )
    print(f"survivor_histograms={histograms}")

    cells, classification, congruence, digest = prime_power_check()
    print(
        f"prime_power_universe=7<=p<{PRIME_BOUND},p_prime,"
        f"1<=k<={EXPONENT_BOUND} cells={cells} "
        f"classification={classification} congruence_family={congruence}"
    )
    print(f"prime_power_semantic_sha256={digest}")

    controls = ((13, 1), (7, 2), (23, 2), (17, 1))
    rows = []
    for prime, exponent in controls:
        height = prime**exponent
        rows.append(
            (
                prime,
                exponent,
                height,
                bin(height),
                has_adjacent_ones(height),
                surviving_multipliers(prime, exponent),
            )
        )
    print(f"boundary_controls={rows}")
    print(
        "sharp_hostile=p:23,k:2,H:529,binary:0b1000010001,"
        "surviving_t:(2,4),candidate_degrees:(1058,2116),"
        "meaning=compiler_not_actual_factor"
    )


if __name__ == "__main__":
    main()
