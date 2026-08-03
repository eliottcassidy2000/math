#!/usr/bin/env python3
"""Exact checks for the all-n unique-protector and sparse-fragility lemmas.

The mathematical proof is recorded in the matching reflection.  This script
checks the explicit rational cell, the congruence-cardinality formula, the
actual safe-shift consequence on a complete small bank, and named boundary
controls.  It also audits the finite-field moment obstruction at the first
unsupported Hamming layer for odd primes.  It makes no LRC(14) claim.
"""

from __future__ import annotations

import ast
import hashlib
import itertools
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def bin_value(n: int, i: int, alpha: Fraction) -> int:
    x = (i * alpha) % 1
    nx = n * x
    return nx.numerator // nx.denominator


def explicit_alpha(n: int, j: int) -> Fraction:
    require(n >= 3 and 1 <= j < n, "invalid explicit-cell parameters")
    if 2 * j <= n - 1:
        t = Fraction(1, 2) + Fraction(1, 4 * n * n)
    else:
        t = Fraction(1, 2 * n)
    return (1 - t / n) / j


def explicit_pattern(n: int, j: int) -> tuple[int, ...]:
    alpha = explicit_alpha(n, j)
    pattern = tuple(bin_value(n, i, alpha) for i in range(1, n))
    for i in range(1, n):
        x = (i * alpha) % 1
        require((n * x).denominator != 1, f"explicit alpha is a breakpoint: {(n, j, i)}")
    endpoint_set = tuple(i for i, b in enumerate(pattern, 1) if b in (0, n - 1))
    require(endpoint_set == (j,), f"explicit singleton failed: {(n, j, endpoint_set)}")
    return pattern


def atomic_patterns(n: int) -> tuple[tuple[int, ...], ...]:
    points = sorted(
        {Fraction(a, n * i) for i in range(1, n) for a in range(n * i + 1)}
    )
    patterns = []
    for left, right in zip(points, points[1:]):
        alpha = (left + right) / 2
        patterns.append(tuple(bin_value(n, i, alpha) for i in range(1, n)))
    require(len(patterns) == len(set(patterns)), f"duplicate atomic patterns at n={n}")
    return tuple(patterns)


def shift_cost(n: int, residue: int) -> int:
    residue %= n
    require(residue != 0, "shift cost is only for nonzero deviations")
    divisor = math.gcd(residue, n)
    return 2 if divisor == 1 else divisor


def actual_bad_shifts(n: int, residue: int, bin_index: int) -> tuple[int, ...]:
    return tuple(
        s
        for s in range(n)
        if (s * residue + bin_index) % n in (0, n - 1)
    )


def safe_shift_for_pattern(
    n: int, vector: tuple[int, ...], pattern: tuple[int, ...]
) -> int | None:
    for s in range(n):
        if all(
            (s * vector[i - 1] + pattern[i - 1]) % n not in (0, n - 1)
            for i in range(1, n)
        ):
            return s
    return None


def vector_cost(n: int, vector: tuple[int, ...]) -> int:
    return sum(shift_cost(n, value) for value in vector if value % n)


def certified_class_count(n: int) -> tuple[int, int, tuple[tuple[int, int], ...]]:
    free_coordinates = n - 2
    residue_weights = Counter(shift_cost(n, value) for value in range(1, n))
    counts = [0] * n
    counts[0] = 1
    for _ in range(free_coordinates):
        next_counts = [0] * n
        for total, count in enumerate(counts):
            if count == 0:
                continue
            next_counts[total] += count
            for weight, multiplicity in residue_weights.items():
                if total + weight < n:
                    next_counts[total + weight] += count * multiplicity
        counts = next_counts
    certified_nonzero = sum(counts) - 1
    all_nonzero = n ** free_coordinates - 1
    return certified_nonzero, all_nonzero, tuple(sorted(residue_weights.items()))


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, math.isqrt(value) + 1))


def critical_polynomial_value(prime: int, value: int) -> int:
    layer = (prime + 1) // 2
    return (
        pow(value, layer, prime)
        - pow(value, layer - 1, prime)
        + layer
        - 1
    ) % prime


def main() -> None:
    print("Micro-staircase unique-protector and sparse-fragility audit")
    print("scope=PROVED_elementary_lemma_plus_FINITE-EXACT_controls;not_LRC")

    singleton_rows = []
    for n in range(3, 31):
        patterns = atomic_patterns(n)
        singleton_counts = Counter()
        for pattern in patterns:
            endpoint_set = tuple(
                i for i, b in enumerate(pattern, 1) if b in (0, n - 1)
            )
            if len(endpoint_set) == 1:
                singleton_counts[endpoint_set[0]] += 1
        explicit = tuple(explicit_pattern(n, j) for j in range(1, n))
        require(len(explicit) == n - 1, "explicit singleton bank has wrong size")
        require(all(singleton_counts[j] > 0 for j in range(1, n)), "missing singleton")
        singleton_rows.append(
            (n, len(patterns), min(singleton_counts.values()), max(singleton_counts.values()))
        )
    print("explicit_singleton_bank=n3..30:PASS")
    print("singleton_rows=(n,patterns,min_per_coordinate,max_per_coordinate)")
    print(tuple(singleton_rows))

    congruence_cases = 0
    for n in range(3, 61):
        for residue in range(1, n):
            cap = shift_cost(n, residue)
            for b in range(n):
                actual = len(actual_bad_shifts(n, residue, b))
                require(actual <= cap, f"shift cap failed: {(n, residue, b, actual, cap)}")
                divisor = math.gcd(residue, n)
                if divisor == 1:
                    require(actual == 2, "unit residue must have two bad shifts")
                else:
                    require(actual in (0, divisor), "nonunit residue cardinality mismatch")
                congruence_cases += 1
    print(f"congruence_cardinality_cases={congruence_cases};status=PASS")

    exhaustive_rows = []
    for n in range(3, 9):
        patterns = tuple(explicit_pattern(n, j) for j in range(1, n))
        certified = 0
        for free in itertools.product(range(n), repeat=n - 2):
            vector = (0,) + free
            if not any(vector):
                continue
            if vector_cost(n, vector) >= n:
                continue
            support_j = next(i for i, value in enumerate(vector, 1) if value)
            safe = safe_shift_for_pattern(n, vector, patterns[support_j - 1])
            require(safe is not None, f"certified vector has no safe shift: {(n, vector)}")
            certified += 1
        zero_patterns = atomic_patterns(n)
        require(
            all(
                any(b in (0, n - 1) for b in pattern)
                for pattern in zero_patterns
            ),
            f"zero-ramp positive control failed at n={n}",
        )
        exhaustive_rows.append((n, certified))
    print(f"exhaustive_normalized_consequence_n3..8={tuple(exhaustive_rows)}")

    for n in (13, 14, 15):
        certified, total, weights = certified_class_count(n)
        print(
            f"normalized_n={n};residue_weight_multiplicities={weights};"
            f"certified_nonzero_classes={certified};all_nonzero_classes={total}"
        )

    prime_cases = 0
    maximum_distinct_roots = 0
    maximum_root_multiplicity = 0
    for prime in (value for value in range(7, 98) if is_prime(value)):
        layer = (prime + 1) // 2
        roots = tuple(
            value
            for value in range(1, prime)
            if critical_polynomial_value(prime, value) == 0
        )
        predicted = tuple(
            value
            for value in (layer % prime, (2 - layer) % prime)
            if value != 0
            and pow(value, (prime - 1) // 2, prime)
            == (-1 if value == layer % prime else 1) % prime
        )
        require(roots == tuple(dict.fromkeys(predicted)), f"Euler root split failed at p={prime}")
        multiplicities = []
        for root in roots:
            derivative = (
                layer * pow(root, layer - 1, prime)
                - (layer - 1) * pow(root, layer - 2, prime)
            ) % prime
            multiplicities.append(2 if derivative == 0 else 1)
        require(sum(multiplicities) < layer, f"critical polynomial could split at p={prime}")
        maximum_distinct_roots = max(maximum_distinct_roots, len(roots))
        maximum_root_multiplicity = max(
            maximum_root_multiplicity, sum(multiplicities)
        )
        prime_cases += 1

    moment_multisets = 0
    for residues in itertools.combinations_with_replacement(range(1, 13), 7):
        first_power = sum(residues) % 13
        if all(
            sum(pow(residue, exponent, 13) for residue in residues) % 13
            == pow(first_power, exponent, 13)
            for exponent in range(2, 7)
        ):
            moment_multisets += 1
    require(moment_multisets == 0, "p=13 critical moment multiset survived")
    forced_roots_13 = tuple(
        value for value in range(1, 13) if critical_polynomial_value(13, value) == 0
    )
    require(forced_roots_13 == (7,), "unexpected p=13 forced roots")
    print(
        "prime_critical_layer=p7..97:PASS;"
        f"prime_cases={prime_cases};max_distinct_roots={maximum_distinct_roots};"
        f"max_total_root_multiplicity={maximum_root_multiplicity};"
        f"p13_moment_multisets={moment_multisets};p13_forced_roots={forced_roots_13}"
    )

    hostile_n = 14
    hostile = [0] * (hostile_n - 1)
    hostile[1] = 7
    hostile[2] = 7
    hostile_vector = tuple(hostile)
    require(vector_cost(hostile_n, hostile_vector) == hostile_n, "hostile is not boundary-cost")
    hostile_patterns = atomic_patterns(hostile_n)
    hostile_witness = None
    for cell_index, pattern in enumerate(hostile_patterns):
        safe = safe_shift_for_pattern(hostile_n, hostile_vector, pattern)
        if safe is not None:
            hostile_witness = (cell_index, safe, pattern)
            break
    require(hostile_witness is not None, "boundary hostile unexpectedly blocks")
    print(
        "strictness_hostile=n14_two_halfturns;cost=14;"
        f"safe_witness={hostile_witness};criterion_makes_no_claim"
    )

    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assert_nodes == 0, "assert node present")
    require(float_literals == 0, "float literal present")
    print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")
    print(f"source_sha256={hashlib.sha256(source.encode('utf-8')).hexdigest()}")
    print("ALL_CHECKS_PASSED")


if __name__ == "__main__":
    main()
