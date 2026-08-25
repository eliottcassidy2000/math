#!/usr/bin/env python3
"""Exact primary audit for THM-4080.

The load-bearing path is Smith elimination over the DVR Z_(p).  It verifies
the closed single-scale two-jet partition, translation/unit-coordinate
invariance, the consecutive-node quadratic-range corollary, and the sharp
cluster-size boundary.
"""

from fractions import Fraction
from hashlib import sha256
import json
from math import comb
import sys


sys.stdout.reconfigure(newline="\n")

EXPECTED_SEMANTIC_SHA256 = "b7cf25da07f12edd250892e42d477dcb929a8f6a081bc46671e2b3ee765e1447"
GATES = 0


def require(condition, label):
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(label)


def vp_integer(value, prime):
    require(value != 0, "valuation needs a nonzero integer")
    value = abs(value)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def vp_fraction(value, prime):
    value = Fraction(value)
    require(value != 0, "valuation needs a nonzero fraction")
    return vp_integer(value.numerator, prime) - vp_integer(value.denominator, prime)


def twojet_matrix(nodes):
    size = 2 * len(nodes)
    return [
        [comb(degree, order) * node ** (degree - order)
         if degree >= order else 0
         for degree in range(size)]
        for node in nodes
        for order in (0, 1)
    ]


def p_smith_exponents(nodes, prime):
    work = [[Fraction(entry) for entry in row] for row in twojet_matrix(nodes)]
    size = len(work)
    exponents = []
    for pivot_index in range(size):
        row, column = min(
            ((row, column)
             for row in range(pivot_index, size)
             for column in range(pivot_index, size)
             if work[row][column]),
            key=lambda pair: vp_fraction(work[pair[0]][pair[1]], prime),
        )
        work[pivot_index], work[row] = work[row], work[pivot_index]
        for current_row in work:
            current_row[pivot_index], current_row[column] = (
                current_row[column], current_row[pivot_index]
            )
        pivot = work[pivot_index][pivot_index]
        exponents.append(vp_fraction(pivot, prime))
        for row in range(pivot_index + 1, size):
            if work[row][pivot_index] == 0:
                continue
            multiplier = work[row][pivot_index] / pivot
            require(vp_fraction(multiplier, prime) >= 0,
                    "DVR-integral row multiplier")
            work[row] = [left - multiplier * right
                         for left, right in zip(work[row], work[pivot_index])]
            require(work[row][pivot_index] == 0, "row clearing")
        for column in range(pivot_index + 1, size):
            if work[pivot_index][column] == 0:
                continue
            multiplier = work[pivot_index][column] / pivot
            require(vp_fraction(multiplier, prime) >= 0,
                    "DVR-integral column multiplier")
            for row in range(pivot_index, size):
                work[row][column] -= multiplier * work[row][pivot_index]
            require(work[pivot_index][column] == 0, "column clearing")
    require(exponents == sorted(exponents), "Smith exponent order")
    expected_total = 4 * sum(
        vp_integer(nodes[right] - nodes[left], prime)
        for left in range(len(nodes))
        for right in range(left + 1, len(nodes))
    )
    require(sum(exponents) == expected_total, "confluent determinant valuation")
    return tuple(exponents)


def single_scale_profile(size, scale_exponent):
    require(size >= 1 and scale_exponent >= 1, "profile domain")
    return tuple(sorted(
        (0, 0)
        + tuple(scale_exponent * value for value in range(1, size))
        + tuple(scale_exponent * value for value in range(size + 1, 2 * size))
    ))


def consecutive_cluster_sizes(node_count, prime):
    quotient, remainder = divmod(node_count, prime)
    return tuple(
        quotient + (1 if residue < remainder else 0)
        for residue in range(prime)
        if quotient + (1 if residue < remainder else 0) > 0
    )


def consecutive_profile(node_count, prime):
    answer = []
    for size in consecutive_cluster_sizes(node_count, prime):
        require(size < prime, "single-scale corollary range")
        answer.extend(single_scale_profile(size, 1))
    return tuple(sorted(answer))


def main():
    require(sys.argv[1:] == [], "no command-line arguments")
    primes = (3, 5, 7, 11)
    single_scale_cases = []
    for prime in primes:
        for scale_exponent in (1, 2, 3):
            scale = prime ** scale_exponent
            for size in range(1, prime):
                unit_sets = (
                    tuple(range(size)),
                    tuple((2 * index + 1) % prime for index in range(size)),
                )
                for variant, units in enumerate(unit_sets):
                    require(len(set(value % prime for value in units)) == size,
                            "unit-coordinate distinctness")
                    translate = 0 if variant == 0 else -3 * prime + 2
                    nodes = tuple(translate + scale * value for value in units)
                    actual = p_smith_exponents(nodes, prime)
                    expected = single_scale_profile(size, scale_exponent)
                    require(actual == expected,
                            ("single-scale profile", prime, scale_exponent,
                             size, variant, actual, expected))
                    single_scale_cases.append(
                        (prime, scale_exponent, size, variant, actual)
                    )

    direct_consecutive_cases = []
    for prime, node_counts in (
        (3, tuple(range(1, 7))),
        (5, tuple(range(1, 21))),
        (7, (1, 6, 7, 8, 13, 14, 15, 20)),
    ):
        for node_count in node_counts:
            require(node_count <= prime * (prime - 1),
                    "quadratic-range node bound")
            actual = p_smith_exponents(tuple(range(node_count)), prime)
            expected = consecutive_profile(node_count, prime)
            require(actual == expected,
                    ("consecutive corollary", prime, node_count, actual, expected))
            direct_consecutive_cases.append((prime, node_count, actual))

    formula_consecutive_cases = []
    for prime in primes:
        for node_count in range(1, prime * (prime - 1) + 1):
            profile = consecutive_profile(node_count, prime)
            determinant_valuation = 4 * sum(
                vp_integer(distance, prime)
                for distance in range(1, node_count)
                for _ in range(node_count - distance)
            )
            require(len(profile) == 2 * node_count,
                    "consecutive profile length")
            require(sum(profile) == determinant_valuation,
                    "consecutive profile determinant")
            formula_consecutive_cases.append(
                (prime, node_count, consecutive_cluster_sizes(node_count, prime), profile)
            )

    hostile_profiles = []
    for prime in (2, 3, 5):
        actual = p_smith_exponents(tuple(prime * index for index in range(prime)), prime)
        false_extension = single_scale_profile(prime, 1)
        require(actual != false_extension, "size-equals-prime hostile must differ")
        hostile_profiles.append((prime, actual, false_extension))

    semantic = {
        "single_scale": [
            [prime, exponent, size, variant, list(profile)]
            for prime, exponent, size, variant, profile in single_scale_cases
        ],
        "direct_consecutive": [
            [prime, count, list(profile)]
            for prime, count, profile in direct_consecutive_cases
        ],
        "formula_consecutive": [
            [prime, count, list(sizes), list(profile)]
            for prime, count, sizes, profile in formula_consecutive_cases
        ],
        "hostiles": [
            [prime, list(actual), list(false_extension)]
            for prime, actual, false_extension in hostile_profiles
        ],
    }
    semantic_digest = sha256(json.dumps(
        semantic, sort_keys=True, separators=(",", ":")
    ).encode("ascii")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    print("CONFLUENT_TWOJET_SINGLE_SCALE_SMITH_THM4080")
    print("status=PROVED_FORMULA;VERIFIED_EXACT")
    print("single_scale_profile=e*(0,0,1,...,s-1,s+1,...,2s-1);scope=1<=s<p")
    print("mechanism=minor_weight_minus_derivative_rows;unit_confluent_matroid_witness")
    print("corollary=consecutive_n<=p*(p-1);partition=union_over_mod_p_cluster_profiles")
    print(f"single_scale_DVR_cases={len(single_scale_cases)}")
    print(f"direct_consecutive_DVR_cases={len(direct_consecutive_cases)}")
    print(f"formula_consecutive_cases={len(formula_consecutive_cases)}")
    for prime, actual, false_extension in hostile_profiles:
        print(f"hostile_s_eq_p(p={prime})=actual:{actual};false_extension:{false_extension}")
    print(f"semantic_sha256={semantic_digest}")
    print(f"gates={GATES}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
