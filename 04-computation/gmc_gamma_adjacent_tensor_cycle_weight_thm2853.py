#!/usr/bin/env python3
"""Exact companion for THM-2853.

Checks the cycle-weighted inclusion-exclusion polynomial, its explicit lower
family and equality boundary, literal permutation models, positive rational
Gamma shapes, and the first-order zero boundary.  Truth gates remain active
under ``python -O``.
"""

from fractions import Fraction
from itertools import permutations, product
from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return tuple(poly)


def add(left, right):
    out = [0] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        out[index] += coefficient
    for index, coefficient in enumerate(right):
        out[index] += coefficient
    return trim(out)


def scale(poly, coefficient):
    return trim(tuple(coefficient * entry for entry in poly))


def multiply(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)


def linear(n):
    """Coefficient tuple of alpha+n."""
    return (n, 1)


def rising(length):
    out = (1,)
    for n in range(length):
        out = multiply(out, linear(n))
    return out


def evaluate(poly, alpha):
    out = Fraction(0)
    power = Fraction(1)
    for coefficient in poly:
        out += coefficient * power
        power *= alpha
    return out


def tensor_inclusion_exclusion(indices):
    k = len(indices)
    total = sum(n + 1 for n in indices)
    out = (0,)
    for mask in range(1 << k):
        term = rising(total - mask.bit_count())
        for i, n in enumerate(indices):
            if (mask >> i) & 1:
                term = multiply(term, linear(n))
        out = add(
            out,
            scale(term, -1 if mask.bit_count() & 1 else 1),
        )
    return trim(out)


def cycle_count(permutation):
    seen = [False] * len(permutation)
    cycles = 0
    for start in range(len(permutation)):
        if seen[start]:
            continue
        cycles += 1
        vertex = start
        while not seen[vertex]:
            seen[vertex] = True
            vertex = permutation[vertex]
    return cycles


def avoiding_cycle_polynomial(indices):
    k = len(indices)
    private_sets = []
    cursor = k
    for n in indices:
        private_sets.append(set(range(cursor, cursor + n)))
        cursor += n
    total = cursor
    out = [0] * (total + 1)
    for permutation in permutations(range(total)):
        predecessor = [0] * total
        for source, target in enumerate(permutation):
            predecessor[target] = source
        good = True
        for i in range(k):
            singleton = permutation[i] == i
            if singleton or predecessor[i] in private_sets[i]:
                good = False
                break
        if good:
            out[cycle_count(permutation)] += 1
    return trim(out)


def lower_polynomial(indices):
    k = len(indices)
    ordinary = sum(indices)
    return scale(
        multiply((0, 1), rising(ordinary)),
        factorial(k - 1),
    )


def gamma_direct(alpha, indices):
    k = len(indices)
    total = sum(n + 1 for n in indices)

    def rising_value(length):
        out = Fraction(1)
        for h in range(length):
            out *= alpha + h
        return out

    answer = Fraction(0)
    for mask in range(1 << k):
        term = rising_value(total - mask.bit_count())
        for i, n in enumerate(indices):
            if (mask >> i) & 1:
                term *= alpha + n
        answer += (-1 if mask.bit_count() & 1 else 1) * term
    return answer


def main():
    permutation_profiles = 0
    coefficient_profiles = 0
    equality_profiles = []

    for k in range(2, 6):
        for indices in product(range(3), repeat=k):
            total = sum(n + 1 for n in indices)
            poly = tensor_inclusion_exclusion(indices)
            lower = lower_polynomial(indices)
            require(
                all(coefficient >= 0 for coefficient in poly),
                f"negative cycle coefficient at {indices}",
            )
            length = max(len(poly), len(lower))
            padded_poly = poly + (0,) * (length - len(poly))
            padded_lower = lower + (0,) * (length - len(lower))
            require(
                all(a >= b for a, b in zip(padded_poly, padded_lower)),
                f"lower-family failure at {indices}",
            )
            expected_equality = k in (2, 3) and all(
                n == 0 for n in indices
            )
            require(
                (padded_poly == padded_lower) == expected_equality,
                f"equality classification failure at {indices}",
            )
            if expected_equality:
                equality_profiles.append((k, indices))
            coefficient_profiles += 1
            if total <= 8:
                require(
                    poly == avoiding_cycle_polynomial(indices),
                    f"literal permutation mismatch at {indices}",
                )
                permutation_profiles += 1

    alphas = (
        Fraction(1, 10),
        Fraction(1, 2),
        Fraction(1),
        Fraction(2),
        Fraction(7, 2),
        Fraction(10),
    )
    numeric_profiles = 0
    for alpha in alphas:
        for k in range(2, 7):
            for indices in product(range(3), repeat=k):
                poly = tensor_inclusion_exclusion(indices)
                raw = gamma_direct(alpha, indices)
                require(
                    raw == evaluate(poly, alpha) and raw > 0,
                    f"Gamma positivity mismatch at {alpha}, {indices}",
                )
                denominator = Fraction(1)
                for n in indices:
                    for h in range(n + 1):
                        denominator *= alpha + h
                require(
                    raw / denominator > 0,
                    f"normalized Gamma positivity failure at {alpha}, {indices}",
                )
                numeric_profiles += 1

    first_order = 0
    for alpha in alphas:
        for n in range(15):
            require(
                gamma_direct(alpha, (n,)) == 0,
                f"first-order boundary failure at {alpha}, {n}",
            )
            first_order += 1

    print("GAMMA ADJACENT TENSOR CYCLE-WEIGHT AUDIT")
    print(f"coefficient_profiles={coefficient_profiles}")
    print(f"literal_permutation_profiles={permutation_profiles}")
    print(f"numeric_positive_profiles={numeric_profiles}")
    print(f"first_order_zero_controls={first_order}")
    print(f"sharp_lower_equalities={tuple(equality_profiles)}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
