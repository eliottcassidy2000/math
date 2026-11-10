#!/usr/bin/env python3
"""Exact low bits of normalized odd A000568, never constructing the full count."""
from collections import Counter
from math import factorial, gcd, prod
from time import perf_counter
import argparse
import json


def v2(n):
    if n <= 0:
        raise ValueError('v2 requires a positive integer')
    return (n & -n).bit_length() - 1


def cutoff_bound(r):
    if r < 1 or r % 2 != 1:
        raise ValueError('odd positive cycle count required')
    return (r - 1) ** 2 // 2 - (r - r.bit_count())


def admissible_counts(bits):
    if bits < 1:
        raise ValueError('positive precision required')
    result = []
    r = 1
    while cutoff_bound(r) < bits:
        result.append(r)
        r += 2
    return result


def partitions_of_length(n, r, minimum=1):
    """Nondecreasing odd partitions, exactly r parts, no hidden filters."""
    if r == 1:
        if n >= minimum and n % 2:
            yield (n,)
        return
    for first in range(minimum, n // r + 1, 2):
        for tail in partitions_of_length(n - first, r - 1, first):
            yield (first,) + tail


def normalized_low_bits(n, bits, statistics=False):
    if n < 1 or n % 2 != 1:
        raise ValueError('positive odd n required')
    modulus = 1 << bits
    total = 0
    visited = retained = 0
    counts = [r for r in admissible_counts(bits) if r <= n]
    for r in counts:
        for parts in partitions_of_length(n, r):
            visited += 1
            cross = sum(gcd(a, b) for i, a in enumerate(parts) for b in parts[:i])
            multiplicities = Counter(parts)
            denominator = prod(part ** m * factorial(m) for part, m in multiplicities.items())
            d = v2(denominator)
            delta = cross - (r - 1) // 2 - d
            if delta < 0:
                raise ArithmeticError('violated nonnegative normalized valuation')
            if delta < bits:
                retained += 1
                odd_denominator = denominator >> d
                total += (1 << delta) * pow(odd_denominator, -1, modulus)
    total %= modulus
    if statistics:
        return total, {'visited_partitions': visited, 'retained_partitions': retained,
                       'cycle_counts': counts, 'precision_bits': bits}
    return total


def coprime_triple_count(n):
    return sum(1 for a, b, c in partitions_of_length(n, 3)
               if a < b < c and gcd(a, b) == gcd(a, c) == gcd(b, c) == 1)


def mod8_formula(n):
    if n == 1 or n == 3:
        return 1
    if n < 1 or n % 2 != 1:
        raise ValueError('positive odd n required')
    return (pow(n, -1, 8) + 2 * pow(n - 2, -1, 8) + 4 * coprime_triple_count(n)) % 8


def davis_integer(n):
    """Independent integer baseline: all cycle types, n! denominator."""
    numerator = 0
    fact = factorial(n)
    for r in range(n % 2 or 2, n + 1, 2):
        for parts in partitions_of_length(n, r):
            exponent = sum((p - 1) // 2 for p in parts)
            exponent += sum(gcd(a, b) for i, a in enumerate(parts) for b in parts[:i])
            denominator = prod(p ** m * factorial(m) for p, m in Counter(parts).items())
            if fact % denominator:
                raise ArithmeticError('noninteger conjugacy class size')
            numerator += (fact // denominator) << exponent
    answer, remainder = divmod(numerator, fact)
    if remainder:
        raise ArithmeticError('Burnside average is not integral')
    return answer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--through', type=int, default=51)
    parser.add_argument('--large-n', type=int, default=1001)
    args = parser.parse_args()
    checks = 0
    examples = []
    t0 = perf_counter()
    for n in range(1, args.through + 1, 2):
        a = davis_integer(n)
        u, rem = divmod(a, 1 << ((n - 1) // 2))
        if rem or not u % 2:
            raise ArithmeticError(('valuation', n))
        for b in range(1, 13):
            actual = normalized_low_bits(n, b)
            if actual != u % (1 << b):
                raise ArithmeticError(('precision mismatch', n, b, actual, u))
            checks += 1
        if n >= 3 and u % 4 != -n % 4:
            raise ArithmeticError(('mod4', n))
        if u % 8 != mod8_formula(n):
            raise ArithmeticError(('mod8', n))
        examples.append({'n': n, 'u_mod_8': u % 8, 'v2_a': v2(a)})
    duration = perf_counter() - t0
    large = []
    for n in (101, 301, args.large_n):
        t0 = perf_counter()
        result, stats = normalized_low_bits(n, 3, True)
        seconds = perf_counter() - t0
        if result != mod8_formula(n):
            raise ArithmeticError(('large independent formula', n))
        large.append({'n': n, 'u_mod_8': result, 'seconds': seconds, **stats})
    negative_checks = 0
    for bad_n, bad_b in [(0, 3), (2, 3), (-1, 3), (3, 0)]:
        try:
            normalized_low_bits(bad_n, bad_b)
        except ValueError:
            negative_checks += 1
        else:
            raise ArithmeticError(('missing domain rejection', bad_n, bad_b))
    n3_naive = (pow(3, -1, 8) + 2 * pow(1, -1, 8)) % 8
    if n3_naive == normalized_low_bits(3, 3):
        raise ArithmeticError('missing n3 hostile')
    huge_n = 10**100 + 1
    output = {'scope': 'odd n; U_n=A000568(n)/2^((n-1)/2); exact residues, not full counts',
              'full_count_crosschecks_through': args.through, 'precision_comparisons': checks,
              'max_checked_precision_bits': 12, 'full_check_seconds': duration,
              'cutoff_examples': {str(b): admissible_counts(b) for b in (1, 2, 3, 5, 6, 14, 15, 32, 64)},
              'examples': examples, 'large_exact_residues': large,
              'huge_formula_example': {'n': str(huge_n), 'u_mod_4': (-huge_n) % 4},
              'hostile_controls': {'n1_formula_fails': normalized_low_bits(1, 2) != -1 % 4,
                                   'n3_naive_mod8': n3_naive,
                                   'n3_actual_mod8': normalized_low_bits(3, 3),
                                   'domain_rejections': negative_checks}}
    print(json.dumps(output, indent=2))

if __name__ == '__main__':
    main()
