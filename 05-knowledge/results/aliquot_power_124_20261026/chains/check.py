#!/usr/bin/env python3
"""Exact finite controls for RESULT.md; standard library, no assert gates."""
import hashlib
import json
from collections import Counter
from functools import reduce
from math import gcd, isqrt
from pathlib import Path


GATES = Counter()


def require(value, group, detail=None):
    GATES[group] += 1
    if not value:
        raise RuntimeError((group, detail))


def factor(n):
    if n < 1:
        raise ValueError("factor domain is positive integers")
    answer = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            answer[p] = answer.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        answer[n] = answer.get(n, 0) + 1
    return answer


def two_address(n):
    if n < 1:
        raise ValueError("address domain is positive integers")
    a = 0
    while n % 2 == 0:
        a += 1
        n //= 2
    return a, n


def iroot(n, k):
    if n < 0 or k < 1:
        raise ValueError("integer root requires n>=0, k>=1")
    lo, hi = 0, n + 1
    while hi - lo > 1:
        middle = (lo + hi) // 2
        if middle ** k <= n:
            lo = middle
        else:
            hi = middle
    return lo


def odd_count(t):
    return (t + 1) // 2


def delayed_even_formula(bound, k):
    result = 0
    for a in range(1, bound.bit_length()):
        ceiling = bound // (2 ** a)
        if ceiling:
            result += odd_count(ceiling)
            result -= odd_count(min(ceiling, iroot(2 ** a, k - 1)))
    return result


def sigma_factored(n):
    result = 1
    for p, e in factor(n).items():
        result *= (p ** (e + 1) - 1) // (p - 1)
    return result


def aliquot(n):
    if n < 1:
        raise ValueError("zero is terminal, not a divisor-sum input")
    return sigma_factored(n) - n


def rad(n):
    result = 1
    for p in factor(n):
        result *= p
    return result


def primitive_edge(n):
    if n < 2:
        raise ValueError("positive abc edge requires n>=2")
    successor = aliquot(n)
    common = gcd(n, successor)
    return (n // common, successor // common, (n + successor) // common), common


def main():
    bound = 20000
    double_count = [0] * (bound + 1)
    double_first = [None] * (bound + 1)
    power_count = [0] * (bound + 1)
    power_first = [None] * (bound + 1)
    combined_count = [0] * (bound + 1)
    combined_first = [None] * (bound + 1)
    proper_power_bases = [set() for _ in range(bound + 1)]

    # Independent producer: literal multiplication, no factorization oracle.
    for r in range(1, bound + 1):
        value = r
        while value <= bound:
            double_count[value] += 1
            if double_first[value] is None:
                double_first[value] = r
            value *= 2
    for r in range(2, bound + 1):
        value, k = r, 1
        while value <= bound:
            power_count[value] += 1
            if power_first[value] is None:
                power_first[value] = r
            if k >= 2:
                proper_power_bases[value].add(r)
            doubled = value
            while doubled <= bound:
                combined_count[doubled] += 1
                if combined_first[doubled] is None:
                    combined_first[doubled] = r
                doubled *= 2
            k += 1
            value *= r

    for n in range(1, bound + 1):
        a, u = two_address(n)
        require(double_count[n] == a + 1 and double_first[n] == u,
                "doubling_fibres", n)
        if n == 1:
            continue
        exponents = factor(n)
        g_all = reduce(gcd, exponents.values())
        allowed = [k for k in range(1, g_all + 1) if g_all % k == 0]
        require(power_count[n] == len(allowed)
                and power_first[n] == iroot(n, g_all), "power_fibres", n)
        if u == 1:
            predicted = sum(a // k for k in range(1, a + 1))
            first = 2
        else:
            g_odd = reduce(gcd, factor(u).values())
            predicted = sum(a // k + 1 for k in range(1, g_odd + 1)
                            if g_odd % k == 0)
            first = iroot(u, g_odd)
        require(combined_count[n] == predicted and combined_first[n] == first,
                "combined_fibres", n)

    require(next(n for n in range(2, bound + 1)
                 if len(proper_power_bases[n]) > 1) == 16,
            "first_occurrence_controls")
    require(double_first[36] == 9 and combined_first[36] == 3,
            "first_occurrence_controls")
    require(double_count[1] == 1 and power_count[1] == combined_count[1] == 0,
            "first_occurrence_controls")

    delayed_receipts = []
    for k in range(2, 10):
        delayed_prefix = [0] * (bound + 1)
        delayed_even_prefix = [0] * (bound + 1)
        for n in range(1, bound + 1):
            a, u = two_address(n)
            # Independent: factor n^k by literal repeated division by two.
            direct = two_address(n ** k)[1] > n
            require(direct == (u ** (k - 1) > 2 ** a), "delay_classifier", (n, k))
            delayed_prefix[n] = delayed_prefix[n - 1] + direct
            delayed_even_prefix[n] = delayed_even_prefix[n - 1] + (direct and n % 2 == 0)
        for ceiling in list(range(1, 257)) + [1000, 4096, 10000, bound]:
            predicted = delayed_even_formula(ceiling, k)
            require(predicted == delayed_even_prefix[ceiling],
                    "delay_census", (ceiling, k))
        require(next(n for n in range(1, bound + 1)
                     if delayed_prefix[n] > delayed_prefix[n - 1]) == 3,
                "first_occurrence_controls", k)
        require(next(n for n in range(2, bound + 1, 2)
                     if delayed_even_prefix[n] > delayed_even_prefix[n - 1]) == 6,
                "first_occurrence_controls", k)
        delayed_receipts.append({"k": k, "bound": bound,
                                 "delayed_even": delayed_even_prefix[bound],
                                 "nondelayed_all": bound - delayed_prefix[bound]})

    for m in range(1, 11):
        ceiling = 4 ** m
        counted = ceiling // 2 - delayed_even_formula(ceiling, 2)
        require(counted == 3 * 2 ** (m - 1) - 1,
                "exact_square_boundary", m)
    for n in range(2, 4097):
        a, u = two_address(n)
        if u == 1:
            require(all(u ** (k - 1) <= 2 ** a for k in range(2, 25)),
                    "power_of_two_hostile", n)
        else:
            k, magnitude = 2, u
            while magnitude <= 2 ** a:
                magnitude *= u
                k += 1
            require(two_address(n ** k)[1] > n
                    and (k == 2 or two_address(n ** (k - 1))[1] <= n),
                    "first_delayed_exponent", (n, k))
    for n in range(1, 65):
        for j in range(0, 6):
            for k in range(1, 7):
                require((2 ** j * n) ** k == 2 ** (j * k) * n ** k,
                        "semidirect_operation", (n, j, k))

    # Independent divisor-sum path: add each proper divisor to its multiples.
    sieve = [0] * (bound + 1)
    for d in range(1, bound // 2 + 1):
        for n in range(2 * d, bound + 1, d):
            sieve[n] += d
    for n in range(1, bound + 1):
        require(aliquot(n) == sieve[n], "divisor_sum_independent", n)
        if n >= 2:
            triple, common = primitive_edge(n)
            a, b, c = triple
            require(a + b == c and min(triple) > 0
                    and gcd(a, b) == gcd(a, c) == gcd(b, c) == 1
                    and common * a == n and common * b == sieve[n]
                    and ((a == b) == (sieve[n] == n)),
                    "primitive_abc_edges", n)
    perfects = [n for n in range(2, bound + 1) if sieve[n] == n]
    require(perfects == [6, 28, 496, 8128], "perfect_finite_controls")
    require(all(primitive_edge(n)[0] == (1, 1, 2) for n in perfects),
            "primitive_scale_hostile")

    diagrams = []
    for a in range(1, 13):
        mersenne = 2 ** (a + 1) - 1
        if factor(mersenne) != {mersenne: 1}:
            continue
        n = 2 ** a * mersenne
        chain = [mersenne ** 2]
        for _ in range(4):
            chain.append(aliquot(chain[-1]))
        require(chain == [mersenne ** 2, 2 ** (a + 1), mersenne, 1, 0]
                and aliquot(n) == n and mersenne ** 2 > n
                and two_address(n ** 2) == (2 * a, mersenne ** 2),
                "perfect_square_diagrams", (a, chain))
        diagrams.append({"a": a, "N": n, "chain": chain})
    require(aliquot(225) == 178 != 16 and aliquot(120) != 120,
            "composite_mersenne_hostile")
    require(two_address(aliquot(6))[1] == 3 != aliquot(two_address(6)[1]),
            "no_aliquot_quotient_conjugacy")
    require(two_address(10 ** 2)[1] > 10 and aliquot(10) != 10,
            "delay_not_perfection")
    require((2 * 3) ** 2 == 36 != 2 * (3 ** 2), "noncommutation_hostile")

    for p in range(2, 32):
        if factor(p) != {p: 1}:
            continue
        for k in range(1, 7):
            n = p ** k
            successor = aliquot(n)
            b = (p - 1) * successor
            require(successor == sum(p ** j for j in range(k))
                    and 1 + b == n and gcd(b, n) == 1,
                    "prime_power_abc_identity", (p, k))
            require(rad(b * n) == p * rad(b),
                    "prime_power_abc_radical", (p, k))
    require(primitive_edge(9)[0] == (9, 4, 13)
            and 1 + (3 - 1) * aliquot(9) == 9, "distinct_abc_maps")
    for bad in (0, 1):
        try:
            primitive_edge(bad)
        except ValueError:
            require(True, "invalid_abc_domain", bad)
        else:
            require(False, "invalid_abc_domain", bad)

    print(json.dumps({
        "status": "FINITE-EXACT; all-parameter conclusions use RESULT.md proofs",
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "universe": {"table_value_bound": bound, "delay_exponents": [2, 9],
                     "census_ceilings": "1..256,1000,4096,10000,20000",
                     "square_boundary_m": [1, 10], "first_exponent_bases": [2, 4096],
                     "semidirect": "n1..64,j0..5,k1..6",
                     "abc_edges": [2, bound], "mersenne_a": [1, 12],
                     "prime_power_abc": "p prime<=31,k1..6"},
        "gates": dict(sorted(GATES.items())), "total_gates": sum(GATES.values()),
        "first_rows_36": {"doubling": double_first[36], "combined": combined_first[36]},
        "delayed_census": delayed_receipts, "perfect_diagrams": diagrams,
        "conditional_boundary": "No numerical ABC inequality or effective constant is tested."
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
