#!/usr/bin/env python3
"""Exact finite controls for the two-clock source-address theorem.

Standard library only; python3 -B and python3 -O -B produce identical files.
No floating-point cancellation or integer convergence claim is used.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import json

CHECKS = Counter()


def check(condition, group):
    CHECKS[group] += 1
    if not condition:
        raise RuntimeError(f"{group}: check {CHECKS[group]} failed")


def compositions(cost, length):
    for cuts in combinations(range(1, cost), length - 1):
        endpoints = (0,) + cuts + (cost,)
        yield tuple(b - a for a, b in zip(endpoints, endpoints[1:]))


def source(word, q):
    p, Q, B = 1, 1, 0
    for a in word:
        p, Q, B = q * p, Q * 2 ** a, q * B + Q
    return -B * pow(p, -1, Q) % Q


def parity_replay(n, q, cost):
    bits = []
    for _ in range(cost):
        bits.append(n % 2)
        n = (q * n + 1) // 2 if n % 2 else n // 2
    return bits


def exact_character_sum(terms, u, depth):
    """Represent in Q[z]/(z^(2^(depth-1))+1), z=e(1/2^depth)."""
    Q = 1 << depth
    half = Q // 2
    coefficients = [F(0)] * half
    for address, mass in terms:
        exponent = int(address * Q) * u % Q
        coefficients[exponent % half] += mass if exponent < half else -mass
    return coefficients


def audit():
    total_words = 0
    data = {}
    for q in (3, 5, 7, 9, 11):
        by_cost = {}
        all_characters = set()
        gram = defaultdict(F)
        for cost in range(1, 13):
            Q = 1 << cost
            addresses = {}
            for m in range(1, cost + 1):
                count = 0
                for word in compositions(cost, m):
                    r = source(word, q)
                    check(r % 2 == 1 and r not in addresses, "binary_partition_injective")
                    expected = [b for a in word for b in [1] + [0] * (a - 1)]
                    check(parity_replay(r, q, cost) == expected, "independent_parity_replay")
                    check(F(r, Q) not in all_characters, "global_character_injective")
                    all_characters.add(F(r, Q))
                    addresses[r] = (m, word)
                    count += 1
                check(count == comb(cost - 1, m - 1), "composition_binomial_count")
            check(set(addresses) == set(range(1, Q, 2)), "binary_partition_surjective")
            by_cost[cost] = addresses
            total_words += len(addresses)
            # Haar on odd seeds: only a difference 0 or 1/2 has nonzero mean.
            for r, (m, word) in addresses.items():
                gram[m, m] += F(1, Q * Q)
                if cost > 1:
                    s = r ^ (Q // 2)
                    n, partner = addresses[s]
                    check(abs(m - n) == 1, "half_turn_changes_one_odd_step")
                    gram[m, n] -= F(1, Q * Q)
                    if word[-1] > 1:
                        expected = word[:-1] + (word[-1] - 1, 1)
                        check(partner == expected, "half_turn_splits_final_block")
            # Independent closed finite-cutoff covariance formula.
            edges = {m: sum((F(comb(A - 2, m - 1), 4 ** A)
                             for A in range(m + 1, cost + 1)), F(0))
                     for m in range(1, cost + 1)}
            for m in range(1, cost + 1):
                for n in range(1, cost + 1):
                    target = F(0)
                    if m == n:
                        target = (F(1, 4) if m == 1 else edges[m - 1]) + edges[m]
                    elif abs(m - n) == 1:
                        target = -edges[min(m, n)]
                    check(gram[m, n] == target, "odd_Haar_tridiagonal_Gram")
            for cutoff in range(1, cost + 1):
                variance = sum((gram[m, n] for m in range(1, cutoff + 1)
                                for n in range(1, cutoff + 1)), F(0)) - F(1, 4)
                check(variance == edges[cutoff], "finite_odd_step_boundary_energy")
            # A nonconstant rational coefficient sequence probes the entire energy.
            a = {m: F((-1) ** m * (m + 1), m + 2) for m in range(1, cost + 1)}
            a[cost + 1] = F(0)
            quadratic = sum((a[m] * a[n] * gram[m, n]
                             for m in range(1, cost + 1)
                             for n in range(1, cost + 1)), F(0))
            energy = a[1] ** 2 / 4 + sum((edges[m] * (a[m] - a[m + 1]) ** 2
                                        for m in range(1, cost + 1)), F(0))
            check(quadratic == energy, "finite_weighted_difference_energy")
        # Independent exact roots-of-unity sums; no trigonometric rounding.
        for depth in range(1, 9):
            for t in (F(1, 2), F(2, 3), F(1)):
                terms = [(F(r, 2 ** A), (t / 2) ** A)
                         for A in range(1, depth + 1) for r in by_cost[A]]
                for u in (0, 1, 2, 3, 4, 8, 11, 16, 19, 63, -1, -5, -17):
                    expected = F(0)
                    for A in range(1, depth + 1):
                        if u % 2 ** A == 0:
                            expected += t ** A / 2
                        elif u % 2 ** (A - 1) == 0:
                            expected -= t ** A / 2
                    polynomial = exact_character_sum(terms, u, depth)
                    check(polynomial == [expected] + [F(0)] * (len(polynomial) - 1),
                          "cost_clock_Ramanujan_identity")
        data[str(q)] = {"words_cost_at_most_12": len(all_characters),
                        "cost_3_word_to_source": {str(w): r for r, (_, w) in by_cost[3].items()}}
    # A cross-length half-turn collision is deliberately retained.
    left, right = F(source((3,), 3), 8), F(source((2, 1), 3), 8)
    check((left, right) == (F(5, 8), F(1, 8)), "odd_Haar_cross_length_hostile")
    check((left - right) % 1 == F(1, 2), "odd_Haar_cross_length_hostile")
    return {"scope": "FINITE-EXACT controls; all-depth proofs in the note; fixed-seed H1 OPEN",
            "multipliers": data, "total_word_cases": total_words,
            "proved_full_Haar_norm_squared": "3^(-m)",
            "proved_odd_Haar_adjacent_inner_product": "-1/(4*3^m)",
            "proved_odd_Haar_partial_sum_error_squared": "1/(4*3^M)",
            "hostile_cross_length_addresses": [str(left), str(right)],
            "checks": dict(sorted(CHECKS.items())), "total_checks": sum(CHECKS.values())}


def main():
    data = audit()
    output = ("Exact source-address two-clock audit\n"
              f"Composition/parity cases: {data['total_word_cases']}\n"
              "q in {3,5,7,9,11}; all costs through 12\n"
              "Odd-Haar Gram matrices: every cost cutoff 1..12\n"
              "Exact cyclotomic cost sums: cutoffs 1..8, 3 damping values, 13 signed/zero seeds\n"
              f"Total explicit checks: {data['total_checks']}\n"
              "All checks passed. Fixed integer estimates and Collatz remain open.\n")
    root = Path(__file__).resolve().parents[2] / "05-knowledge" / "results"
    stem = "source_address_two_clocks_20261004"
    (root / (stem + ".json")).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    (root / (stem + ".out")).write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
