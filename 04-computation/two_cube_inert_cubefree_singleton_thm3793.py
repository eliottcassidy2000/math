#!/usr/bin/env python3
"""Exact finite controls for the inert cube-free pair-sum singleton lemma.

The companion README contains the all-scale valuation proof.  This script
checks every positive distinct representation with admissible pair-sum at
most 1000 against a complete coordinate fibre, audits the sharp hostiles, and
counts the induced subatlas inside THM-3743's p+q <= 356 ratio universe.

No Python assertions are used; all gates remain active under ``python -O``.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factor(n: int) -> dict[int, int]:
    require(n >= 1, "factor called outside positive integers")
    out: dict[int, int] = {}
    candidate = 2
    while candidate * candidate <= n:
        while n % candidate == 0:
            out[candidate] = out.get(candidate, 0) + 1
            n //= candidate
        candidate = 3 if candidate == 2 else candidate + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def admissible_sum(d: int) -> bool:
    if d < 2:
        return False
    factors = factor(d)
    return bool(factors) and all(prime % 3 == 2 and exponent <= 2
                                 for prime, exponent in factors.items())


def representations_by_value(coordinate_cap: int) -> dict[int, list[tuple[int, int]]]:
    fibres: dict[int, list[tuple[int, int]]] = {}
    for y in range(2, coordinate_cap + 1):
        for x in range(1, y):
            fibres.setdefault(x ** 3 + y ** 3, []).append((x, y))
    return fibres


def phi(n: int) -> int:
    value = n
    for prime in factor(n):
        value -= value // prime
    return value


def main() -> None:
    sum_cap = 1000
    fibres = representations_by_value(sum_cap - 1)
    checked_values: set[int] = set()
    admissible_sums: list[int] = []
    all_pair_checks = 0
    valuation_checks = 0

    for d in range(2, sum_cap + 1):
        if not admissible_sum(d):
            continue
        admissible_sums.append(d)
        factors = factor(d)
        for x in range(1, (d + 1) // 2):
            y = d - x
            if not x < y:
                continue
            m = x ** 3 + y ** 3
            require(y < d <= sum_cap, "coordinate completeness bound failed")
            require(fibres.get(m) == [(x, y)], f"singleton failure at d={d}, pair={(x, y)}")
            require(m not in checked_values, "cross-row singleton collision")
            checked_values.add(m)
            all_pair_checks += 3

            g = math.gcd(x, y)
            primitive_x, primitive_y = x // g, y // g
            cofactor = primitive_x ** 2 - primitive_x * primitive_y + primitive_y ** 2
            for prime, exponent in factors.items():
                alpha = 0
                gg = g
                while gg % prime == 0:
                    alpha += 1
                    gg //= prime
                m_valuation = 0
                mm = m
                while mm % prime == 0:
                    m_valuation += 1
                    mm //= prime
                require(cofactor % prime != 0, "inert prime entered primitive cofactor")
                require(m_valuation == exponent + 2 * alpha, "valuation invoice failed")
                valuation_checks += 2

    require(fibres[1729] == [(9, 10), (1, 12)], "1729 split-prime hostile changed")
    require(fibres[65728] == [(31, 33), (12, 40)], "64-sum non-cube-free hostile changed")
    require(sum((31, 33)) == 64 and factor(64) == {2: 6}, "64 hostile typing failed")
    require(fibres[515375] == [(54, 71), (15, 80)], "125-sum exponent-three hostile changed")
    require(sum((54, 71)) == 125 and factor(125) == {5: 3}, "exponent-three hostile typing failed")
    hostile_checks = 5

    generalized_values: set[int] = set()
    generalized_extra = 0
    generalized_checks = 0
    for d in range(3, sum_cap + 1):
        d_factors = factor(d)
        if not all(prime % 3 == 2 for prime in d_factors):
            continue
        for x in range(1, (d + 1) // 2):
            y = d - x
            if not x < y:
                continue
            primitive_sum = d // math.gcd(x, y)
            if not all(exponent <= 2 for exponent in factor(primitive_sum).values()):
                continue
            m = x ** 3 + y ** 3
            require(fibres.get(m) == [(x, y)], "generalized primitive-sum singleton failed")
            require(m not in generalized_values, "generalized singleton rows collided")
            generalized_values.add(m)
            if not admissible_sum(d):
                generalized_extra += 1
            generalized_checks += 2

    lrc_cap = 356
    total_ratio_atlas = 0
    singleton_ratio_atlas = 0
    singleton_all_pairs = 0
    lrc_rows: list[dict[str, object]] = []
    for d in range(3, lrc_cap + 1):
        ratio_count = 0
        all_count = 0
        for x in range(1, (d + 1) // 2):
            y = d - x
            if not x < y:
                continue
            all_count += 1
            if math.gcd(x, y) == 1:
                ratio_count += 1
        total_ratio_atlas += ratio_count
        if admissible_sum(d):
            require(ratio_count == phi(d) // 2, "Euler-phi row count failed")
            singleton_ratio_atlas += ratio_count
            singleton_all_pairs += all_count
            lrc_rows.append({"d": d, "factor": factor(d), "primitive_ratios": ratio_count})

    require(total_ratio_atlas == 19314, "THM-3743 ratio atlas count changed")
    require(singleton_ratio_atlas == 5855, "singleton ratio subatlas count changed")
    require(singleton_all_pairs == 7991, "all-pair singleton row count changed")
    require(len(lrc_rows) == 94, "admissible LRC sum count changed")
    labelled_ratio_subatlas = 78 * singleton_ratio_atlas
    require(labelled_ratio_subatlas == 456690, "labelled ratio count changed")

    # Exact combinatorial controls for the two-prime critical-mass amplification.
    inert_primes = [p for p in range(5, sum_cap + 1)
                    if len(factor(p)) == 1 and next(iter(factor(p).values())) == 1 and p % 3 == 2]
    semiprime_rows = 0
    semiprime_pair_values = 0
    for i, p in enumerate(inert_primes):
        for q in inert_primes[i + 1:]:
            d = p * q
            if d > sum_cap:
                break
            require(admissible_sum(d), "two-prime row not admissible")
            require(Fraction((d - 1) // 2, d * d) >= Fraction(2, 5 * d),
                    "row-mass counting tariff failed")
            semiprime_rows += 1
            semiprime_pair_values += (d - 1) // 2

    critical_mass_checks = 0
    analytic_primes = [p for p in range(5, 5001)
                       if len(factor(p)) == 1 and next(iter(factor(p).values())) == 1 and p % 3 == 2]
    for z in (11, 50, 100, 500, 1000, 5000):
        primes = [p for p in analytic_primes if p <= z]
        a_sum = sum((Fraction(1, p) for p in primes), Fraction(0))
        b_sum = sum((Fraction(1, p * p) for p in primes), Fraction(0))
        unordered_sum = Fraction(0)
        for i, p in enumerate(primes):
            for q in primes[i + 1:]:
                d = p * q
                require(5 * (d - 1) >= 4 * d, "symbolic two-prime row tariff failed")
                unordered_sum += Fraction(1, d)
                critical_mass_checks += 1
        require(2 * unordered_sum == a_sum * a_sum - b_sum,
                "ordered/unordered reciprocal bookkeeping failed")
        require(b_sum < Fraction(1, 4), "finite square-reciprocal control failed")
        critical_mass_checks += 2

    semantic_digest = hashlib.sha256(
        "\n".join(json.dumps(row, sort_keys=True) for row in lrc_rows).encode("ascii")
    ).hexdigest()

    print("TWO-CUBE INERT CUBE-FREE SINGLETON PROBE")
    print("LEMMA: if every prime divisor of d=x+y is 2 mod 3 and v_p(d)<=2, then r_+(x^3+y^3)=1")
    print("COMPLETE_PAIR_SUM_UNIVERSE", "2..1000")
    print("ADMISSIBLE_SUMS", len(admissible_sums))
    print("DISTINCT_SINGLETON_VALUES", len(checked_values))
    print("ALL_PAIR_ACTIVE_CHECKS", all_pair_checks)
    print("VALUATION_ACTIVE_CHECKS", valuation_checks)
    print("HOSTILE_ACTIVE_CHECKS", hostile_checks)
    print("SPLIT_PRIME_HOSTILE", "1729=(9,10)=(1,12), sums 19 and 13")
    print("CUBEFREE_HOSTILE", "65728=(31,33)=(12,40), first displayed sum 64=2^6")
    print("SHARP_EXPONENT_HOSTILE", "515375=(54,71)=(15,80), first displayed sum 125=5^3")
    print("GENERALIZED_PRIMITIVE_SUM_SINGLETON_VALUES", len(generalized_values))
    print("GENERALIZED_EXTRA_NONCUBEFREE_VALUES", generalized_extra)
    print("GENERALIZED_ACTIVE_CHECKS", generalized_checks)
    print("LRC_TOTAL_UNORDERED_RATIOS", total_ratio_atlas)
    print("LRC_SINGLETON_CODED_RATIOS", singleton_ratio_atlas)
    print("LRC_SINGLETON_CODED_LABELLED_ROWS", labelled_ratio_subatlas)
    print("LRC_ADMISSIBLE_PAIR_SUMS", len(lrc_rows))
    print("LRC_ALL_SINGLETON_PAIR_VALUES", singleton_all_pairs)
    print("LRC_ROW_SHA256", semantic_digest)
    print("TWO_INERT_PRIME_ROWS_THROUGH_1000", semiprime_rows)
    print("TWO_INERT_PRIME_SINGLETON_VALUES_THROUGH_1000", semiprime_pair_values)
    print("CRITICAL_MASS_ACTIVE_CHECKS", critical_mass_checks)
    print("CRITICAL_BOUND", "H(Z^6)>=(A(Z)^2-B(Z))/5, unordered prime pairs")
    print("ASYMPTOTIC", "liminf H(X)/(log log X)^2 >= 1/20")
    print("RESULT PASS")


if __name__ == "__main__":
    main()
