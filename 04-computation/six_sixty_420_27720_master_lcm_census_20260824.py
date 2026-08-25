#!/usr/bin/env python3
"""Exact structural audit for the 6, 60, 420, 27720 repository pattern.

The default report checks the master initial-segment LCM mechanism and the
main derived clocks found in the canon.  ``--scan-corpus`` additionally gives
a deliberately semantics-free token census over ``04-computation`` and
``05-knowledge/results``.  The latter is a contamination diagnostic, not
evidence that every matching token has the same mathematical meaning.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import gcd, lcm, prod
from pathlib import Path
import re


TARGETS = (6, 60, 420, 27720)
SCRIPT_REL = Path("04-computation/six_sixty_420_27720_master_lcm_census_20260824.py")
OUTPUT_REL = Path("05-knowledge/results/six_sixty_420_27720_master_lcm_census_20260824.out")


def primes_upto(bound: int) -> list[int]:
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, int(bound**0.5) + 1):
        if sieve[p]:
            sieve[p * p : bound + 1 : p] = b"\x00" * (((bound - p * p) // p) + 1)
    return [n for n in range(2, bound + 1) if sieve[n]]


def valuation(n: int, p: int) -> int:
    answer = 0
    while n % p == 0:
        n //= p
        answer += 1
    return answer


def radical(n: int) -> int:
    answer = 1
    d = 2
    while d * d <= n:
        if n % d == 0:
            answer *= d
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    return answer * n if n > 1 else answer


def lcm_upto(bound: int) -> int:
    answer = 1
    for n in range(1, bound + 1):
        answer = lcm(answer, n)
    return answer


def mobius(n: int) -> int:
    answer = 1
    for p in primes_upto(int(n**0.5) + 1):
        if n % p:
            continue
        n //= p
        answer = -answer
        if n % p == 0:
            return 0
        while n % p == 0:
            n //= p
    if n > 1:
        answer = -answer
    return answer


def vgrid_moment(L: int, v: int) -> Fraction:
    mertens = [0] * (L + 1)
    for n in range(1, L + 1):
        mertens[n] = mertens[n - 1] + mobius(n)
    weights = [0] + [mertens[L // d] for d in range(1, L + 1)]
    answer = Fraction()
    for d in range(1, L + 1):
        d_prime = d // gcd(d, v)
        for e in range(1, L + 1):
            e_prime = e // gcd(e, v)
            answer += Fraction(weights[d] * weights[e], lcm(d_prime, e_prime))
    return answer


def thm4042_product_formula(P: int) -> int:
    """Equation (7b) of THM-4042, independently from the gcd identity."""
    answer = 1
    for ell in primes_upto(P - 2):
        power = ell
        exponent = 1
        while power * ell <= P - 2:
            power *= ell
            exponent += 1
        answer *= ell ** min(exponent, valuation(P - 1, ell) + 1)
    return answer


def subset_lcm_census(max_size: int | None) -> tuple[int, Counter[int]]:
    counts: Counter[int] = Counter()
    upper = 12 if max_size is None else max_size
    total = 0
    for size in range(1, upper + 1):
        for subset in combinations(range(1, 13), size):
            counts[lcm(*subset)] += 1
            total += 1
    return total, counts


def rank_by_frequency(counts: Counter[int], target: int) -> int:
    return 1 + sum(multiplicity > counts[target] for multiplicity in counts.values())


def structural_report() -> None:
    print("STRUCTURAL EXACT AUDIT")
    print("targets:", TARGETS)
    indices = (3, 5, 7, 11)
    ladder_values = tuple(lcm_upto(index) for index in indices)
    assert ladder_values == TARGETS
    assert all(b % a == 0 for a, b in zip(TARGETS, TARGETS[1:]))
    print("initial-segment LCM identity:", tuple(zip(indices, ladder_values)))
    print("successive sampled ratios:", tuple(b // a for a, b in zip(TARGETS, TARGETS[1:])))
    print("full intervening LCM staircase n=3..11:", tuple(lcm_upto(n) for n in range(3, 12)))
    print("THM-4044 observer tower: 6|60|420|27720; root sets grow and kernels shrink")

    prime_bound = 10_000
    prime_set = set(primes_upto(prime_bound))
    equality_primes: list[int] = []
    M = 1
    for m in range(1, prime_bound - 1):
        M = lcm(M, m)
        P = m + 2
        if P not in prime_set or P < 5:
            continue
        Pi_gcd = gcd(M, (P - 1) * radical(M))
        Pi_product = thm4042_product_formula(P)
        assert Pi_gcd == Pi_product
        if Pi_gcd == M:
            equality_primes.append(P)
    assert equality_primes == [5, 7, 13]
    print("THM-4042 master identity verified for every prime P<=10000")
    print("M_P=L_(P-2); Pi_P = gcd(M_P, (P-1) rad(M_P)) for prime P>=5")
    print("P=3 scope boundary: exact Pi_3=2, proposed gcd right side=1")
    print("Pi_P=L_(P-2) primes P<=10000:", equality_primes)
    previous_prime = {5: 3, 7: 5, 11: 7, 13: 11, 17: 13}
    for P in previous_prime:
        M = lcm_upto(P - 2)
        Pi = gcd(M, (P - 1) * radical(M))
        previous_lcm = lcm_upto(previous_prime[P])
        print(
            f"P={P}: Pi={Pi}, L_(P-2)={M}, "
            f"L_previous_prime={previous_lcm}, previous-prime equality={Pi == previous_lcm}"
        )

    print("THM-1420 L_floor(n/2) target bands:")
    print("  6: n=6..7; 60: n=10..13; 420: n=14..15; 27720: n=22..25")
    print("THM-1105 blocker speeds L_Q at Q=3,5,7,11:", ladder_values)
    for L, Q in zip(indices, ladder_values):
        assert vgrid_moment(L, 0) == 1
        for p in primes_upto(L):
            assert vgrid_moment(L, Q // p) == Fraction(1, p)
    assert vgrid_moment(5, 2) == Fraction(47, 30)
    assert vgrid_moment(5, 4) == Fraction(3, 5)
    print("THM-879 exact minimal v-periods at L=3,5,7,11:", ladder_values)
    print("THM-879 valuation control at L=5:", "S_2=47/30; S_4=3/5")
    P_121 = prod(primes_upto(121))
    L_121 = lcm_upto(121)
    U_121 = 84 * P_121
    V_121 = 84 * L_121
    assert U_121 % 121 == 55 and V_121 % 121 == 0
    assert all(V_121 % q == 0 for q in range(1, 122))
    print("THM-3848 valuation-depth hostile: 84*primorial(121)=55 mod 121; 84*L_121=0 mod every q<=121")

    B_420 = (1, 2, 3, 4, 5, 6)
    B_27720 = (3, 5, 6, 8, 9, 11)
    base_420 = lcm(*B_420)
    base_27720 = lcm(*B_27720)
    assert (base_420, 7 * base_420) == (60, 420)
    assert (base_27720, 7 * base_27720) == (3960, 27720)
    print("THM-563 period clocks 7*lcm(B):", (base_420, 7 * base_420), (base_27720, 7 * base_27720))

    R = (1, 2, 3, 5, 7, 11, 13)
    differences = sorted({abs(a - b) for a, b in combinations(R, 2)})
    assert differences == [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12]
    difference_clock = lcm(*differences)
    assert difference_clock == 3960
    assert lcm(difference_clock, 14) == 27720
    assert lcm(60, 14) == 420
    print("THM-1226 difference set:", differences)
    print("parallel modulus joins:", "lcm(60,14)=420; lcm(3960,14)=27720")
    print("commuting scale relation:", "3960=66*60; 27720=66*420")

    moduli_4000 = tuple(range(6, 12))
    assert lcm(*moduli_4000) == 27720
    assert prod(moduli_4000) == 12 * 27720
    print("THM-4000 CRT observer:", "lcm(6,7,8,9,10,11)=27720; product=12*27720")
    assert lcm(*range(2, 13)) == 27720
    assert 13 * lcm(*range(2, 13)) == 360360
    print("THM-823 replay continuation:", "lcm(2,...,12)*13=27720*13=360360")

    endpoint_moduli = (33, 3, 44, 55, 165, 198)
    assert 14 * lcm(*endpoint_moduli) == 27720
    print("THM-2258 endpoint common mesh:", "14*lcm(33,3,44,55,165,198)=27720")

    ruler = (1, 9, 10, 11, 12, 14)
    ruler_lcm = lcm(*ruler)
    assert ruler_lcm == 13860
    assert 14 * ruler_lcm == 194040
    assert (14 * ruler_lcm) // 7 == 27720
    print("THM-3106 derived grid:", "lcm(E)=13860; 14*lcm(E)=194040; /7=27720")

    limited_total, limited_counts = subset_lcm_census(max_size=6)
    full_total, full_counts = subset_lcm_census(max_size=None)
    assert limited_total == 2509 and full_total == 4095
    expected_limited = {6: (10, 53), 60: (123, 2), 420: (97, 5), 27720: (15, 44)}
    expected_full = {6: 10, 60: 132, 420: 132, 27720: 192}
    print("THM-1415 subset-LCM census on {1,...,12}, subset size<=6:")
    for target in TARGETS:
        pair = (limited_counts[target], rank_by_frequency(limited_counts, target))
        assert pair == expected_limited[target]
        print(f"  {target}: count={pair[0]}/{limited_total}, frequency-rank={pair[1]}")
    print("all nonempty subsets of {1,...,12}:")
    for target in TARGETS:
        assert full_counts[target] == expected_full[target]
        print(f"  {target}: count={full_counts[target]}/{full_total}")

    odd_primes = (3, 5, 7, 11, 13)
    sun_scales = tuple(24 * prod(odd_primes[:k]) for k in range(1, len(odd_primes) + 1))
    assert sun_scales == (72, 360, 2520, 27720, 360360)
    print("THM-4036 Sun hostile scales (partial overlap only):", sun_scales)
    trimode_mass = sum((Fraction(1, d) for d in (1, 2, 3, 4, 5, 7)), Fraction())
    assert trimode_mass == Fraction(1019, 420)
    print("Fibonacci/Farey finite root-mass overlap:", "1+1/2+1/3+1/4+1/5+1/7=1019/420")
    assert all(n in set(primes_upto(500)) for n in (7, 61, 421))
    assert 27721 == 19 * 1459
    print("hostile L_n+1 control:", "7,61,421 prime; 27721=19*1459 composite")
    print("all structural assertions passed")


def corpus_scan(repo: Path) -> None:
    roots = (repo / "04-computation", repo / "05-knowledge" / "results")
    excluded = {repo / SCRIPT_REL, repo / OUTPUT_REL}
    patterns = {target: re.compile(rb"(?<![0-9])" + str(target).encode() + rb"(?![0-9])") for target in TARGETS}
    token_counts = Counter()
    line_counts = Counter()
    file_counts = Counter()
    files_scanned = 0
    bytes_scanned = 0
    all_four_files = 0
    denominator_27720 = 0
    for root in roots:
        for path in sorted(root.rglob("*")):
            if not path.is_file() or path in excluded:
                continue
            data = path.read_bytes()
            files_scanned += 1
            bytes_scanned += len(data)
            present = 0
            for target, pattern in patterns.items():
                matches = pattern.findall(data)
                token_counts[target] += len(matches)
                if matches:
                    file_counts[target] += 1
                    present += 1
                    line_counts[target] += sum(bool(pattern.search(line)) for line in data.splitlines())
            if present == len(TARGETS):
                all_four_files += 1
            denominator_27720 += len(re.findall(rb"/27720(?![0-9])", data))
    print("RAW TOKEN CENSUS (SEMANTICS-FREE CONTAMINATION DIAGNOSTIC)")
    print("universe: regular files below 04-computation and 05-knowledge/results")
    print("excluded self/output:", SCRIPT_REL.as_posix(), OUTPUT_REL.as_posix())
    print("files scanned:", files_scanned)
    print("bytes scanned:", bytes_scanned)
    for target in TARGETS:
        print(
            f"{target}: exact tokens={token_counts[target]}, "
            f"matching lines={line_counts[target]}, matching files={file_counts[target]}"
        )
    print("files containing all four exact tokens:", all_four_files)
    print("literal /27720 denominator tokens:", denominator_27720)
    if token_counts[27720]:
        print("/27720 share of 27720 tokens: {:.6%}".format(denominator_27720 / token_counts[27720]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-corpus", action="store_true")
    args = parser.parse_args()
    structural_report()
    if args.scan_corpus:
        repo = Path(__file__).resolve().parents[1]
        corpus_scan(repo)


if __name__ == "__main__":
    main()
