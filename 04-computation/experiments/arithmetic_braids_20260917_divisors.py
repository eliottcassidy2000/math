"""Exact divisor-balance audit and ordered almost-prime sandwiches.

No third-party dependencies. Run from the repository root:
  python 04-computation/experiments/arithmetic_braids_20260917_divisors.py

Default universe: integers 2..10000 for direct divisor audit, exponent
profiles of support 1..9 with exponents 1..8, and centers 6k, 1<=k<=10^6.
All verification uses explicit exceptions and survives Python -O.
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from fractions import Fraction
from itertools import combinations_with_replacement
import json
from math import prod
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor_trial(n):
    result = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def divisor_counts_direct(n):
    factors = factor_trial(n)
    divisors = [1]
    for p, a in factors.items():
        divisors = [d * p**b for d in divisors for b in range(a + 1)]
    proper = [d for d in divisors if 1 < d < n]
    squarefree = [d for d in proper if all(d % (p * p) for p in factors)]
    primes = [d for d in proper if d in factors]
    return len(proper), len(squarefree), len(primes)


def profile_counts(profile):
    r = len(profile)
    squarefree = all(a == 1 for a in profile)
    prime = profile == (1,)
    return prod(a + 1 for a in profile) - 2, (1 << r) - 1 - squarefree, r - prime


def is_balance_profile(profile):
    return profile in ((1,), (3,), (1, 1, 2))


def defect_expected(profile):
    f, s, u = profile_counts(profile)
    r = len(profile)
    if profile == (1,):
        return 0
    if max(profile) == 1:
        return -r
    return prod(a + 1 for a in profile) - (1 << r) - r - 1


def omega_sieve(limit):
    omega = bytearray(limit + 1)
    for p in range(2, limit + 1):
        if omega[p] != 0:
            continue
        power = p
        while power <= limit:
            for multiple in range(power, limit + 1, power):
                omega[multiple] += 1
            power *= p
    return omega


def multiplicative_partitions(total, length, minimum=2):
    """Ordered factor profiles b1<=...<=br, bi>=2, product=total."""
    if length == 1:
        if total >= minimum:
            yield (total,)
        return
    b = minimum
    while b**length <= total:
        if total % b == 0:
            for rest in multiplicative_partitions(total // b, length - 1, b):
                yield (b,) + rest
        b += 1


def all_profiles_at_defect(d):
    """Complete support/exponent classification at fixed defect, no size cap."""
    results = []
    if d == 0:
        results.append((1,))
    elif d <= -2:
        results.append((1,) * (-d))
    r = 1
    while True:
        minimum_defect = (1 << (r - 1)) - r - 1
        if r >= 2 and minimum_defect > d:
            break
        target = (1 << r) + r + 1 + d
        if minimum_defect <= d:
            for factors in multiplicative_partitions(target, r):
                profile = tuple(b - 1 for b in factors)
                if max(profile) > 1:
                    require(defect_expected(profile) == d, ("inverse defect", d, profile))
                    results.append(profile)
        r += 1
    return sorted(results)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--centers", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.centers >= 1000, "Need at least 1000 centers for controls")
    outpath = args.output or Path(__file__).with_suffix(".json")

    count = 0
    zeros = set()
    sharp_minima = {}
    for r in range(1, 10):
        minimum = None
        minimizers = []
        for profile in combinations_with_replacement(range(1, 9), r):
            f, s, u = profile_counts(profile)
            d = f - s - u
            require(s >= u, ("S>=U", profile))
            require(d == defect_expected(profile), ("defect", profile))
            require((d == 0) == is_balance_profile(profile), ("classification", profile))
            if d == 0:
                zeros.add(profile)
            if max(profile) > 1:
                if minimum is None or d < minimum:
                    minimum, minimizers = d, [profile]
                elif d == minimum:
                    minimizers.append(profile)
            count += 1
        require(minimum == (1 << (r - 1)) - r - 1, ("minimum", r))
        require(minimizers == [(1,) * (r - 1) + (2,)], ("minimizer", r))
        sharp_minima[str(r)] = {"defect": minimum, "profile": minimizers[0]}

    equality_values = []
    for n in range(2, 10_001):
        profile = tuple(sorted(factor_trial(n).values()))
        direct = divisor_counts_direct(n)
        require(direct == profile_counts(profile), ("direct divisors", n))
        f, s, u = direct
        require((f == s + u) == is_balance_profile(profile), ("integer classification", n))
        if f == s + u and len(equality_values) < 30:
            equality_values.append(n)
    controls = {str(n): divisor_counts_direct(n) for n in (2, 4, 8, 12, 16, 30, 60, 210, 420)}
    require(controls["2"] == (0, 0, 0), "prime positive")
    require(controls["8"] == (2, 1, 1), "cube positive")
    require(controls["60"] == (10, 7, 3), "p2qr positive")
    require(controls["12"] == (4, 3, 2), "p2q hostile")
    complete_small_defects = {str(d): all_profiles_at_defect(d) for d in range(-6, 13)}
    require(all_profiles_at_defect(0) == [(1,), (1, 1, 2), (3,)], "inverse equality profiles")
    require(all_profiles_at_defect(1) == [(1, 3), (4,)], "defect one profiles")
    require(all_profiles_at_defect(2) == [(2, 2), (5,)], "defect two profiles")

    omega = omega_sieve(6 * args.centers + 1)
    for n in range(2, 10_001):
        require(omega[n] == sum(factor_trial(n).values()), ("independent omega", n))

    matrix = [[0] * 4 for _ in range(4)]
    raw = Counter()
    first = {}
    checkpoints = {}
    scales = sorted({100, 1000, 10_000, 100_000, 1_000_000, args.centers})
    scales = {k for k in scales if k <= args.centers}
    first_orientation_inequality = None
    for k in range(1, args.centers + 1):
        left, right = 6 * k - 1, 6 * k + 1
        a, b = omega[left], omega[right]
        matrix[min(a, 4) - 1][min(b, 4) - 1] += 1
        raw[(a, b)] += 1
        key = f"{a},{b}"
        if a <= 3 and b <= 3 and key not in first:
            first[key] = {"k": k, "center": 6 * k, "left": left, "right": right,
                          "left_factors": factor_trial(left), "right_factors": factor_trial(right)}
        if first_orientation_inequality is None and raw[(1, 2)] != raw[(2, 1)]:
            first_orientation_inequality = {"k": k, "PS": raw[(1, 2)], "SP": raw[(2, 1)]}
        if k <= 10_000:
            lf, rf = factor_trial(left), factor_trial(right)
            left_negative = sum(e for p, e in lf.items() if p % 6 == 5)
            right_negative = sum(e for p, e in rf.items() if p % 6 == 5)
            require(left_negative % 2 == 1, ("left mod6 parity", k))
            require(right_negative % 2 == 0, ("right mod6 parity", k))
        if k in scales:
            checkpoints[str(k)] = {
                "max_center": 6 * k,
                "matrix_omega_1_2_3_ge4": [row[:] for row in matrix],
                "within_omega_le3": sum(matrix[i][j] for i in range(3) for j in range(3)),
                "PS_minus_SP": raw[(1, 2)] - raw[(2, 1)],
            }
    require(sum(map(sum, matrix)) == args.centers, "matrix exhaustion")
    require(len(first) == 9, "all nine states must have positive controls")
    require(first_orientation_inequality == {"k": 4, "PS": 1, "SP": 0}, "orientation hostile")

    # A second implementation of the full 4x4 matrix at k<=1000.
    independent = [[0] * 4 for _ in range(4)]
    for k in range(1, 1001):
        a = sum(factor_trial(6 * k - 1).values())
        b = sum(factor_trial(6 * k + 1).values())
        independent[min(a, 4) - 1][min(b, 4) - 1] += 1
    require(independent == checkpoints["1000"]["matrix_omega_1_2_3_ge4"], "independent matrix")

    # CRT-local reflection k -> -k exchanges left/right divisibility words.
    local_reflection = []
    for modulus, primes in ((5, (5,)), (35, (5, 7)), (385, (5, 7, 11))):
        hist = Counter()
        for k in range(modulus):
            left_hits = sum((6 * k - 1) % p == 0 for p in primes)
            right_hits = sum((6 * k + 1) % p == 0 for p in primes)
            hist[(left_hits, right_hits)] += 1
        require(all(hist[(a, b)] == hist[(b, a)] for a, b in hist), "CRT reflection")
        polynomial = Counter({(0, 0): 1})
        for p in primes:
            next_polynomial = Counter()
            for (a, b), coefficient in polynomial.items():
                next_polynomial[(a, b)] += (p - 2) * coefficient
                next_polynomial[(a + 1, b)] += coefficient
                next_polynomial[(a, b + 1)] += coefficient
            polynomial = next_polynomial
        require(hist == polynomial, "CRT three-branch generating polynomial")
        mean = sum(Fraction(1, p) for p in primes)
        covariance = sum(Fraction(a * b * v, modulus) for (a, b), v in hist.items()) - mean**2
        require(covariance == -sum(Fraction(1, p * p) for p in primes), "CRT covariance")
        local_reflection.append({"modulus": modulus, "primes": primes,
                                 "mean_each": str(mean), "covariance": str(covariance),
                                 "matrix": {f"{a},{b}": v for (a, b), v in sorted(hist.items())}})

    result = {
        "status": "FINITE-EXACT; universal proofs are in companion markdown",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "universe": {"centers_k": [1, args.centers], "centers_W": [6, 6 * args.centers],
                     "endpoint_convention": "ordered (W-1,W+1), Omega counted with multiplicity",
                     "no_filter_on_centers": True, "matrix_tail": "Omega>=4",
                     "divisor_integer_audit": [2, 10000], "profile_support": [1, 9], "profile_exponents": [1, 8]},
        "divisor_profiles_checked": count,
        "equality_profiles": sorted(zeros),
        "equality_values_first30": equality_values,
        "sharp_minimum_nonsquarefree_defect": sharp_minima,
        "complete_defect_profiles_minus6_to12": complete_small_defects,
        "divisor_controls_F_S_U": controls,
        "sandwich_checkpoints": checkpoints,
        "first_sandwich_witnesses": first,
        "first_PS_SP_inequality": first_orientation_inequality,
        "all_omega_pairs": {f"{a},{b}": v for (a, b), v in sorted(raw.items())},
        "mod6_factor_parity_audit_k": [1, 10000],
        "CRT_local_reflection": local_reflection,
        "independent_checks": {"trial_factor_omega": [2, 10000], "trial_factor_matrix_k": [1, 1000]},
        "all_checks_passed": True,
    }
    outpath.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"output": str(outpath), "all_checks_passed": True,
                      "profiles": count, "matrix": matrix, "checkpoints": checkpoints}, indent=2))


if __name__ == "__main__":
    main()
