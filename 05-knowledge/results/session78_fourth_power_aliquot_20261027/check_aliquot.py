#!/usr/bin/env python3
"""Exact, dependency-free Session 78 aliquot certificate.

No assertion statements: every verification remains active under python -O.
The global inverse-fiber completeness proof is in ALIQUOT.md; the search below
implements its proved bound and explicit prime-square exception.
"""

from __future__ import annotations

import argparse
import json
from math import isqrt
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_positive(n: int) -> None:
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("the aliquot map has positive-integer domain")


def proper_divisors(n: int) -> list[int]:
    require_positive(n)
    result = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            if d < n:
                result.add(d)
            if n // d < n:
                result.add(n // d)
    return sorted(result)


def aliquot(n: int) -> int:
    return sum(proper_divisors(n))


def factorization(n: int) -> list[tuple[int, int]]:
    require_positive(n)
    out = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def sigma_factorized(n: int) -> int:
    result = 1
    for p, e in factorization(n):
        result *= (p ** (e + 1) - 1) // (p - 1)
    return result


def is_prime(n: int) -> bool:
    return n >= 2 and all(n % p for p in range(2, isqrt(n) + 1))


def odd_part(n: int) -> int:
    require_positive(n)
    while n % 2 == 0:
        n //= 2
    return n


def divisor_sieve(limit: int) -> list[int]:
    """Independent all-proper-multiples computation; no factorization."""
    values = [0] * (limit + 1)
    for d in range(1, limit // 2 + 1):
        for multiple in range(2 * d, limit + 1, d):
            values[multiple] += d
    return values


def reconstruct_from_doubling_pair(x: int, y: int) -> int:
    """Invert (s(n), s(2n)); reject every pair outside the exact image."""
    if any(not isinstance(v, int) or isinstance(v, bool) for v in (x, y)):
        raise ValueError("observations must be exact integers")
    c = y - 2 * x
    d = y - x
    if c <= 0 or d <= 0:
        raise ValueError("positive carry and positive difference required")
    power_of_two = d & -d
    odd_core = d // power_of_two
    m = 2 * c - odd_core
    if m <= 0 or sigma_factorized(m) != c:
        raise ValueError("recovered odd part fails the divisor-sum guard")
    return power_of_two * m


def reconstruct_from_prime_dilation_pair(p: int, x: int, y: int) -> int:
    """Invert (s(n),s(p*n)) for known prime p; validate the exact image."""
    if not isinstance(p, int) or isinstance(p, bool) or not is_prime(p):
        raise ValueError("the dilation parameter must be a known prime")
    if any(not isinstance(v, int) or isinstance(v, bool) for v in (x, y)):
        raise ValueError("observations must be exact integers")
    c = y - p * x
    d = y - x
    if c <= 0 or d <= 0:
        raise ValueError("positive carry and positive difference required")
    power_of_p = 1
    p_free_d = d
    while p_free_d % p == 0:
        p_free_d //= p
        power_of_p *= p
    numerator = p * c - p_free_d
    if numerator % (p - 1):
        raise ValueError("recovered p-free part is not integral")
    m = numerator // (p - 1)
    if m <= 0 or m % p == 0 or sigma_factorized(m) != c:
        raise ValueError("recovered p-free part fails the divisor-sum guard")
    return power_of_p * m


def complete_fiber(target: int) -> dict:
    if target < 2:
        raise ValueError("finite-fiber theorem here requires target >= 2")
    bound = (target - 1) ** 2 // 4
    candidates = list(range(1, bound + 1))
    prime_square = (target - 1) ** 2 if is_prime(target - 1) else None
    if prime_square is not None and prime_square > bound:
        candidates.append(prime_square)
    fiber = [n for n in candidates if aliquot(n) == target]
    return {
        "target": target,
        "ordinary_source_bound": bound,
        "possible_prime_square_exception": prime_square,
        "complete_fiber": fiber,
        "proper_divisor_witnesses": {
            str(n): proper_divisors(n) for n in fiber
        },
    }


def chain(seed: int) -> list[int]:
    result = []
    n = seed
    while n and n not in result:
        result.append(n)
        n = aliquot(n)
    result.append(n)
    return result


def main() -> dict:
    limit = 20000
    sieve = divisor_sieve(limit)
    for n in range(1, limit + 1):
        direct = aliquot(n)
        require(direct == sieve[n], f"direct/sieve mismatch at {n}")
        require(direct == sigma_factorized(n) - n,
                f"direct/factorized mismatch at {n}")
        factors = factorization(n)
        if direct > 1:
            prime_square = len(factors) == 1 and factors[0][1] == 2
            if prime_square:
                require(n == (direct - 1) ** 2,
                        f"prime-square exception failure at {n}")
            else:
                require(4 * n <= (direct - 1) ** 2,
                        f"ordinary source-bound failure at {n}")

    expected_fibers = {
        5: [], 10: [14], 11: [21], 12: [121], 13: [27, 35],
        25: [95, 119, 143], 37: [155, 203, 299, 323],
        50: [40, 94], 100: [124, 194],
    }
    fibers = []
    for target, expected in expected_fibers.items():
        certificate = complete_fiber(target)
        require(certificate["complete_fiber"] == expected,
                f"fiber mismatch for {target}")
        # Independent cross-check through the weaker, also proved universal
        # bound (target - 1)^2. This includes every prime-square exception.
        universal_bound = (target - 1) ** 2
        broad_fiber = [n for n in range(1, universal_bound + 1)
                       if sieve[n] == target]
        require(broad_fiber == expected,
                f"independent universal-bound fiber mismatch for {target}")
        certificate["universal_bound_crosscheck"] = universal_bound
        fibers.append(certificate)

    doubling_checks = 10000
    for n in range(1, doubling_checks + 1):
        carry = sigma_factorized(odd_part(n))
        require(sieve[2 * n] == 2 * sieve[n] + carry,
                f"doubling law failure at {n}")
        require(reconstruct_from_doubling_pair(sieve[n], sieve[2 * n]) == n,
                f"two-observation reconstruction failed at {n}")

    # Exhaustively check the pair-image guard on a rectangle. For this
    # rectangle d<=150 and c<=150, so every possible reconstructed source
    # satisfies n=2^(a+1)c-d < 2dc <= 45000. The independently sieved source
    # image therefore gives a complete reference, not a sampled reference.
    image_limit = 45000
    image_sieve = divisor_sieve(11 * image_limit)
    source_pairs = {}
    for n in range(1, image_limit + 1):
        pair = (image_sieve[n], image_sieve[2 * n])
        require(pair not in source_pairs,
                f"two-observation injectivity failed at {n}")
        source_pairs[pair] = n
    pair_guard_checks = 0
    pair_guard_accepted = 0
    for x in range(51):
        for y in range(1, 151):
            reference = source_pairs.get((x, y))
            try:
                reconstructed = reconstruct_from_doubling_pair(x, y)
            except ValueError:
                require(reference is None,
                        f"valid pair wrongly rejected: {(x, y)}")
            else:
                require(reconstructed == reference,
                        f"invalid pair wrongly accepted: {(x, y)}")
                pair_guard_accepted += 1
            pair_guard_checks += 1

    prime_dilation_reports = []
    for p in (2, 3, 5, 7, 11):
        for n in range(1, 5001):
            m = n
            while m % p == 0:
                m //= p
            x, y = image_sieve[n], image_sieve[p * n]
            require(y == p * x + sigma_factorized(m),
                    f"prime-dilation carry law failed at {(p, n)}")
            require(reconstruct_from_prime_dilation_pair(p, x, y) == n,
                    f"prime-dilation reconstruction failed at {(p, n)}")
        # For x in [0,30], y in [1,100], positive c,d are <=100.
        # Any actual source has n=(p^(a+1)c-d)/(p-1)
        # < p*d*c/(p-1) <= 20000. Thus the 45000-source sieve
        # below is a complete reference for arbitrary-pair image checking.
        reference = {}
        for n in range(1, image_limit + 1):
            pair = (image_sieve[n], image_sieve[p * n])
            require(pair not in reference,
                    f"prime-dilation injectivity failed at {(p, n)}")
            reference[pair] = n
        accepted = 0
        for x in range(31):
            for y in range(1, 101):
                source = reference.get((x, y))
                try:
                    recovered = reconstruct_from_prime_dilation_pair(p, x, y)
                except ValueError:
                    require(source is None,
                            f"prime-dilation valid pair rejected: {(p, x, y)}")
                else:
                    require(recovered == source,
                            f"prime-dilation invalid pair accepted: {(p, x, y)}")
                    accepted += 1
        prime_dilation_reports.append({
            "prime": p, "carry_and_reconstruction_checks": 5000,
            "injectivity_checks_through": image_limit,
            "pair_guard_checks": 3100, "pair_guard_accepted": accepted,
            "pair_guard_rectangle": {"x": [0, 30], "y": [1, 100]},
        })
    try:
        reconstruct_from_prime_dilation_pair(4, 1, 16)
    except ValueError:
        pass
    else:
        raise RuntimeError("composite dilation p=4 accepted")
    require(aliquot(3) == 1 and aliquot(12) == 16
            and (4 * (16 - 4) - (16 - 1)) // 3 == 11
            and 11 != 3,
            "composite-dilation formal inverse hostile control failed")

    tower_checks = 0
    for m in range(1, 202, 2):
        c = sigma_factorized(m)
        for a in range(9):
            n = 2 ** a * m
            expected = 2 ** a * (2 * c - m) - c
            require(aliquot(n) == expected, f"tower formula failed at {n}")
            tower_checks += 1

    expected_chains = {
        5: [5, 1, 0],
        10: [10, 8, 7, 1, 0],
        25: [25, 6, 6],
        50: [50, 43, 1, 0],
        100: [100, 117, 65, 19, 1, 0],
        74: [74, 40, 50, 43, 1, 0],
        70: [70, 74, 40, 50, 43, 1, 0],
    }
    for seed, expected in expected_chains.items():
        require(chain(seed) == expected, f"chain mismatch from {seed}")

    require(aliquot(27) == aliquot(35) == 13, "scalar collision lost")
    require(aliquot(54) == 66 and aliquot(70) == 74,
            "doubling collision control failed")
    require(sigma_factorized(27) == 40 and sigma_factorized(35) == 48,
            "carry values failed")
    require(aliquot(54) != aliquot(70),
            "incorrect scalar-only doubling quotient was accepted")
    require(aliquot(50) != 2 * aliquot(25),
            "incorrect bare doubling law was accepted")
    require(complete_fiber(12)["complete_fiber"] == [121]
            and 121 > complete_fiber(12)["ordinary_source_bound"],
            "prime-square exception was dropped")
    for invalid in (0, -1):
        try:
            aliquot(invalid)
        except ValueError:
            pass
        else:
            raise RuntimeError(f"invalid domain input {invalid} accepted")
    for invalid in (0, 1):
        try:
            complete_fiber(invalid)
        except ValueError:
            pass
        else:
            raise RuntimeError(f"invalid finite-fiber target {invalid} accepted")

    semiprime_pairs = [
        [p, 36 - p] for p in range(3, 18)
        if is_prime(p) and is_prime(36 - p) and p < 36 - p
    ]
    require(semiprime_pairs == [[5, 31], [7, 29], [13, 23], [17, 19]],
            "37 semiprime pair list mismatch")
    require([p * q for p, q in semiprime_pairs] == expected_fibers[37],
            "37 complete fiber is not exactly those semiprimes")
    require(sigma_factorized(25) == 31 and is_prime(31),
            "25 divisor-sum prime bridge failed")
    require(aliquot(5 * sigma_factorized(25)) == 37,
            "5,25,37 factorization bridge failed")
    require(aliquot(100) == 9 * 13 and sigma_factorized(9) == 13,
            "100,13 factorization bridge failed")
    prime_square_bridges = []
    for p in (3, 5, 7, 11, 13, 17, 19):
        q = p * p + p + 1
        slope = 2 * q - p * p
        actual = aliquot(p * q)
        defect = (p + 1) * (aliquot(q) - 1)
        require(sigma_factorized(p * p) == aliquot(p ** 3) == q,
                f"prime-square carry identity failed at {p}")
        require(actual - slope == defect,
                f"prime-square bridge defect failed at {p}")
        require((actual == slope) == is_prime(q),
                f"prime-square bridge iff guard failed at {p}")
        prime_square_bridges.append({
            "p": p, "q": q, "q_is_prime": is_prime(q),
            "slope": slope, "source_pq": p * q,
            "actual_aliquot": actual, "defect": defect,
        })
    require(aliquot(39) == 17 and aliquot(155) == 37
            and aliquot(399) == 241 and aliquot(399) != 65,
            "prime-square bridge examples or hostile control failed")
    for a in range(21):
        value_3 = sigma_factorized(2 ** a * 3) - 2 ** a * 3
        value_7 = sigma_factorized(2 ** a * 7) - 2 ** a * 7
        require(value_7 - value_3 == 4 * (2 ** a - 1),
                f"modular-trace counterexample failed at {a}")
        require(value_7 % 4 == value_3 % 4,
                f"modular-trace collision failed at {a}")
        value_27 = sigma_factorized(2 ** a * 27) - 2 ** a * 27
        value_35 = sigma_factorized(2 ** a * 35) - 2 ** a * 35
        require(value_35 - value_27 == 8 * (2 ** a - 1),
                f"fiber-13 modular-trace formula failed at {a}")
        require(value_35 % 8 == value_27 % 8,
                f"fiber-13 modulo-8 collision failed at {a}")
    try:
        reconstruct_from_doubling_pair(1, 4)
    except ValueError:
        pass
    else:
        raise RuntimeError("fake pair (1,4) accepted without sigma guard")

    return {
        "status": "PASS",
        "scope": "self-contained proofs plus exact finite integer checks",
        "three_independent_aliquot_methods_checked_through": limit,
        "source_bound_checks_through": limit,
        "doubling_checks": doubling_checks,
        "two_observation_reconstruction_checks": doubling_checks,
        "pair_image_guard_rectangle": {"x": [0, 50], "y": [1, 150]},
        "pair_image_guard_checks": pair_guard_checks,
        "pair_image_guard_accepted": pair_guard_accepted,
        "pair_image_guard_reference_source_bound": image_limit,
        "two_observation_injectivity_checks_through": image_limit,
        "prime_dilation_checks": prime_dilation_reports,
        "modulo_4_trace_collision_checks": 21,
        "modulo_8_fiber_13_trace_collision_checks": 21,
        "tower_formula_checks": tower_checks,
        "complete_inverse_fibers": fibers,
        "doubling_towers": {
            "odd_part_5": {"carry": 6, "inputs": [5, 10, 20, 40],
                           "outputs": [aliquot(n) for n in [5, 10, 20, 40]]},
            "odd_part_25": {"carry": 31, "inputs": [25, 50, 100],
                            "outputs": [aliquot(n) for n in [25, 50, 100]]},
        },
        "forward_chains": {str(seed): values
                           for seed, values in expected_chains.items()},
        "scalar_collision": {
            "sources": [27, 35], "common_aliquot": 13,
            "carries": [40, 48], "doubled_aliquots": [66, 74],
        },
        "37_semiprime_factor_pairs": semiprime_pairs,
        "guarded_prime_square_carry_bridges": prime_square_bridges,
        "hostile_controls_passed": [
            "bare doubling s(2n)=2s(n) rejected",
            "scalar-only deterministic doubling update rejected",
            "prime-square exception retained",
            "nonpositive aliquot inputs rejected",
            "targets 0 and 1 excluded from finite-fiber routine",
            "fake observation pair (1,4) rejected by sigma guard",
            "modulo 4 traces do not determine the source",
            "modulo 8 traces cannot separate the complete 13-fiber",
            "composite carry q=57 rejects the prime-only bridge at p=7",
            "composite dilation p=4 rejected; formal inverse gives wrong n=11",
        ],
        "nonclaims": [
            "no definition of s(0)",
            "no global classification of untouchable numbers",
            "no general termination claim for aliquot sequences",
            "no fourth-power representation conclusion from numeric coincidence",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    payload = json.dumps(main(), indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
