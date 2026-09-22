"""Exact standard-library certificates for append-digit prime patterns.

Universal mechanisms are proved in the companion note. Finite universes
and complete primality checks are explicit. No probable-prime tests used.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from math import gcd, isqrt, lcm, prod
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factor_trial(n: int) -> list[list[int]]:
    require(n >= 1, "positive integer required")
    result = []
    divisor = 2
    while divisor * divisor <= n:
        exponent = 0
        while n % divisor == 0:
            n //= divisor
            exponent += 1
        if exponent:
            result.append([divisor, exponent])
        divisor += 1
    if n > 1:
        result.append([n, 1])
    return result


def prime_trial(n: int) -> bool:
    return n >= 2 and all(n % d != 0 for d in range(2, isqrt(n) + 1))


def verify_factorization(n: int, factors: list[list[int]]) -> None:
    product = 1
    for prime, exponent in factors:
        require(prime_trial(prime), "factor not independently trial-prime")
        require(exponent >= 1, "positive exponent required")
        product *= prime**exponent
    require(product == n, "factorization product mismatch")


def phi(n: int) -> int:
    result = n
    for prime, _ in factor_trial(n):
        result = result // prime * (prime - 1)
    return result


def order(base: int, modulus: int) -> int:
    require(gcd(base, modulus) == 1, "unit required")
    candidate = phi(modulus)
    for prime, _ in factor_trial(candidate):
        while candidate % prime == 0 and pow(base, candidate // prime, modulus) == 1:
            candidate //= prime
    require(pow(base, candidate, modulus) == 1, "order return failed")
    for prime, _ in factor_trial(candidate):
        require(pow(base, candidate // prime, modulus) != 1, "order minimality failed")
    return candidate


def decimal_term(k: int) -> int:
    return (10 ** (k + 1) - 7) // 3


def append_term(base: int, digit: int, k: int) -> int:
    return 1 + digit * base * ((base**k - 1) // (base - 1))


def first_composites() -> dict:
    repeated_rows = []
    for k in range(1, 9):
        value = decimal_term(k)
        require(value == int("3" * k + "1"), "decimal digit encoding failure")
        factors = factor_trial(value)
        verify_factorization(value, factors)
        is_prime = factors == [[value, 1]]
        require(is_prime == (k <= 7), "first decimal composite boundary failed")
        repeated_rows.append({"k": k, "digits": k + 1, "value": value,
                              "factorization": factors, "prime": is_prime})
    require(repeated_rows[-1]["factorization"] == [[17, 1], [19607843, 1]],
            "first decimal composite factors failed")
    fermat_rows = []
    for n in range(6):
        value = 2 ** (2**n) + 1
        factors = factor_trial(value)
        verify_factorization(value, factors)
        is_prime = factors == [[value, 1]]
        require(is_prime == (n <= 4), "first Fermat composite boundary failed")
        fermat_rows.append({"n": n, "value": value, "factorization": factors,
                            "prime": is_prime})
    require(fermat_rows[-1]["factorization"] == [[641, 1], [6700417, 1]],
            "Fermat F5 factorization failed")
    return {"decimal_k_1_through_8": repeated_rows, "Fermat_n_0_through_5": fermat_rows,
            "primality_method": "complete trial division; factors independently checked through integer square root"}


def general_gcd_controls() -> dict:
    checks = 0
    bases_digits = 0
    for base in range(2, 17):
        for digit in range(1, base):
            bases_digits += 1
            c = (digit - 1) * base + 1
            terms = [append_term(base, digit, k) for k in range(26)]
            for k in range(25):
                require(terms[k + 1] == base * terms[k] + c, "append recurrence failed")
                require(gcd(terms[k], base * c) == 1, "append unit invariant failed")
                require(terms[k] == 1 + sum(digit * base**i for i in range(1, k + 1)),
                        "independent positional encoding failed")
            for m in range(13):
                for gap in range(1, 13):
                    repunit = (base**gap - 1) // (base - 1)
                    require(gcd(terms[m], terms[m + gap]) == gcd(terms[m], repunit),
                            "general gcd law failed")
                    checks += 1
    return {"base_range": [2, 16], "digits": "1<=d<b", "m_range": [0, 12],
            "gap_range": [1, 12], "base_digit_pairs": bases_digits, "gcd_checks": checks}


def window_certificate() -> dict:
    rows = []
    candidates = set()
    expected = [3, 3, 9, 3, 3, 9, 3, 3, 999, 33, 3, 117, 3, 3]
    for gap in range(1, 15):
        common = gcd(10**gap - 1, 7**gap - 1)
        require(common == expected[gap - 1], "window gcd table failure")
        factors = factor_trial(common)
        verify_factorization(common, factors)
        candidates.update(prime for prime, _ in factors)
        rows.append({"gap": gap, "gcd": common, "factorization": factors})
    require(candidates == {3, 11, 13, 37}, "candidate prime set failure")
    groups = []
    for prime in [11, 13, 37]:
        powers = []
        value = 1
        while value not in powers:
            powers.append(value)
            value = value * 10 % prime
        require(value == 1 and 7 not in powers, "forbidden target group test failed")
        groups.append({"prime": prime, "powers_of_ten": powers})
    terms = [decimal_term(k) for k in range(116)]
    checked_pairs = 0
    for m in range(101):
        for gap in range(1, 15):
            require(gcd(terms[m], terms[m + gap]) == 1, "independent initial window check failed")
            checked_pairs += 1
    require(gcd(decimal_term(1), decimal_term(16)) == 31, "sharp window witness failed")
    require(order(10, 31) == 15 and pow(10, 2, 31) == 7, "sharp order/target failed")
    return {"difference_certificate": rows, "candidate_primes": sorted(candidates),
            "excluding_subgroups": groups, "independent_pairs_checked": checked_pairs,
            "independent_m_range": [0, 100], "sharp_indices": [1, 16], "sharp_gcd": 31}


def prime_power_clocks() -> list[dict]:
    result = []
    expected = {17: [(8, 16), (248, 272), (2152, 4624)],
                31: [(1, 15), (226, 465), (3481, 14415)]}
    for prime in [17, 31]:
        h = order(10, prime)
        require((pow(10, h, prime**2) - 1) % prime**2 != 0,
                "ordinary prime-power lift hypothesis failed")
        for exponent in [1, 2, 3]:
            modulus = prime**exponent
            period = order(10, modulus)
            require(period == h * prime ** (exponent - 1), "prime-power period failure")
            # Independent recurrence traversal records all zeros in one cycle.
            value = 1
            hits = []
            for k in range(period):
                if value == 0:
                    hits.append(k)
                value = (10 * value + 21) % modulus
            require(value == 1 and len(hits) == 1, "affine cycle root uniqueness failed")
            require((hits[0], period) == expected[prime][exponent - 1], "prime-power first hit failure")
            require(pow(10, hits[0] + 1, modulus) == 7, "independent target congruence failed")
            result.append({"prime": prime, "exponent": exponent, "modulus": modulus,
                           "first_index": hits[0], "period": period})
    return result


def strip_product(value: int, earlier_product: int) -> int:
    while True:
        common = gcd(value, earlier_product)
        if common == 1:
            return value
        value //= common


def strip_individually(value: int, earlier_values: list[int]) -> int:
    for previous in earlier_values:
        while True:
            common = gcd(value, previous)
            if common == 1:
                break
            value //= common
    return value


def primitive_part_controls() -> dict:
    earlier_product = 1
    earlier = []
    rows = []
    no_new = []
    for k in range(1, 201):
        value = decimal_term(k)
        residual = strip_product(value, earlier_product)
        require(residual == strip_individually(value, earlier), "independent primitive stripping mismatch")
        if residual == 1:
            no_new.append(k)
        rows.append({"k": k, "primitive_residual": residual, "old_prime_part": value // residual})
        earlier_product *= value
        earlier.append(value)
    require(no_new == [], "finite decimal primitive-prime census changed")
    require(all(row["old_prime_part"] == 1 for row in rows[:15]), "first fifteen supports not disjoint")
    mersenne = []
    earlier_product = 1
    exceptions = []
    for n in range(1, 41):
        value = 2**n - 1
        residual = strip_product(value, earlier_product)
        if residual == 1:
            exceptions.append(n)
        mersenne.append({"n": n, "value": value, "primitive_residual": residual})
        earlier_product *= value
    require(exceptions == [1, 6], "finite Mersenne exceptions failed")
    fermat = [2 ** (2**n) + 1 for n in range(9)]
    for n in range(1, 9):
        require(fermat[n] - 2 == prod(fermat[:n]), "Fermat product identity failed")
        for m in range(n):
            require(gcd(fermat[m], fermat[n]) == 1, "Fermat coprimality failed")
    require(order(2, 3) == 2 and order(2, 7) == 3 and order(2, 9) == 6,
            "63 order-lift control failed")
    require(order(2, 641) == order(2, 6700417) == 64, "F5 factor orders failed")
    return {"decimal_k_range": [1, 200], "decimal_rows": rows,
            "decimal_without_new_prime": no_new, "Mersenne_n_range": [1, 40],
            "Mersenne_rows": mersenne, "Mersenne_exceptions": exceptions,
            "Fermat_pairwise_coprime_n_range": [0, 8], "Fermat_values": fermat,
            "order_lift_63": {"mod3": 2, "mod7": 3, "mod9": 6},
            "F5_prime_factor_orders": {"641": 64, "6700417": 64}}


def finite_prime_pool_control() -> dict:
    primes = [17, 31]
    period = lcm(*(order(10, p) for p in primes))
    values = [decimal_term(t * period) for t in [1, 2, 3]]
    require(all(value % p == 1 for value in values for p in primes), "finite prime pool avoidance failed")
    return {"prime_pool": primes, "return_period": period, "index_multipliers": [1, 2, 3],
            "all_residues": 1, "does_not_certify_prime_terms": True}


def fermat_prime_decimal_clocks() -> dict:
    rows = []
    for prime, expected_hit in [(17, 8), (257, 226), (65537, 29252)]:
        require(prime_trial(prime) and prime % 40 == 17, "Fermat clock prime hypothesis failed")
        period = order(10, prime)
        require(period == prime - 1, "Fermat prime generator theorem failed")
        value = 1
        hits = []
        for k in range(period):
            if value == 0:
                hits.append(k)
            value = (10 * value + 21) % prime
        require(value == 1 and hits == [expected_hit], "Fermat-prime occurrence clock failed")
        require(pow(10, expected_hit + 1, prime) == 7, "Fermat target check failed")
        negative_count = sum(2 * (10 * j % prime) > prime for j in range(1, (prime + 1) // 2))
        t = (prime - 17) // 40
        require(negative_count == 10 * t + 5, "Gauss pairing count failed")
        require(pow(10, (prime - 1) // 2, prime) == prime - 1, "Gauss sign failed")
        rows.append({"prime": prime, "order_of_ten": period,
                     "first_decimal_index": expected_hit,
                     "negative_residue_count": negative_count,
                     "seven_in_generated_group": True})
    for k in range(201):
        require(decimal_term(k) % 3 == decimal_term(k) % 5 == 1,
                "denominator/nonunit hostile failed")
    return {"nonexceptional_rows": rows,
            "prime3": {"order_of_ten": 1, "divides_any_decimal_term": False,
                       "reason": "R_k=1 mod3; numerator divisibility cannot survive division by3 blindly"},
            "prime5": {"order_of_ten": None, "divides_any_decimal_term": False,
                       "reason": "ten is not a unit; R_k=1 mod5"}}


def run() -> dict:
    return {"status": "FINITE-EXACT certificates; scoped universal proofs in the companion note",
            "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
            "first_composites": first_composites(), "general_gcd": general_gcd_controls(),
            "sharp_fifteen_window": window_certificate(), "prime_power_clocks": prime_power_clocks(),
            "primitive_parts": primitive_part_controls(), "finite_prime_pool": finite_prime_pool_control(),
            "Fermat_prime_decimal_clocks": fermat_prime_decimal_clocks(),
            "checks": "PASS; explicit RuntimeError checks remain active under -O",
            "limitations": ["No probable-prime calls", "No all-index primitive-prime theorem for R_k",
                            "No prime-term infinitude theorem", "No Collatz descent consequence"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_suffix(".json"))
    args = parser.parse_args()
    output = run()
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8", newline="\n")
    print(json.dumps({"checks": output["checks"], "output": str(args.output),
                      "general_gcd_checks": output["general_gcd"]["gcd_checks"],
                      "decimal_primitive_indices_checked": 200,
                      "sharp_pairwise_coprime_window": 15}, sort_keys=True))
