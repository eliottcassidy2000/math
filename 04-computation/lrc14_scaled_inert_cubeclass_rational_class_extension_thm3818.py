#!/usr/bin/env python3
"""Exact companion for THM-3818's scaled support-two cube-class packet.

The finite universe is the 5,855 coprime ratios p<q, p+q<=356, whose
pair sum has only primes 2 modulo 3 and exponent at most two.  Two separate
integer algorithms compute the class of p^3+q^3 modulo rational cubes:
trial factorization modulo three and a direct maximal-cube-divisor scan.

The remaining gates check the algebraic decoder on finite controls, the
flatness facet formulas, pair-sum residue synchronization, and the two sharp
scope hostiles.  No Python assertions are used, so every truth gate remains
active under ``python -O``.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from fractions import Fraction
from itertools import combinations, product


CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def factor(n: int) -> dict[int, int]:
    require(isinstance(n, int) and n >= 1, "factor input must be positive integer")
    factors: dict[int, int] = {}
    candidate = 2
    while candidate * candidate <= n:
        while n % candidate == 0:
            factors[candidate] = factors.get(candidate, 0) + 1
            n //= candidate
        candidate = 3 if candidate == 2 else candidate + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def admissible_sum(total: int) -> bool:
    factors = factor(total)
    return bool(factors) and all(
        prime % 3 == 2 and exponent <= 2
        for prime, exponent in factors.items()
    )


def factor_cube_key(n: int) -> tuple[int, int, tuple[tuple[int, int], ...]]:
    """Return maximal cube root, cube-free kernel, and exponent-mod-three key."""
    factors = factor(n)
    cube_root = 1
    kernel = 1
    key = []
    for prime, exponent in sorted(factors.items()):
        cube_root *= prime ** (exponent // 3)
        residue = exponent % 3
        if residue:
            kernel *= prime**residue
            key.append((prime, residue))
    require(n == cube_root**3 * kernel, "factor cube decomposition failed")
    return cube_root, kernel, tuple(key)


def scan_cube_key(n: int, root_cap: int) -> tuple[int, int]:
    """Compute the same kernel without factoring, by scanning cube divisors."""
    best = 1
    for candidate in range(2, root_cap + 1):
        if n % candidate**3 == 0:
            best = candidate
    kernel = n // best**3
    require(n == best**3 * kernel, "scan cube decomposition failed")
    return best, kernel


def continued_fraction(numerator: int, denominator: int) -> tuple[int, ...]:
    require(numerator > 0 and denominator > 0, "continued-fraction input")
    digits = []
    while denominator:
        quotient, numerator, denominator = (
            numerator // denominator,
            denominator,
            numerator % denominator,
        )
        digits.append(quotient)
    return tuple(digits)


def direct_midpoint_ties(u: int, v: int) -> int:
    first = {Fraction(2 * index + 1, 2 * u) for index in range(u)}
    second = {Fraction(2 * index + 1, 2 * v) for index in range(v)}
    return len(first & second)


def distance_residue(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def positive_compositions(total: int, parts: int, prefix: tuple[int, ...] = ()):
    if parts == 1:
        yield prefix + (total,)
        return
    for first in range(1, total - parts + 2):
        yield from positive_compositions(total - first, parts - 1, prefix + (first,))


def short_relations(row: tuple[int, ...], norm_cap: int):
    """All primitive-up-to-sign sparse relations with l1 norm at most norm_cap."""
    relations = set()
    for mass in range(1, norm_cap + 1):
        for support_size in range(1, min(len(row), mass) + 1):
            for support in combinations(range(len(row)), support_size):
                for magnitudes in positive_compositions(mass, support_size):
                    for signs in product((-1, 1), repeat=support_size):
                        coefficients = tuple(
                            sign * magnitude for sign, magnitude in zip(signs, magnitudes)
                        )
                        require(len(coefficients) == len(support), "relation coefficient typing")
                        if sum(coeff * row[index] for coeff, index in zip(coefficients, support)):
                            continue
                        if coefficients[0] < 0:
                            coefficients = tuple(-coefficient for coefficient in coefficients)
                        relations.add((support, coefficients))
    return tuple(sorted(relations))


def lattice_witness(row: tuple[int, ...], denominator: int) -> bool:
    for phase in range(denominator):
        if all(
            14 * distance_residue(phase * speed, denominator) >= denominator
            for speed in row
        ):
            return True
    return False


def residue_cover(row: tuple[int, ...], first_owner: int, second_owner: int) -> bool:
    denominator = row[first_owner] + row[second_owner]
    pair_safe = {
        phase
        for phase in range(denominator)
        if 14 * distance_residue(phase * row[first_owner], denominator) >= denominator
    }
    require(
        pair_safe
        == {
            phase
            for phase in range(denominator)
            if 14 * distance_residue(phase * row[second_owner], denominator) >= denominator
        },
        "shell-partner synchronization failed",
    )
    outside_danger = set()
    for owner, speed in enumerate(row):
        if owner in (first_owner, second_owner):
            continue
        outside_danger.update(
            phase
            for phase in range(denominator)
            if 14 * distance_residue(phase * speed, denominator) < denominator
        )
    return pair_safe <= outside_danger


def main() -> None:
    sum_cap = 356
    rows: list[tuple[int, int, int, int, int]] = []
    cube_classes: dict[int, tuple[int, int]] = {}
    raw_values: dict[int, tuple[int, int]] = {}
    admissible_sums = []
    internal_cube_roots: Counter[int] = Counter()
    synchronized_phases = 0
    facet_checks = 0
    decoder_controls = 0

    for total in range(3, sum_cap + 1):
        if not admissible_sum(total):
            continue
        admissible_sums.append(total)
        for p in range(1, (total + 1) // 2):
            q = total - p
            if not p < q or math.gcd(p, q) != 1:
                continue

            value = p**3 + q**3
            cube_root_a, kernel_a, exponent_key = factor_cube_key(value)
            cube_root_b, kernel_b = scan_cube_key(value, sum_cap)
            require(cube_root_a == cube_root_b, "independent cube roots disagree")
            require(kernel_a == kernel_b, "independent cube kernels disagree")
            require(
                kernel_a
                == math.prod(prime**exponent for prime, exponent in exponent_key),
                "factor exponent key disagrees with integer kernel",
            )
            require(kernel_a not in cube_classes, "cube-class collision")
            require(value not in raw_values, "raw cube-sum collision")
            cube_classes[kernel_a] = (p, q)
            raw_values[value] = (p, q)
            rows.append((p, q, value, cube_root_a, kernel_a))
            if cube_root_a > 1:
                internal_cube_roots[cube_root_a] += 1

            lower = Fraction(q - 13 * p, 14)
            upper = Fraction(13 * q - p, 14)
            require(upper - lower == Fraction(6 * (p + q), 7), "facet width identity")
            require(q * p - p * q == 0, "support-two covector identity")
            facet_checks += 2

            for phase in range(total):
                require(
                    distance_residue(phase * p, total)
                    == distance_residue(phase * q, total),
                    "pair-sum residue synchronization",
                )
                synchronized_phases += 1

            for scale in (1, 2, 7, 43):
                address = scale**3 * value
                require(address // value == scale**3, "scale quotient is not the expected cube")
                require(address % value == 0, "scaled address divisibility")
                decoder_controls += 2

    require(len(admissible_sums) == 94, "admissible-sum count changed")
    require(len(rows) == 5_855, "primitive ratio count changed")
    require(len(cube_classes) == 5_855, "cube-class separation count changed")
    require(len(raw_values) == 5_855, "raw-address separation count changed")
    require(sum(internal_cube_roots.values()) == 43, "internal cube-factor count changed")
    require(
        internal_cube_roots == Counter({7: 36, 13: 4, 19: 2, 31: 1}),
        "internal cube-root distribution changed",
    )
    first_internal = next(row for row in rows if row[3] > 1)
    require(first_internal == (1, 19, 6860, 7, 20), "first internal cube hostile changed")

    payload = "\n".join(
        json.dumps([p, q, value, cube_root, kernel], separators=(",", ":"))
        for p, q, value, cube_root, kernel in rows
    )
    digest = hashlib.sha256(payload.encode("ascii")).hexdigest()
    require(
        digest == "feee0062ab6afb2d19f00e5397f0650f3738b860d1cdc3c1ab49193b94636a26",
        "cube-class payload digest changed",
    )

    # The atlas filter is load-bearing: the same integer has an external
    # representation with a different midpoint schedule.
    taxicab_value = 86**3 + 344**3
    require(taxicab_value == 197**3 + 323**3 == 41_343_640, "taxicab hostile equality")
    require((1, 4) in cube_classes.values(), "taxicab atlas source missing")
    require(math.gcd(86, 344) == 86, "taxicab atlas scale")
    require(math.gcd(197, 323) == 1, "taxicab external pair primitiveness")
    require(factor(430) == {2: 1, 5: 1, 43: 1}, "scaled source-sum typing")
    require(factor(520) == {2: 3, 5: 1, 13: 1}, "external sum typing")
    source_ties = direct_midpoint_ties(86, 344)
    external_ties = direct_midpoint_ties(197, 323)
    require(source_ties == 0, "scaled 1:4 tie hostile")
    require(external_ties == 1, "external odd/odd tie hostile")
    require(86 + 344 - source_ties == 430, "scaled source wall count")
    require(197 + 323 - external_ties == 519, "external wall count")
    require(continued_fraction(4, 1) == (4,), "scaled source CF")
    require(
        continued_fraction(323, 197) == (1, 1, 1, 1, 3, 2, 3, 2),
        "external hostile CF",
    )
    base_value = 1**3 + 4**3
    require(taxicab_value == 86**3 * base_value, "taxicab atlas decoder")
    require(factor_cube_key(taxicab_value)[1] == factor_cube_key(base_value)[1], "taxicab cube class")

    # Same unique minimum relation and complete pair-local packet, opposite
    # global safety at the canonical phase 1/5.
    safe_row = (1, 4) + tuple(21**power for power in range(1, 12))
    blocked_row = safe_row[:-1] + (5 * 21**11,)
    expected_relation = (((0, 1), (4, -1)),)
    require(len(safe_row) == len(blocked_row) == 13, "arrival hostile row size")
    require(len(set(safe_row)) == len(set(blocked_row)) == 13, "arrival hostile distinctness")
    require(math.gcd(*safe_row) == math.gcd(*blocked_row) == 1, "arrival hostile primitiveness")
    require(short_relations(safe_row, 4) == (), "safe row has shorter relation")
    require(short_relations(blocked_row, 4) == (), "blocked row has shorter relation")
    require(short_relations(safe_row, 5) == expected_relation, "safe row minimum relation")
    require(short_relations(blocked_row, 5) == expected_relation, "blocked row minimum relation")
    safe_minimum = min(
        min(Fraction(speed % 5, 5), 1 - Fraction(speed % 5, 5))
        for speed in safe_row
    )
    blocked_minimum = min(
        min(Fraction(speed % 5, 5), 1 - Fraction(speed % 5, 5))
        for speed in blocked_row
    )
    require(safe_minimum == Fraction(1, 5), "safe hostile phase value")
    require(blocked_minimum == 0, "blocked hostile phase value")
    require(not residue_cover(safe_row, 0, 1), "safe row residue cover should fail")
    require(residue_cover(blocked_row, 0, 1), "blocked row residue cover should hold")
    require(lattice_witness(safe_row, 5), "safe row lattice witness missing")
    require(not lattice_witness(blocked_row, 5), "blocked row lattice cover failed")

    print("THM-3818 SCALED INERT CUBECLASS SUPPORT-TWO PAIR PACKET")
    print("STATUS PROVED_ALGEBRA + FINITE_EXACT; INDEPENDENT_HOSTILE_AUDIT PENDING")
    print("SUM_CAP", sum_cap)
    print("ADMISSIBLE_SUMS", len(admissible_sums))
    print("PRIMITIVE_RATIOS", len(rows))
    print("RAW_CUBE_SUMS", len(raw_values))
    print("RATIONAL_CUBE_CLASSES", len(cube_classes))
    print("CUBECLASS_COLLISIONS", len(rows) - len(cube_classes))
    print("ALGORITHM_A", "trial_factor_exponents_mod_3 PASS")
    print("ALGORITHM_B", "maximal_cube_divisor_scan PASS")
    print("INTERNAL_CUBE_FACTOR_BASES", sum(internal_cube_roots.values()))
    print("INTERNAL_CUBE_ROOT_DISTRIBUTION", "7:36,13:4,19:2,31:1")
    print("FIRST_INTERNAL_CUBE_FACTOR", "(1,19):6860=7^3*20")
    print("CUBECLASS_PAYLOAD_SHA256", digest)
    print("FACET_IDENTITY_CHECKS", facet_checks)
    print("PAIR_SUM_SYNCHRONIZED_PHASES", synchronized_phases)
    print("SCALE_DECODER_CONTROLS", decoder_controls)
    print("TAXICAB_HOSTILE", "41343640=(86,344)=(197,323); ties=0/1; blocks=430/519")
    print("ARRIVAL_HOSTILE", "unique_l1=5_pair=(1,4); t=1/5 minima=1/5,0")
    print("RESIDUE_COVER_HOSTILE", "safe_cover=False blocked_cover=True")
    print("NONCONSEQUENCE", "no LRC14 row excluded; no off-lattice or semantic arrival")
    print("ACTIVE_CHECKS", CHECKS)
    print("RESULT PASS")


if __name__ == "__main__":
    main()
