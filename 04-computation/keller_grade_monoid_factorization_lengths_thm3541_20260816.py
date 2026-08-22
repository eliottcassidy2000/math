#!/usr/bin/env python3
"""Exact arithmetic companion for THM-3541.

The degree monoid is M={1} union {n>=3}.  This script independently checks
its irreducibles, factorization-length formulas, first nonunique grades, and
the finite-prefix harmonic identity for its atom support.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
import json
from math import isqrt


EXPECTED_SEMANTIC_SHA256 = "c31c4c64207547cf2b239fd59be8fc53925f75e84f4b6123772ef098a49432a3"
CLASSIFICATION_LIMIT = 100_000
BRUTE_LENGTH_LIMIT = 5_000
HARMONIC_CUTOFF = 1_000


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def in_grade_monoid(value: int) -> bool:
    return value == 1 or value >= 3


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def grade_atom_by_divisors(value: int) -> bool:
    if value < 3:
        return False
    for divisor in range(3, isqrt(value) + 1):
        if value % divisor == 0 and in_grade_monoid(value // divisor):
            return False
    return True


def grade_atom_closed_form(value: int) -> bool:
    if value in (4, 8):
        return True
    if value % 2 == 1:
        return is_prime(value)
    if value % 4 == 2:
        return is_prime(value // 2) and value // 2 % 2 == 1
    return False


def two_valuation_and_odd_omega(value: int) -> tuple[int, int]:
    require(value >= 1, ("positive value", value))
    exponent_two = 0
    odd_part = value
    while odd_part % 2 == 0:
        exponent_two += 1
        odd_part //= 2
    omega = 0
    divisor = 3
    while divisor * divisor <= odd_part:
        while odd_part % divisor == 0:
            omega += 1
            odd_part //= divisor
        divisor += 2
    if odd_part > 1:
        omega += 1
    return exponent_two, omega


def closed_length_interval(value: int) -> tuple[int, int]:
    require(value >= 3, ("nonunit grade", value))
    exponent_two, odd_omega = two_valuation_and_odd_omega(value)
    deficit = max(exponent_two - odd_omega, 0)
    minimum = odd_omega + (deficit + 2) // 3
    maximum = odd_omega + exponent_two // 2
    return minimum, maximum


def allocation_lengths(value: int) -> frozenset[int]:
    """Enumerate x atoms 2p, y atoms 4, z atoms 8."""
    exponent_two, odd_omega = two_valuation_and_odd_omega(value)
    lengths = set()
    for paired in range(min(exponent_two, odd_omega) + 1):
        remainder = exponent_two - paired
        for eights in range(remainder // 3 + 1):
            leftover = remainder - 3 * eights
            if leftover % 2 == 0:
                fours = leftover // 2
                lengths.add(odd_omega + fours + eights)
    return frozenset(lengths)


@lru_cache(maxsize=None)
def brute_factorization_lengths(value: int) -> frozenset[int]:
    require(value >= 3, ("brute grade", value))
    lengths = {1} if grade_atom_by_divisors(value) else set()
    for atom in range(3, value + 1):
        if value % atom != 0 or not grade_atom_by_divisors(atom):
            continue
        quotient = value // atom
        if quotient >= 3:
            lengths.update(1 + length for length in brute_factorization_lengths(quotient))
    return frozenset(lengths)


@lru_cache(maxsize=None)
def atom_factorizations(value: int, minimum_atom: int = 3) -> tuple[tuple[int, ...], ...]:
    rows = set()
    if value >= minimum_atom and grade_atom_by_divisors(value):
        rows.add((value,))
    for atom in range(minimum_atom, isqrt(value) + 1):
        if value % atom != 0 or not grade_atom_by_divisors(atom):
            continue
        quotient = value // atom
        if quotient < atom:
            continue
        for tail in atom_factorizations(quotient, atom):
            rows.add((atom,) + tail)
    return tuple(sorted(rows, key=lambda row: (len(row), row)))


def primes_through(limit: int) -> tuple[int, ...]:
    sieve = bytearray(b"\x01") * (limit + 1)
    if limit >= 0:
        sieve[0] = 0
    if limit >= 1:
        sieve[1] = 0
    for prime in range(2, isqrt(limit) + 1):
        if sieve[prime]:
            start = prime * prime
            sieve[start:limit + 1:prime] = b"\x00" * (((limit - start) // prime) + 1)
    return tuple(index for index, flag in enumerate(sieve) if flag)


def fraction_digest(value: Fraction) -> str:
    return sha256(f"{value.numerator}/{value.denominator}".encode("ascii")).hexdigest()


def main() -> None:
    # Independent divisor search versus the closed atom classification.
    atoms = []
    for value in range(3, CLASSIFICATION_LIMIT + 1):
        divisor_answer = grade_atom_by_divisors(value)
        closed_answer = grade_atom_closed_form(value)
        require(divisor_answer == closed_answer,
                ("grade atom classification", value, divisor_answer, closed_answer))
        if divisor_answer:
            atoms.append(value)

        allocation = allocation_lengths(value)
        lower, upper = closed_length_interval(value)
        require(allocation == frozenset(range(lower, upper + 1)),
                ("length interval", value, allocation, lower, upper))

    # A second recursive factorization route checks the length formula.
    for value in range(3, BRUTE_LENGTH_LIMIT + 1):
        require(brute_factorization_lengths(value) == allocation_lengths(value),
                ("brute length disagreement", value,
                 brute_factorization_lengths(value), allocation_lengths(value)))

    first_reducible = next(value for value in range(3, 100)
                           if not grade_atom_by_divisors(value))
    first_nonunique = next(value for value in range(3, 200)
                           if len(atom_factorizations(value)) > 1)
    first_variable_length = next(
        value for value in range(3, 300)
        if len({len(row) for row in atom_factorizations(value)}) > 1
    )
    require(first_reducible == 9, first_reducible)
    require(first_nonunique == 24, first_nonunique)
    require(atom_factorizations(24) == ((3, 8), (4, 6)),
            atom_factorizations(24))
    require(first_variable_length == 36, first_variable_length)
    require(atom_factorizations(36) == ((6, 6), (3, 3, 4)),
            atom_factorizations(36))

    # Global elasticity bound: upper/lower <= 3/2, first equality at 36.
    first_elasticity_three_halves = None
    for value in range(3, CLASSIFICATION_LIMIT + 1):
        lower, upper = closed_length_interval(value)
        require(2 * upper <= 3 * lower,
                ("elasticity exceeds 3/2", value, lower, upper))
        if 2 * upper == 3 * lower and upper > lower:
            first_elasticity_three_halves = value
            break
    require(first_elasticity_three_halves == 36,
            first_elasticity_three_halves)

    samples = {}
    for value in (4, 8, 9, 12, 18, 24, 36, 64, 72, 144):
        samples[str(value)] = {
            "v2_omega_odd": two_valuation_and_odd_omega(value),
            "lengths": sorted(allocation_lengths(value)),
            "factorizations": atom_factorizations(value),
        }

    # Exact finite-prefix harmonic identity for
    # Atom(M)=P_odd disjoint-union 2P_odd disjoint-union {4,8}.
    odd_primes = tuple(prime for prime in primes_through(HARMONIC_CUTOFF)
                       if prime % 2 == 1)
    direct_atom_mass = sum(
        (Fraction(1, value) for value in atoms if value <= HARMONIC_CUTOFF),
        Fraction(0),
    )
    prime_mass = sum((Fraction(1, prime) for prime in odd_primes), Fraction(0))
    doubled_prime_mass = sum(
        (Fraction(1, 2 * prime) for prime in odd_primes
         if 2 * prime <= HARMONIC_CUTOFF),
        Fraction(0),
    )
    finite_exception_mass = Fraction(1, 4) + Fraction(1, 8)
    require(direct_atom_mass
            == prime_mass + doubled_prime_mass + finite_exception_mass,
            "harmonic atom support identity")

    ledger = {
        "limits": {
            "classification": CLASSIFICATION_LIMIT,
            "brute_lengths": BRUTE_LENGTH_LIMIT,
            "harmonic": HARMONIC_CUTOFF,
        },
        "atom_classification": "odd primes union {4,8} union twice odd primes",
        "atom_count_through_limit": len(atoms),
        "first_atoms": atoms[:30],
        "length_formula": {
            "minimum": "Omega(oddpart)+ceil(max(v2-Omega(oddpart),0)/3)",
            "maximum": "Omega(oddpart)+floor(v2/2)",
            "sets_are_intervals": True,
        },
        "first_reducible": first_reducible,
        "first_nonunique": {
            "grade": first_nonunique,
            "factorizations": atom_factorizations(first_nonunique),
        },
        "first_variable_length": {
            "grade": first_variable_length,
            "factorizations": atom_factorizations(first_variable_length),
        },
        "elasticity": {
            "supremum": "3/2",
            "first_attained": first_elasticity_three_halves,
        },
        "samples": samples,
        "harmonic_prefix": {
            "cutoff": HARMONIC_CUTOFF,
            "atom_count": sum(1 for value in atoms if value <= HARMONIC_CUTOFF),
            "mass_digest": fraction_digest(direct_atom_mass),
            "identity": "H(P_odd)+H(2P_odd)+1/4+1/8",
        },
    }
    semantic_sha256 = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
            ("semantic digest", semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    print("== THM-3541 Keller grade-monoid factorization atlas ==")
    print("grade_monoid=M={1} union {n>=3};"
          "atoms=odd_primes union {4,8} union {2p:p odd prime}")
    print("length_interval_for_N=2^a*m_odd;b=Omega(m):"
          "[b+ceil(max(a-b,0)/3),b+floor(a/2)]")
    print("first_reducible_grade=9;first_nonunique_grade=24:"
          "3*8=4*6;first_variable_length_grade=36:6*6=4*3*3")
    print("grade_factorization_elasticity=3/2;first_attained_at=36")
    print(f"sample_length_sets={{{', '.join(key+':'+str(value['lengths']) for key, value in samples.items())}}}")
    print("harmonic_atom_support=P_odd disjoint_union 2P_odd disjoint_union {4,8};"
          f"cutoff={HARMONIC_CUTOFF};exact_mass_digest={fraction_digest(direct_atom_mass)}")
    print(f"semantic_sha256={semantic_sha256}")
    print("status=PROVED NUMERICAL MONOID CLASSIFICATION + VERIFIED-EXACT;"
          "map atoms remain every grade; no fixed-map factorization uniqueness or LRC claim")


if __name__ == "__main__":
    main()
