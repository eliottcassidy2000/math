#!/usr/bin/env python3
"""Exact companion for THM-3706's checksum-depth no-go theorem."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
P = 13
N = 14
K = N * P

PINS = {
    "04-computation/lrc_padic_scalar_seam_total_speed_checksum_thm3679.py":
        "cb9e80f688254bf6c83927c20306e33a3e716ddd16e1e70febf577cc96baecde",
    "05-knowledge/results/lrc_padic_scalar_seam_total_speed_checksum_thm3679.out":
        "c357dc4643e251b42807d1c086b27928d2c3b1f2ffa9b1c3c9875f402d4c2c8d",
    "04-computation/lrc14_disc_bernoulli_dedekind_kps_S128.py":
        "58b430c961e4ad952a4843e32eced5175dd62dce1b87a81bfd6c7cdfbc730087",
    "05-knowledge/results/lrc14_disc_bernoulli_dedekind_kps_S128.out":
        "00b3aad54c779de8e63a1823671f31e1f6cbe2db752bf8e5bb9364dbb422e70e",
}


def require(condition, payload):
    if condition is not True:
        raise RuntimeError(payload)


def valuation(integer: int, prime: int) -> int:
    require(integer != 0, "valuation of zero")
    value = abs(integer)
    exponent = 0
    while value % prime == 0:
        value //= prime
        exponent += 1
    return exponent


def rational_valuation(value: Fraction, prime: int) -> int:
    return valuation(value.numerator, prime) - valuation(value.denominator, prime)


def sawtooth(value: Fraction) -> Fraction:
    fractional = value - value.numerator // value.denominator
    return Fraction() if fractional == 0 else fractional - Fraction(1, 2)


def dedekind(first: int, second: int) -> Fraction:
    require(gcd(first, second) == 1 and second > 0, (first, second))
    return sum(
        sawtooth(Fraction(residue, second))
        * sawtooth(Fraction(first * residue, second))
        for residue in range(1, second)
    )


def reciprocity_ray_value(level: int, quotient: int) -> Fraction:
    """s(level,level*quotient+1), simplified using reciprocity."""
    denominator = level * quotient + 1
    return Fraction(
        level * quotient * (quotient - level),
        12 * denominator,
    )


def exact_linear_lift(constant: int, slope: int, prime: int, depth: int) -> int:
    """Construct positive m with v_p(constant+slope*m)=depth."""
    require(gcd(slope, prime) == 1, (slope, prime))
    if depth == 0:
        candidate = 1
        while valuation(constant + slope * candidate, prime) != 0:
            candidate += 1
        return candidate
    modulus = prime**depth
    root = (-constant * pow(slope, -1, modulus)) % modulus
    require(root % prime != 0, (constant, slope, prime, depth, root))
    candidates = tuple(root + digit * modulus for digit in range(prime))
    exact = tuple(
        candidate
        for candidate in candidates
        if candidate > 0 and valuation(constant + slope * candidate, prime) == depth
    )
    deeper = tuple(
        candidate
        for candidate in candidates
        if candidate > 0 and valuation(constant + slope * candidate, prime) > depth
    )
    require(len(exact) == prime - 1 and len(deeper) == 1, (exact, deeper))
    return min(exact)


def typed_row(multiplier: int):
    return (1, 2, 3, 4, 5, 11, 13, 26, K * multiplier)


def strict_typed_row(multiplier: int):
    return (
        1,
        2,
        3,
        4,
        5,
        11,
        11 * P,
        12 * P**2,
        P**2 * K * multiplier,
    )


def is_primitive_distinct_positive(row):
    return min(row) > 0 and len(set(row)) == len(row) and gcd(*row) == 1


def is_cover(level: int, row) -> bool:
    return all(any(speed % modulus == 0 for speed in row) for modulus in range(2, level + 1))


def main():
    for relative, expected in PINS.items():
        actual = sha256((ROOT / relative).read_bytes()).hexdigest()
        require(actual == expected, ("parent hash drift", relative, actual))

    # Direct finite-sum controls for the reciprocity identity.  The theorem's
    # all-m proof is the standard one-step reciprocity calculation.
    reciprocity_controls = ((5, 5), (8, 7), (12, 11), (14, 13), (14, 91))
    for level, quotient in reciprocity_controls:
        require(
            dedekind(level, level * quotient + 1)
            == reciprocity_ray_value(level, quotient),
            (level, quotient),
        )

    base_dedekind = dedekind(14, 183)
    require(base_dedekind == Fraction(-91, 1098), base_dedekind)
    require(-12 * 183 * base_dedekind == 182, base_dedekind)

    typed_examples = []
    for total_depth in range(1, 7):
        inner_depth = total_depth - 1
        multiplier = exact_linear_lift(5, 14, P, inner_depth)
        row = typed_row(multiplier)
        require(is_primitive_distinct_positive(row), row)
        require(all(value % P for value in row[:6]), row)
        require(all(value % P == 0 for value in row[6:]), row)
        require(sum(row) == P * (5 + 14 * multiplier), row)
        require(valuation(sum(row), P) == total_depth, row)
        typed_examples.append((total_depth, multiplier, sum(row), total_depth + 1))

    strict_examples = []
    for total_depth in range(3, 9):
        inner_depth = total_depth - 3
        multiplier = exact_linear_lift(1, 14, P, inner_depth)
        row = strict_typed_row(multiplier)
        require(multiplier % P != 0, multiplier)
        require(is_primitive_distinct_positive(row), row)
        require(tuple(valuation(value, P) for value in row[6:]) == (1, 2, 3), row)
        require(sum(row) == P**3 * (1 + 14 * multiplier), row)
        require(valuation(sum(row), P) == total_depth, row)
        strict_examples.append((total_depth, multiplier, sum(row), total_depth + 1))

    cover_examples = []
    for total_depth in range(1, 7):
        inner_depth = total_depth - 1
        multiplier = exact_linear_lift(6, 14, P, inner_depth)
        row = tuple(range(1, 13)) + (K * multiplier,)
        require(multiplier % P != 0, multiplier)
        require(is_primitive_distinct_positive(row), row)
        require(is_cover(14, row), row)
        require(sum(row) == P * (6 + 14 * multiplier), row)
        require(valuation(sum(row), P) == total_depth, row)
        denominator = K * multiplier + 1
        fixed_value = reciprocity_ray_value(14, P * multiplier)
        binding_value = reciprocity_ray_value(14 * multiplier, P)
        require(denominator == 14 * (P * multiplier) + 1, denominator)
        require(denominator == (14 * multiplier) * P + 1, denominator)
        require(rational_valuation(fixed_value, P) == 1, fixed_value)
        require(rational_valuation(binding_value, P) == 1, binding_value)
        if multiplier in (1, 7):
            require(dedekind(14, denominator) == fixed_value, multiplier)
            require(dedekind(14 * multiplier, denominator) == binding_value, multiplier)
        cover_examples.append(
            (
                total_depth,
                multiplier,
                sum(row),
                rational_valuation(fixed_value, P),
                rational_valuation(binding_value, P),
            )
        )

    # Uniform arithmetic analogue for n=p+1.  These finite controls test the
    # algebra at several primes; the proof for every p>=5 is symbolic.
    general_controls = 0
    primes = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
    for prime in primes:
        level = prime + 1
        far_factor = prime * level
        constant = (prime - 1) // 2
        for inner_depth in range(5):
            multiplier = exact_linear_lift(constant, level, prime, inner_depth)
            row = tuple(range(1, prime)) + (far_factor * multiplier,)
            require(multiplier % prime != 0, (prime, inner_depth, multiplier))
            require(is_primitive_distinct_positive(row), row)
            require(is_cover(level, row), (level, row))
            require(
                sum(row) == prime * (constant + level * multiplier),
                (prime, multiplier, sum(row)),
            )
            require(valuation(sum(row), prime) == inner_depth + 1, (prime, multiplier))
            denominator = far_factor * multiplier + 1
            fixed_value = reciprocity_ray_value(level, prime * multiplier)
            binding_value = reciprocity_ray_value(level * multiplier, prime)
            require(denominator == level * (prime * multiplier) + 1, denominator)
            require(denominator == (level * multiplier) * prime + 1, denominator)
            require(
                rational_valuation(fixed_value, prime) == 1,
                (prime, multiplier, fixed_value),
            )
            require(
                rational_valuation(binding_value, prime) == 1,
                (prime, multiplier, binding_value),
            )
            general_controls += 1

    print(f"s(14,183)={base_dedekind}")
    print(f"-12*183*s(14,183)={-12 * 183 * base_dedekind}")
    print(f"reciprocity_direct_controls={len(reciprocity_controls)}")
    print("typed_depth_examples=" + repr(typed_examples))
    print("strict_valuation_123_examples=" + repr(strict_examples))
    print("cover14_deep_well_examples=" + repr(cover_examples))
    print("cover14_all_rows_primitive_distinct_cover=True")
    print("cover14_fixed_and_binding_dedekind_v13=1_at_every_listed_depth")
    print(f"general_prime_controls={general_controls}")
    print("general_rule=n=p+1,total_vp_unbounded,fixed_and_binding_dedekind_vp=1")
    print("PASS")


if __name__ == "__main__":
    main()
