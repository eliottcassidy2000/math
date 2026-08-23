#!/usr/bin/env python3
"""Bounded controls for THM-3730's quantitative inert-prime mass.

All theorem inequalities are proved symbolically in THM-3730.
Floating evaluations here are hostile controls, not proof dependencies.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import fsum, gamma, isqrt, log


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def primes_up_to(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0:2] = b"\x00\x00"
    for value in range(2, isqrt(limit) + 1):
        if flags[value]:
            flags[value * value:limit + 1:value] = b"\x00" * (
                (limit - value * value) // value + 1
            )
    return tuple(index for index, flag in enumerate(flags) if flag)


def normalized_mass(prime):
    terms = (
        (left**3 + (prime - left)**3) ** (-2.0 / 3.0)
        for left in range(1, (prime - 1) // 2 + 1)
    )
    return prime * fsum(terms)


def taylor_lower(prime):
    numerator = (prime - 1) * (
        201 * prime**5 + 75 * prime**4 + 33 * prime**3
        + 33 * prime**2 + 40 * prime + 40
    )
    return Fraction(numerator, 252 * prime**6)


def stable_digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def main():
    prime_limit = 20_000
    collision_limit = 5_000
    primes = primes_up_to(prime_limit)
    inert = tuple(p for p in primes if p >= 5 and p % 3 == 2)

    kappa = gamma(4.0 / 3.0) ** 2 / (2.0 * gamma(5.0 / 3.0))
    integral = 2.0 * kappa
    endpoint = 2.0 ** (1.0 / 3.0)
    uniform = Fraction(6478, 9375)

    require(not tuple(range(1, (2 - 1) // 2 + 1)), "p=2 must be empty")
    require(taylor_lower(5) == uniform, "uniform endpoint identity")

    minimum_mass = (float("inf"), None)
    minimum_taylor_margin = (float("inf"), None)
    maximum_scaled_error = (0.0, None, None)
    sample_rows = []
    for prime in inert:
        mass = normalized_mass(prime)
        lower = taylor_lower(prime)
        require(lower >= uniform, ("Taylor monotonicity", prime, lower))
        require(mass + 2e-13 >= float(lower),
                ("Taylor mass lower", prime, mass, lower))
        riemann_lower = integral - endpoint / prime
        riemann_upper = integral + (2.0 ** (4.0 / 3.0) - 1.0) / prime
        require(mass + 2e-13 >= riemann_lower,
                ("Riemann lower", prime, mass, riemann_lower))
        require(mass <= riemann_upper + 2e-13,
                ("Riemann upper", prime, mass, riemann_upper))

        if mass < minimum_mass[0]:
            minimum_mass = (mass, prime)
        margin = mass - float(lower)
        if margin < minimum_taylor_margin[0]:
            minimum_taylor_margin = (margin, prime)
        scaled_error = prime * (mass - integral)
        if abs(scaled_error) > maximum_scaled_error[0]:
            maximum_scaled_error = (abs(scaled_error), prime, scaled_error)
        if prime in (5, 11, 71, 101, 10007, 19997):
            sample_rows.append((
                prime,
                format(mass, ".15f"),
                format(float(lower), ".15f"),
                format(scaled_error, ".15f"),
            ))

    constructed = {}
    collision_primes = tuple(
        p for p in primes if 5 <= p <= collision_limit and p % 3 == 2
    )
    for prime in collision_primes:
        row = []
        for left in range(1, (prime - 1) // 2 + 1):
            right = prime - left
            value = left**3 + right**3
            row.append(value)
            previous = constructed.setdefault(value, (prime, left, right))
            require(previous == (prime, left, right),
                    ("cross-prime collision", value, previous,
                     (prime, left, right)))
        require(all(a > b for a, b in zip(row, row[1:])),
                ("within-prime monotonicity", prime))

    prime_mass_rows = []
    for cutoff in (100, 1000, 5000, 10_000, 20_000):
        partial = fsum(
            normalized_mass(p) / p
            for p in inert if p <= cutoff
        )
        prime_mass_rows.append((
            cutoff,
            format(partial, ".15f"),
            format(partial / log(log(cutoff)), ".15f"),
        ))

    semantic = (
        prime_limit,
        len(inert),
        tuple(sample_rows),
        (format(minimum_mass[0], ".15f"), minimum_mass[1]),
        (format(minimum_taylor_margin[0], ".15f"),
         minimum_taylor_margin[1]),
        (format(maximum_scaled_error[0], ".15f"),
         maximum_scaled_error[1],
         format(maximum_scaled_error[2], ".15f")),
        collision_limit,
        len(constructed),
        tuple(prime_mass_rows),
    )
    print("quantitative two-cube inert-prime critical-mass probe")
    print("prime_limit=" + repr(prime_limit))
    print("inert_prime_count=" + repr(len(inert)))
    print("kappa=" + format(kappa, ".15f"))
    print("integral_2kappa=" + format(integral, ".15f"))
    print("uniform_constant=" + repr(uniform))
    print("sample_rows=" + repr(tuple(sample_rows)))
    print("minimum_normalized_mass=" + repr((
        format(minimum_mass[0], ".15f"), minimum_mass[1]
    )))
    print("minimum_taylor_margin=" + repr((
        format(minimum_taylor_margin[0], ".15f"),
        minimum_taylor_margin[1],
    )))
    print("maximum_abs_scaled_riemann_error=" + repr((
        format(maximum_scaled_error[0], ".15f"),
        maximum_scaled_error[1],
        format(maximum_scaled_error[2], ".15f"),
    )))
    print("cross_prime_constructed_values=" + repr(len(constructed)))
    print("prime_mass_rows=" + repr(tuple(prime_mass_rows)))
    print("semantic_sha256=" + stable_digest(semantic))
    print("scope=finite_controls_only;all_scale_proof_is_in_THM-3730")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
