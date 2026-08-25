#!/usr/bin/env python3
"""Exact companion for THM-4027 (Sun 2-4-6-8 modular solubility).

All binomial coefficients are evaluated as integers with ``math.comb``.  In
particular, this script never divides by a factorial modulo a prime; see the
MISTAKE-363 guardrail.
"""

from __future__ import annotations

import hashlib
import math
from fractions import Fraction


DEGREES = (2, 4, 6, 8)
SMALL_ODD_PRIMES = (
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def prime_power_factors(n: int) -> list[tuple[int, int]]:
    factors = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            exponent = 0
            while n % p == 0:
                n //= p
                exponent += 1
            factors.append((p, exponent))
        p += 1
    if n > 1:
        factors.append((n, 1))
    return factors


def floor_log_p(r: int, p: int) -> int:
    exponent = 0
    power = 1
    while power * p <= r:
        power *= p
        exponent += 1
    return exponent


def binomial_period_prime_power(r: int, p: int, exponent: int) -> int:
    """A Vandermonde-certified period modulo p**exponent."""
    return p ** (exponent + floor_log_p(r, p))


def binomial_period(r: int, modulus: int) -> int:
    period = 1
    for p, exponent in prime_power_factors(modulus):
        period = math.lcm(period, binomial_period_prime_power(r, p, exponent))
    return period


def value_set(r: int, modulus: int) -> set[int]:
    period = binomial_period(r, modulus)
    values = {math.comb(k, r) % modulus for k in range(period)}
    # Hostile boundary check for the stated period, including k<r.
    for k in range(period + r + 1):
        require(
            math.comb(k + period, r) % modulus == math.comb(k, r) % modulus,
            "binomial period failure",
        )
    return values


def sumsets(sets: list[set[int]], modulus: int) -> set[int]:
    current = {0}
    for values in sets:
        current = {(a + b) % modulus for a in current for b in values}
    return current


def regular_triangular_values(p: int) -> tuple[set[int], dict[int, int]]:
    """Values C(w,2) having 2w-1 nonzero mod the odd prime p."""
    witnesses: dict[int, int] = {}
    for w in range(p):
        if (2 * w - 1) % p != 0:
            witnesses.setdefault(math.comb(w, 2) % p, w)
    return set(witnesses), witnesses


def finite_prime_certificate(p: int) -> tuple[list[int], list[int], str]:
    require(p % 2 == 1, "odd prime required")
    regular, regular_witness = regular_triangular_values(p)
    higher_sets = [value_set(r, p) for r in (4, 6, 8)]
    periods = [binomial_period(r, p) for r in DEGREES]
    sizes = [len(regular), *(len(values) for values in higher_sets)]

    # Store one regular-triangular tuple for every residue.  This is stronger
    # than mere surjectivity and is the exact sidecar needed for Hensel lift.
    higher_witnesses = []
    for r, period in zip((4, 6, 8), periods[1:]):
        by_value: dict[int, int] = {}
        for k in range(period):
            by_value.setdefault(math.comb(k, r) % p, k)
        higher_witnesses.append(by_value)
    triple_witness: dict[int, tuple[int, int, int]] = {}
    for x_value, x in higher_witnesses[0].items():
        for y_value, y in higher_witnesses[1].items():
            for z_value, z in higher_witnesses[2].items():
                triple_witness.setdefault((x_value + y_value + z_value) % p, (x, y, z))

    witness_rows: list[tuple[int, int, int, int, int]] = []
    for target in range(p):
        witness = None
        for w_value in sorted(regular):
            w = regular_witness[w_value]
            higher = triple_witness.get((target - w_value) % p)
            if higher is not None:
                witness = (target, w, *higher)
                break
        require(witness is not None, "missing regular finite-prime witness")
        if witness is None:
            raise RuntimeError("unreachable witness guard")
        _, w, x, y, z = witness
        require((2 * w - 1) % p != 0, "singular triangular witness")
        require(
            sum(math.comb(k, r) for k, r in zip((w, x, y, z), DEGREES)) % p
            == target,
            "finite-prime witness square-back failure",
        )
        witness_rows.append(witness)

    encoded = ";".join(",".join(map(str, row)) for row in witness_rows).encode("ascii")
    digest = hashlib.sha256(encoded).hexdigest()
    require(sumsets([regular, *higher_sets], p) == set(range(p)), "sumset gap")
    return periods, sizes, digest


def composite_control(limit: int) -> str:
    encoded_rows = []
    for modulus in range(2, limit + 1):
        sets = [value_set(r, modulus) for r in DEGREES]
        represented = sumsets(sets, modulus)
        require(represented == set(range(modulus)), "composite control gap")
        encoded_rows.append(
            f"{modulus}:"
            + ",".join(str(binomial_period(r, modulus)) for r in DEGREES)
            + ":"
            + ",".join(str(len(values)) for values in sets)
        )
    return hashlib.sha256(";".join(encoded_rows).encode("ascii")).hexdigest()


def main() -> None:
    reciprocal_sum = sum(Fraction(1, r) for r in DEGREES)
    require(reciprocal_sum == Fraction(25, 24), "reciprocal-sum regression")
    print(f"degrees={DEGREES}")
    print(f"reciprocal_sum={reciprocal_sum} excess={reciprocal_sum - 1}")
    print(
        "large_prime_bound="
        "(p-1)/2+ceil(p/4)+ceil(p/6)+ceil(p/8)-3>=p for p>=89"
    )
    for p in SMALL_ODD_PRIMES:
        periods, sizes, digest = finite_prime_certificate(p)
        print(
            f"p={p:2d} periods={periods} regular-plus-images={sizes} "
            f"sumset={p} witness_sha256={digest}"
        )

    # The two-adic lift changes exactly the next bit and no lower bit.
    for exponent in range(1, 13):
        step = 1 << (exponent + 1)
        for w in range(0, 257):
            difference = math.comb(w + step, 2) - math.comb(w, 2)
            require(
                difference == (1 << exponent) * (2 * w + step - 1),
                "two-adic identity failure",
            )
            require(
                difference % (1 << (exponent + 1)) == (1 << exponent),
                "two-adic digit toggle failure",
            )
    print("two_adic_lift=PASS exponents=1..12 w=0..256")

    digest = composite_control(300)
    print(f"composite_control=PASS moduli=2..300 sha256={digest}")
    print("universal_modular_solubility=PROVED_WITH_EXACT_FINITE_SIDECAR")


if __name__ == "__main__":
    main()
