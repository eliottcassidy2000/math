#!/usr/bin/env python3
"""Exact finite control for the THM-4026 Eisenstein-norm sidecar.

For a putative representation, with S=C(y,6)+C(z,8), the identities in
THM-4026 give

    6*(N-S)+1 = Norm_Z[omega]((A+B)/2 + B*omega).

Every rational prime p=2 (mod 3) has even valuation in an Eisenstein norm.
This companion applies that necessary row-level condition to the complete
(z,y) box.  It is a search diagnostic, not another nonrepresentation proof:
surviving rows still need the thin x-image and triangular coordinates.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb


TARGET = 896_315_812_331_399
CUTOFFS = (10, 100, 1_000)
EXPECTED_ALIVE = {10: 187_434, 100: 74_911, 1_000: 61_188}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def primes_up_to(bound: int) -> list[int]:
    sieve = bytearray(b"\x01") * (bound + 1)
    sieve[0:2] = b"\x00\x00"
    for divisor in range(2, int(bound**0.5) + 1):
        if sieve[divisor]:
            start = divisor * divisor
            count = (bound - start) // divisor + 1
            sieve[start : bound + 1 : divisor] = b"\x00" * count
    return [p for p in range(2, bound + 1) if sieve[p]]


def row_norm_targets() -> list[int]:
    targets: list[int] = []
    z = 7
    while comb(z, 8) <= TARGET:
        atom8 = comb(z, 8)
        y = 5
        while atom8 + comb(y, 6) <= TARGET - 1:
            targets.append(6 * (TARGET - atom8 - comb(y, 6)) + 1)
            y += 1
        z += 1
    return targets


def check_norm_witness(n: int, witness: tuple[int, int, int, int]) -> None:
    w, x, y, z = witness
    require(
        sum(comb(t, k) for t, k in zip(witness, (2, 4, 6, 8))) == n,
        "witness equality",
    )
    a = x * x - 3 * x + 1
    b = 2 * w - 1
    require(a % 2 == 1 and b % 2 == 1, "norm parity")
    u = (a + b) // 2
    norm = u * u - u * b + b * b
    higher = comb(y, 6) + comb(z, 8)
    require(norm == 6 * (n - higher) + 1, "Eisenstein norm identity")
    require(4 * a + 5 == (2 * x - 3) ** 2, "thin quartic image")


def has_odd_valuation(value: int, prime: int) -> bool:
    parity = False
    while value % prime == 0:
        value //= prime
        parity = not parity
    return parity


def main() -> None:
    check_norm_witness(TARGET - 1, (33_663_667, 9_433, 16, 9))
    check_norm_witness(TARGET + 1, (40_920_205, 6_138, 22, 13))
    combinadic = ((281, 8), (279, 7), (234, 6), (212, 5),
                  (188, 4), (136, 3), (43, 2), (15, 1))
    require(sum(comb(n, k) for n, k in combinadic) == TARGET,
            "rank-eight combinadic identity")

    norm_targets = row_norm_targets()
    require(len(norm_targets) == 248_160, "complete admissible (z,y) rows")
    inert_primes = [p for p in primes_up_to(CUTOFFS[-1]) if p % 3 == 2]
    alive = list(range(len(norm_targets)))
    prime_cursor = 0

    print("THM4026_EISENSTEIN_INERT_ROW_SIEVE")
    print(f"target={TARGET}")
    print(f"admissible_yz_rows={len(norm_targets)}")
    print("norm_identity_controls=N-1,N+1 PASS")
    print("combinadic_control=PASS")

    for cutoff in CUTOFFS:
        before = len(alive)
        while prime_cursor < len(inert_primes) and inert_primes[prime_cursor] <= cutoff:
            prime = inert_primes[prime_cursor]
            alive = [
                index
                for index in alive
                if not has_odd_valuation(norm_targets[index], prime)
            ]
            prime_cursor += 1
        require(len(alive) == EXPECTED_ALIVE[cutoff], f"cutoff {cutoff} survivor count")
        print(
            f"cutoff={cutoff} inert_primes_used={prime_cursor} "
            f"newly_killed={before-len(alive)} alive={len(alive)}"
        )

    print(f"final_alive_fraction={Fraction(len(alive), len(norm_targets))}")
    print("scope=FINITE-EXACT necessary-row-sieve; not a representation certificate")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
