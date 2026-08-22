#!/usr/bin/env python3
"""Exact companion for THM-3645: local-prime support and mod-24 cones."""

from __future__ import annotations

import ast
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path


CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def digest(value: object) -> str:
    return sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()


def primes_up_to(limit: int) -> tuple[int, ...]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                (limit - p * p) // p + 1
            )
    return tuple(p for p in range(2, limit + 1) if sieve[p])


def prime_divisors(value: int, primes: tuple[int, ...]) -> tuple[int, ...]:
    n = value
    factors: list[int] = []
    for p in primes:
        if p * p > n:
            break
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
    if n > 1:
        factors.append(n)
    return tuple(factors)


def candidates(limit: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (m, n)
        for n in range(3, limit + 1, 2)
        for m in range(n // 2 + 1, n)
        if m % 2 == 0 and gcd(m, n) == 1
    )


def conic_coefficients(m: int, n: int) -> tuple[int, int]:
    return n * n * (4 * m * m - n * n), 4 * (
        2 * m * m + 2 * m * n - n * n
    )


def soluble_mod_bruteforce(m: int, n: int, p: int) -> bool:
    coefficient, constant = conic_coefficients(m, n)
    left = {(3 * w * w) % p for w in range(p)}
    right = {(coefficient * z * z + constant) % p for z in range(p)}
    return not left.isdisjoint(right)


def legendre(a: int, p: int) -> int:
    require(p > 2, ("Legendre prime", p))
    a %= p
    if a == 0:
        return 0
    value = pow(a, (p - 1) // 2, p)
    require(value in (1, p - 1), ("Euler criterion", a, p, value))
    return 1 if value == 1 else -1


def support_data(m: int, n: int) -> tuple[int, int, int]:
    u = n - m
    v = 2 * m - n
    t = 2 * m + n
    require(u > 0 and v > 0, ("positive transform", m, n))
    require(u % 2 == v % 2 == t % 2 == 1, ("odd transform", m, n))
    require(gcd(u, v) == 1, ("primitive transform", m, n))
    require(gcd(n, v) == gcd(n, t) == gcd(v, t) == 1,
            ("pairwise support", m, n))
    require(4 * m * m - n * n == v * t, ("coefficient split", m, n))
    return u, v, t


def obstruction_primes(m: int, n: int, factor_primes: tuple[int, ...]) -> tuple[int, ...]:
    """All mod-prime obstructions, using the proved support reduction."""
    u, v, t = support_data(m, n)
    bad: list[int] = []
    if not (m % 3 == n % 3 and n % 3 != 0):
        bad.append(3)
    for p in set(prime_divisors(n, factor_primes) + prime_divisors(v, factor_primes)):
        if p not in (2, 3) and legendre(6, p) == -1:
            bad.append(p)
    for p in prime_divisors(t, factor_primes):
        if p not in (2, 3) and legendre(-2, p) == -1:
            bad.append(p)
    return tuple(sorted(set(bad)))


def survives_screen(m: int, n: int, screen_limit: int,
                    factor_primes: tuple[int, ...]) -> bool:
    return not any(p <= screen_limit for p in obstruction_primes(m, n, factor_primes))


def first_brute_obstruction(m: int, n: int,
                            primes: tuple[int, ...]) -> int | None:
    return next((p for p in primes if not soluble_mod_bruteforce(m, n, p)), None)


def main() -> None:
    factor_primes = primes_up_to(2000)
    screen_997 = tuple(p for p in factor_primes if p <= 997)
    require(len(screen_997) == 168 and screen_997[-1] == 997, "prime screen")

    # Independent small-universe brute force compares the least structural
    # obstruction with the least brute-force obstruction.  The theorem's
    # hostile audit separately checked every candidate-prime pair here.
    hostile_primes = primes_up_to(97)
    for m, n in candidates(31):
        structural = next(
            (p for p in obstruction_primes(m, n, factor_primes) if p <= 97),
            None,
        )
        brute = first_brute_obstruction(m, n, hostile_primes)
        require(structural == brute, ("small brute disagreement", m, n,
                                      structural, brute))

    slopes_401 = candidates(401)
    survivors_401 = tuple(
        slope for slope in slopes_401
        if not obstruction_primes(*slope, factor_primes)
    )
    require(len(slopes_401) == 8195, "401 candidate count")
    require(len(survivors_401) == 104, "401 survivor count")
    require(digest(survivors_401)
            == "2b1e24ba10e59deac190e1b4d20a0afe3792c138ecc84e06fa22e1671448ff8e",
            "THM-3640 survivor digest")

    slopes_997 = candidates(997)
    survivors_997 = tuple(
        slope for slope in slopes_997
        if not obstruction_primes(*slope, factor_primes)
    )
    require(len(slopes_997) == 50_507, "997 candidate count")
    require(len(survivors_997) == 481, "997 survivor count")
    residue_histogram = {
        residue: sum(1 for m, n in survivors_997 if (2 * m - n) % 24 == residue)
        for residue in (5, 23)
    }
    require(residue_histogram == {5: 238, 23: 243}, "997 residue histogram")
    residue_pairs = tuple(sorted({(n % 24, (2 * m - n) % 24)
                                  for m, n in survivors_997}))
    require(residue_pairs == ((5, 23), (23, 5)), "two residue cones")

    denominators_after_401 = tuple(sorted({n for _m, n in survivors_997 if n > 401}))
    require(denominators_after_401[0] == 431, "first new denominator")
    local_431 = tuple(slope for slope in survivors_997 if slope[1] == 431)
    require(len(local_431) == 8, "431 local survivor count")

    # Find the first slope that passes the fixed p<=997 screen but violates
    # the full-screen residue cones.  This is deliberately a different
    # quantifier from the n<=997 theorem above.
    first_fixed_screen_escape = None
    for slope in candidates(1039):
        m, n = slope
        if survives_screen(m, n, 997, factor_primes):
            pair = (n % 24, (2 * m - n) % 24)
            if pair not in ((5, 23), (23, 5)):
                first_fixed_screen_escape = slope
                break
    require(first_fixed_screen_escape == (512, 1019), "first fixed-screen escape")
    escape = first_fixed_screen_escape
    require(first_brute_obstruction(*escape, screen_997) is None,
            "escape brute screen")
    require(first_brute_obstruction(*escape, (1019,)) == 1019,
            "hidden prime obstruction")
    u, v, t = support_data(*escape)
    require((u, v, t) == (507, 5, 2043), "escape transform")
    require(legendre(6, 1019) == -1,
            "escape mechanism")

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source raw LF")
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(source_bytes.decode("utf-8")))),
            "Python assert node present")

    semantic = {
        "survivors_401": survivors_401,
        "survivors_997": survivors_997,
        "residue_pairs": residue_pairs,
        "local_431": local_431,
        "fixed_screen_escape": escape,
        "escape_uvt": (u, v, t),
    }
    print("== THM-3645 two-cube local-prime support and mod-24 cones ==")
    print("support=(U=n-m,V=2m-n,T=2m+n);4m^2-n^2=V*T")
    print("odd_prime_conditions=p|nV:(6/p)=1;p|T:(-2/p)=1;off_support:automatic")
    print("p3_condition=m=n!=0(mod3);p2_condition=automatic")
    print("complete_screen_rule=screen_limit>=n sees n,V,T/3 prime support")
    print(f"n401=(candidates:{len(slopes_401)},survivors:{len(survivors_401)},digest:{digest(survivors_401)})")
    print(f"n997=(candidates:{len(slopes_997)},survivors:{len(survivors_997)},digest:{digest(survivors_997)})")
    print(f"residue_pairs={residue_pairs};V_histogram={residue_histogram}")
    print(f"first_new_denominator={denominators_after_401[0]};local_431={local_431}")
    print(f"fixed_screen_997_first_escape={escape};UVT={(u,v,t)};hidden_obstruction=1019")
    print(f"semantic_sha256={digest(semantic)}")
    print(f"CHECKS={CHECKS}")
    print("status=PROVED LOCAL-ALGEBRA + FINITE-EXACT;mod-prime screen only")
    print("scope=no p-adic sufficiency,no Pell-orbit sufficiency,no infinite admissible-slope family")


if __name__ == "__main__":
    main()
