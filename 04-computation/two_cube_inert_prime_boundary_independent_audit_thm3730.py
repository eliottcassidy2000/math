#!/usr/bin/env python3
"""Independent hostile probe for the inert-prime two-cube boundary proof.

This deliberately does not import the candidate probe.  It checks the coarse
cube-free sieve, root lifting, the original singleton conclusion, and the
stronger observation that cube-freeness is unnecessary.
"""

from __future__ import annotations

from hashlib import sha256
from math import isqrt


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
    return tuple(i for i, flag in enumerate(flags) if flag)


def cofactor(prime, left):
    return 3 * left * left - 3 * prime * left + prime * prime


def cube_free(value, primes):
    for ell in primes:
        if ell**3 > value:
            return True
        if value % (ell**3) == 0:
            return False
    raise RuntimeError(("prime table too short", value, primes[-1]))


def digest(value):
    return sha256(repr(value).encode("ascii")).hexdigest()


def main():
    prime_limit = 5000
    primes = primes_up_to(prime_limit)
    inert = tuple(p for p in primes if p > 3 and p % 3 == 2)

    # Independent replay of the original coarse cube-free count.
    cube_free_census = []
    for p in inert:
        good = 0
        for x in range(1, (p + 1) // 2):
            q = cofactor(p, x)
            require(p * p // 4 < q < p * p, ("range", p, x, q))
            require(q % p != 0 and q % 3 == 1, ("local exclusions", p, x, q))
            good += cube_free(q, primes)
        if p >= 1728:
            require(12 * good > p - 3, ("strict coarse bound", p, good))
            require(13 * good >= p, ("one-thirteenth bound", p, good))
        cube_free_census.append((p, good, (p - 1) // 2))

    # Direct roots modulo ell^3: every root is the unique lift of a simple
    # root modulo ell, including the characteristic-two edge case.
    lifting_checks = 0
    for p in (5, 11, 71, 101):
        for ell in primes:
            if ell in (3, p) or ell**3 >= p * p:
                continue
            roots_ell = tuple(
                r for r in range(ell) if cofactor(p, r) % ell == 0
            )
            roots_cube = tuple(
                r for r in range(ell**3)
                if cofactor(p, r) % (ell**3) == 0
            )
            require(len(roots_ell) <= 2, ("too many base roots", p, ell))
            require(
                len(roots_cube) == len(roots_ell),
                ("nonunique or missing lift", p, ell, roots_ell, roots_cube),
            )
            require(
                tuple(sorted(r % ell for r in roots_cube)) == roots_ell,
                ("wrong lifted residue", p, ell, roots_ell, roots_cube),
            )
            lifting_checks += 1

    # Stronger finite hostile: use every x, not just cube-free cofactors.
    # If p <= coordinate_limit then m < p^3, so every positive summand in an
    # alternative representation is < p <= coordinate_limit.  The box below
    # is therefore complete for every target it checks.
    coordinate_limit = 800
    targets = {}
    for p in inert:
        if p > coordinate_limit:
            break
        for x in range(1, (p + 1) // 2):
            y = p - x
            value = x**3 + y**3
            previous = targets.setdefault(value, (x, y, p))
            require(previous == (x, y, p), ("source collision", value, previous))

    fibres = {value: [] for value in targets}
    for u in range(1, coordinate_limit):
        for v in range(u + 1, coordinate_limit):
            value = u**3 + v**3
            if value in fibres:
                fibres[value].append((u, v, u + v))
    for value, source in targets.items():
        x, y, p = source
        require(
            fibres[value] == [(x, y, p)],
            ("strong singleton failed", value, source, fibres[value]),
        )

    # Load-bearing hostile for the original presentation: the first
    # non-cube-free cofactor is still a singleton, showing that cube-freeness
    # is sufficient but not necessary.
    hostile_p, hostile_x = 71, 16
    hostile_q = cofactor(hostile_p, hostile_x)
    hostile_m = hostile_p * hostile_q
    require(hostile_q == 7**4, ("hostile cofactor", hostile_q))
    require(not cube_free(hostile_q, primes), "hostile unexpectedly cube-free")
    require(
        fibres[hostile_m] == [(16, 55, 71)],
        ("hostile singleton", fibres[hostile_m]),
    )

    # Cross-prime collision check for the stronger all-x family on a larger
    # bounded universe.  This tests the consequence object directly.
    all_x_seen = {}
    for p in inert:
        for x in range(1, (p + 1) // 2):
            y = p - x
            value = x**3 + y**3
            previous = all_x_seen.setdefault(value, (p, x, y))
            require(previous == (p, x, y), ("cross-prime collision", value))

    tail = tuple(row for row in cube_free_census if row[0] >= 1728)
    semantic = (
        prime_limit,
        len(inert),
        len(cube_free_census),
        min(12 * good - (p - 3) for p, good, _ in tail),
        lifting_checks,
        coordinate_limit,
        len(targets),
        hostile_p,
        hostile_x,
        hostile_q,
        len(all_x_seen),
    )
    print("independent two-cube boundary referee probe")
    print("prime_limit=" + repr(prime_limit))
    print("inert_prime_count=" + repr(len(inert)))
    print("strict_cube_free_minimum_margin=" + repr(semantic[3]))
    print("direct_root_lifting_checks=" + repr(lifting_checks))
    print("complete_singleton_targets=" + repr(len(targets)))
    print("first_non_cube_free_singleton=(71,16,55,2401)")
    print("all_x_cross_prime_values=" + repr(len(all_x_seen)))
    print("semantic_sha256=" + digest(semantic))
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
