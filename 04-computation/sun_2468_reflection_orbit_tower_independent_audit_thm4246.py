#!/usr/bin/env python3
"""Independent tuple-orbit audit for the Sun reflection tower."""

from __future__ import annotations

from hashlib import sha256
from itertools import product
from math import comb


RANKS = (2, 4, 6, 8)


def factor(n: int) -> list[tuple[int, int]]:
    out = []
    p = 2
    while p * p <= n:
        if n % p:
            p = 3 if p == 2 else p + 2
            continue
        a = 0
        while n % p == 0:
            n //= p
            a += 1
        out.append((p, a))
        p = 3 if p == 2 else p + 2
    if n > 1:
        out.append((n, 1))
    return out


def period(m: int, k: int) -> int:
    ans = 1
    for p, a in factor(m):
        extra = 0
        q = k
        while q >= p:
            q //= p
            extra += 1
        ans *= p ** (a + extra)
    return ans


def value(n: int, k: int, m: int) -> int:
    return (comb(n, k) if n >= k else 0) % m


def audit_modulus(m: int) -> tuple[list[list[int]], list[int], int]:
    periods = [period(m, k) for k in RANKS]
    beta = [[0] * m for _ in range(5)]
    direct = [0] * m
    tuples = 0
    for ns in product(*(range(p) for p in periods)):
        tuples += 1
        target = sum(value(n, k, m) for n, k in zip(ns, RANKS)) % m
        direct[target] += 1
        active = 0
        canonical = True
        for n, k, p in zip(ns, RANKS, periods):
            mate = (k - 1 - n) % p
            assert value(n, k, m) == value(mate, k, m)
            if mate == n:
                continue
            active += 1
            if n > mate:
                canonical = False
                break
        if canonical:
            beta[active][target] += 1
    for t in range(m):
        assert direct[t] == sum((1 << d) * beta[d][t] for d in range(5))
    return beta, direct, tuples


def main() -> None:
    print("SUN_2468_REFLECTION_ORBIT_TOWER_INDEPENDENT")
    digest = sha256()
    total_tuples = 0
    for m in (3, 5, 7, 9, 11, 13):
        beta, direct, tuples = audit_modulus(m)
        total_tuples += tuples
        for t in range(m):
            row = f"{m},{t}," + ",".join(str(beta[d][t]) for d in range(5)) + f",{direct[t]}\n"
            digest.update(row.encode("ascii"))
        print(f"m={m} tuples={tuples} fibre_mass={sum(direct)}")
    print(f"total_tuples={total_tuples}")
    print(f"audit_sha256={digest.hexdigest()}")
    print("verdict=direct_tuple_orbits_match_weighted_tower")


if __name__ == "__main__":
    main()
