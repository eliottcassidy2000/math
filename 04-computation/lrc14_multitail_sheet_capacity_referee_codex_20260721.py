#!/usr/bin/env python3
"""Exact integer referee for THM-2064."""

from __future__ import annotations

from fractions import Fraction
from math import gcd, lcm
import random


def require(condition: bool) -> None:
    if not condition:
        raise AssertionError("exact referee check failed")


def abs_residue(k: int, modulus: int) -> int:
    r = k % modulus
    return min(r, modulus - r)


def dangerous(k: int, modulus: int) -> bool:
    return 14 * abs_residue(k, modulus) < modulus


def danger_cap(t: int) -> int:
    return (t + 6) // 7


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def tail_data(N: int, a: int, w: int) -> tuple[int, int, int, int]:
    Q = N * a
    g = gcd(w, Q)
    h = Q // g
    u = w // g
    d = gcd(N, h)
    return h, u, d, h // d


def check_open_arc_bound(limit: int) -> int:
    cases = 0
    for h in range(1, limit + 1):
        for d in divisors(h):
            t = h // d
            for residue_class in range(d):
                count = sum(dangerous(residue_class + d * j, h) for j in range(t))
                require(count <= danger_cap(t))
                cases += 1
    return cases


def check_random_fibers(trials: int) -> tuple[int, int]:
    rng = random.Random(2064)
    certified = 0
    for _ in range(trials):
        N = rng.randint(1, 60)
        a = rng.randint(2, 20)
        Q = N * a
        r = rng.randrange(N)
        tails = [rng.randint(1, 3 * Q) for _ in range(rng.randint(1, 5))]
        data = [tail_data(N, a, w) for w in tails]
        L = lcm(N, *(h for h, _, _, _ in data))
        T = L // N
        fiber = [r + N * j for j in range(T)]

        capacity = Fraction(0, 1)
        for h, unit, _d, t in data:
            count = sum(dangerous(unit * k, h) for k in fiber)
            require(count * t <= T * danger_cap(t))
            capacity += Fraction(danger_cap(t), t)

        if capacity < 1:
            require(any(all(not dangerous(unit * k, h) for h, unit, _d, _t in data) for k in fiber))
            certified += 1
    return trials, certified


def check_two_tail_exception(limit: int) -> int:
    exceptional = []
    for tx in range(2, limit + 1):
        for ty in range(2, limit + 1):
            capacity = Fraction(danger_cap(tx), tx) + Fraction(danger_cap(ty), ty)
            if capacity >= 1:
                exceptional.append((tx, ty))
    require(exceptional == [(2, 2)])
    return (limit - 1) ** 2


def check_primitive_dyadic_box() -> tuple[int, int]:
    antecedents = 0
    primitive_antecedents = 0
    for a in range(2, 25):
        for N in range(1, 37):
            two_sheet = []
            for w in range(1, 73):
                _h, _u, _d, t = tail_data(N, a, w)
                if t == 2:
                    two_sheet.append(w)
            for x in two_sheet:
                for y in two_sheet:
                    antecedents += 1
                    require(a % 2 == 0)
                    require(x % (a // 2) == 0 and y % (a // 2) == 0)
                    require((x // (a // 2)) % 2 == 1 and (y // (a // 2)) % 2 == 1)
                    if gcd(gcd(a, x), y) == 1:
                        primitive_antecedents += 1
                        require(a == 2 and x % 2 == 1 and y % 2 == 1)
    return antecedents, primitive_antecedents


def main() -> None:
    print("LRC14 MULTI-TAIL SHEET-CAPACITY AUDIT -- exact integer fibers")
    print(f"open-arc cosets through period 180: {check_open_arc_bound(180)} PASS")
    trials, certified = check_random_fibers(12000)
    print(f"multi-tail fibers: {trials} checked, {certified} capacity-certified PASS")
    print(f"two-tail sheet-order pairs through 240: {check_two_tail_exception(240)} PASS")
    antecedents, primitive = check_primitive_dyadic_box()
    print(f"dyadic antecedents: {antecedents}; primitive antecedents: {primitive} PASS")
    print("PASS")


if __name__ == "__main__":
    main()
