#!/usr/bin/env python3
"""Exact controls for the scaled inert cube-class support-two packet.

The mathematical proof is in THM-3818.  This companion exhausts its bounded
5,855-ratio/78-placement atlas, tests nonprimitive scales through divisor
decoding, reconstructs pair-sum-grid safety from ambient residues, and keeps
split-prime, exponent-three, and off-grid hostiles active.

No Python assertions are used, so every gate survives ``python -O``.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction


CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def factor(n: int) -> dict[int, int]:
    require(n >= 1, "factor input must be positive")
    out: dict[int, int] = {}
    candidate = 2
    while candidate * candidate <= n:
        while n % candidate == 0:
            out[candidate] = out.get(candidate, 0) + 1
            n //= candidate
        candidate = 3 if candidate == 2 else candidate + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def admissible_primitive_sum(d: int) -> bool:
    fs = factor(d)
    return bool(fs) and all(p % 3 == 2 and e <= 2 for p, e in fs.items())


def divisors(fs: dict[int, int]) -> list[int]:
    values = [1]
    for prime, exponent in fs.items():
        old = tuple(values)
        power = 1
        for _ in range(exponent):
            power *= prime
            values.extend(v * power for v in old)
    return sorted(values)


def representations_from_sum_divisors(m: int) -> tuple[tuple[int, int], ...]:
    """Recover every positive distinct x<y with x^3+y^3=m."""
    found: set[tuple[int, int]] = set()
    for d in divisors(factor(m)):
        numerator = d * d * d - m
        denominator = 3 * d
        if numerator <= 0 or numerator % denominator:
            continue
        product = numerator // denominator
        disc = d * d - 4 * product
        if disc <= 0:
            continue
        root = math.isqrt(disc)
        if root * root != disc or (d - root) % 2:
            continue
        x, y = (d - root) // 2, (d + root) // 2
        if 0 < x < y and x ** 3 + y ** 3 == m:
            found.add((x, y))
    return tuple(sorted(found))


def placed_covector(p: int, q: int, i: int, j: int) -> tuple[int, ...]:
    row = [0] * 13
    row[i] = q
    row[j] = -p
    return tuple(row)


def decode_covector(a: tuple[int, ...]) -> tuple[int, int, int, int]:
    positive = [(i, v) for i, v in enumerate(a) if v > 0]
    negative = [(i, -v) for i, v in enumerate(a) if v < 0]
    require(len(positive) == len(negative) == 1, "not a placed support-two covector")
    i, q = positive[0]
    j, p = negative[0]
    require(math.gcd(p, q) == 1, "covector is not primitive")
    return p, q, i, j


def residue_schedule(residues: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    safe: list[int] = []
    for k in range(modulus):
        distances = []
        for rho in residues:
            r = (rho * k) % modulus
            distances.append(min(r, modulus - r))
        if all(14 * distance >= modulus for distance in distances):
            safe.append(k)
    return tuple(safe)


def direct_schedule(speeds: tuple[int, ...], modulus: int) -> tuple[int, ...]:
    safe: list[int] = []
    for k in range(modulus):
        good = True
        for speed in speeds:
            value = Fraction(speed * k, modulus)
            frac = value - value.numerator // value.denominator
            if min(frac, 1 - frac) < Fraction(1, 14):
                good = False
                break
        if good:
            safe.append(k)
    return tuple(safe)


def main() -> None:
    ratios: list[tuple[int, int]] = []
    values: dict[int, tuple[int, int]] = {}
    complete_fibres: dict[int, list[tuple[int, int]]] = {}
    for y in range(2, 356):
        for x in range(1, y):
            complete_fibres.setdefault(x ** 3 + y ** 3, []).append((x, y))

    placement_checks = 0
    for d in range(3, 357):
        if not admissible_primitive_sum(d):
            continue
        for p in range(1, (d + 1) // 2):
            q = d - p
            if not p < q or math.gcd(p, q) != 1:
                continue
            ratios.append((p, q))
            m = p ** 3 + q ** 3
            require(m not in values, "cube address collision inside admissible atlas")
            values[m] = (p, q)
            require(complete_fibres[m] == [(p, q)], "bounded cube fibre is not singleton")
            require(representations_from_sum_divisors(m) == ((p, q),),
                    "divisor decoder disagrees with complete coordinate fibre")
            for i in range(13):
                for j in range(i + 1, 13):
                    a = placed_covector(p, q, i, j)
                    pp, qq, ii, jj = decode_covector(a)
                    require((pp, qq, ii, jj) == (p, q, i, j),
                            "placed covector decoder failed")
                    require(sum(abs(v) for v in a) == d, "covector l1 mass changed")
                    require(a[i] * p + a[j] * q == 0, "covector relation failed")
                    require(({k for k, v in enumerate(a) if v > 0},
                             {k for k, v in enumerate(a) if v < 0}) == ({i}, {j}),
                            "exposed face signs lost placement")
                    placement_checks += 4

    require(len(ratios) == 5855, "THM-3793 subatlas count changed")
    require(len(values) == 5855, "cube-class count changed")
    require(placement_checks == 5855 * 78 * 4, "labelled placement count changed")

    scale_cases: list[tuple[int, int, int]] = []
    for p, q in ((1, 4), (2, 9), (5, 6), (1, 24)):
        require(admissible_primitive_sum(p + q), "scale control left admissible ratios")
        for g in (1, 2, 4, 5, 25, 125, 256, 2000):
            require(all(prime % 3 == 2 for prime in factor(g)),
                    "scale control acquired a split prime")
            m = g ** 3 * (p ** 3 + q ** 3)
            reps = representations_from_sum_divisors(m)
            require(reps == ((g * p, g * q),), "all-scale singleton decoder failed")
            x, y = reps[0]
            decoded_g = math.gcd(x, y)
            require((decoded_g, x // decoded_g, y // decoded_g) == (g, p, q),
                    "scale/primitive-ratio reconstruction failed")
            scale_cases.append((g, p, q))

    require(representations_from_sum_divisors(1729) == ((1, 12), (9, 10)),
            "split-prime 1729 hostile changed")
    require(representations_from_sum_divisors(515375) == ((15, 80), (54, 71)),
            "inert exponent-three hostile changed")

    positive_row = tuple(list(range(1, 11)) + [12, 13, 14])
    positive_modulus = 11
    positive_residues = tuple(v % positive_modulus for v in positive_row)
    positive_schedule = residue_schedule(positive_residues, positive_modulus)
    require(positive_schedule == tuple(range(1, 11)), "positive grid control changed")
    require(positive_schedule == direct_schedule(positive_row, positive_modulus),
            "residue/direct grid schedules disagree")

    lifted = tuple(v if i in (0, 9) else v + positive_modulus * (i + 1)
                   for i, v in enumerate(positive_row))
    require(tuple(v % positive_modulus for v in lifted) == positive_residues,
            "ambient residue lift changed")
    require(direct_schedule(lifted, positive_modulus) == positive_schedule,
            "pair-sum schedule was not residue invariant")

    ap_row = tuple(range(1, 14))
    ap_schedule = residue_schedule(tuple(v % 5 for v in ap_row), 5)
    require(ap_schedule == (), "AP pair-sum grid should be empty")
    require(direct_schedule(ap_row, 5) == (), "AP direct pair-sum grid changed")
    require(all(min(Fraction(v, 14) % 1, 1 - (Fraction(v, 14) % 1))
                >= Fraction(1, 14) for v in ap_row),
            "AP off-grid loneliness hostile changed")

    semantic = {
        "ratios": len(ratios),
        "placements": len(ratios) * 78,
        "scale_cases": scale_cases,
        "positive_schedule": positive_schedule,
        "ap_schedule": ap_schedule,
        "first": ratios[0],
        "last": ratios[-1],
    }
    semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")

    print("THM3818 SCALED INERT CUBECLASS PAIR PACKET")
    print("RATIO_UNIVERSE=p<q;gcd=1;p+q<=356;inert_cube_free_sum")
    print(f"RATIOS={len(ratios)}")
    print(f"LABELLED_PLACEMENTS={len(ratios) * 78}")
    print(f"PLACEMENT_ACTIVE_CHECKS={placement_checks}")
    print(f"ALL_SCALE_CONTROLS={len(scale_cases)}")
    print("DECODER=M_to_unique_(gp,gq)_to_(g,p,q)")
    print("GRID_DECODER=(ambient_speeds_mod_g(p+q))_to_exact_1/14_safety_mask")
    print(f"POSITIVE_GRID_CONTROL={','.join(map(str, positive_schedule))}")
    print("OFF_GRID_HOSTILE=AP13_pair_(1,4):grid_empty_but_t=1/14_safe")
    print("COLLISION_HOSTILES=1729_split;515375_inert_exponent_three")
    print(f"CHECKS={CHECKS}")
    print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
