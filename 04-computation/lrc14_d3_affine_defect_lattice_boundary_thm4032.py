#!/usr/bin/env python3
"""Exact audit for THM-4032's d=3 affine defect boundary."""

from fractions import Fraction as F
from itertools import combinations
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")
CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def distance(speed, x):
    y = (speed * x) % 1
    return min(y, 1 - y)


def bad_sets(y, speeds):
    return tuple(
        tuple(j for j in range(3) if distance(w, (y + j) / 3) < F(1, 14))
        for w in speeds
    )


def fully_spoiled(y, speeds):
    return set().union(*map(set, bad_sets(y, speeds))) == set(range(3))


def direct_phase(speeds):
    """Independent wall-cell semantics on y in [0,1]."""
    walls = {F(0), F(1)}
    for w in speeds:
        for j in range(3):
            for integer in range(-1, 2 * w + 2):
                for sign in (-1, 1):
                    wall = F(3, w) * (F(integer) + sign * F(1, 14)) - j
                    if 0 <= wall <= 1:
                        walls.add(wall)
    walls = sorted(walls)
    probes = set(walls)
    probes.update((x + y) / 2 for x, y in zip(walls, walls[1:]))
    return next((y for y in sorted(probes) if fully_spoiled(y, speeds)), None)


def defects(a, b, c, A, B, C):
    return (
        3 * b * A - 3 * a * B + a * b,
        3 * c * A - 3 * a * C + 2 * a * c,
        3 * c * B - 3 * b * C + b * c,
    )


def interval_point(a, b, c, A, B, C):
    centres = (F(A, a), F(B, b) - F(1, 3), F(C, c) - F(2, 3))
    radii = (F(1, 14 * a), F(1, 14 * b), F(1, 14 * c))
    left = max(x - r for x, r in zip(centres, radii))
    right = min(x + r for x, r in zip(centres, radii))
    return (left + right) / 2 if left < right else None


def oriented_lattice(a, b, c):
    for A in range(a):
        for B in range(2 * b + 1):
            for C in range(2 * c + 1):
                n_ab, n_ac, n_bc = defects(a, b, c, A, B, C)
                if not 14 * abs(n_ab) < 3 * (a + b):
                    continue
                if not 14 * abs(n_ac) < 3 * (a + c):
                    continue
                if not 14 * abs(n_bc) < 3 * (b + c):
                    continue
                point = interval_point(a, b, c, A, B, C)
                require(point is not None, ("Helly", a, b, c, A, B, C))
                require(c * n_ab - b * n_ac + a * n_bc == 0, "circuit")
                require(n_ab % 3 == (a * b) % 3, "ab mod3")
                require(n_ac % 3 == (2 * a * c) % 3, "ac mod3")
                require(n_bc % 3 == (b * c) % 3, "bc mod3")
                y = (3 * point) % 1
                require(fully_spoiled(y, (a, b, c)), "lattice reconstructs phase")
                return A, B, C, n_ab, n_ac, n_bc, point, y
    return None


def lattice(sp):
    a, b, c = sp
    return oriented_lattice(a, b, c) or oriented_lattice(a, c, b)


def defect_bound(x, y):
    return (3 * (x + y) - 1) // 14


def oriented_defect(a, b, c):
    for n_ab in range(-defect_bound(a, b), defect_bound(a, b) + 1):
        if n_ab % 3 != (a * b) % 3 or n_ab % gcd(a, b):
            continue
        for n_ac in range(-defect_bound(a, c), defect_bound(a, c) + 1):
            if n_ac % 3 != (2 * a * c) % 3 or n_ac % gcd(a, c):
                continue
            numerator = b * n_ac - c * n_ab
            if numerator % a:
                continue
            n_bc = numerator // a
            if abs(n_bc) > defect_bound(b, c):
                continue
            if n_bc % 3 != (b * c) % 3 or n_bc % gcd(b, c):
                continue
            k_ab = (n_ab - a * b) // 3
            k_ac = (n_ac - 2 * a * c) // 3
            k_bc = (n_bc - b * c) // 3
            for A in range(a):
                if (b * A - k_ab) % a or (c * A - k_ac) % a:
                    continue
                B = (b * A - k_ab) // a
                C = (c * A - k_ac) // a
                require(c * B - b * C == k_bc, "CRT third equation")
                require(
                    defects(a, b, c, A, B, C) == (n_ab, n_ac, n_bc),
                    "realize defects",
                )
                return n_ab, n_ac, n_bc, A, B, C
            raise AssertionError(("unrealized", a, b, c, n_ab, n_ac, n_bc))
    return None


def defect_certificate(sp):
    a, b, c = sp
    return oriented_defect(a, b, c) or oriented_defect(a, c, b)


def gcd_gates(sp):
    return all(14 * gcd(x, y) < 3 * (x + y) for x, y in combinations(sp, 2))


def safe_phase(speeds):
    walls = {F(0), F(1)}
    for w in speeds:
        for integer in range(w + 1):
            for sign in (-1, 1):
                x = (F(integer) + sign * F(1, 14)) / w
                if 0 <= x <= 1:
                    walls.add(x)
    walls = sorted(walls)
    probes = set(walls)
    probes.update((x + y) / 2 for x, y in zip(walls, walls[1:]))
    return next(
        (x for x in sorted(probes) if min(distance(w, x) for w in speeds) >= F(1, 14)),
        None,
    )


def main():
    values = tuple(x for x in range(1, 24) if x % 3)
    profiles = 0
    first_false_gate = None
    for sp in combinations(values, 3):
        direct = direct_phase(sp)
        affine = lattice(sp)
        defect = defect_certificate(sp)
        require((direct is not None) == (affine is not None), ("direct/lattice", sp))
        require((direct is not None) == (defect is not None), ("direct/defect", sp))
        if gcd_gates(sp) and direct is None and first_false_gate is None:
            first_false_gate = sp
        profiles += 1

    minimal = None
    for maximum in range(1, 30):
        vals = tuple(x for x in range(1, maximum + 1) if x % 3)
        for sp in combinations(vals, 3):
            y = direct_phase(sp)
            if y is not None:
                minimal = (maximum, sp, y, bad_sets(y, sp), lattice(sp))
                break
        if minimal:
            break

    # The divided pack H={1,...,10} has safe phases k/11.  Pair entries
    # 1,10 and the eight body entries 2,...,9 lift back by a factor of three.
    pack_hostile = None
    for maximum in range(1, 30):
        vals = tuple(x for x in range(1, maximum + 1) if x % 3)
        for sp in combinations(vals, 3):
            for k in range(1, 11):
                y = F(k, 11)
                if fully_spoiled(y, sp):
                    pack_hostile = (maximum, sp, y, bad_sets(y, sp))
                    break
            if pack_hostile:
                break
        if pack_hostile:
            break

    require(minimal is not None, "minimal exists")
    require(pack_hostile is not None, "pack hostile exists")
    a, b, c = pack_hostile[1]
    body = tuple(3 * x for x in range(2, 10)) + (a, b, c)
    pair = (3, 30)
    row = body + pair
    require(len(set(row)) == 13, ("typed row distinct", row))
    require(gcd(*body) == 1, "typed body primitive")
    lonely = safe_phase(row)
    require(lonely is not None, "typed row has safe phase")

    print("LRC14_D3_AFFINE_DEFECT_LATTICE_BOUNDARY_THM4032")
    print("iff=direct_phase=affine_centres=mod3_defect_circuit")
    print(f"profiles={profiles};first_false_gcd_gate={first_false_gate}")
    print(f"minimal_distinct_geometric={minimal}")
    print(f"minimal_H10_pack_hostile={pack_hostile}")
    print(
        f"typed_row={row};safe_phase={lonely};"
        f"clearance={min(distance(w, lonely) for w in row)}"
    )
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
