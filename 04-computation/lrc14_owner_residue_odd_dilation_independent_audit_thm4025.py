#!/usr/bin/env python3
"""Independent integer-cross-product audit of THM-4025.

Unlike the primary companion, this engine stores gap records as integer
numerator/denominator pairs and decides the gate by one cleared inequality.
It independently recovers the minimal ordinary-order revival and checks the
divisibility implication.  It does not import the primary script.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def residue_record(t: int, owner: int, second: bool) -> tuple[int, int]:
    divisor = gcd(t, owner)
    if second:
        numerator = (9 * t + 16 * owner) % (126 * divisor)
        denominator = 126 * owner
    else:
        numerator = (3 * t - 4 * owner) % (42 * divisor)
        denominator = 42 * owner
    require(numerator != 0, f"nonzero residue t={t},owner={owner},second={second}")
    return numerator, denominator


def ratio_less(first: tuple[int, int], second: tuple[int, int]) -> bool:
    return first[0] * second[1] < second[0] * first[1]


def closed_integer(
    t: int, U: int, first: tuple[int, int], second: tuple[int, int]
) -> bool:
    """Decide D+epsilon_1+epsilon_2>2/63 after clearing denominators."""
    dn, dd = t * (2 * U - 1), 84 * U * (U - 1)
    an, ad = first
    bn, bd = second
    left = 63 * (dn * ad * bd + an * dd * bd + bn * dd * ad)
    right = 2 * dd * ad * bd
    return left > right


def scan(bound: int):
    rows = []
    for t in range(11, bound + 1, 2):
        best1: tuple[int, int] | None = None
        best2: tuple[int, int] | None = None
        for U in range(1, t + 1):
            first = residue_record(t, U, False)
            second = residue_record(t, U, True)
            if best1 is None or ratio_less(first, best1):
                best1 = first
            if best2 is None or ratio_less(second, best2):
                best2 = second
            if U < max(11, 3 * t // 4 + 1):
                continue
            require(best1 is not None and best2 is not None, f"records t={t},U={U}")
            scale = gcd(t, U)
            rows.append(
                (
                    t,
                    U,
                    (t // scale, U // scale),
                    scale,
                    closed_integer(t, U, best1, best2),
                    best1,
                    best2,
                )
            )
    return rows


def main() -> None:
    # Raw numerator/modulus scaling, before any Fraction construction.
    owner_scaling_checks = 0
    for t in range(3, 202, 2):
        for owner in range(1, t + 1):
            for multiplier in (1, 3, 5, 7, 9):
                for second in (False, True):
                    base = residue_record(t, owner, second)
                    scaled = residue_record(multiplier * t, multiplier * owner, second)
                    require(
                        scaled == (multiplier * base[0], multiplier * base[1]),
                        f"raw scaling {t},{owner},{multiplier},{second}",
                    )
                    owner_scaling_checks += 1

    rows = scan(501)
    by_ray: dict[tuple[int, int], list[tuple]] = {}
    for row in rows:
        by_ray.setdefault(row[2], []).append(row)

    revivals = []
    divisor_checks = 0
    for ray, ray_rows in by_ray.items():
        ray_rows.sort(key=lambda row: row[3])
        seen_survivor = None
        states = {row[3]: row[4] for row in ray_rows}
        for row in ray_rows:
            if not row[4] and seen_survivor is None:
                seen_survivor = row
            if row[4] and seen_survivor is not None:
                revivals.append((row[0], ray, seen_survivor[:5], row[:5]))
                break
        for scale, closed in states.items():
            if closed:
                continue
            for multiplier in range(1, max(states) // scale + 1, 2):
                target = multiplier * scale
                if target in states:
                    require(not states[target], f"survivor multiple {ray},{scale},{target}")
                    divisor_checks += 1
    revivals.sort()
    require(revivals and revivals[0][0:2] == (55, (5, 4)), "minimal revival")

    lookup = {(row[0], row[1]): row for row in rows}
    row45 = lookup[(45, 36)]
    row55 = lookup[(55, 44)]
    require(not row45[4] and row55[4], "hostile states")
    require(row45[5:] == ((1, 966), (1, 2772)), "row45 records")
    require(row55[5:] == ((1, 1722), (17, 5166)), "row55 records")

    # Independent exact evaluation of the two cleared margins.
    D45 = F(45 * 71, 84 * 36 * 35)
    B45 = F(2, 63) - F(1, 966) - F(1, 2772)
    D55 = F(55 * 87, 84 * 44 * 43)
    B55 = F(2, 63) - F(1, 1722) - F(17, 5166)
    require(D45 - B45 == F(-97, 595056), "row45 margin")
    require(D55 - B55 == F(63, 28208), "row55 margin")

    # Direct semigroup transport controls.
    semigroup_checks = 0
    for first in range(1, 80):
        for second in range(1, 80):
            star = 2 * first * second - first - second + 1
            require(2 * star - 1 == (2 * first - 1) * (2 * second - 1), "star")
            semigroup_checks += 1

    print("THM4025_OWNER_RESIDUE_ODD_DILATION_INDEPENDENT_AUDIT")
    print("engine=integer_records_and_one_cleared_gate_inequality")
    print("minimal_ordinary_revival=ray_(5,4);(45,36)=S;(55,44)=C")
    print("records_(45,36)=1/966,1/2772;records_(55,44)=1/1722,17/5166")
    print("margins=-97/595056,63/28208")
    print(f"scan_cells_through_501={len(rows)};rays={len(by_ray)};revivals={len(revivals)}")
    print(f"owner_scaling_checks={owner_scaling_checks}")
    print(f"divisibility_survival_checks={divisor_checks};semigroup_checks={semigroup_checks}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
