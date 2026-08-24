#!/usr/bin/env python3
"""Independent exact audit for THM-4004.

This path is deliberately separate from the primary certificate.  It checks
the open-arc orbit count by endpoint perturbation, branch multiplicities,
prime thresholds, hostile scope, and the crossing and role-swap arithmetic.
"""

from __future__ import annotations

from fractions import Fraction as F
from functools import lru_cache
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")


def require(condition: bool, label: str, emit: bool = True) -> None:
    if not condition:
        raise RuntimeError(label)
    if emit:
        print(f"PASS {label}")


def distance(speed: int, phase: F) -> F:
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def in_open_danger(phase: F) -> bool:
    residue = phase % 1
    return min(residue, 1 - residue) < F(1, 14)


def orbit_formula(order: int) -> int:
    return (order + 6) // 7


def branch_formula(modulus: int, delta: int) -> int:
    common = gcd(modulus, delta)
    return common * orbit_formula(modulus // common)


@lru_cache(maxsize=None)
def direct_open_arc_max(order: int) -> int:
    """Maximize an order-m orbit in the strict danger arc by exact walls."""
    epsilon = F(1, 100 * 7 * order)
    points = tuple(F(k, order) for k in range(order))
    anchors = points if order <= 50 else (points[0], points[order // 2])
    candidates = []
    for point in anchors:
        for shift in (
            -F(1, 14) + epsilon - point,
            F(1, 14) - epsilon - point,
        ):
            candidates.append(sum(in_open_danger(other + shift) for other in points))
    return max(candidates)


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


orbit_checks = 0
for order in range(1, 401):
    require(
        direct_open_arc_max(order) == orbit_formula(order),
        f"open-orbit count m={order}",
        emit=False,
    )
    orbit_checks += 1

branch_checks = 0
for modulus in range(2, 151):
    for delta in range(1, 2 * modulus + 1):
        common = gcd(modulus, delta)
        direct = common * direct_open_arc_max(modulus // common)
        require(direct == branch_formula(modulus, delta), "branch multiplicity", emit=False)
        branch_checks += 1

normalized = {order: F(orbit_formula(order), order) for order in range(2, 401)}
require(normalized[2] == F(1, 2), "two-comb equality at order two")
require(
    all(value < F(1, 2) for order, value in normalized.items() if order >= 3),
    "two-comb strictness for every audited order at least three",
)

primes = tuple(value for value in range(2, 10000) if is_prime(value))
require(
    all(3 * orbit_formula(ell) < ell for ell in primes if ell >= 5),
    "three-comb prime gate for all audited primes ell>=5",
)
require(3 * orbit_formula(3) == 3, "ell=3 exact equality")
require(3 * orbit_formula(2) > 2, "ell=2 exceeds branch budget")

pair = (1, 4)
divisible3 = (6, 9, 15, 18, 21, 24, 27, 30)
detuned3 = (1, 10, 11)
pack3 = tuple(value // 3 for value in divisible3) + pair
require(tuple(sorted(pack3)) == tuple(range(1, 11)), "ell=3 pack H={1,...,10}")
y3 = F(1, 11)
require(min(distance(value, y3) for value in pack3) == F(1, 11), "ell=3 pack witness")
branches3 = tuple((y3 + branch) / 3 for branch in range(3))
require(
    tuple(distance((1, 11, 10)[k], branches3[k]) for k in range(3))
    == (F(1, 33), F(0), F(1, 33)),
    "ell=3 labelled combs saturate selected branches",
)
speeds3 = detuned3 + divisible3 + (3, 12)
require(
    min(distance(value, F(4, 33)) for value in speeds3) == F(1, 11),
    "ell=3 control is lonely at another pack phase",
)

divisible5 = (10, 15, 25, 30, 35, 40, 45, 50)
detuned5 = (1, 2, 3)
pack5 = tuple(value // 5 for value in divisible5) + pair
require(tuple(sorted(pack5)) == tuple(range(1, 11)), "ell=5 pack H={1,...,10}")
branches5 = tuple((F(1, 11) + branch) / 5 for branch in range(5))
speeds5 = detuned5 + divisible5 + (5, 20)
require(
    tuple(min(distance(value, phase) for value in speeds5) for phase in branches5)
    == (F(1, 55), F(1, 11), F(1, 11), F(1, 11), F(1, 11)),
    "ell=5 positive control",
)

z = F(13, 20)
require(
    min(distance(value, z) for value in (3, 7)) == F(1, 20)
    and min(distance(value, z + F(1, 2)) for value in (3, 7)) == F(1, 20),
    "reduced pair (3,7) is a valid two-lift hostile",
)
z_smaller = F(3, 16)
require(
    min(distance(value, z_smaller) for value in (3, 5)) == F(1, 16)
    and min(distance(value, z_smaller + F(1, 2)) for value in (3, 5)) == F(1, 16),
    "(3,7) is valid but not minimal",
)

q_height = 91**6
require(q_height == 567_869_252_041, "crossing height Q")
integer_floor = q_height // 177 + 1
require(integer_floor == 3_208_300_859, "integer crossing floor")
require(177 * (integer_floor - 1) < q_height < 177 * integer_floor, "strict rounding")

role_checks = 0
for maximum in (2, 3, 17, 1000):
    for multiplier in range(1, maximum):
        require(F(6 * multiplier, 7 * maximum) < F(6, 7) < 1, "role swap", emit=False)
        role_checks += 1

print(f"orbit_checks={orbit_checks}")
print(f"branch_checks={branch_checks}")
print(f"prime_checks={sum(ell >= 5 for ell in primes)}")
print(f"role_checks={role_checks}")
print(f"integer_crossing_floor={integer_floor}")
print("STATUS=PASS; (3,7) IS VALID BUT NOT MINIMAL")
