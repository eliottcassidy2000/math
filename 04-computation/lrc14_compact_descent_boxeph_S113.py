#!/usr/bin/env python3
"""Corrected exact audit of the one S113 compact boundary row.

The original version checked only moduli 2..13, used QMAX=250, removed the
maximum rather than the speed exposing the dilated AP, and printed an
unsupported equivalence claim.  This version uses the complete pair-sum
ruler from THM-668/THM-1002 and states only what the displayed row proves.
"""

from fractions import Fraction as F
from math import gcd, isqrt


def divisors(number: int) -> set[int]:
    result: set[int] = set()
    for candidate in range(1, isqrt(number) + 1):
        if number % candidate == 0:
            result.add(candidate)
            result.add(number // candidate)
    return result


def exact_M(values: tuple[int, ...]) -> F:
    rulers: set[int] = set()
    for index, left in enumerate(values):
        for right in values[index:]:
            rulers.update(divisors(left + right))

    best = F(0)
    for denominator in rulers:
        if denominator < 2:
            continue
        for numerator in range(1, denominator):
            if gcd(numerator, denominator) != 1:
                continue
            clearance = min(
                min((numerator * speed) % denominator,
                    denominator - (numerator * speed) % denominator)
                for speed in values
            )
            best = max(best, F(clearance, denominator))
    return best


def covers(values: tuple[int, ...], upper: int) -> bool:
    return all(
        any(speed % modulus == 0 for speed in values)
        for modulus in range(2, upper + 1)
    )


V = (2, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 22, 24)
maximum_deletion = V[:-1]
ap_deletion = tuple(speed for speed in V if speed != 13)
rho = F(V[-1], V[-2])
M_V = exact_M(V)
M_max_deletion = exact_M(maximum_deletion)
M_ap_deletion = exact_M(ap_deletion)
descent_bound = rho * M_max_deletion / (rho + 1)

assert covers(V, 14)
assert ap_deletion == tuple(2 * value for value in range(1, 13))
assert M_V == M_ap_deletion == F(1, 13)
assert M_max_deletion == F(1, 12)
assert descent_bound == F(1, 23) < F(1, 13)

print("CORRECTED S113 SINGLE-ROW AUDIT")
print(f"V={V}")
print(f"primitive={gcd(*V) == 1} Cover14={covers(V, 14)} rho={rho}")
print(f"M(V)={M_V}")
print(f"delete maximum 24: M={M_max_deletion}, descent bound={descent_bound}<1/13")
print(f"delete 13: core=2*[12], M={M_ap_deletion}")
print("scope=one exact boundary row; no 100-row/5-of-15 bank is stored")
print("logic=compact 1/13 floor is sufficient for LRC14, not proved equivalent")
print("remaining=tight-deletion extraction (crown collapse) + n=12 equality rigidity")
