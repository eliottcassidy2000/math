#!/usr/bin/env python3
"""Exact regression for the symbolic fixed-bank theorem THM-2161.

For Q>=14, the theorem uses

    L = lcm(2,...,Q),
    V_Q = {1,...,6} union {L+7,...,L+13}.

Every labelled residue modulo q<=Q agrees with the tight AP {1,...,13}, but
the phase x=(L+8)/(8(L+14)) has minimum distance x>3/41.  This script checks
the theorem's consequence object by exact integer/Fraction arithmetic on
several hostile prefix banks.  The replay is a regression, not the proof of
the all-Q quantifier, which is algebraic in THM-2161.
"""

from collections import Counter
from fractions import Fraction as QF
from math import gcd, lcm

H = QF(3, 41)
AP = tuple(range(1, 14))


def circle_distance(x: QF) -> QF:
    y = x % 1
    return min(y, 1 - y)


def modulus_numerator(speeds: tuple[int, ...], q: int) -> int:
    return max(
        min(min((a * v) % q, q - (a * v) % q) for v in speeds)
        for a in range(1, q)
    )


def construct(bound: int) -> tuple[int, tuple[int, ...], QF]:
    assert bound >= 14
    L = lcm(*range(2, bound + 1))
    speeds = tuple(range(1, 7)) + tuple(L + j for j in range(7, 14))
    phase = QF(L + 8, 8 * (L + 14))
    return L, speeds, phase


def audit(bound: int) -> tuple[int, int, QF, int, QF, int]:
    L, speeds, phase = construct(bound)

    # The exact definition used in THM-2146/THM-2161 is
    # defect(V)=|V\\{1,...,13}|, not distance from the closest affine AP.
    assert len(speeds) == 13
    assert len(set(speeds)) == 13
    assert tuple(sorted(speeds)) == speeds
    assert gcd(*speeds) == 1
    assert len(set(speeds) - set(AP)) == 7
    assert set(speeds) & set(AP) == set(range(1, 7))

    best = QF(0)
    best_q = 0
    for q in range(2, bound + 1):
        assert L % q == 0
        assert all(speeds[i - 1] % q == i % q for i in range(1, 14))
        assert Counter(v % q for v in speeds) == Counter(i % q for i in AP)
        actual = QF(modulus_numerator(speeds, q), q)
        control = QF(modulus_numerator(AP, q), q)
        assert actual == control
        assert control <= QF(1, 14) < H
        if actual > best:
            best, best_q = actual, q
    assert best == QF(1, 14)
    assert best_q == 14

    distances = tuple(circle_distance(v * phase) for v in speeds)
    minimum = min(distances)
    assert minimum == phase
    assert distances[0] == phase
    assert distances[-1] == phase
    assert minimum > H

    # The phase denominator is adaptive: after reduction it exceeds the bank.
    M = L // 8
    assert L == 8 * M
    assert phase == QF(M + 1, 8 * M + 14)
    common = gcd(M + 1, 8 * M + 14)
    assert common == gcd(M + 1, 6) <= 6
    reduced_denominator = phase.denominator
    assert reduced_denominator >= (L + 14) // 6
    assert L % (bound * (bound - 1)) == 0
    assert reduced_denominator > bound

    return L, speeds[-1], minimum, reduced_denominator, best, best_q


def main() -> None:
    bounds = (14, 15, 16, 17, 23, 41, 47, 64, 97)
    print("THM-2161 FIXED-PREFIX MODULUS BLINDNESS REGRESSION")
    print("V_Q={1,...,6,L+7,...,L+13}; L=lcm(2,...,Q); threshold h=3/41")
    print("all arithmetic is exact; the all-Q proof is symbolic in THM-2161")
    print()
    print(" Q  L digits  max-speed digits  bank best (q)  phase minimum  phase-denominator digits")
    for bound in bounds:
        L, maximum, minimum, denominator, best, best_q = audit(bound)
        print(
            f"{bound:>2}  {len(str(L)):>8}  {len(str(maximum)):>16}"
            f"  {str(best):>9} ({best_q:>2})  {str(minimum):>13}"
            f"  {len(str(denominator)):>24}"
        )
    print()
    print("controls: Q=14 equality sensor; Q=47/97 large banks; every phase denominator exceeds Q")
    print("result: every tested prefix is AP-identical to all its sensors but has a strict 3/41 phase")
    print("all assertions passed")


if __name__ == "__main__":
    main()
