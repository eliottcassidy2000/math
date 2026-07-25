#!/usr/bin/env python3
"""Exact regression for the symbolic fixed-bank theorem THM-2161.

For Q>=14, the theorem uses

    L = lcm(2,...,Q),
    V_Q = {1,...,6} union {L+7,...,L+13}.
    W_Q = {1,...,6,L} union {L+8,...,L+13}.

Every labelled residue modulo q<=Q agrees with the tight AP {1,...,13}, but
the phase x=(L+8)/(8(L+14)) has minimum distance x>3/41.  The covering row
W_Q contains L, so every scalar certificate m_q/q is zero throughout the bank,
but the same phase has the same minimum.  This script checks both consequence
objects by exact integer/Fraction arithmetic on several hostile prefix banks.
The replay is a regression, not the proof of the all-Q quantifier, which is
algebraic in THM-2161.
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


def construct(bound: int) -> tuple[int, tuple[int, ...], tuple[int, ...], QF]:
    assert bound >= 14
    L = lcm(*range(2, bound + 1))
    mimic = tuple(range(1, 7)) + tuple(L + j for j in range(7, 14))
    covering = tuple(range(1, 7)) + (L,) + tuple(L + j for j in range(8, 14))
    phase = QF(L + 8, 8 * (L + 14))
    return L, mimic, covering, phase


def audit(bound: int) -> tuple[int, int, QF, QF, int, QF, int]:
    L, mimic, covering, phase = construct(bound)

    # The exact definition used in THM-2146/THM-2161 is
    # defect(V)=|V\\{1,...,13}|, not distance from the closest affine AP.
    for speeds in (mimic, covering):
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
        assert all(mimic[i - 1] % q == i % q for i in range(1, 14))
        assert Counter(v % q for v in mimic) == Counter(i % q for i in AP)
        actual = QF(modulus_numerator(mimic, q), q)
        control = QF(modulus_numerator(AP, q), q)
        assert actual == control
        assert control <= QF(1, 14) < H
        assert modulus_numerator(covering, q) == 0
        if actual > best:
            best, best_q = actual, q
    assert best == QF(1, 14)
    assert best_q == 14

    mimic_distances = tuple(circle_distance(v * phase) for v in mimic)
    covering_distances = tuple(circle_distance(v * phase) for v in covering)
    mimic_minimum = min(mimic_distances)
    covering_minimum = min(covering_distances)
    assert mimic_minimum == covering_minimum == phase
    assert mimic_distances[0] == mimic_distances[-1] == phase
    assert covering_distances[0] == covering_distances[-1] == phase
    assert circle_distance(L * phase) > QF(1, 4) > phase > H

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

    return (
        L,
        mimic[-1],
        mimic_minimum,
        covering_minimum,
        reduced_denominator,
        best,
        best_q,
    )


def main() -> None:
    bounds = (14, 15, 16, 17, 23, 41, 47, 64, 97)
    print("THM-2161 FIXED-PREFIX MODULUS BLINDNESS REGRESSION")
    print("V_Q={1,...,6,L+7,...,L+13}; L=lcm(2,...,Q); threshold h=3/41")
    print("W_Q={1,...,6,L,L+8,...,L+13}; every scalar bank certificate is zero")
    print("all arithmetic is exact; the all-Q proof is symbolic in THM-2161")
    print()
    print(" Q  L digits  max-speed digits  V bank best (q)  V/W phase minimum  phase-denominator digits")
    for bound in bounds:
        L, maximum, minimum, cover_minimum, denominator, best, best_q = audit(bound)
        assert minimum == cover_minimum
        print(
            f"{bound:>2}  {len(str(L)):>8}  {len(str(maximum)):>16}"
            f"  {str(best):>9} ({best_q:>2})  {str(minimum):>13}"
            f"  {len(str(denominator)):>24}"
        )
    print()
    print("controls: Q=14 equality sensor; Q=47/97 large banks; every phase denominator exceeds Q")
    print("result: V is AP-identical; W is covering with zero scalar certificates; both escape 3/41")
    print("all assertions passed")


if __name__ == "__main__":
    main()
