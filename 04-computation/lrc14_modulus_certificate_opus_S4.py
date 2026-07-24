#!/usr/bin/env python3
"""
lrc14_modulus_certificate_opus_S4.py    opus-2026-07-24-S4

A STRUCTURAL (non-counting) toolkit for the defect >= 7 wall of HYP-9024, where klein's
covering lemma dies (2*k*h >= 1) and every measure/counting argument is vacuous.

MODULUS CERTIFICATE.  For any modulus q and any a, tau = a/q gives ||v tau|| = dist(va mod q)/q,
hence
        gap(V)  >=  m_q(V)/q,      m_q(V) := max_{a=1..q-1} min_{v in V} dist(va mod q).
This is elementary and EXACT, and -- crucially -- m_q(V) depends only on the speeds MOD q.
So it is completely immune to the far speeds being unbounded, which is exactly what defeats
the counting route at d >= 7.

DIVISIBILITY COROLLARY (q | no speed => m_q >= 1):
        if some q divides NO speed of V, then gap(V) >= 1/q.
Since 1/13 > 3/41, a near-tight config (gap <= 3/41) must have EVERY q in {2,...,13} dividing
some speed; equivalently it must contain multiples of 7, 8, 9, 10, 11, 12 and 13.
(Verified: AP, GW and {1..11,13,36} all satisfy this; it is a genuine necessary condition.)
"""
import numpy as np
from fractions import Fraction as Fr
TH = Fr(3, 41)

def m_q(V, q):
    best = 0
    for a in range(1, q):
        m = min(min((v * a) % q, q - (v * a) % q) for v in V)
        if m > best: best = m
    return best

def modulus_bound(V, QMAX=60):
    """best certified lower bound on gap(V) from moduli q <= QMAX, and the witnessing q."""
    best, bq = Fr(0), None
    for q in range(2, QMAX + 1):
        val = Fr(m_q(V, q), q)
        if val > best: best, bq = val, q
    return best, bq

def missing_divisors(V):
    return [q for q in range(2, 14) if not any(v % q == 0 for v in V)]

if __name__ == "__main__":
    print("MODULUS CERTIFICATE  gap(V) >= m_q(V)/q   (depends only on V mod q)")
    print("=" * 74)
    print("near-tight configs must NOT be certified by any q (their gap really is <= 3/41):")
    for nm, V in [("AP {1..13}", list(range(1, 14))),
                  ("GW {1..11,13,24}", list(range(1, 12)) + [13, 24]),
                  ("{1..11,13,36}", list(range(1, 12)) + [13, 36])]:
        b, q = modulus_bound(V)
        print(f"   {nm:20s} best bound {str(b):>8} (q={q})  <= 3/41 ? {b <= TH}   "
              f"missing divisors {missing_divisors(V)}")
    print()
    print("DIVISIBILITY COROLLARY: gap<=3/41 forces multiples of 7,8,9,10,11,12,13 to be present.")
    print("  AP has them as themselves; GW needs 24 for q=12; {1..11,13,36} needs 36 for q=12.")
    print()
    print("defect>=7 witnesses (all certified by a SMALL q; none comes near 3/41):")
    for V in [[1,2,3,4,5,6,14,16,18,20,22,24,26],
              [1,2,3,4,5,6,21,24,27,30,33,36,39],
              [1,4,6,7,10,11,17,21,32,56,80,103,171]]:
        b, q = modulus_bound(V)
        print(f"   {str(V):46s} bound {str(b):>8} (q={q})  > 3/41 ? {b > TH}")
