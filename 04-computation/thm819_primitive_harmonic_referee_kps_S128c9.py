#!/usr/bin/env python3
"""THM-819 referee: m({1..k}; 1/(k+2)) = (2/((k+1)(k+2))) * sum_{u<=k, gcd(u,k+1)=1} 1/u, k=1..30 exact.
kind-pasteur-2026-07-15-S128 (cont.9). Output: 05-knowledge/results/thm819_primitive_harmonic_referee_kps_S128c9.out"""
from fractions import Fraction as F
from math import gcd

def is_prime(q):
    return q >= 2 and all(q % d for d in range(2, int(q ** 0.5) + 1))

def good_measure(speeds, lam):
    pieces = []
    for w in speeds:
        r = lam / w
        for k in range(w):
            c = F(k, w); lo, hi = c - r, c + r
            if lo < 0: pieces.append((F(0), hi)); pieces.append((lo + 1, F(1)))
            elif hi > 1: pieces.append((lo, F(1))); pieces.append((F(0), hi - 1))
            else: pieces.append((lo, hi))
    pieces.sort(); tot = F(0); cur = F(0)
    for lo, hi in pieces:
        if lo > cur: tot += lo - cur
        cur = max(cur, hi)
    if cur < 1: tot += 1 - cur
    return tot
allok = True
for k in range(1, 31):
    q = k + 1
    m = good_measure(list(range(1, k + 1)), F(1, k + 2))
    pred = F(2, q * (k + 2)) * sum(F(1, u) for u in range(1, k + 1) if gcd(u, q) == 1)
    ok = m == pred; allok &= ok
    kind = 'PRIME' if is_prime(q) else 'comp'
    print('k=%2d q=%2d %-5s  %s' % (k, q, kind, 'EXACT' if ok else '** FAIL'))
print('PRIMITIVE HARMONIC LAW k=1..30: ' + ('ALL EXACT' if allok else 'FAILURES'))
