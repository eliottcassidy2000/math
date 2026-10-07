#!/usr/bin/env python3
"""Exhaustive check of the local lemmas of the coalescence proof on all states |k| <= 10, |N| <= 4000, both bits.
f = e / 3^max(k,0) (|f| = |e| for k <= 0).  c = beta (k >= 0), beta xor (e mod 2) (k < 0).
 L1: |f'| <= A(c)|f| + 1/2, A(0) = 1/2, A(1) = 3/2      (all steps with k != 0, and k = 0 runs)
 L2: flips at |k| >= 1 move |k| toward 0 iff c = 1
 L3: departures from k = 0 (e odd): |f'| <= |f|/2 + 1/6
 L4: runs (e even) at level h >= 1: new v2 = w-1 if w < m; = m-1 if w > m (c=1), w-1 (c=0); >= m if w = m and c = 1; m = v2(3^h - 1)
 L5: k = 0 runs: e -> e/2 or 3e/2 exactly
"""
from fractions import Fraction as Fr
from coal_check import step
def v2(x): x = abs(x); return (x & -x).bit_length() - 1 if x else 10**9
bad = {f"L{i}": 0 for i in range(1, 6)}; n = 0
for k in range(-10, 11):
    for N in range(-4000, 4001):
        if k == 0 and N == 0: continue
        e = Fr(N, 3**max(0, -k)); f = abs(e) / 3**max(k, 0)
        sigma = N & 1
        for beta in (0, 1):
            n += 1
            k2, N2 = step(k, N, beta)
            e2 = Fr(N2, 3**max(0, -k2)); f2 = abs(e2) / 3**max(k2, 0)
            c = beta if k >= 0 else beta ^ sigma
            A = Fr(1, 2) if c == 0 else Fr(3, 2)
            if sigma == 1 and k == 0:
                if not f2 <= f / 2 + Fr(1, 6): bad["L3"] += 1
            else:
                if not f2 <= A * f + Fr(1, 2): bad["L1"] += 1
            if sigma == 1 and k != 0:
                toward = abs(k2) < abs(k)
                if toward != (c == 1): bad["L2"] += 1
            if sigma == 0 and k != 0:
                h = abs(k); m = v2(3**h - 1); w = v2(N)   # mirror: v2 of N equals v2 of e (3-power units)
                w2 = v2(N2)
                if N == 0: ok = (w2 == (m - 1 if c == 1 else w))   # e = 0 (w = inf): stays 0 on c = 0
                elif w < m: ok = (w2 == w - 1)
                elif w > m: ok = (w2 == (m - 1 if c == 1 else w - 1))
                else: ok = (w2 == m - 1) if c == 0 else (w2 >= m)
                if not ok: bad["L4"] += 1
            if sigma == 0 and k == 0:
                if not (k2 == 0 and (e2 == e / 2 if beta == 0 else e2 == 3 * e / 2)): bad["L5"] += 1
print(f"checked {n} (state, bit) pairs; failures: {bad}")
