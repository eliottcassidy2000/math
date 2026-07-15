#!/usr/bin/env python3
"""
tutte_staircase_referee_kps_S128c8.py
=====================================
kind-pasteur-2026-07-15-S128 (cont.8).  REFEREE for THM-805 + the LRC harmonic-measure bridge.

(1) closed form T(N_n;x,y) = prod_{k=1}^{n-2} (x-1+[k]_y)  vs brute-force subset expansion, n=4..7
(2) specialization table
(3) THE HARMONIC LADDER: m({1..k}; lambda=1/(k+2)) ?= H_k / C(k+2,2)  EXACT, k=2..12
(4) Sudler landscape sanity: |T(N_n;1,e(2pi i theta))| minima vs the good set of {1..n-2} (float)
"""
import sys
from fractions import Fraction as F
from math import comb, sin, pi, log
from collections import defaultdict

sys.stdout.reconfigure(line_buffering=True)

def poly_mul(A, B):
    C = defaultdict(int)
    for (i, j), a in A.items():
        for (k, l), b in B.items():
            C[(i + k, j + l)] += a * b
    return dict(C)

def poly_add(A, B):
    C = defaultdict(int, A)
    for k, v in B.items():
        C[k] += v
    return {k: v for k, v in C.items() if v != 0}

def closed_form(n):
    T = {(0, 0): 1}
    for k in range(1, n - 1):
        fac = {(1, 0): 1, (0, 0): -1}
        for i in range(k):
            fac = poly_add(fac, {(0, i): 1})
        T = poly_mul(T, fac)
    return {k: v for k, v in T.items() if v != 0}

def brute_tutte(n):
    """subset expansion on N_n: bundles k=1..n-2 between consecutive path nodes."""
    bundles = list(range(1, n - 1))
    E = []
    for bi, k in enumerate(bundles):
        E += [bi] * k
    m = len(E)
    r_full = len(bundles)
    T = defaultdict(int)
    for S in range(1 << m):
        cnt = defaultdict(int)
        na = 0
        for e in range(m):
            if (S >> e) & 1:
                cnt[E[e]] += 1
                na += 1
        rA = sum(1 for b in cnt if cnt[b] > 0)
        T[(r_full - rA, na - rA)] += 1
    # T(x,y) = sum (x-1)^(r-rA) (y-1)^(|A|-rA): expand into monomials
    P = {}
    for (a, b), c in T.items():
        # (x-1)^a (y-1)^b * c
        term = {(0, 0): c}
        xa = {(1, 0): 1, (0, 0): -1}
        yb = {(0, 1): 1, (0, 0): -1}
        for _ in range(a):
            term = poly_mul(term, xa)
        for _ in range(b):
            term = poly_mul(term, yb)
        P = poly_add(P, term)
    return {k: v for k, v in P.items() if v != 0}

def peval(P, x, y):
    return sum(c * x**i * y**j for (i, j), c in P.items())

print("=" * 92)
print("(1)+(2) closed form vs brute force; specializations")
print("=" * 92)
for n in range(4, 8):
    A = closed_form(n)
    B = brute_tutte(n)
    ok = A == B
    kap = peval(A, 1, 1)
    forests = peval(A, 2, 1)
    connsub = peval(A, 1, 2)
    acyc = peval(A, 2, 0)
    # chromatic P(q) = (-1)^r q^c T(1-q,0) at q=3 as spot: (n-1 nodes path -> 3*2^(n-2))
    q = 3
    chrom = (-1) ** (n - 2) * q * peval(A, 1 - q, 0)
    print("n=%d closed==brute: %s | kappa=%d=(n-2)!:%s forests=%d=(n-1)!:%s conn=%d acyc=%d=2^(n-2):%s chrom(3)=%d=3*2^%d:%s"
          % (n, ok, kap, kap == 1 if n == 3 else "OK" if kap == __import__('math').factorial(n - 2) else "FAIL",
             forests, "OK" if forests == __import__('math').factorial(n - 1) else "FAIL",
             connsub, acyc, "OK" if acyc == 2 ** (n - 2) else "FAIL",
             chrom, n - 2, "OK" if chrom == 3 * 2 ** (n - 2) else "FAIL"))
    assert ok and acyc == 2 ** (n - 2)

print()
print("=" * 92)
print("(3) THE HARMONIC LADDER: m({1..k}; 1/(k+2)) vs H_k / C(k+2,2), exact")
print("=" * 92)

def good_measure(speeds, lam):
    pieces = []
    for w in speeds:
        r = lam / w
        for k in range(w):
            c = F(k, w)
            lo, hi = c - r, c + r
            if lo < 0:
                pieces.append((F(0), hi)); pieces.append((lo + 1, F(1)))
            elif hi > 1:
                pieces.append((lo, F(1))); pieces.append((F(0), hi - 1))
            else:
                pieces.append((lo, hi))
    pieces.sort()
    tot = F(0)
    cur = F(0)
    for lo, hi in pieces:
        if lo > cur:
            tot += lo - cur
        cur = max(cur, hi)
    if cur < 1:
        tot += 1 - cur
    return tot

all_ok = True
for k in range(2, 13):
    lam = F(1, k + 2)
    m_exact = good_measure(list(range(1, k + 1)), lam)
    Hk = sum(F(1, j) for j in range(1, k + 1))
    pred = Hk / comb(k + 2, 2)
    ok = m_exact == pred
    all_ok &= ok
    print("k=%2d  m({1..k};1/%d) = %-18s  H_k/C(k+2,2) = %-18s  %s"
          % (k, k + 2, m_exact, pred, "EQUAL" if ok else "** DIFFER ratio=%.6f" % float(m_exact / pred)))
print("HARMONIC LADDER: %s" % ("CONFIRMED EXACTLY k=2..12" if all_ok else "NOT the identity -- see ratios"))

print()
print("=" * 92)
print("(4) Sudler landscape (float sanity): argmax of |T(N_n;1,e)| vs good set of {1..n-2}")
print("=" * 92)
for n in [7, 14]:
    K = n - 2
    lam = 1.0 / n
    best = []
    for j in range(1, 5000):
        th = j / 10000.0   # (0, 1/2)
        s = 0.0
        okg = True
        for k in range(1, K + 1):
            v = abs(sin(pi * k * th))
            s += log(max(v, 1e-300))
            d = abs(k * th - round(k * th))
            if d < lam:
                okg = False
        s -= K * log(abs(sin(pi * th)))
        best.append((s, th, okg))
    best.sort()
    lows = best[:6]
    goods = [b for b in best if b[2]]
    print("n=%2d: |T| SMALLEST at theta=%s (good:%s); good-set thetas have |T|-rank mean %.2f of %d"
          % (n, ["%.4f" % b[1] for b in lows[:3]], [b[2] for b in lows[:3]],
             (sum(sorted(range(len(best)), key=lambda i: best[i][0]).index(best.index(g)) for g in goods[:200]) / max(1, len(goods[:200]))) if goods else -1, len(best)))
print("\nREFEREE DONE")
