#!/usr/bin/env python3
"""Coloured twisted Burnside: k-arc-coloured tournaments (THM-4720's T_{2k}(n)) versus even (2k)-coloured complete graphs.
A (2k)-colouring c: pairs -> Z/2k; an automorphism g is odd iff the number of edges reversed by g (relative to the order
orientation) that carry an ODD colour is odd. Brute force over types for n <= 4; Burnside formula for comparison."""
import itertools
from math import gcd, factorial
from fractions import Fraction
from collections import Counter
def partitions(n, m=None):
    if m is None: m = n
    if n == 0: yield (); return
    for k in range(min(n, m), 0, -1):
        for r in partitions(n - k, k): yield (k,) + r
def zlam(p):
    c = Counter(p); z = 1
    for L, m in c.items(): z *= L ** m * factorial(m)
    return z
def T_formula(K, n):
    tot = Fraction(0)
    for lam in partitions(n):
        if any(L % 2 == 0 for L in lam): continue
        c2 = sum((L - 1) // 2 for L in lam) + sum(gcd(lam[i], lam[j]) for i in range(len(lam)) for j in range(i + 1, len(lam)))
        tot += Fraction(K ** c2, zlam(lam))
    assert tot.denominator == 1
    return int(tot)
def brute(K, n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    idx = {p: t for t, p in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    seen = set(); count = 0
    for col in itertools.product(range(K), repeat=len(pairs)):
        if col in seen: continue
        orbit = set(); even = True
        for g in perms:
            img = [0] * len(pairs); rev_odd = 0
            for t, (u, v) in enumerate(pairs):
                a, b = g[u], g[v]
                q = (min(a, b), max(a, b)); img[idx[q]] = col[t]
            img = tuple(img); orbit.add(img)
            if img == col:  # automorphism: count reversed edges with odd colour
                s = 0
                for t, (u, v) in enumerate(pairs):
                    if g[u] > g[v] and col[t] % 2 == 1: s ^= 1
                if s: even = False
        seen |= orbit
        if even: count += 1
    return count
for K in [2, 4, 6]:
    for n in range(2, 5):
        b = brute(K, n); f = T_formula(K, n)
        print(f"K = {K} colours (k = {K//2} arc colours), n = {n}: even {K}-coloured complete graphs {b}, T_{K}({n}) = {f} {'OK' if b == f else 'MISMATCH'}")
