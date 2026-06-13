#!/usr/bin/env python3
"""THM-479 (claudebox-2026-06-11-S2): level law + branch-split Burnside for switching
classes of tournaments (A049313). Pure python, exact F2 bitmask linear algebra.
Part A: direct verification (solvability == level law; Burnside totals) for n <= NMAX_DIRECT.
Part B: closed form to n = 16 with odd/even branch split; graph control (A002854)."""
import sys
from itertools import combinations
from math import factorial, gcd
from collections import Counter
from fractions import Fraction

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield []; return
    for p in range(min(n, mx), 0, -1):
        for rest in partitions(n - p, p): yield [p] + rest

def zmu(mu):
    z = 1
    for l, m in Counter(mu).items(): z *= (l ** m) * factorial(m)
    return z

def v2(x):
    r = 0
    while x % 2 == 0: x //= 2; r += 1
    return r

def is_level(mu):
    odd = all(l % 2 for l in mu)
    return odd or (all(l % 2 == 0 for l in mu) and len({v2(l) for l in mu}) == 1)

# ---------- Part A: direct F2 linear algebra ----------
def analyze(n, mu):
    pairs = list(combinations(range(n), 2)); idx = {p: i for i, p in enumerate(pairs)}
    m = len(pairs)
    pi = []; start = 0
    for l in mu:
        pi += [start + (i + 1) % l for i in range(l)]; start += l
    P = [0] * m; t = 0
    for k, (i, j) in enumerate(pairs):
        a, b = pi[i], pi[j]
        e2 = idx[(min(a, b), max(a, b))]
        P[k] = e2
        if a > b: t |= (1 << e2)
    cols = [(1 << P[k]) ^ (1 << k) for k in range(m)]
    def span(vecs):
        basis = []
        for v in vecs:
            for b in basis: v = min(v, v ^ b)
            if v: basis.append(v); basis.sort(reverse=True)
        return basis
    imb = span(cols); rank = len(imb); null = m - rank
    cuts = []
    for v in range(n):
        c = 0
        for k, (i, j) in enumerate(pairs):
            if i == v or j == v: c |= (1 << k)
        cuts.append(c)
    cutb = span(cuts)
    big = span(list(imb) + list(cutb))
    tt = t
    for b in big:
        if tt ^ b < tt: tt ^= b
    solvable = (tt == 0)
    dim_cap = rank + len(cutb) - len(big)
    fix = (1 << (null + dim_cap - (n - 1))) if solvable else 0
    return solvable, fix

# ---------- Part B: closed form ----------
def o2(mu):
    return sum(l // 2 for l in mu) + sum(gcd(mu[i], mu[j]) for i in range(len(mu)) for j in range(i + 1, len(mu)))

def closed(n):
    t_odd = Fraction(0); t_lev = Fraction(0); t_graph = Fraction(0)
    for mu in partitions(n):
        k = len(mu); o = o2(mu); z = zmu(mu)
        if all(l % 2 for l in mu): t_odd += Fraction(2 ** (o - k + 1), z)
        elif all(l % 2 == 0 for l in mu) and len({v2(l) for l in mu}) == 1:
            t_lev += Fraction(2 ** (o - k), z)
        ko = sum(1 for l in mu if l % 2)
        t_graph += Fraction(2 ** (o - k + (1 if ko else 0)), z)
    return t_odd, t_lev, t_graph

A049313 = [1, 1, 1, 2, 2, 6, 12, 79, 792, 19576, 886288, 75369960, 11856006240,
           3467430423264, 1893448825054528, 1938818712501985736]   # OEIS, n=1..16
A002854 = [1, 1, 2, 3, 7, 16, 54, 243, 2038, 33120, 1182004, 87723296, 12886193064,
           3633057074584, 1944000150734320, 1967881448329407496]

if __name__ == '__main__':
    NMAX_DIRECT = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    print("== Part A: direct F2 verification (law + Burnside) ==")
    for n in range(1, NMAX_DIRECT + 1):
        total = 0; viol = 0
        for mu in partitions(n):
            s, fix = analyze(n, mu)
            if s != is_level(mu): viol += 1
            total += fix * (factorial(n) // zmu(mu))
        cnt = total // factorial(n)
        print(f"n={n}: count={cnt} (OEIS {A049313[n-1]}) law-violations={viol}")
        assert viol == 0 and cnt == A049313[n - 1]
    print("\n== Part B: closed form, branch split, graph control ==")
    print("n : total | N_odd | N_lev | graphs(A002854)")
    for n in range(1, 17):
        a, b, g = closed(n)
        T = a + b
        assert T == A049313[n - 1], (n, T)
        assert g == A002854[n - 1], (n, g)
        if n >= 3: assert a.denominator == 1 and b.denominator == 1, (n, a, b)
        print(f"{n} : {T} | {a} | {b} | {g}")
    print("\nall checks passed (A049313 and A002854 reproduced; branches integral for n>=3)")
