#!/usr/bin/env python3
"""wheels3.py: base-2 Zsigmondy wheels, three questions per period N (klein-2026-11-16-two3).
 (1) exact: does one class per prime p | 2^N - 1 (class modulo ord_p(2)) cover Z/N?
 (2) nested-dyadic criterion: write N = 2^a m, m odd.  Use the dyadic primes (orders 2^i, i <= a; always available: every
     prime factor of 2^(2^(i-1)) + 1 has order exactly 2^i) as nested classes covering all but one class modulo 2^a, then
     ask whether Z/m is covered by one class modulo oddpart(ord_p(2)) for each remaining prime.  Sufficient for (1).
 (3) Brier-compatible: cover of Z/N using only primes with EVEN order (a prime with odd order can serve one sign only).
usage: wheels3.py NMAX"""
import sys
from sympy import factorint, n_order, divisors
NMAX = int(sys.argv[1])
def wheel(N):
    out = []
    for d in divisors(N):
        if d == 1: continue
        for p in factorint(2**d - 1):
            if n_order(2, p) == d: out.append((p, d))
    return sorted(out)
def cover(N, mods, limit=20_000_000):
    full = (1 << N) - 1
    cm = {}
    for d in set(mods):
        for r in range(d):
            m = 0
            for x in range(r, N, d): m |= 1 << x
            cm[(d, r)] = m
    used = [False] * len(mods); nodes = [0]
    def rec(covered):
        nodes[0] += 1
        if covered == full: return True
        if nodes[0] > limit: raise TimeoutError
        x = (~covered & full); x = (x & -x).bit_length() - 1
        tried = set()
        for i, d in enumerate(mods):
            if used[i] or d in tried: continue
            tried.add(d); used[i] = True
            if rec(covered | cm[(d, x % d)]): return True
            used[i] = False
        return False
    if sum(1 / d for d in mods) < 1 - 1e-12: return False
    try: return rec(0)
    except TimeoutError: return None
def oddpart(d):
    while d % 2 == 0: d //= 2
    return d
print("N | 2^a m | exact cover | nested-dyadic criterion (odd wheel on Z/m) | Brier (even orders only) | odd-part multiset")
for N in range(2, NMAX + 1):
    W = wheel(N); mods = [d for _, d in W]
    a = 0; m = N
    while m % 2 == 0: m //= 2; a += 1
    ex = cover(N, mods)
    odd_mods = [oddpart(d) for d in mods if oddpart(d) > 1]
    nested = cover(m, odd_mods) if m > 1 else (a >= 1 and all(any(d == 2**i for d in mods) for i in range(1, a + 1)) and False)
    # m = 1 (N a power of two): the leftover class mod 2^a can only be covered by a second prime of some dyadic order
    if m == 1:
        nested = sum(1 for d in mods) > a   # more dyadic primes than levels (e.g. two primes of order 64)
        nested = cover(N, mods) if nested else False
    even_mods = [d for d in mods if d % 2 == 0]
    brier = cover(N, even_mods)
    mult = {}
    for d in odd_mods: mult[d] = mult.get(d, 0) + 1
    fmt = lambda v: 'YES' if v else ('timeout' if v is None else 'no')
    if ex or nested or brier or sum(1 / d for d in mods) >= 1:
        print(f"{N} | 2^{a}*{m} | {fmt(ex)} | {fmt(nested)} | {fmt(brier)} | {dict(sorted(mult.items()))}", flush=True)
