#!/usr/bin/env python3
"""brier_sheets.py: self-dual covers.  A prime p with even order d = ord_p(2) divides k*2^n - 1 for n = r (mod d) and then
k*2^n + 1 for n = r + d/2 (mod d).  A single wheel serves both signs of the same k iff residues r_p (even-order primes
only) can be chosen so that both sheets {r_p mod d_p} and {r_p + d_p/2 mod d_p} cover Z/N.  Exact search over N <= NMAX.
usage: brier_sheets.py NMAX"""
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
def two_sheet(N, mods):
    full = (1 << N) - 1
    def cls(d, r):
        m = 0
        for x in range(r % d, N, d): m |= 1 << x
        return m
    used = [False] * len(mods); nodes = [0]
    def rec(cm, cp):
        nodes[0] += 1
        if cm == full and cp == full: return True
        if nodes[0] > 5_000_000: raise TimeoutError
        if cm != full:
            x = (~cm & full); x = (x & -x).bit_length() - 1; sheet = 0
        else:
            x = (~cp & full); x = (x & -x).bit_length() - 1; sheet = 1
        tried = set()
        for i, d in enumerate(mods):
            if used[i] or d in tried: continue
            tried.add(d); used[i] = True
            r = x if sheet == 0 else x - d // 2
            if rec(cm | cls(d, r), cp | cls(d, r + d // 2)): return True
            used[i] = False
        return False
    if sum(1 / d for d in mods) < 1 - 1e-12: return False
    try: return rec(0, 0)
    except TimeoutError: return None
print("N | even-order moduli | sum 1/d | self-dual two-sheet cover?")
for N in range(2, NMAX + 1):
    W = wheel(N); even = [d for _, d in W if d % 2 == 0]
    dens = sum(1 / d for d in even)
    if dens < 1 - 1e-12: continue
    res = two_sheet(N, even)
    print(f"{N} | {sorted(even)} | {dens:.3f} | {'YES' if res else ('timeout' if res is None else 'no')}", flush=True)
