#!/usr/bin/env python3
"""wheels2.py: Zsigmondy wheels in base b.  W_b(N) = the primes p not dividing b with ord_p(b) | N (the primes of b^N - 1
other than those dividing b - 1 with order 1... all primes p | b^N - 1, p not | b).  A covering of Z/N by one residue class
per prime, the class of p taken modulo ord_p(b), is a covering set of period N for k*b^n +- 1.  Exact search: cover the
smallest uncovered residue with an unused prime (branch over primes), so non-existence is decided exactly.
usage: wheels2.py base NMAX   -> one line per N: wheel size, missing moduli among divisors of N (Zsigmondy gaps),
density, cover? (YES/no), a witness (prime: residue) when it exists."""
import sys, time
from sympy import factorint, n_order, divisors
b = int(sys.argv[1]); NMAX = int(sys.argv[2])
def wheel(N):
    out = []
    for d in divisors(N):
        if d == 1: continue
        for p in factorint(b**d - 1):
            if b % p == 0: continue
            if n_order(b, p) == d: out.append((p, d))
    return sorted(out)
def cover(N, mods):
    """exact: does one class per modulus (mods, with multiplicity) cover Z/N?  returns a witness list or None."""
    full = (1 << N) - 1
    classmask = {}
    for d in set(mods):
        for r in range(d):
            m = 0
            for x in range(r, N, d): m |= 1 << x
            classmask[(d, r)] = m
    used = [False] * len(mods); choice = [None] * len(mods)
    sys.setrecursionlimit(10000)
    nodes = [0]
    def rec(covered):
        nodes[0] += 1
        if covered == full: return True
        if nodes[0] > 5_000_000: raise TimeoutError
        # smallest uncovered residue
        x = (~covered & full); x = (x & -x).bit_length() - 1
        tried = set()
        for i, d in enumerate(mods):
            if used[i] or d in tried: continue
            tried.add(d)                       # same modulus, same class: symmetric
            used[i] = True; choice[i] = x % d
            if rec(covered | classmask[(d, x % d)]): return True
            used[i] = False; choice[i] = None
        return False
    try:
        ok = rec(0)
    except TimeoutError:
        return 'timeout'
    return [(mods[i], choice[i]) for i in range(len(mods)) if used[i]] if ok else None
print(f"base {b}: N | #primes | moduli (order:multiplicity) | missing divisors of N (no prime of that order) | sum 1/ord | cover | witness (order:class)")
for N in range(2, NMAX + 1):
    t = time.time()
    try:
        W = wheel(N)
    except Exception as e:
        print(f"{N} | factorisation failed: {e}"); continue
    mods = [d for _, d in W]
    missing = [d for d in divisors(N) if d > 1 and d not in mods]
    dens = sum(1 / d for d in mods)
    mult = {}
    for d in mods: mult[d] = mult.get(d, 0) + 1
    if dens < 1 - 1e-12:
        res = 'no(density)'; wit = ''
    else:
        w = cover(N, mods)
        res = 'timeout' if w == 'timeout' else ('YES' if w else 'no')
        wit = '' if not isinstance(w, list) else ' '.join(f"{d}:{r}" for d, r in sorted(w))
    print(f"{N} | {len(W)} | {' '.join(f'{d}:{m}' for d, m in sorted(mult.items()))} | {missing} | {dens:.3f} | {res} | {wit}  ({time.time()-t:.1f}s)", flush=True)
