#!/usr/bin/env python3
"""EXTENSION of S576 Lemma G verification (monad-compute).
The S576 script verified  G(v) := mu(safe set of S=S'u{v}) == Phi(C)  exactly
for EVEN n = 6,8,10,12,14.  This script extends the exact check to:
  - ODD n (7,9,11,13,15)  -- previously untested
  - LARGER n (16,18,20)   -- previously untested
to confirm the circuit-to-gap functional identity holds for all n, not just
small even n.  Same definitions as lrc_circuit_to_gap_functional_s576.py.
monad-compute-2026-06-03."""
from fractions import Fraction as F
from math import gcd, floor, ceil
import random

def dist(x):
    x %= 1
    return min(x, 1-x)

def G_components(Sp, n):
    THR = F(1, n); eps = set([F(0)])
    for u in Sp:
        for k in range(u+1):
            for s in (-1, 1):
                eps.add(F(k*n+s, n*u) % 1)
    pts = sorted(eps); comps = []; L = len(pts)
    for i in range(L):
        a = pts[i]; b = pts[(i+1) % L]
        ln = (b-a) if b > a else (b-a+1); mid = (a+ln/2) % 1
        if all(dist(u*mid) > THR for u in Sp):
            comps.append((a, ln))
    return comps

def safe_measure(V, n):
    THR = F(1, n); eps = set([F(0)])
    for v in V:
        for k in range(v+1):
            for s in (-1, 1):
                eps.add(F(k*n+s, n*v) % 1)
    pts = sorted(eps); meas = F(0); L = len(pts)
    for i in range(L):
        a = pts[i]; b = pts[(i+1) % L]
        ln = (b-a) if b > a else (b-a+1); mid = (a+ln/2) % 1
        if all(dist(v*mid) >= THR for v in V):
            meas += ln
    return meas

def interval_band_overlap(lo, hi, n):
    inv = F(1, n); tot = F(0)
    k0 = floor(lo-inv); k1 = ceil(hi+inv)
    for k in range(k0, k1+1):
        l = max(lo, k-inv); r = min(hi, k+inv)
        if r > l:
            tot += r-l
    return tot

def Phi(Sp, v, n):
    tot = F(0)
    for (lo, ln) in G_components(Sp, n):
        a = lo; b = lo+ln
        ov = interval_band_overlap(v*a, v*b, n)
        tot += v*ln - ov
    return tot/v

def prim(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v//g for v in V))

def main():
    rng = random.Random(7676)
    print("EXTENSION of S576 Lemma G: verify  G(v)=mu(safe) == Phi(C)  exactly")
    print("over multiple-of-n configs, for ODD n and LARGER n.")
    print()
    grand_ok = grand_tot = 0
    for n in [7, 9, 11, 13, 15, 16, 18, 20]:
        m = n-1; tot = 0; ok = 0; maxerr = F(0)
        attempts = 0
        target = 600 if n <= 15 else 300
        while tot < target and attempts < 60000:
            attempts += 1
            others = set(rng.sample([x for x in range(1, n+8) if x % n != 0], m-1))
            w = rng.randint(1, 3); v = n*w
            if v in others:
                continue
            V = prim(tuple(sorted(others | {v})))
            if len(V) != m or not any(x % n == 0 for x in V):
                continue
            mults = [x for x in V if x % n == 0]; vv = mults[0]
            Sp = tuple(x for x in V if x != vv)
            if not G_components(Sp, n):
                continue
            tot += 1
            direct = safe_measure(V, n); phi = Phi(Sp, vv, n)
            if direct == phi:
                ok += 1
            else:
                maxerr = max(maxerr, abs(direct-phi))
        parity = "odd " if n % 2 else "even"
        print(f"  n={n:2d} ({parity}): {tot} configs; Phi==mu(safe) exactly: {ok}/{tot}; max|err|={float(maxerr):.2e}", flush=True)
        grand_ok += ok; grand_tot += tot
    print()
    print(f"TOTAL: {grand_ok}/{grand_tot} exact matches "
          f"({'ALL EXACT' if grand_ok == grand_tot else 'MISMATCH FOUND'})")

if __name__ == '__main__':
    main()
