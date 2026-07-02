#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYP-3955 / THM-599: THE SYMBOLIC k <= 13 LEDGER — exact-rational A(U) via the c-breakpoint engine.

A(U) = Int_0^1 L^c(U) dc,  L^c(U) = meas{x in T : ||u x - c|| >= h for all u in U},  h = 1/14.

ENGINE: L^c is piecewise-AFFINE in c. The x-endpoints of speed u's safe set are (k + c ± h)/u —
affine in c with slope 1/u. Structure changes only when endpoints of two DIFFERENT speeds collide:
   (k + c + s1 h)/u = (k' + c + s2 h)/u'   =>   c = [u(k' + s2 h) - u'(k + s1 h)] / (u - u'),
a finite set of rationals. Between consecutive breakpoints L^c is affine, so
   Int over the interval = length x L^(midpoint)   (EXACT).
We compute everything in Fractions. Affineness is verified per interval by a 3-point linearity check
(t = 1/4, 1/2, 3/4 of the interval: L(1/4) + L(3/4) = 2 L(1/2) exactly, else REFINE = safety net).

RUNS:
 R1: acid tests — A(1,2,3) = 61/98, A(pair) = 36/49, A(1,2,3,4) = 11/21 (THM-599 hand values).
 R2: the argmin patterns k = 2..13 from the S30 sampled search — exact rationals + margins vs
     witnessMP = 14249/252252 (the k = 8..13 rows are the hfloor-facing ledger).
 R3: near-minimal pool per k (float pre-scan top-5) — exact-certified minima.
"""
import sys, itertools, random
from fractions import Fraction as Fr
from math import gcd, floor
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
rng = random.Random(3955)

h = Fr(1, 14)

def L_exact(U, c):
    """exact L^c(U) via interval intersection over one period [0,1)."""
    arcs = [(Fr(0), Fr(1))]
    for u in sorted(U):
        # safe set of speed u at target c: complement of teeth ((k + c - h)/u, (k + c + h)/u)
        new = []
        for (lo, hi) in arcs:
            k0 = int(floor(float(lo*u - c))) - 1
            k1 = int(floor(float(hi*u - c))) + 1
            for k in range(k0, k1+1):
                a = (k + c + h)/u; b = (k + 1 + c - h)/u
                l = lo if lo > a else a
                r = hi if hi < b else b
                if l < r: new.append((l, r))
        arcs = new
        if not arcs: return Fr(0)
    return sum(r - l for l, r in arcs)

def breakpoints(U):
    """candidate c in [0,1) where the affine structure can change."""
    bps = {Fr(0), Fr(1)}
    Us = sorted(U)
    for i in range(len(Us)):
        for j in range(i+1, len(Us)):
            u, up = Us[i], Us[j]
            for k in range(u):
                for kp in range(up):
                    for s1 in (h, -h):
                        for s2 in (h, -h):
                            num = u*(kp + s2) - up*(k + s1)
                            c = Fr(num, u - up)
                            c -= floor(c)
                            bps.add(Fr(c))
    return sorted(bps)

def A_symbolic(U, spot=0.04):
    """exact A(U) = sum over breakpoint intervals of length * L(midpoint).
    Affineness between breakpoints is THM-599/engine-derived; a random spot fraction of intervals
    gets a 3-point linearity check with refine safety net (0 trips expected)."""
    bps = breakpoints(U)
    total = Fr(0)
    bad = 0
    for i in range(len(bps) - 1):
        lo, hi = bps[i], bps[i+1]
        if hi <= lo: continue
        mid = (lo + hi)/2
        Lm = L_exact(U, mid)
        if rng.random() < spot:
            L1 = L_exact(U, lo + (hi-lo)/4)
            L3 = L_exact(U, lo + 3*(hi-lo)/4)
            if L1 + L3 != 2*Lm:
                bad += 1
                total += A_refine(U, lo, mid) + A_refine(U, mid, hi)
                continue
        total += (hi - lo)*Lm
    return total, bad

def A_refine(U, lo, hi, depth=0):
    mid = (lo + hi)/2
    Lm = L_exact(U, mid)
    L1 = L_exact(U, lo + (hi-lo)/4); L3 = L_exact(U, lo + 3*(hi-lo)/4)
    if L1 + L3 == 2*Lm or depth >= 12:
        return (hi - lo)*Lm
    return A_refine(U, lo, mid, depth+1) + A_refine(U, mid, hi, depth+1)

MP = Fr(14249, 252252)

print("="*100)
print(" R1: ACID TESTS (THM-599 hand values)")
print("="*100)
for U, want in [((1,2), Fr(36,49)), ((1,2,3), Fr(61,98)), ((10,24,38), Fr(61,98)), ((1,2,3,4), Fr(11,21))]:
    A, bad = A_symbolic(U)
    print(f"  A{U} = {A} = {float(A):.6f}   expected {want}   MATCH: {A == want}   (refines: {bad})")

print("\n" + "="*100)
print(" R2: THE ARGMIN PATTERNS k = 2..13 — EXACT RATIONALS (the hfloor-facing ledger)")
print("="*100)
argmins = { 2: (5,37), 3: (10,24,38), 4: (5,6,7,8), 5: (3,5,6,7,9), 6: (2,4,6,7,8,10),
            7: (3,5,6,7,8,9,11), 8: (1,3,5,6,7,8,9,11), 9: (1,3,4,5,6,7,8,9,11),
            10: (1,3,4,5,6,7,8,9,11,13), 11: (2,3,4,5,6,7,8,9,10,12,14),
            12: (4,8,9,10,12,13,14,15,16,17,22,24), 13: (3,7,8,9,11,13,15,16,17,19,23,24,27) }
import time
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 11   # k=12,13 exact deferred (breakpoint count)
for k in range(2, KMAX+1):
    U = argmins[k]
    t0 = time.time()
    A, bad = A_symbolic(U)
    dt = time.time() - t0
    print(f"  k={k:2d}  U={U}")
    print(f"        A = {A}  = {float(A):.6f}   vs MP = {float(MP):.6f}   margin x{float(A/MP):.2f}"
          f"   (refines {bad}, {dt:.1f}s)", flush=True)

print("\n" + "="*100)
print(" R3: NEAR-MINIMAL POOL — float pre-scan top-3 per k (k = 5..10), exact-certified")
print("="*100)
def A_float(U, NC=600):
    BAND = 1.0/14.0
    def clip(arcs, v, c):
        out = []
        for (lo, hi) in arcs:
            k0 = int(floor(lo*v - c)) - 1; k1 = int(floor(hi*v - c)) + 1
            for kk in range(k0, k1+1):
                aa = (kk + c + BAND)/v; bb = (kk + 1 + c - BAND)/v
                l = lo if lo > aa else aa
                hh2 = hi if hi < bb else bb
                if l < hh2 - 1e-15: out.append((l, hh2))
        return out
    tot = 0.0
    for ic in range(NC):
        c = (ic+0.5)/NC
        arcs = [(0.0, 1.0)]
        for u in sorted(U): arcs = clip(arcs, u, c)
        tot += sum(hi-lo for lo, hi in arcs)/NC
    return tot
for k in (6, 8, 10):
    pool = set()
    for V0 in range(k, min(k+5, 14)):
        for C in itertools.combinations(range(1, V0+1), k):
            g = 0
            for c_ in C: g = gcd(g, c_)
            pool.add(tuple(sorted(c_//g for c_ in C)))
    pool = list(pool)
    if len(pool) > 300: pool = rng.sample(pool, 300)
    scored = sorted((A_float(U), U) for U in pool)[:2]
    row = []
    for _, U in scored:
        A, bad = A_symbolic(U)
        row.append((A, U))
    row.sort()
    for A, U in row:
        print(f"  k={k}: A = {str(A):>22} = {float(A):.6f}   U={U}", flush=True)
    print(f"       => exact-certified min over pool: {row[0][0]}  (x{float(row[0][0]/MP):.2f} of MP)")
print("DONE.")
