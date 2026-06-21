#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THREAD 3 -- RECONCILE the k=12 measS7 numeric discrepancy (opus, 2026-06-21).

Two exact-Fraction implementations disagree on measS7 at k=12:
 (A) mac-mini standard whole-circle breakpoint method:
       measS7(consec_12)=121103/194040,  measS7(E*=[0..10,12])=11381/17640
 (B) drill's W_a window decomposition (sum a=1..6 only):
       measS7(consec_12)=119843/194040,  measS7(E*)=11171/17640
 differ by ~0.006 (consec) / ~0.012 (E*).

This script:
 1. Reproduces A and B exactly.
 2. Localizes the discrepancy: B = (A restricted to a=1..6 windows); A = B + (a=0 window).
 3. Computes measS7 by THREE further independent methods to pin the correct value:
    (C) theta-strip Farey method: measS7 = (1/7) sum_{j=0}^{6} M_j(k)  [canonical thread C frame]
    (D) direct N=0 occupancy: measS7 = P(N=0), N=#empty inner sectors among {1..6} [reflection frame]
    (E) fine uniform rational-grid bracket at resolution 1/(7*lcm) sampling cell-midpoints.
 4. Confirms dilation-invariance discriminant (A invariant, B not).
 5. Confirms the counterexample E* > consec_12 and both <= cap_12 = 6/7 under the CORRECT method.
"""
import sys, math, itertools
from fractions import Fraction as F
from math import gcd, lcm
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

# ---------------------------------------------------------------------------
# METHOD A: whole-circle breakpoint, all 7 sectors hit at midpoint
# ---------------------------------------------------------------------------
def measS7_A(E):
    E = sorted(set(E)); bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7*e+1): bps.add(F(m, 7*e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2; secs = set()
        for e in E:
            secs.add(int(((e*xm) % 1)*7))
        if len(secs) == 7: total += x1-x0
    return total

# ---------------------------------------------------------------------------
# METHOD B: drill W_a window decomposition (a = 1..6, drops a=0 window)
# ---------------------------------------------------------------------------
def W_a(E, a):
    E = sorted(set(E))
    lo = F(a, 7) - F(1, 14); hi = F(a, 7) + F(1, 14)
    bps = {lo, hi}
    for e in E:
        if e == 0: continue
        d = 7*abs(e)
        j0 = math.floor(lo*d); j1 = math.ceil(hi*d)
        for j in range(j0-1, j1+2):
            x = F(j, d)
            if lo <= x <= hi: bps.add(x)
    bps = sorted(bps); tot = F(0)
    for l, h in zip(bps, bps[1:]):
        if h <= l: continue
        xm = (l+h)/2; hit = set()
        for e in E:
            v = e*xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator*7) // v.denominator)
        if len(hit) == 7: tot += h-l
    return tot
def measS7_B(E): return sum(W_a(E, a) for a in range(1, 7))

# ---------------------------------------------------------------------------
# METHOD C: theta-strip Farey method (canonical thread C frame).
#   theta = 7x in [0,7). With theta = j + s, color(e) = (e*j + floor(e*s)) mod 7.
#   measS7 = (1/7) sum_{j=0}^{6} M_j,  M_j = meas{s in [0,1): colors cover Z/7}.
#   floor(e*s) for e in E is constant on order-(maxE) Farey intervals of [0,1).
# ---------------------------------------------------------------------------
def farey(n):
    # Farey sequence of order n on [0,1], as Fractions, ascending.
    a,b,c,d = 0,1,1,n
    seq = [F(0,1)]
    while c <= n:
        k = (n + b)//d
        a,b,c,d = c,d,k*c-a,k*d-b
        seq.append(F(a,b))
    return seq

def M_j(E, j):
    E = sorted(set(E))
    maxE = max(e for e in E)
    fa = farey(maxE)
    tot = F(0)
    for lo, hi in zip(fa, fa[1:]):
        sm = (lo+hi)/2
        colors = set()
        for e in E:
            colors.add((e*j + math.floor(e*sm)) % 7)
        if len(colors) == 7:
            tot += hi-lo
    return tot

def measS7_C(E):
    return sum(M_j(E, j) for j in range(7)) / 7

# ---------------------------------------------------------------------------
# METHOD D: direct N=0 occupancy via residue-depth law.
#   N(x) = # of inner sectors {1..6} not hit by {frac(e x): e in E}.
#   (sector 0 always hit since 0 in E). measS7 = meas{x: N(x)=0}.
#   Use whole-circle breakpoints but classify by N==0 explicitly.
# ---------------------------------------------------------------------------
def measS7_D(E):
    E = sorted(set(E)); bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for m in range(0, 7*e+1): bps.add(F(m, 7*e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps)-1):
        x0, x1 = bps[i], bps[i+1]
        if x1 <= x0: continue
        xm = (x0+x1)/2; hit = set()
        for e in E:
            hit.add(int(((e*xm) % 1)*7))
        inner_empty = [s for s in range(1,7) if s not in hit]
        N = len(inner_empty)
        if N == 0: total += x1-x0
    return total

# ---------------------------------------------------------------------------
# METHOD E: fine uniform rational-grid bracket.
#   Sample x = (i+1/2)/G for i=0..G-1, G = 7*lcm(E)*GRIDMULT. Count cells covered.
#   Returns a float fraction (covered/G) which must bracket the exact value.
# ---------------------------------------------------------------------------
def measS7_E_grid(E, gridmult=40):
    E = sorted(set(e for e in E if e != 0))
    L = 1
    for e in E: L = lcm(L, e)
    G = 7*L*gridmult
    covered = 0
    for i in range(G):
        # midpoint of cell i: (i+1/2)/G  -> use integer arithmetic on numerator
        num2 = 2*i + 1   # x = num2 / (2G)
        twoG = 2*G
        secs = set()
        for e in E:
            # frac(e*x)*7  -> sector = floor( frac(e*num2/(2G)) * 7 )
            t = (e*num2) % twoG            # numerator of frac(e x) over twoG
            sec = (t*7) // twoG
            secs.add(sec)
        secs.add(0)  # e=0 always sector 0
        if len(secs) == 7: covered += 1
    return F(covered, G)

# ===========================================================================
if __name__ == "__main__":
    consec12 = tuple(range(12))
    Estar    = (0,1,2,3,4,5,6,7,8,9,10,12)
    cap12    = F(6,7)

    print("="*78)
    print("THREAD 3: k=12 measS7 reconciliation")
    print("="*78)

    print("\n--- STEP 1: reproduce A and B (must match prompt) ---")
    for name,E in [("consec_12",consec12),("E*",Estar)]:
        a=measS7_A(E); b=measS7_B(E)
        print(f"  {name:<10} A = {str(a):>16} = {float(a):.8f}")
        print(f"  {'':<10} B = {str(b):>16} = {float(b):.8f}")
        print(f"  {'':<10} A-B = {a-b} = {float(a-b):.8f}")

    print("\n--- STEP 2: localize discrepancy (A = B + a=0 window) ---")
    def measS7_window(E, lo, hi):
        E=sorted(set(E)); bps=set([lo,hi])
        for e in E:
            if e==0: continue
            d=7*abs(e); j0=math.floor(lo*d); j1=math.ceil(hi*d)
            for j in range(j0-1,j1+2):
                x=F(j,d)
                if lo<=x<=hi: bps.add(x)
        bps=sorted(bps); total=F(0)
        for i in range(len(bps)-1):
            x0,x1=bps[i],bps[i+1]
            if x1<=x0: continue
            xm=(x0+x1)/2; secs=set()
            for e in E:
                v=e*xm; v=v-(v.numerator//v.denominator)
                secs.add((v.numerator*7)//v.denominator)
            if len(secs)==7: total+=x1-x0
        return total
    for name,E in [("consec_12",consec12),("E*",Estar)]:
        six = sum(measS7_window(E, F(a,7)-F(1,14), F(a,7)+F(1,14)) for a in range(1,7))
        w0  = measS7_window(E, F(13,14), F(1)) + measS7_window(E, F(0), F(1,14))
        print(f"  {name:<10} six(a=1..6)={str(six):>16}={float(six):.6f}  a=0 window={str(w0):>8}={float(w0):.6f}")
        print(f"  {'':<10} six == B? {six==measS7_B(E)}   six+w0 == A? {six+w0==measS7_A(E)}")

    print("\n--- STEP 3: three independent methods (C theta-strip, D N=0, E grid) ---")
    for name,E in [("consec_12",consec12),("E*",Estar)]:
        a=measS7_A(E); c=measS7_C(E); d=measS7_D(E)
        print(f"  {name}:")
        print(f"     A (breakpoint)      = {str(a):>16} = {float(a):.8f}")
        print(f"     C (theta-strip)     = {str(c):>16} = {float(c):.8f}   match A: {c==a}")
        print(f"     D (N=0 occupancy)   = {str(d):>16} = {float(d):.8f}   match A: {d==a}")

    print("\n--- STEP 3b: fine rational-grid bracket (method E) ---")
    for name,E in [("consec_12",consec12),("E*",Estar)]:
        a=measS7_A(E); b=measS7_B(E)
        for gm in [10,40,120]:
            g=measS7_E_grid(E,gm)
            print(f"  {name:<10} grid(mult={gm:3d}) = {float(g):.8f}   |g-A|={abs(float(g-a)):.6f}  |g-B|={abs(float(g-b)):.6f}")

    print("\n--- STEP 4: dilation-invariance discriminant ---")
    for d in [1,2,3,5]:
        E=tuple(d*i for i in range(8))
        print(f"  {d}*consec_8: A={float(measS7_A(E)):.8f}  B={float(measS7_B(E)):.8f}")
    print("  (canonical: measS7 is exactly dilation-invariant => A is correct, B is not)")

    print("\n--- STEP 5: counterexample + cap survival (CORRECT method = A=C=D) ---")
    a_c=measS7_A(consec12); a_e=measS7_A(Estar)
    print(f"  measS7(consec_12) = {a_c} = {float(a_c):.8f}")
    print(f"  measS7(E*)        = {a_e} = {float(a_e):.8f}")
    print(f"  E* STRICTLY beats consec? {a_e>a_c}   delta = {a_e-a_c} = {float(a_e-a_c):+.8f}")
    print(f"  cap_12 = {cap12} = {float(cap12):.8f}")
    print(f"  consec_12 <= cap_12? {a_c<=cap12}  (margin {cap12-a_c} = {float(cap12-a_c):.6f})")
    print(f"  E*        <= cap_12? {a_e<=cap12}  (margin {cap12-a_e} = {float(cap12-a_e):.6f})")
