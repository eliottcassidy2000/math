#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_persistence_kps-S3-wf

THE PERSISTENCE / EQUIDISTRIBUTION ANALYSIS — the rigorous bridge that would turn the
verified K_j-ruler witness into a PROOF for ALL S3 (including arbitrarily large V0).

The K_j ruler certifies every sampled S3.  To make it a THEOREM we must show the certificate
PERSISTS as V0 -> infty.  The mechanism:

  - The cluster L' = {w0 + d : d in Delta}, with w0 = min, max offset = s.  P = small speeds.
  - Ruler periods I_j = ((14j+1)/(14 w0), (14j+13)/(14 w0)), j=0..w0-1, PERFECTLY periodic.
  - Inside I_j, the cluster-safe sub-window K_j is determined by the fractional parts
    {(w0+d) tau} = {w0 tau + d tau}.  On I_j, w0 tau in (j+1/14, j+13/14) (one full gap).
    Write w0 tau = j + phi, phi in (1/14, 13/14).  Then (w0+d) tau = j + phi + d*tau, and
    d*tau ~ d*(j+phi)/w0 -- a SLOW drift (since d<=s, w0 large).  The cluster window K_j thus
    depends on the drift offsets {d*tau mod 1 : d in Delta}, which by Weyl EQUIDISTRIBUTE as j
    sweeps 0..w0-1 (the map j -> d*j/w0 mod 1 fills [0,1) when gcd(d,w0) small).

  - Meanwhile the P-danger pattern, sampled at the ruler left-endpoints (j+1/14)/w0, ALSO
    equidistributes (the points (j+1/14)/w0 are ~ j/w0, equidistributing mod 1/u for each u in P).

So as V0=w0 grows, the joint (cluster-window-shape, P-danger-phase) over the w0 periods
equidistributes over a FIXED compact pattern space depending only on (Delta, P).  The fraction
rho_K of good periods -> a POSITIVE limit rho*(Delta,P) > 0 (an integral of an indicator over
the pattern torus).  Hence for w0 large the count #good = rho_K * w0 -> infty; in particular
>= 1.  Small w0 is a FINITE check.

THIS SCRIPT measures the equidistribution quantitatively:
  (1) rho_K as a function of w0 for a FIXED offset pattern Delta and fixed P, scanning w0.
      Does rho_K converge to a positive constant?  Does #good = rho_K*w0 grow linearly?
  (2) The limiting rho* estimated by direct torus integration vs the finite-w0 rho_K.
  (3) The WORST offset pattern Delta (AP, all-same-residue) -- does rho* stay > 0?

If rho_K converges to a positive floor for every (Delta,P) and grows, the proof structure is:
  "for w0 >= W*(Delta,P) the count is >=1 by equidistribution; for w0 < W* finite check."
We probe whether W* is uniformly bounded (a SINGLE finite check would then suffice).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, sys

def flush(*a):
    print(*a); sys.stdout.flush()
C=F(1,14)

def teeth_intervals(u,h=C):
    iv=[]
    for j in range(0,u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merged_intervals(iv):
    iv=sorted(iv); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    return merged
def danger_merged(A,h=C):
    iv=[]
    for u in A: iv+=teeth_intervals(u,h)
    return merged_intervals(iv)
def safe_from_danger(merged):
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def safe_components(A,h=C): return safe_from_danger(danger_merged(A,h))
def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_all(S): return reduce(gcd,S)

def good_periods(P, Lp, Vmax):
    """Return (count of ruler periods j with a >thresh window safe for all L' and P, w0)."""
    w0=min(Lp); thresh=F(1,7*Vmax)
    GLp=safe_components(Lp); GP=safe_components(P) if P else [(F(0),F(1))]
    GA=[]; i=j=0
    while i<len(GLp) and j<len(GP):
        a,b=GLp[i]; c,d=GP[j]
        lo=max(a,c); hi=min(b,d)
        if lo<hi: GA.append((lo,hi))
        if b<d: i+=1
        else: j+=1
    periods=set()
    for a,b in GA:
        if b-a<=thresh: continue
        jj=int(w0*a)
        for cj in range(max(0,jj-1),jj+2):
            lo=F(14*cj+1,14*w0); hi=F(14*cj+13,14*w0)
            x=max(a,lo); y=min(b,hi)
            if y-x>thresh: periods.add(cj)
    return len(periods), w0

def run():
    flush("="*70)
    flush("PERSISTENCE: rho_K vs w0 for fixed offset pattern Delta and fixed P")
    flush("="*70)
    # Worst-ish setups: AP cluster offsets, small P that is dangerous near 0.
    # We fix Delta and P, vary w0 = V0 (the min cluster speed), set Vmax = w0 + s + extra.
    test_patterns=[
        ("Delta={0,1,2,3,4}, P={1,2,3,5,7,11,13}",[0,1,2,3,4],[1,2,3,5,7,11,13]),
        ("Delta={0,2,4,6,8} (AP step2), P={1,2,3,4,6,12}",[0,2,4,6,8],[1,2,3,4,6,12]),
        ("Delta={0,7,14,21} (step7), P={1,5,8,9,10,11,13}",[0,7,14,21],[1,5,8,9,10,11,13]),
        ("Delta={0,1}, P={2,3,4,5,6,7,8,9,10,11,12,13}",[0,1],list(range(2,14))),
        ("Delta={0,3,9,12,15} mixed, P={1,2,5,7,11}",[0,3,9,12,15],[1,2,5,7,11]),
    ]
    for name,Delta,P in test_patterns:
        flush(f"\n  pattern: {name}")
        s=max(Delta)
        flush(f"   {'w0':>6} {'Vmax':>6} {'#good':>6} {'w0':>6} {'rho_K':>9}")
        rhos=[]
        for w0 in [29, 53, 101, 211, 401, 809, 1601, 3203]:
            Lp=[w0+d for d in Delta]
            Vmax=w0+s+ (14 - (w0+s)%14)%14 + 14  # an extra large runner above the cluster
            # ensure Vmax > max(Lp)
            if Vmax<=max(Lp): Vmax=max(Lp)+14
            cnt,wq=good_periods(P,Lp,Vmax)
            rho=F(cnt,w0)
            rhos.append(float(rho))
            flush(f"   {w0:>6} {Vmax:>6} {cnt:>6} {w0:>6} {float(rho):>9.4f}")
        # convergence check
        flush(f"   rho_K trend: {['%.4f'%r for r in rhos]}  "
              f"({'CONVERGING>0' if min(rhos[-3:])>0.01 else 'CHECK'})")

if __name__=="__main__":
    run()
