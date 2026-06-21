#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C (mac-mini) -- pin down the pgf crossover.

Finding from lrc14_threadC_Nlaw_orders: consec MAXIMIZES the prob. gen. fn
  pgf_E(z) = E[z^{N_E}] = sum_n P(N_E=n) z^n
for SMALL z, and LOSES near z=1 (right-tail dominated).  Since
  measS7(E) = P(N_E=0) = pgf_E(0),
the relevant regime is z near 0.  Define the (closed) consec-domination interval
  I(E) = { z in [0,1] : pgf_consec(z) >= pgf_E(z) }.
We want z* = inf over competitors E of sup{ z : [0,z] subset I(E) }, i.e. the
largest z* such that consec maximizes pgf on [0,z*] over ALL E.

HYP-C1 (falsifiable): there is a UNIVERSAL z* in (0,1) (the same across k=8,9,10)
such that consec maximizes pgf(z) for all z in [0,z*].  If z*>0 cleanly, this is
a NEW one-sided certificate that DOMINATES P(N=0)-max (the z=0 endpoint).

We compute EXACTLY (Fractions): for each competitor E, the polynomial
  D_E(z) = pgf_consec(z) - pgf_E(z) = sum_n (P_consec(N=n)-P_E(N=n)) z^n.
D_E(0) = p0_consec - p0_E >= 0 (VERIFIED).  We find the SMALLEST positive root of
D_E in (0,1] -- the first z where E catches up to consec.  z*(k) = min over E of
that first root.  (If D_E>0 on all (0,1], first root = 1 or none -> never caught.)
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce

def occupancy_pi(E):
    E=sorted(set(E))
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1)
    pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)]+=hi-lo
    return pi

def N_law(E):
    pi=occupancy_pi(E); nl=[F(0)]*8
    for h in range(8): nl[7-h]+=pi[h]
    return nl

def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

def poly_val(coeffs, z):  # coeffs[n] for z^n, z rational
    return sum(c*(z**n) for n,c in enumerate(coeffs))

def smallest_pos_root_upper(coeffs):
    """Return the smallest z in (0,1] where poly D(z) becomes <=0, via exact
    bisection bracketing on a fine rational grid then exact Sturm-free bisection.
    Returns (z_lo, z_hi) bracketing the first sign change, or None if D>0 on (0,1].
    coeffs: D_E(z), D_E(0)>=0 assumed (consec wins at 0)."""
    # if D>0 at all grid points incl 1, declare never-caught (return ('>1',))
    # else bracket first sign change then bisect to tight rational interval.
    N=2000
    prev_z=F(0); prev=poly_val(coeffs,prev_z)
    # ensure D(0)>=0
    if prev<0: return ('neg_at_0', prev)
    for i in range(1,N+1):
        z=F(i,N); v=poly_val(coeffs,z)
        if v<=0:
            lo,hi=prev_z,z
            for _ in range(60):
                mid=(lo+hi)/2
                if poly_val(coeffs,mid)>0: lo=mid
                else: hi=mid
            return ('root',(lo,hi))
        prev_z,prev=z,v
    return ('never',)  # D>0 on (0,1]

def run(k,W):
    C=consec(k); nlC=N_law(C)
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
    bank=[E for E in bank if primitive(E) and list((0,)+tuple(sorted(set(E)))[1:])!=None]
    bank=[E for E in bank if list(E)!=C]
    zstar=F(1); zstarE=None
    never=0; caught=0
    for E in bank:
        nl=N_law(list(E))
        D=[nlC[n]-nl[n] for n in range(8)]
        r=smallest_pos_root_upper(D)
        if r[0]=='never': never+=1; continue
        if r[0]=='neg_at_0':
            return ('P0FAIL',E)
        caught+=1
        lo,hi=r[1]
        if hi<zstar:
            zstar=hi; zstarE=list(E); zstar_lo=lo
    return dict(k=k,W=W,n=len(bank),never=never,caught=caught,
                zstar_hi=zstar,zstar_lo=(zstar_lo if zstarE else None),zstarE=zstarE)

if __name__=="__main__":
    print("="*78)
    print("THREAD C: pgf crossover z* -- largest z s.t. consec maxes E[z^N] on [0,z*]")
    print("="*78)
    results=[]
    for k,W in [(8,13),(9,13),(10,12),(8,15),(9,15)]:
        res=run(k,W)
        if isinstance(res,tuple):
            print(f"k={k}: P(N=0) FAIL at {res[1]} -- abort"); continue
        results.append(res)
        print(f"\n--- k={k}, span<={W}, competitors={res['n']} ---")
        print(f"  never-caught (consec wins pgf on ALL of (0,1]): {res['never']}")
        print(f"  caught (E overtakes somewhere in (0,1]):        {res['caught']}")
        zl,zh=res['zstar_lo'],res['zstar_hi']
        print(f"  z*  (binding crossover) in [{float(zl):.6f}, {float(zh):.6f}]")
        print(f"  binding competitor E = {res['zstarE']}")
    print("\n--- SUMMARY: binding z* across k ---")
    for r in results:
        print(f"  k={r['k']} W={r['W']}: z* ~ {float(r['zstar_hi']):.6f}   (E={r['zstarE']})")
    if results:
        gmin=min(float(r['zstar_hi']) for r in results)
        print(f"\n  GLOBAL min z* over tested (k,W) = {gmin:.6f}")
        print(f"  => HYP-C1 candidate: consec maximizes E[z^N] for all z in [0, {gmin:.4f}].")
    print("DONE.")
