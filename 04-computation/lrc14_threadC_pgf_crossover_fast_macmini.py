#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD C (mac-mini) -- pgf crossover z*, FAST (float pre-screen + exact confirm).

For each competitor E define D_E(z) = pgf_consec(z) - pgf_E(z), a degree<=7 poly
with D_E(0) = p0_consec - p0_E >= 0.  The crossover root r(E) = smallest z in (0,1]
with D_E(z)=0.  z*(k) = min_E r(E).  consec maximizes E[z^N] on [0, z*(k)].

Fast: evaluate D_E on a coarse float grid to find competitors that drop earliest;
then take the global-min candidate and locate its first root EXACTLY by bisection
in Fractions (only on the few binding competitors).
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

def first_root_float(Df):
    """Df: float coeffs of D(z). Return smallest z in (0,1] where D<=0, else 2.0."""
    G=400
    prev=Df[0]
    if prev<-1e-15: return 0.0
    for i in range(1,G+1):
        z=i/G
        v=sum(c*z**n for n,c in enumerate(Df))
        if v<=0: return z
    return 2.0

def first_root_exact(D):
    """D: Fraction coeffs, D(0)>=0. Bracket on grid then bisect exactly."""
    def val(z): return sum(c*(z**n) for n,c in enumerate(D))
    G=400; prev_z=F(0)
    if val(F(0))<0: return None
    for i in range(1,G+1):
        z=F(i,G); v=val(z)
        if v<=0:
            lo,hi=prev_z,z
            for _ in range(80):
                mid=(lo+hi)/2
                if val(mid)>0: lo=mid
                else: hi=mid
            return (lo,hi)
        prev_z=z
    return None  # never caught on (0,1]

def run(k,W,topm=25):
    C=consec(k); nlC=N_law(C)
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),k-1)]
    bank=[E for E in bank if primitive(E) and list(E)!=C]
    cand=[]; never=0
    p0fail=None
    for E in bank:
        nl=N_law(list(E))
        D=[nlC[n]-nl[n] for n in range(8)]
        if D[0]<0: p0fail=E; break
        Df=[float(c) for c in D]
        r=first_root_float(Df)
        if r>1.0: never+=1
        else: cand.append((r,E,D))
    if p0fail: return ('P0FAIL',p0fail)
    cand.sort(key=lambda t:t[0])
    # exact-confirm the smallest few float roots
    best=None
    for r,E,D in cand[:topm]:
        ex=first_root_exact(D)
        if ex is None: continue
        lo,hi=ex
        if best is None or hi<best[0]:
            best=(hi,lo,E)
    return dict(k=k,W=W,n=len(bank),never=never,caught=len(cand),
                zstar_hi=(best[0] if best else None),
                zstar_lo=(best[1] if best else None),
                zstarE=(best[2] if best else None))

if __name__=="__main__":
    print("="*78)
    print("THREAD C: pgf crossover z* (FAST exact-confirmed)")
    print("  consec maximizes E[z^N] for all z in [0, z*]; measS7=E[z^N]|_{z=0}")
    print("="*78)
    results=[]
    for k,W in [(8,13),(9,13),(10,12),(8,15),(9,15),(10,15)]:
        res=run(k,W)
        if isinstance(res,tuple):
            print(f"k={k}: P(N=0) FAIL at {res[1]}"); continue
        results.append(res)
        zl,zh=res['zstar_lo'],res['zstar_hi']
        print(f"\n--- k={k}, span<={W}, competitors={res['n']} ---")
        print(f"  never-caught: {res['never']}   caught: {res['caught']}")
        if zh is not None:
            print(f"  z* in [{float(zl):.6f}, {float(zh):.6f}]   binding E={res['zstarE']}")
        else:
            print(f"  consec maxes pgf on ALL (0,1]!")
    print("\n--- SUMMARY ---")
    for r in results:
        zh=r['zstar_hi']
        print(f"  k={r['k']} W={r['W']}: z*~{(float(zh) if zh else 1.0):.6f}  E={r['zstarE']}")
    vals=[float(r['zstar_hi']) for r in results if r['zstar_hi'] is not None]
    if vals:
        print(f"\n  GLOBAL min z* = {min(vals):.6f}")
    print("DONE.")
