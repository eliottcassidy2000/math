#!/usr/bin/env python3
"""
lrc14_fourier_z7_consec_vs_rival_macmini_0620.py  (mac-mini-2026-06-20)

ANGLE B, decisive comparison at k=8.

We established (companion script lrc14_fourier_z7_resonance):
  measS7(E) = sum_{a in (Z/7)^k, 7 | sum a_e e}  Khat(a) J(E,a)         (*)
  - the support {a : 7 | sum a_e e} (relation lattice mod 7) depends ONLY on the
    residue-multiplicity vector rho(E)=(#{e: e=r mod7})_r;
  - Khat(a) depends ONLY on the multiset of a-values;
  - measS7 is NOT a function of rho alone -> the integer spacing enters through J(E,a).
  - consec_k uniquely has the MOST BALANCED rho (max-min<=1).

This script DECOMPOSES the gap measS7(consec8) - measS7(rival8) inside (*):
  (1) split the signed sum by Fourier weight t = #nonzero coords of a;
  (2) split by the residue-set touched by the support of a;
  (3) identify exactly which a-modes carry consec's advantage.

GOAL: locate a clean SIGNED sub-identity where consec dominates, or prove (honestly)
that the advantage is irreducibly aggregate (the C1 dead-end) even in Fourier space.

stdlib only. Exact Fractions for measS7; complex for the Fourier split (cross-checked
against the exact value).
"""
import sys, itertools, cmath, math
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)
w = cmath.exp(2j*math.pi/7)

def measS7_geom(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E: secs.add(int((e*xm % 1)*7))
        if len(secs)==7: total+=x1-x0
    return total

def J_exact(E, avec):
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=0j
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; val=1+0j
        for e,a in zip(E,avec):
            if a==0: continue
            val *= w**(a*int((e*xm % 1)*7))
        total += val*float(x1-x0)
    return total

def cM(M,a): return sum(w**(-a*j) for j in range(7) if j not in M)/7
_kc={}
def Khat(avec):
    key=tuple(sorted(avec))
    if key in _kc: return _kc[key]
    tot=0j
    for r in range(8):
        for M in itertools.combinations(range(7),r):
            p=1+0j
            for a in avec:
                p*=cM(M,a)
                if abs(p)<1e-18: break
            tot+=((-1)**r)*p
    _kc[key]=tot
    return tot

def fourier_split(E):
    """Return (total, by_weight dict, by_supportsize dict) of the signed sum (*)."""
    k=len(E)
    byw=defaultdict(float); bys=defaultdict(float); tot=0.0
    for r in range(0,k+1):
        for S in itertools.combinations(range(k), r):
            if r==0:
                a=[0]*k; c=(Khat(a)*J_exact(E,a)).real
                byw[0]+=c; bys[0]+=c; tot+=c; continue
            for vals in itertools.product(range(1,7), repeat=r):
                a=[0]*k
                for idx,v in zip(S,vals): a[idx]=v
                if sum(x*e for x,e in zip(a,E))%7!=0: continue
                c=(Khat(a)*J_exact(E,a)).real
                if abs(c)<1e-15: continue
                byw[r]+=c; bys[r]+=c; tot+=c
    return tot, byw

if __name__=="__main__":
    print("#"*84)
    print("# Fourier-on-Z/7: consec8 vs strongest rival -- where does consec win?")
    print("#"*84)

    consec8=list(range(8))
    # strongest known non-consec competitor at k=8 (from crux: residual non-consec max 0.2736)
    # plus the runner-up consec-like shapes; pick a few rivals.
    rivals=[
        ("consec8 ", consec8),
        ("drop1_8 ", [0,1,2,3,4,5,6,8]),     # near-consec, one gap
        ("drop2_8 ", [0,1,2,3,4,5,7,8]),
        ("AP2step  ", [0,1,2,3,4,5,6,7]),    # ==consec, sanity dup
        ("spread8  ", [0,1,2,3,4,5,6,9]),
        ("twoblk8  ", [0,1,2,3,7,8,9,10]),
    ]
    print("\nmeasS7 exact (rivals):")
    for lab,E in rivals:
        print(f"  {lab}: {measS7_geom(E)} = {float(measS7_geom(E)):.6f}")

    print("\n"+"="*84)
    print("Fourier weight split of (*): contrib_t for t = #nonzero a-coords")
    print("="*84)
    splits={}
    for lab,E in [("consec8 ",consec8),("drop1_8 ",[0,1,2,3,4,5,6,8]),
                  ("spread8 ",[0,1,2,3,4,5,6,9])]:
        tot,byw=fourier_split(E)
        splits[lab]=byw
        s7=float(measS7_geom(E))
        prof=" ".join(f"t{t}:{byw.get(t,0.0):+.5f}" for t in range(9) if abs(byw.get(t,0.0))>1e-9)
        print(f"  {lab} measS7={s7:.6f} recon={tot:.6f}")
        print(f"      [{prof}]")

    print("\n"+"="*84)
    print("Per-weight ADVANTAGE of consec8 over rivals (consec - rival):")
    print("="*84)
    cons=splits["consec8 "]
    for lab in ["drop1_8 ","spread8 "]:
        riv=splits[lab]
        print(f"  consec - {lab}:")
        for t in range(9):
            d=cons.get(t,0.0)-riv.get(t,0.0)
            if abs(d)>1e-9:
                print(f"      t={t}: {d:+.6f}")
        tot_d=sum(cons.get(t,0.0)-riv.get(t,0.0) for t in range(9))
        print(f"      TOTAL: {tot_d:+.6f}  (sign of every weight level? -> per-weight extremality?)")
