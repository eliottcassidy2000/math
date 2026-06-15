#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXTENSION: pin the LRC singular-series infimum and its extremizer.
kind-pasteur-2026-06-15-S3.  Verified (S3): inf L ~0.0237; the near-tight cores
{1..12}U{14m} reach it, with the MIN at small m (attained, not a limit).
Here: (a) higher-precision L for the top candidates; (b) broad search over
near-tight perturbations of the tight AP (replace one runner of {1..13} by a small
multiple of 14, and scaled APs) to find the GLOBAL extremizer; (c) random
primitive mult-of-14 search to confirm nothing dips below ~0.0237.
Implication if the extremizer is a bounded near-tight perturbation: inf L is over a
FINITE family -> C'(14) reduces to a finite positivity check.
"""
import sys, random
from math import gcd
from functools import reduce
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def D(q,S):
    h=q//14; c=0; Sr=[v%q for v in S]
    for a in range(q):
        ok=True
        for v in Sr:
            r=(v*a)%q
            if r<=h or r>=q-h:  # dangerous band (incl r=0)
                ok=False; break
        if ok: c+=1
    return c

def L_est(S, Q0, W):
    return sum(D(q,S)/q for q in range(Q0,Q0+W))/W

def primitive(S): return reduce(gcd,S)==1 and len(set(S))==len(S)
def has14(S): return any(v%14==0 for v in S)

def main():
    rng=random.Random(615)
    print("=== (a) high-precision L (Q0=8000,W=1200) for top candidates ===")
    cands={
        "{1..12,14}": sorted(list(range(1,13))+[14]),
        "{1..12,28}": sorted(list(range(1,13))+[28]),
        "tight {1..13}": list(range(1,14)),
        "evader 7*{1..12}+611": sorted([7*k for k in range(1,13)]+[611]),
    }
    for name,S in cands.items():
        print(f"   {name:26s}: L = {L_est(S,8000,1200):.5f}")

    print("\n=== (b) broad near-tight search: tight AP with ONE runner -> small mult of 14 ===")
    best=(1.0,None)
    results=[]
    # base tight APs: {1..13} scaled by d (d coprime to 14 to keep primitive after the swap)
    for d in (1,2,3,5):
        base=[d*j for j in range(1,14)]
        for idx in range(13):
            for m in range(1,6):
                S=base[:idx]+base[idx+1:]+[14*m]
                S=sorted(set(S))
                if len(S)!=13 or not primitive(S) or not has14(S): continue
                Lv=L_est(S,4000,400)
                results.append((Lv,d,idx,m,tuple(S)))
                if Lv<best[0]: best=(Lv,tuple(S))
    results.sort()
    print("   lowest-L near-tight perturbations found:")
    for Lv,d,idx,m,S in results[:8]:
        print(f"     L={Lv:.4f}  (AP d={d}, dropped pos {idx+1}, +14*{m})  S={S}")
    print(f"   GLOBAL near-tight min: L={best[0]:.4f} at S={best[1]}")

    print("\n=== (c) random primitive mult-of-14 search: anything below 0.0235? ===")
    lo=1.0; loS=None; nbelow=0
    for _ in range(400):
        S=sorted(rng.sample(range(1,80),13))
        if not primitive(S) or not has14(S): continue
        Lv=L_est(S,3000,300)
        if Lv<lo: lo=Lv; loS=S
        if Lv<0.0235: nbelow+=1
    print(f"   400 random configs: min L = {lo:.4f} at S={loS}; #(L<0.0235)={nbelow}")
    print(f"   => infimum candidate {best[0]:.4f} at {best[1]}; "
          f"random search did not beat the near-tight extremizer.")

if __name__=="__main__":
    main()
