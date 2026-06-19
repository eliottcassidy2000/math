#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) crux: is "consecutive minimizes the 1/7-empty-window functional" provable by a
GAP-CONTRACTION monotonicity (shrink a sorted-offset gap by 1 ⟹ functional non-increasing)?
This is a ROUTE-VALIDITY test (the agents verified consec-extremal exhaustively; what's missing
is a PROOF, and the gap-contraction lemma THM-530-D is the sharpest candidate). If contraction
monotonically lowers the functional toward consec, that IS the proof; if it ever raises it, dead.
kind-pasteur-2026-06-18-S8.

EWLB_A(E) = meas{ x : ∃ window W_j=(j/14,(j+2)/14) (j=0..6) with NO orbit point frac(e x), e∈E }
  — exact rational lower bound on μ_{1/7}(E) (THM-530-E). consec_k minimizes it (verified).
We compute EWLB_A EXACTLY (staggered-window sweep) and test the contraction lemma.
"""
import sys; sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
from fractions import Fraction as Fr
from math import gcd
from functools import reduce
import itertools, random

def EWLB(E):
    """exact meas{x: some window (j/14,(j+2)/14), j=0..6, is empty of {frac(e x):e in E}}."""
    E=[e for e in E if e!=0]  # e=0 sits at 0; it occupies windows containing 0 (j=0: (0,1/7) yes; j=6:(6/14,8/14)=(3/7,4/7) no). Keep e=0 effect:
    # e=0 -> frac(0)=0, in window j iff 0 in (j/14,(j+2)/14) -> j=13? no; (0,2/14): 0 is the OPEN left end, NOT inside. So e=0 occupies NO open window. Drop it.
    windows=[(Fr(j,14),Fr(j+2,14)) for j in range(7)]
    # breakpoints: x where some frac(e x) crosses a window boundary b in {j/14,(j+2)/14}
    bps=set([Fr(0),Fr(1)])
    bset=set()
    for j in range(7):
        bset.add(Fr(j,14)); bset.add(Fr(j+2,14))
    for e in E:
        for b in bset:
            # e x ≡ b (mod 1) -> x = (m + b)/e, m=0..e-1
            for m in range(e):
                bps.add(Fr(m,1)/e + b/e if False else (Fr(m)+b)/e)
    bps={p for p in bps if 0<=p<=1}
    pts=sorted(bps)
    total=Fr(0)
    def occ(xmid):
        # which windows are occupied by some frac(e*xmid)?
        occupied=[False]*7
        for e in E:
            f=(e*xmid)%1
            for j in range(7):
                lo,hi=windows[j]
                # window wraps if hi>1 (it doesn't: max (8/14)<1). open interval
                if lo< f <hi: occupied[j]=True
        return occupied
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        mid=(a+b)/2
        oc=occ(mid)
        if not all(oc):  # some window empty
            total+=b-a
    return total

def primitive(E): return reduce(gcd,[e for e in E if e!=0])==1 if any(E) else False

def gaps(E):
    s=sorted(E); return [s[i]-s[i-1] for i in range(1,len(s))]

def contract_gap(E, idx):
    """reduce the idx-th sorted gap by 1 (shift the upper part down by 1); return sorted set or None."""
    s=sorted(E)
    if s[idx]-s[idx-1] <= 1: return None  # gap already 1, can't contract
    t=s[:]
    for i in range(idx, len(t)): t[i]-=1
    if len(set(t))!=len(t) or min(t)<0: return None
    return tuple(t)

if __name__=="__main__":
    k=8
    consec=tuple(range(k))
    Ec=EWLB(consec)
    print(f"EWLB(consec_{k}={consec}) = {Ec} = {float(Ec):.6f}  (memory: 407/588={float(Fr(407,588)):.6f})")

    print("\n=== TEST 1: does consec MINIMIZE EWLB? (reproduce, small sweep) ===")
    rng=random.Random(8); below=0; tested=0; mn=(Fr(2),None)
    for _ in range(300):
        sp=rng.randint(k-1, 16)
        S=tuple(sorted({0}|set(rng.sample(range(1,sp+1), k-1))))
        if len(S)!=k or not primitive(S): continue
        v=EWLB(S); tested+=1
        if v<mn[0]: mn=(v,S)
        if v<Ec-Fr(1,10**9): below+=1
    print(f"  tested {tested}; min EWLB={float(mn[0]):.6f} at {mn[1]}; #below consec={below} (want 0)")

    print("\n=== TEST 2 (THE ROUTE-VALIDITY CHECK): gap-contraction monotonicity ===")
    print("  for random E, contract each contractible gap; does EWLB NOT increase (move toward consec)?")
    rng=random.Random(99); viol=0; chk=0; examples=[]
    for _ in range(120):
        sp=rng.randint(k+1, 18)
        S=tuple(sorted({0}|set(rng.sample(range(1,sp+1), k-1))))
        if len(S)!=k or not primitive(S): continue
        vS=EWLB(S)
        for idx in range(1,k):
            Sc=contract_gap(S, idx)
            if Sc is None or not primitive(Sc): continue
            vC=EWLB(Sc); chk+=1
            if vC > vS + Fr(1,10**9):  # contraction INCREASED EWLB -> lemma FAILS
                viol+=1
                if len(examples)<5: examples.append((float(vS),S,float(vC),Sc,idx))
    print(f"  checked {chk} contractions; VIOLATIONS (EWLB increased) = {viol}")
    if examples:
        print("  contraction-INCREASES-EWLB examples (would KILL the monotone route):")
        for vS,S,vC,Sc,idx in examples: print(f"    {S} EWLB={vS:.5f} --contract gap {idx}--> {Sc} EWLB={vC:.5f}")
    else:
        print("  NO violations found — gap-contraction monotonicity HOLDS on this sample (route is viable).")
