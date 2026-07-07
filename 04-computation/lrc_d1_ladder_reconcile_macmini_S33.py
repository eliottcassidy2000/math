#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S33 (HYP-4642) -- FINISH the crux: (a) nail the d=1 ladder bound
(M({1..11,x}) >= 2/25 unless the 12-AP), (b) reconcile the defect-stratification (opus-S123)
with the transversal dichotomy (my S32b / concurrent pair-blocking): is every d>=3 family a
NON-transversal (so opus's 'd>=3 GREEN via mod-25' holds), or are there d>=3 transversals
(which mod-25 canNOT clear -- a hole to patch with a non-mod-25 witness)?
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
GAP_LO,GAP_HI=F(1,13),F(2,25)

def Mfast(S):
    S=sorted(set(S)); Q=set()
    for v in S: Q.add(2*v)
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0)
    for q in Q:
        for a in range(1,q):
            mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
            if F(mn,q)>best: best=F(mn,q)
    return best

def transversal(S):
    return len(set(frozenset({v%25,(-v)%25}) for v in S if v%5!=0))==10

def longest_ap(S):
    S=sorted(set(S)); idx=set(S); best=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            g=S[j]-S[i]; L=2; nx=S[j]+g
            while nx in idx: L+=1; nx+=g
            best=max(best,L)
    return best

if __name__=="__main__":
    # (a) d=1 base ladder: {1..11} u {x}
    print("(a) d=1 ladder  M({1..11,x}):  gap = (1/13,2/25)")
    ingap=0; vals={}
    for x in range(12,80):
        S=[i for i in range(1,12)]+[x]
        if len(set(S))<12: continue
        M=Mfast(S); vals[x]=M
        if GAP_LO<M<GAP_HI: ingap+=1
        tag = "=1/13 (AP)" if M==F(1,13) else ("** IN GAP **" if GAP_LO<M<GAP_HI else ("=2/25 wall" if M==GAP_HI else ""))
        if x<=30 or tag: print(f"    x={x:>3}: M={str(M):>7}  {tag}")
    above=[M for M in vals.values() if M>F(1,13)]
    print(f"  min M above 1/13 (x=12..79): {min(above)}  (target 2/25={F(2,25)});  # IN-GAP: {ingap}")

    # (b) reconcile defect vs transversal
    print("\n(b) cross-tab defect d=12-longestAP  vs  transversal (pair-blocker) mod 25:")
    cnt={'d>=3 & transversal (mod-25 FAILS)':0,'d>=3 & non-transversal (mod-25 clears)':0,
         'd<3 & transversal':0,'d<3 & non-transversal':0}
    ingap_d3_tr=[]
    tested=0
    for g in range(1,5):
      for c in range(1,4):
        base=[c*(1+g*j) for j in range(9)]
        span=max(base); pool=[x for x in range(1,span+6) if x not in base]
        for defs in combinations(pool[:15],3):
            S=sorted(base+list(defs))
            if len(set(S))!=12: continue
            S=[v//reduce(gcd,S) for v in S]
            if len(set(S))!=12: continue
            d=12-longest_ap(S); tr=transversal(S)
            key=f"d{'>=3' if d>=3 else '<3'} & {'transversal' if tr else 'non-transversal'}"
            if d>=3 and tr:
                cnt['d>=3 & transversal (mod-25 FAILS)']+=1
                M=Mfast(S)
                if GAP_LO<M<GAP_HI: ingap_d3_tr.append((tuple(S),M))
            elif d>=3: cnt['d>=3 & non-transversal (mod-25 clears)']+=1
            elif tr: cnt['d<3 & transversal']+=1
            else: cnt['d<3 & non-transversal']+=1
            tested+=1
            if tested>50000: break
        if tested>50000: break
      if tested>50000: break
    print(f"    tested {tested}:")
    for k,v in cnt.items(): print(f"      {k}: {v}")
    print(f"    d>=3 & transversal families IN GAP: {len(ingap_d3_tr)}")
    print("    => 'd>=3 & transversal'>0 means opus's 'd>=3 GREEN via mod-25' is INCOMPLETE for those")
    print("       (mod-25 rotation fails on transversals); but if none are IN GAP they're still loose")
    print("       via another witness -- the transversal dichotomy is the cleaner PARTITION.")
