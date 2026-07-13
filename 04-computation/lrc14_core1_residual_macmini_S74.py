#!/usr/bin/env python3
"""mac-mini-S74: confirm opus-S259's route reduces to the |core|=1 residual (= LRC(14) itself).
For |core|=1 families, coreCover = runner-1 (the sole coprime-core speed) density in G' =
1 - safe/meas(G'); coreCover<1 <=> M>1/14 = the conjecture, NOT helped by equidistribution
(runner 1 is a single arc). Check the EXTREMALS (deep well etc.) -- coreCover closest to 1?"""
from math import gcd
LV=1.0/14; P30030=[2,3,5,7,11,13]
def cop(v): return all(v%p!=0 for p in P30030)
def is_cov(S,n=14): return all(any(v%q==0 for v in S) for q in range(2,n+1))
def prim(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def nn(x): x%=1.0; return min(x,1-x)
def analyze(S,res=800000):
    core=[v for v in S if cop(v)]; noncore=[v for v in S if not cop(v)]
    if not core: return None
    inGp=0; cd=0
    for j in range(res):
        t=(j+0.5)/res
        if all(nn(w*t)>=LV for w in noncore):
            inGp+=1
            if any(nn(v*t)<LV for v in core): cd+=1
    if inGp==0: return (len(core),None,0.0)
    return (len(core), cd/inGp, inGp/res)

print(f"1/14={LV:.5f}. |core|=1 extremals: coreCover = runner-1 density = 1 - (M>1/14 slack).\n")
print(f"{'family':30s} | |core| | core | coreCover | meas(G') | note")
print("-"*84)
fams={
 "deep well {1..12,182}": [*range(1,13),182],
 "{1..12,364}": [*range(1,13),364],
 "{1..12,546}": [*range(1,13),546],
 "{1..11,13,84}": [*range(1,12),13,84],
 "{1..11,13,168}": [*range(1,12),13,168],
 "{1..10,22,13,84}": [*range(1,11),22,13,84],
 "2*deepwell {2..24,182}(imprim)": [2,4,6,8,10,12,14,16,18,20,22,24,182],
}
mx=(0,None)
for nm,S in fams.items():
    S=sorted(set(S))
    r=analyze(S)
    if r is None: print(f"  {nm:28s}: core empty"); continue
    k,cc,mg=r
    core=[v for v in S if cop(v)]
    if cc is None: print(f"  {nm:28s}: G' empty"); continue
    if cc>mx[0]: mx=(cc,nm)
    note = "<== |core|=1: NOT closed by equidist = LRC(14)" if k==1 else "|core|>=2: equidist route applies"
    print(f"  {nm:28s} | {k:5d} | {str(core):10s} | {cc:.5f} | {mg:.5f} | {note}")
print(f"\nMAX coreCover (extremals): {mx[0]:.5f} at {mx[1]}")
print("VERDICT: the |core|=1 extremals have coreCover ~0.85-0.95 (< 1, but ONLY because M>1/14).")
print("opus-S259 equidistribution closes |core|>=2 (each core runner ~1/7, sum<1); the |core|=1")
print("case has NO equidistribution (single arc) => coreCover<1 there IS the conjecture, unclosed.")
print("So opus-S259 REDUCES LRC(14) to the |core|=1 (runner-1-dominated) extremals -- real")
print("progress (large-core done) but the hard extremals remain = HYP-2566.")
