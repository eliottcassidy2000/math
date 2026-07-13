#!/usr/bin/env python3
"""mac-mini-S73: pursue opus-S257's dual certificate via the EXPLICIT test measure nu = Leb|_G'.
c=14/183. Core = runners coprime to 30030 (=2*3*5*7*11*13). G' = {t: all NON-core runners safe
at level c}. Dual cert nu=Leb|_G' gives INT W dnu = (core-danger measure in G')/meas(G').
UNION BOUND within G': core-danger <= |core|*2c, so safe-in-G' >= meas(G') - |core|*2c > 0
whenever meas(G') > |core|*2c => M(v) >= c. Test: does this CLOSE loose families and FAIL only
for the deep well (tight, small G' = the knife-edge)?"""
from fractions import Fraction as F
from math import gcd
c=F(14,183); twoc=2*c
P30030=[2,3,5,7,11,13]
def coprime30030(v): return all(v%p!=0 for p in P30030)

def safe_measure(runners, level):
    """meas{t in [0,1): ||v t||>=level for all v in runners}, exact via breakpoints."""
    if not runners: return F(1)
    pts=set([F(0),F(1)])
    for v in runners:
        for m in range(0,v+1):
            for s in (level,-level):
                t=(m+s)/v
                if 0<=t<=1: pts.add(t)
    Pl=sorted(pts); tot=F(0)
    for a,b in zip(Pl,Pl[1:]):
        mid=(a+b)/2
        if all(min((v*mid)%1,1-((v*mid)%1))>=level for v in runners): tot+=b-a
    return tot

def safe_measure_intersect(coreset, noncoreset, level):
    """meas(G' AND core-safe) = meas{t: ALL runners safe} = safe_measure(all)."""
    return safe_measure(coreset+noncoreset, level)

def Mval(S,Qmax):
    best=F(0)
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=min(min((a*v)%q, q-((a*v)%q)) for v in S)
            if mind>0 and F(mind,q)>best: best=F(mind,q)
    return best

print(f"c=14/183={float(c):.5f}, 2c={float(twoc):.5f}. Dual cert nu=Leb|_G' closes if meas(G')>|core|*2c.\n")
print(f"{'family':34s} | |core| | meas(G') | |core|*2c | closes? | safe-in-G' | M")
print("-"*104)
fams={
 "deep well {1..12,182}": [*range(1,13),182],
 "{1..11,13,84}": [*range(1,12),13,84],
 "{1..10,22,13,84}": [*range(1,11),22,13,84],
 "{2,4,6,8,10,12,14,16,18,20,22,24,182}": [2,4,6,8,10,12,14,16,18,20,22,24,182], # imprim (2*..)
 "{1..12,364}": [*range(1,13),364],
 "{1,2,3,4,5,6,7,8,9,10,11,12,910}":[*range(1,13),910],
}
for nm,S in fams.items():
    S=sorted(set(S))
    if len(S)!=13: continue
    core=[v for v in S if coprime30030(v)]
    noncore=[v for v in S if not coprime30030(v)]
    mGp=safe_measure(noncore,c)
    bound=len(core)*twoc
    closes = mGp>bound
    safeall=safe_measure(S,c)  # meas(G' AND core-safe) = full safe set
    M=Mval(S,min(2*max(S),260))
    print(f"{nm:34s} | {len(core):5d} | {float(mGp):.5f} | {float(bound):.5f} | {str(closes):5s} | {float(safeall):.6f} | {float(M):.5f}")
print()
print("READING: if the criterion CLOSES loose families (meas(G')>|core|*2c => safe pt) and FAILS")
print("only for the deep well (tight, meas(G')<=|core|*2c, safe measure 0 = knife-edge), then the")
print("dual cert nu=Leb|_G' PROVES the loose stratum, deep well by S255 rigidity. Check the split.")
