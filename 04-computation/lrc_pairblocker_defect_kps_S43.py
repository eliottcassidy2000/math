#!/usr/bin/env python3
"""
kps-S43 (creative work on the crux): the pair-blocking rigidity (mac-mini S32 case 2) is the
WHOLE remaining crux.  Two checks that whittle it:

(A) VERIFY opus S123's 'd>=3 GREEN via kps mod-25': it relies on 'd>=3 => not a pair-blocker
    (=> mod-25-clearable => M>=2/25 by LRCMod25Floor)'.  So NO pair-blocker should have d>=3.
    Check the defect distribution of pair-blockers -- if all are d<=2, opus's d>=3 GREEN is SOUND
    (the pair-blocker residual lives entirely in d<=2, matching the d=1,d=2 open strata).

(B) The non-AP blockers have M>=1/12 (mac-mini); a small-modulus covering system leaves a residual
    (S43: 1182/~10k uncleared by {7,9,10,11,12,13,14,24,26}).  So the rigidity is NOT a simple
    finite-modulus covering -- it needs the AP-minimality (mod-13 prime tight-locus, S12).  Confirm
    the residual blockers are comfortably loose (M>=1/12) and characterize them.

Pair-blocker = 12 speeds, no multiple of 25, residues mod 25 hit all 10 unit +-pairs.
"""
from fractions import Fraction
import numpy as np, random
from math import gcd
from functools import reduce

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        a=np.arange(1,q); r=np.outer(va,a)%q; d=np.minimum(r,q-r); col=d.min(axis=0); qb=int(col.max())
        if qb*bd>bn*q: bn,bd,bc=qb,q,int(a[col.argmax()])
    return Fraction(bn,bd)
def lap(v):
    s=set(v); vs=sorted(v); b=1
    for i in range(len(vs)):
        for j in range(i+1,len(vs)):
            dd=vs[j]-vs[i]; L=2; x=vs[j]+dd
            while x in s: L+=1; x+=dd
            b=max(b,L)
    return b
PAIRS=[{u,25-u} for u in [1,2,3,4,6,7,8,9,11,12]]
def is_blocker(v):
    if any(x%25==0 for x in v): return False
    R=set(x%25 for x in v); return all(p&R for p in PAIRS)

random.seed(2)
base=list(range(1,13)); blockers=[]; seen=set()
for _ in range(200000):
    v=base[:]
    for _ in range(random.randint(1,4)):
        i=random.randrange(12)
        v[i]=random.choice([v[i]+25,v[i]+50,(25-v[i]%25) if v[i]%25 else v[i],random.randint(1,60)])
    v=tuple(sorted(set(x for x in v if x>0)))
    if len(v)!=12 or v in seen: continue
    seen.add(v)
    if reduce(gcd,v)!=1: continue
    if is_blocker(v): blockers.append(v)
print(f"pair-blockers found: {len(blockers)}\n", flush=True)

print("=== (A) defect distribution of pair-blockers (opus d>=3 GREEN needs: no blocker has d>=3) ===", flush=True)
defdist={}; d3=[]
for v in blockers:
    d=12-lap(v); defdist[d]=defdist.get(d,0)+1
    if d>=3: d3.append(v)
print(f"  defect(d=12-longestAP) distribution: {dict(sorted(defdist.items()))}", flush=True)
print(f"  pair-blockers with d>=3: {len(d3)}  {'<-- would be a HOLE in d>=3 GREEN!' if d3 else '(none: d>=3 => not a blocker => mod-25-clearable => opus d>=3 GREEN is SOUND)'}", flush=True)
if d3:
    for v in d3[:5]: print(f"    d>=3 blocker {v}: M={Mw(v)}", flush=True)

print("\n=== (B) M-spectrum of pair-blockers: only AP < 2/25; others >= 1/12 ===", flush=True)
tight=[v for v in blockers if Mw(v)<Fraction(2,25)]
below_1_12=[v for v in blockers if Mw(v)<Fraction(1,12) and list(v)!=list(range(1,13))]
print(f"  blockers with M<2/25: {len(tight)} = {[list(v) for v in tight]}", flush=True)
print(f"  NON-AP blockers with M<1/12: {len(below_1_12)} (mac-mini: should be 0) {below_1_12[:3]}", flush=True)
mmin=min((Mw(v) for v in blockers if list(v)!=list(range(1,13))), default=None)
print(f"  min M over NON-AP blockers: {mmin} (~{float(mmin):.4f}) -- the second-smallest blocker value (>2/25 => (G))", flush=True)

print("\n=== READING ===", flush=True)
print("  If (A) shows no d>=3 blocker: opus's 'd>=3 GREEN via kps mod-25' is SOUND -- the pair-blocker", flush=True)
print("  residual = exactly the d<=2 strata (d=1,d=2 open, d=0=AP). The whole crux is the pair-blocking", flush=True)
print("  rigidity restricted to d<=2, and (B) says its M-spectrum jumps 1/13 -> [second value] > 2/25.", flush=True)
print("  The second-smallest blocker M is the number (G) turns on.", flush=True)
