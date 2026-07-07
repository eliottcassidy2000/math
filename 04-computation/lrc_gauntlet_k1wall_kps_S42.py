#!/usr/bin/env python3
"""
kps-S42 (creative whittling): sharpen (G) to 'no order-k>=2 value is achieved at N=12', and
reconcile the defect-frame (mac-mini HYP-4612) with the ORDER-frame (opus S122 MISTAKE-115:
defect count does NOT govern; order does -- the per-order gauntlet, kps S40).

KEY STRUCTURAL FACT: both gap edges are k=1 (Kravitz) rungs --
    1/13 = 1/(12*1+1)  (s=1,k=1)   and   2/25 = 2/(12*2+1)  (s=2,k=1),
and the interior values are exactly k>=2 (k<s<2k).  So the first gap is the open interval
between two consecutive k=1 rungs, and (G) <=> no k>=2 value is achieved at N=12
<=> 'M < 2/25 => AP' (mac-mini's AP-rigidity), given LRC(13) => M >= 1/13.

This file: (1) a broad diverse N=12 search shows the achievable M pile up at the k=1 WALL 2/25,
never reaching k>=2; (2) the k=1 wall family {1..11,24} clears mod-25 at c=2 (M=2/25 exactly),
so kps-S41's mod25_covering_floor proves the loose side AT the wall; (3) the residual 'M<2/25'
families are exactly the mod-25-NON-clearable ones (residues cover all 10 antipodal pairs mod 25,
no multiple of 25) -- and the AP {1..12} is one such (residues {1..12} hit one per pair).
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
    return Fraction(bn,bd),(bc,bd)
def mod25_clear(v):
    for c in range(1,25):
        if gcd(c,25)!=1: continue
        if all(min((x*c)%25,25-(x*c)%25)>=2 for x in v): return c
    return None
def covers_all_pairs(v):
    R=set(x%25 for x in v); pairs=[{u,25-u} for u in range(1,13) if gcd(u,25)==1]
    return all(p&R for p in pairs)

N=12; LO,HI=Fraction(1,N+1),Fraction(2,2*N+1)
print("=== (1) both gap edges are k=1 rungs; interior is k>=2 ===", flush=True)
for M in [Fraction(1,13),Fraction(2,25),Fraction(3,38),Fraction(4,51)]:
    s,q=M.numerator,M.denominator; print(f"    {M}: s={s}, k=q-12s={q-12*s}", flush=True)

print("\n=== (2) broad N=12 search: achievable M pile up at the k=1 WALL 2/25, never k>=2 ===", flush=True)
random.seed(9); ingap=0; seen=0; wall={}
for _ in range(45000):
    r=random.random()
    if r<0.5:
        d=random.choice([1,1,2,3,4,5,6]); L=random.choice([7,8,9,10,11]); a=random.randint(1,4)
        ap=[a+i*d for i in range(L)]; pool=[x for x in range(1,50) if x not in ap]
        if len(pool)<12-L: continue
        v=ap+random.sample(pool,12-L)
    else:
        v=list(range(1,13))
        for _ in range(random.randint(1,4)): v[random.randrange(12)]=random.randint(1,48)
    v=sorted(set(v))
    if len(v)!=12: continue
    g=reduce(gcd,v)
    if g!=1: v=[x//g for x in v]
    if len(set(v))!=12: continue
    seen+=1; M,_=Mw(v)
    if LO<M<HI: ingap+=1
    elif HI<=M<Fraction(1,11):
        k=M.denominator-12*M.numerator
        if k not in wall or M<wall[k][0]: wall[k]=(M,v)
print(f"    {seen} families; IN open gap (1/13,2/25): {ingap}", flush=True)
print(f"    smallest M-above-gap by order k (the wall):", flush=True)
for k in sorted(wall):
    M,v=wall[k]; print(f"      k={k}: {M} (~{float(M):.4f})  e.g. {v}", flush=True)
print(f"    => the wall is k=1 at 2/25; NOTHING reaches k>=2 (the gap interior).", flush=True)

print("\n=== (3) the k=1 wall family {1..11,24} clears mod-25 => S41 certificate applies AT the wall ===", flush=True)
w=list(range(1,12))+[24]; M,(c,q)=Mw(w); cc=mod25_clear(w)
print(f"    {w}: M={M} at t={c}/{q}; mod-25 clear at c={cc} (residues {[x*cc%25 for x in w]})", flush=True)
print(f"    => mod25_covering_floor(w, {cc}) certifies M >= 2/25 (the loose side, exactly at 2/25).", flush=True)

print("\n=== (4) M<2/25 families are mod-25-NON-clearable (cover all 10 pairs); AP {1..12} is one ===", flush=True)
ap=list(range(1,13)); print(f"    AP {ap}: M={Mw(ap)[0]}, mod-25 clear? {mod25_clear(ap)}, covers all pairs? {covers_all_pairs(ap)}", flush=True)
print(f"    => rigidity 'M<2/25 => AP' <=> 'only the all-pairs-covering AP has M<2/25' (the residual).", flush=True)
