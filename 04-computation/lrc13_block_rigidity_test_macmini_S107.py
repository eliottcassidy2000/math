#!/usr/bin/env python3
"""mac-mini-S107: SKEPTICAL test of the LRC(13) 12-block tightness rigidity BEFORE proving it.
Claim: every primitive 12-set A with M(A)=1/13 equals {1..12}. But LRC(14) (n=13) has a NON-AP tight
instance (Goddyn-Wong {1..11,13,24}, THM-733/734). Does an analogous doubling make the 12-set claim
FALSE? Test GW-analogs + single-element doublings + a wider structured search."""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
def M_exact(S):
    best=F(0); dens=set()
    for a,b in combinations(sorted(set(S)),2):
        dens.add(a+b)
        if b>a: dens.add(b-a)
    for q in dens:
        for m in range(1,q):
            num=q
            for v in S:
                r=(v*m)%q; d=min(r,q-r)
                if d<num: num=d
                if num*13<q: break
            if F(num,q)>best: best=F(num,q)
    return best
def prim(S): return reduce(gcd,S)==1
base=list(range(1,13))
print(f"base {{1..12}}: M={M_exact(base)} (expect 1/13)")
print("\n(1) GW-analog candidates (drop one, double another) -- do any hit M=1/13?")
cands={
 "{1..10,12,22} (drop11,add22)": [*range(1,11),12,22],
 "{1..11,24} (drop12,add24)":     [*range(1,12),24],
 "{1..11,13,24} = LRC(14) GW (13-set, control)": [*range(1,12),13,24],
 "{2..13} (dilate? no, =shift)":  list(range(2,14)),
}
for j in range(1,13):  # double element j
    S=sorted(set([x for x in base if x!=j]+[2*j]))
    if len(S)==12: cands[f"double {j}->{2*j}"]=S
hits=[]
for name,S in cands.items():
    if len(set(S))!=12: continue
    M=M_exact(S); p=prim(S)
    tag=""
    if M==F(1,13) and p and sorted(S)!=base: tag="  <<< NON-{1..12} TIGHT (rigidity FALSE!)"; hits.append(sorted(S))
    print(f"   {name}: M={M}={float(M):.5f} prim={p}{tag}")
print(f"\n(2) wider structured search: 12-subsets of {{1..26}}, primitive, M=1/13, not {{1..12}}:")
# too many to brute-force C(26,12); sample + target near-tight (contain {1..10} + 2 larger)
import random; rng=random.Random(13)
found=[]; checked=0
for _ in range(200000):
    small=list(range(1,11))  # keep {1..10} (likely needed for tightness), pick 2 more from {11..26}
    extra=rng.sample(range(11,27),2)
    S=sorted(small+extra)
    if len(set(S))!=12 or not prim(S): continue
    checked+=1
    if M_exact(S)==F(1,13) and S!=base: found.append(S)
print(f"   checked {checked} sets ({{1..10}}+2 from {{11..26}}); M=1/13 non-{{1..12}}: {len(sorted(set(map(tuple,found))))}")
for S in sorted(set(map(tuple,found)))[:8]: print(f"     {list(S)}")
print()
if hits or found:
    print("=> RIGIDITY IS FALSE: non-{1..12} primitive 12-sets with M=1/13 exist (GW-type). Do NOT prove it;")
    print("   the tight 12-sets are {1..12} PLUS a GW family, exactly as LRC(14) has {AP, GW-doubling}.")
else:
    print("=> no non-{1..12} tight 12-set found (GW-analogs all M!=1/13); rigidity SURVIVES -- attempt proof.")
