#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Independent per-subset derivation of FULL E[alpha2] at n=8 (and the 3-5 cross
term generally). alpha2 = unordered pairs of vertex-disjoint odd directed cycles.
For n=8: 3-3 disjoint + 3-5 disjoint (5-5 needs 10 vertices). Disjoint subsets have
independent internal arcs => P(both cycles) = P(cyc1)*P(cyc2). Flushed, exact."""
import sys
from itertools import combinations, product, permutations
from fractions import Fraction
if hasattr(sys.stdout,"reconfigure"): sys.stdout.reconfigure(encoding="utf-8")

def prob_dcycle(subset):
    sub=tuple(sorted(subset))
    pairs=list(combinations(sorted(sub,reverse=True),2))
    free=[p for p in pairs if p[0]-p[1]!=1]
    forced=[p for p in pairs if p[0]-p[1]==1]
    k=len(sub);s0=sub[0];rest=sub[1:]
    tot=Fraction(0)
    for bits in product((0,1),repeat=len(free)):
        adj={}
        for p in forced: adj[(p[0],p[1])]=1;adj[(p[1],p[0])]=0
        for idx,p in enumerate(free):
            if bits[idx]==0: adj[(p[0],p[1])]=1;adj[(p[1],p[0])]=0
            else: adj[(p[0],p[1])]=0;adj[(p[1],p[0])]=1
        c=0
        for perm in permutations(rest):
            seq=(s0,)+perm;ok=True
            for i in range(k):
                if adj.get((seq[i],seq[(i+1)%k]),0)==0: ok=False;break
            if ok:c+=1
        tot+=c
    return tot/(1<<len(free))

def a2_components(n):
    trips=list(combinations(range(1,n+1),3))
    fives=list(combinations(range(1,n+1),5))
    pc3={t:prob_dcycle(t) for t in trips}
    pc5={t:prob_dcycle(t) for t in fives}
    # 3-3 disjoint
    a33=Fraction(0)
    for i in range(len(trips)):
        for j in range(i+1,len(trips)):
            if set(trips[i]).isdisjoint(trips[j]):
                a33+=pc3[trips[i]]*pc3[trips[j]]
    # 3-5 disjoint
    a35=Fraction(0)
    for t3 in trips:
        s3=set(t3)
        for t5 in fives:
            if s3.isdisjoint(t5):
                a35+=pc3[t3]*pc5[t5]
    # 5-5 disjoint
    a55=Fraction(0)
    for i in range(len(fives)):
        for j in range(i+1,len(fives)):
            if set(fives[i]).isdisjoint(fives[j]):
                a55+=pc5[fives[i]]*pc5[fives[j]]
    return a33,a35,a55

print("Per-subset FULL alpha2 decomposition (3-3 + 3-5 + 5-5):", flush=True)
for n in range(6,11):
    a33,a35,a55=a2_components(n)
    full=a33+a35+a55
    print(f"n={n}: a33={a33}  a35={a35}  a55={a55}  FULL={full}", flush=True)
    if n==8:
        print(f"   CHECK n=8: a33(tri)={a33} vs claim 173/8={Fraction(173,8)} "
              f"({'OK' if a33==Fraction(173,8) else 'X'})", flush=True)
        print(f"   CHECK n=8: a35(cross)={a35} vs claim 447/32={Fraction(447,32)} "
              f"({'OK' if a35==Fraction(447,32) else 'X'})", flush=True)
        print(f"   CHECK n=8: full={full} vs claim 1139/32={Fraction(1139,32)} "
              f"({'OK' if full==Fraction(1139,32) else 'X'})", flush=True)
