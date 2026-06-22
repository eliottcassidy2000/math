#!/usr/bin/env python3
"""Realizability lens: the THREE-GAP (Steinhaus) clock as the LRC obstruction object.
(1) Only AP clusters have the <=3-gap property for all x (three-gap rigidity is unique to APs).
(2) That rigidity makes sector coverage 'all-or-nothing' => N_E most BIMODAL => AP maximizes L_y (Node 2).
kind-mendel-2026-06-22-S7."""
from fractions import Fraction as F
from math import gcd, floor
from functools import reduce
from itertools import combinations
def gl(xs): return reduce(lambda a,b:a*b//gcd(a,b),[x for x in xs if x],1)

def num_gaps_distribution(E, samples=2000):
    "for x in (0,1), # distinct circular gap-lengths among {frac(e x)}; report max over x"
    mx=0; hist={}
    for i in range(1,samples):
        x=F(i, samples)
        pts=sorted(set((e*x)%1 for e in E))
        if len(pts)<2: continue
        gaps=set()
        for a,b in zip(pts, pts[1:]+[pts[0]+1]):
            gaps.add(b-a)
        mx=max(mx,len(gaps)); hist[len(gaps)]=hist.get(len(gaps),0)+1
    return mx, hist

def NE_distribution(E):
    "exact distribution of N_E(x)=# empty inner 1/7-sectors (sectors 1..6); returns {N: measure}"
    pos=[e for e in E if e]
    D=7*gl(pos); bset=set([0,D])
    for e in pos:
        step=D//(7*e); x=0
        while x<=D: bset.add(x); x+=step
    bps=sorted(bset); dist={}
    for a,b in zip(bps,bps[1:]):
        if b<=a: continue
        mid2=a+b
        hit=set()
        for e in E:
            hit.add((7*((e*mid2)%(2*D)))//(2*D))
        # empty among inner sectors 1..6 (sector 0 always hit by e=0)
        Nempty=len([j for j in range(1,7) if j not in hit])
        dist[Nempty]=dist.get(Nempty,F(0))+F(b-a,D)
    return dist
def bimodality(dist):
    "fraction of mass at the extremes N=0 and N=6 (high => bimodal/all-or-nothing)"
    tot=sum(dist.values())
    return float((dist.get(0,F(0))+dist.get(6,F(0)))/tot)

print("=== (1) THREE-GAP rigidity: max # distinct gaps over x (AP should be <=3) ===")
k=8
tests={'AP/consec {0..7}':list(range(8)),
       'perturb {0,1,2,3,4,5,6,8}':[0,1,2,3,4,5,6,8],
       'spread {0,1,2,4,7,11,16,22}':[0,1,2,4,7,11,16,22],
       'dilated 2*{0..7}? primitive {0,2,4,6,8,10,12,14}':[0,2,4,6,8,10,12,14]}
for name,E in tests.items():
    mx,hist=num_gaps_distribution(E)
    print(f"  {name:42s} max distinct gaps = {mx}  hist={dict(sorted(hist.items()))}")

print("\n=== (2) BIMODALITY of N_E (mass at extremes N=0,6) and L_y proxy; AP should be most bimodal ===")
for name,E in tests.items():
    d=NE_distribution(E); bm=bimodality(d)
    p0=float(d.get(0,F(0)))   # N=0 => all sectors hit = cover
    print(f"  {name:42s} bimodality={bm:.4f}  p0(N=0)={p0:.4f}  N-dist={ {k:round(float(v),3) for k,v in sorted(d.items())} }")
