#!/usr/bin/env python3
"""Observer-blind LRC: the SYMMETRIC trienerment on the n points P={0,v_1,..,v_{n-1}}.
Each pair {i,j} -> state by ||(p_i-p_j)t|| vs 1/n. danger graph D(t); SAFE graph; point
p lonely <=> p isolated in D(t) <=> p universal in safe graph. Explore the symmetric
loneliness (which points lonely when), the geometry. opus-2026-06-03-S583 round 1."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
def dist(x): x%=1; return min(x,1-x)
def crit_times(P):
    n=len(P); cs=set()
    for i in range(n):
        for j in range(n):
            if i==j: continue
            d=abs(P[i]-P[j])
            if d:
                for m in range(0,2*d+1): cs.add(F(m,2*d)%1)  # include boundary crossings
    cs.discard(F(0))
    return sorted(cs)
def danger_edges(P,t,n,delta):
    return [(i,j) for i,j in combinations(range(n),2) if dist((P[i]-P[j])*t)<delta]
def lonely_points(P,t,n,delta):
    de=danger_edges(P,t,n,delta); inc=set()
    for i,j in de: inc.add(i); inc.add(j)
    return [p for p in range(n) if p not in inc]
def main():
    print("ROUND 1: symmetric trienerment / danger graph -- which points get lonely & when")
    # n = total points (incl observer 0); delta=1/n
    examples = {
      "AP n=7 {0..6}": [0,1,2,3,4,5,6],
      "n=7 {0,1,2,3,5,7,11}": [0,1,2,3,5,7,11],
      "n=5 {0,1,2,3,4}": [0,1,2,3,4],
      "n=5 {0,2,3,5,7}": [0,2,3,5,7],
    }
    for name,P in examples.items():
        P=sorted(set(P)); n=len(P); delta=F(1,n)
        ts=crit_times(P)
        # for each point, does it get lonely at some critical-interval midpoint?
        lonely_of={p:[] for p in range(n)}
        pts=ts
        for k in range(len(pts)):
            a=pts[k]; b=pts[(k+1)%len(pts)]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
            for p in lonely_points(P,mid,n,delta): lonely_of[p].append(float(mid))
        allcov=all(lonely_of[p] for p in range(n))
        obs_lonely=bool(lonely_of[0])
        cnt={p:len(lonely_of[p]) for p in range(n)}
        print(f"  {name}: every point lonely at some t? {allcov}; observer(0) lonely? {obs_lonely}; "
              f"#lonely-intervals per point={cnt}")
if __name__=='__main__': main()
