#!/usr/bin/env python3
"""Round 4 (bounce to Lemma A, discrepancy route): does circuit-free => safe MEASURE
bounded below? and is the relevant energy the 3-term count (not 4-term)? Cross-check
the discrepancy premise. opus-2026-06-03-S577."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
def dist(x): x%=1; return min(x,1-x)
def safe_measure(V,n):
    THR=F(1,n); eps=set([F(0)])
    for v in V:
        for k_ in range(v+1):
            for s in(-1,1): eps.add(F(k_*n+s,n*v)%1)
    pts=sorted(eps); meas=F(0); L=len(pts)
    for i in range(L):
        a=pts[i]; b=pts[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(v*mid)>THR for v in V): meas+=ln
    return float(meas)
def n3(V):
    s=set(V); return sum(1 for a,b in combinations(sorted(V),2) if a+b in s)
def n4(V):
    sums={}
    for a,b in combinations(sorted(V),2): sums.setdefault(a+b,0); sums[a+b]+=1
    return sum(c*(c-1)//2 for c in sums.values())
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    rng=random.Random(404)
    print("Lemma A measure route: min safe-measure(level delta) for circuit-free vs 3-term-rich;")
    print("and does measure track 3-term (not 4-term)?")
    for k in [6,8,10]:
        n=k+1; B=2*k+6
        cf_minmu=None; rich_minmu=None
        # 4-term-rich but circuit-free: measure still big?
        cf4_minmu=None
        for _ in range(6000):
            V=prim(tuple(sorted(rng.sample(range(1,B+1),k))))
            if len(V)!=k: continue
            mu=safe_measure(V,n); t3=n3(V)
            if t3==0:
                cf_minmu=mu if cf_minmu is None else min(cf_minmu,mu)
                if n4(V)>=3: cf4_minmu=mu if cf4_minmu is None else min(cf4_minmu,mu)
            elif t3>=3:
                rich_minmu=mu if rich_minmu is None else min(rich_minmu,mu)
        print(f"  k={k:2d}: circuit-free min safe-measure={cf_minmu:.4f} (4-term-rich-but-cf min={('%.4f'%cf4_minmu) if cf4_minmu else 'n/a'}); "
              f"3-term-rich(>=3) min safe-measure={rich_minmu:.4f}")
if __name__=='__main__': main()
