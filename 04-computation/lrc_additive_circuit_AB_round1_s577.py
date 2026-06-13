#!/usr/bin/env python3
"""Lemma A (circuit-free => G>=delta via equidistribution) vs Lemma B (3-term relation
= a fold). Round 1: machinery + de-risk both premises + verify fold=shield.
k = #nonzero speeds, delta = 1/(k+1).  opus-2026-06-03-S577."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
def dist(x): x%=1; return min(x,1-x)
def Mexact(V):
    """M(S)=max_t min_i ||v_i t||, exact via crossing times t=m/(v_i +- v_j)."""
    cands=set()
    for i in range(len(V)):
        for j in range(len(V)):
            if i==j: continue
            for D in (V[i]+V[j], abs(V[i]-V[j])):
                if D==0: continue
                for m in range(1,D): cands.add(F(m,D))
        cands.add(F(1,2*V[i]))  # single-runner peak
    best=F(0)
    for t in cands:
        mn=min(dist(v*t) for v in V)
        if mn>best: best=mn
    return best
def three_terms(V):
    s=set(V); return [(a,b,a+b) for a,b in combinations(sorted(V),2) if a+b in s]
def four_terms(V):  # additive quadruples a+b=c+d (energy), distinct unordered
    s={}; cnt=0
    sums={}
    for a,b in combinations(sorted(V),2):
        sums.setdefault(a+b,[]).append((a,b))
    for val,prs in sums.items():
        if len(prs)>=2: cnt+=len(prs)*(len(prs)-1)//2
    return cnt
def prim(V):
    g=0
    for v in V: g=gcd(g,v)
    return tuple(sorted(v//g for v in V))
def main():
    print("== fold = shield check: for v_c=v_a+v_b, at pinch t=m/v_c runner c is at integer ==")
    for V in [(1,2,3),(1,3,4),(2,3,5),(1,4,5)]:
        tt=three_terms(V); ok=all(dist(c*F(1,c))==0 for (a,b,c) in tt)
        print(f"   V={V}: 3-terms={tt}; c at 0 on its clock (shield): {ok}")
    print()
    print("== Lemma A premise: circuit-free (no 3-term) min M vs delta=1/(k+1) ==")
    rng=random.Random(7)
    for k in [4,5,6,7,8,9,10]:
        delta=F(1,k+1); B=3*k+4
        minM=None; argmin=None; n_cf=0; min4rich=None
        for _ in range(6000):
            V=prim(tuple(sorted(rng.sample(range(1,B+1),k))))
            if len(V)!=k: continue
            if three_terms(V): continue          # circuit-free only
            n_cf+=1
            M=Mexact(V)
            if minM is None or M<minM: minM=M; argmin=V
            if four_terms(V)>=2:                  # 4-term-rich but circuit-free
                if min4rich is None or M<min4rich: min4rich=M
        mar=float(minM-delta) if minM else None
        m4=float(min4rich-delta) if min4rich else None
        print(f"   k={k:2d}: delta={float(delta):.4f}; circuit-free min M={float(minM):.4f} (margin {mar:+.4f}); "
              f"4-term-rich circuit-free margin {('%.4f'%m4) if m4 is not None else 'n/a'}; argmin={argmin}")
if __name__=='__main__': main()
