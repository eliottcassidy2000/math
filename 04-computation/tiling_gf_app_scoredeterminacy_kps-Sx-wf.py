#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Score-determinacy test (flushed): is c3 / c5 / alpha2-tri determined by the
sorted score vector across the tiling space? Claim: c3 YES, c5 NO, alpha2 NO.
Also tests whether (score, A^2-derived) enrichment determines them: claim score+A^2
works at n<=6 but FAILS at n=7."""
import sys
from itertools import product, combinations, permutations
from collections import defaultdict
from math import comb
if hasattr(sys.stdout,"reconfigure"): sys.stdout.reconfigure(encoding="utf-8")

def tiles(n): return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]

def cyc(sub):
    sub=tuple(sorted(sub)); s0=sub[0]; rest=sub[1:]
    return [(s0,)+p for p in permutations(rest)]

def count(seqs,out):
    c=0
    for seq in seqs:
        ok=True;k=len(seq)
        for i in range(k):
            if not (out[seq[i]]>>seq[(i+1)%k])&1: ok=False;break
        if ok:c+=1
    return c

def analyze(n):
    T=tiles(n);F=len(T)
    trips=list(combinations(range(1,n+1),3)); ts={t:cyc(t) for t in trips}
    fives=list(combinations(range(1,n+1),5)); fs={t:cyc(t) for t in fives}
    tpairs=[(i,j) for i in range(len(trips)) for j in range(i+1,len(trips))
            if set(trips[i]).isdisjoint(trips[j])]
    base=[0]*(n+1)
    for k in range(n,1,-1): base[k]|=(1<<(k-1))
    by_score=defaultdict(lambda: {"c3":set(),"c5":set(),"a2":set()})
    by_score_a2mat=defaultdict(lambda: {"c5":set(),"a2":set()})
    for bv in product((0,1),repeat=F):
        out=base[:]
        for (a,b),bit in zip(T,bv):
            if bit==0: out[a]|=(1<<b)
            else: out[b]|=(1<<a)
        sv=[bin(out[i]).count("1") for i in range(1,n+1)]
        key=tuple(sorted(sv))
        c3=comb(n,3)-sum(comb(s,2) for s in sv)
        tri_ind={t:count(ts[t],out) for t in trips}
        c5=sum(count(fs[t],out) for t in fives)
        a2=sum(tri_ind[trips[i]]*tri_ind[trips[j]] for (i,j) in tpairs)
        by_score[key]["c3"].add(c3); by_score[key]["c5"].add(c5); by_score[key]["a2"].add(a2)
        # A^2 enrichment: multiset of (row sums of A^2) = #2-paths per vertex -> sorted
        # build A as 0/1 matrix and compute A^2 row sums
        A=[[ (out[i]>>j)&1 for j in range(1,n+1)] for i in range(1,n+1)]
        a2mat=[[sum(A[i][k]*A[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
        a2feat=tuple(sorted(sum(row) for row in a2mat))
        key2=(key,a2feat)
        by_score_a2mat[key2]["c5"].add(c5); by_score_a2mat[key2]["a2"].add(a2)
    c3det=all(len(v["c3"])==1 for v in by_score.values())
    c5det=all(len(v["c5"])==1 for v in by_score.values())
    a2det=all(len(v["a2"])==1 for v in by_score.values())
    c5det2=all(len(v["c5"])==1 for v in by_score_a2mat.values())
    a2det2=all(len(v["a2"])==1 for v in by_score_a2mat.values())
    print(f"n={n}: SCORE-only determines:  c3={c3det}  c5={c5det}  a2tri={a2det}", flush=True)
    print(f"n={n}: SCORE+A^2-rowsums determines:  c5={c5det2}  a2tri={a2det2}", flush=True)
    if not c5det:
        for k,v in by_score.items():
            if len(v["c5"])>1:
                print(f"    c5 counterexample @score {k}: c5 in {sorted(v['c5'])}", flush=True); break
    if not a2det:
        for k,v in by_score.items():
            if len(v["a2"])>1:
                print(f"    a2 counterexample @score {k}: a2 in {sorted(v['a2'])}", flush=True); break

for n in (5,6,7):
    analyze(n)
