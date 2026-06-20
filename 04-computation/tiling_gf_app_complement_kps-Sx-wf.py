#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Complement-halving (2x lossless) check: c3,c5,alpha2-tri invariant under
tournament complement T->T^op (reverse all arcs). If invariant, the half-tiling
(THM-549) 2x fold is real and lossless for these OCF data. Flushed."""
import sys
from itertools import product, combinations, permutations
from math import comb
if hasattr(sys.stdout,"reconfigure"): sys.stdout.reconfigure(encoding="utf-8")

def tiles(n): return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]
def cyc(sub):
    sub=tuple(sorted(sub));s0=sub[0];rest=sub[1:]
    return [(s0,)+p for p in permutations(rest)]
def count(seqs,out):
    c=0
    for seq in seqs:
        ok=True;k=len(seq)
        for i in range(k):
            if not (out[seq[i]]>>seq[(i+1)%k])&1: ok=False;break
        if ok:c+=1
    return c

def check(n):
    T=tiles(n);F=len(T)
    trips=list(combinations(range(1,n+1),3));ts={t:cyc(t) for t in trips}
    fives=list(combinations(range(1,n+1),5));fs={t:cyc(t) for t in fives}
    tpairs=[(i,j) for i in range(len(trips)) for j in range(i+1,len(trips))
            if set(trips[i]).isdisjoint(trips[j])]
    base=[0]*(n+1)
    for k in range(n,1,-1): base[k]|=(1<<(k-1))
    bad=0;checked=0
    for bv in product((0,1),repeat=F):
        out=base[:]
        for (a,b),bit in zip(T,bv):
            if bit==0: out[a]|=(1<<b)
            else: out[b]|=(1<<a)
        # complement: reverse all arcs
        outc=[0]*(n+1)
        for i in range(1,n+1):
            for j in range(1,n+1):
                if (out[i]>>j)&1: outc[j]|=(1<<i)
        def feats(o):
            ti={t:count(ts[t],o) for t in trips}
            c3=sum(ti.values())
            c5=sum(count(fs[t],o) for t in fives)
            a2=sum(ti[trips[i]]*ti[trips[j]] for (i,j) in tpairs)
            return (c3,c5,a2)
        if feats(out)!=feats(outc): bad+=1
        checked+=1
    print(f"n={n}: {checked} tilings, complement violations={bad} -> "
          f"{'LOSSLESS 2x OK' if bad==0 else 'BROKEN'}", flush=True)

for n in (4,5,6):
    check(n)
