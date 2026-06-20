#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Brute E[alpha2] (disjoint directed-triangle pairs) over ALL tilings, n=5,6,7.
Resolves the discrepancy between per-subset and claimed list 0,0,0,15/16,93/16,173/8."""
import sys
from itertools import combinations, permutations, product
from math import comb
from fractions import Fraction as Fr
if hasattr(sys.stdout,"reconfigure"): sys.stdout.reconfigure(encoding="utf-8")

def tiles(n): return [(a,b) for a in range(3,n+1) for b in range(1,a-1)]
def build_adj(n,T,bv):
    adj=[[0]*(n+1) for _ in range(n+1)]
    for k in range(n,1,-1): adj[k][k-1]=1
    for (a,b),bit in zip(T,bv):
        if bit==0: adj[a][b]=1
        else: adj[b][a]=1
    return adj
def is_dir_tri(adj,s):
    a,b,c=s
    return sorted([adj[a][b]+adj[a][c],adj[b][a]+adj[b][c],adj[c][a]+adj[c][b]])==[1,1,1]
def alpha2_tri(adj,n):
    tris=[frozenset(s) for s in combinations(range(1,n+1),3) if is_dir_tri(adj,s)]
    cnt=0
    for i in range(len(tris)):
        for j in range(i+1,len(tris)):
            if tris[i].isdisjoint(tris[j]): cnt+=1
    return cnt

for n in range(5,8):
    T=tiles(n); F=len(T); tot=1<<F; s=0
    for bv in product((0,1),repeat=F):
        adj=build_adj(n,T,bv); s+=alpha2_tri(adj,n)
    print(f" n={n}: brute E[alpha2_tri] = {Fr(s,tot)}   (claim list says n=6:15/16, n=7:93/16)")
