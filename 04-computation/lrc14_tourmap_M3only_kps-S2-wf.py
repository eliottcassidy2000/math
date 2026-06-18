"""M3-only forbidden-class check (lonely-tau halfplane majority). Fast, flushed."""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product
import sys

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def gmin(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def all_opt_taus(S):
    b=F(0)
    for t in cand(S):
        v=gmin(S,t)
        if v>b: b=v
    return [t for t in cand(S) if gmin(S,t)==b]
def canon_key(adj,m):
    best=None
    for p in permutations(range(m)):
        bits=tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m) if i!=j)
        if best is None or bits<best: best=bits
    return best
def h_count(adj,m): return sum(1 for p in permutations(range(m)) if all(adj[p[k]][p[k+1]] for k in range(m-1)))
def score_seq(adj,m): return tuple(sorted(sum(adj[i][j] for j in range(m) if j!=i) for i in range(m)))
def num_3cycles(adj,m):
    c=0
    for a,b,cc in combinations(range(m),3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c+=1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c+=1
    return c
def primitive(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def all_iso(m):
    seen={}
    pairs=list(combinations(range(m),2))
    for bits in product([0,1],repeat=len(pairs)):
        adj=[[0]*m for _ in range(m)]
        for (i,j),b in zip(pairs,bits):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        kk=canon_key(adj,m)
        if kk not in seen: seen[kk]=(h_count(adj,m),num_3cycles(adj,m),score_seq(adj,m))
    return seen
def method3(S):
    taus=all_opt_taus(S); m=len(S); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            sc=0
            for t in taus:
                fi=(S[i]*t)%1; fj=(S[j]*t)%1
                wi=0 if fi<F(1,2) else 1; wj=0 if fj<F(1,2) else 1
                if wi!=wj: sc+=(1 if wi<wj else -1)
                else: sc+=(1 if fi<fj else -1) if fi!=fj else 0
            if sc>0: adj[i][j]=1
            elif sc<0: adj[j][i]=1
            else: adj[i][j]=1 if S[i]<S[j] else 0; adj[j][i]=1-adj[i][j]
    return adj

free5=all_iso(5); free4=all_iso(4)
print("M3 (lonely-tau halfplane majority) forbidden-class check", flush=True)
for maxsp in [16,20,24]:
    realized=set(); cnt=0
    for S in combinations(range(1,maxsp+1),5):
        if not primitive(S): continue
        cnt+=1
        realized.add(canon_key(method3(S),5))
    forb=set(free5)-realized
    reg=any(free5[k][2]==(2,2,2,2,2) for k in forb)
    print(f"  n=5 speeds 1..{maxsp}: {cnt} sets; realized={len(realized)}/12; forbidden={len(forb)}; regular-forbidden={reg}", flush=True)
    print(f"     forbidden: {sorted(free5[k] for k in forb)}", flush=True)
# n=4
realized=set()
for S in combinations(range(1,26),4):
    if not primitive(S): continue
    realized.add(canon_key(method3(S),4))
print(f"  n=4 speeds 1..25: realized={len(realized)}/4; forbidden={sorted(free4[k] for k in (set(free4)-realized))}", flush=True)
print("DONE", flush=True)
