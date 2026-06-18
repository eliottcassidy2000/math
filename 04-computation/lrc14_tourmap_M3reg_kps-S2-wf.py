from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations
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
def opt(S):
    b=F(0)
    for t in cand(S):
        v=gmin(S,t)
        if v>b: b=v
    return [t for t in cand(S) if gmin(S,t)==b]
def score(adj,m): return tuple(sorted(sum(adj[i][j] for j in range(m) if j!=i) for i in range(m)))
def m3(S):
    taus=opt(S); m=len(S); adj=[[0]*m for _ in range(m)]
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
# search for regular (2,2,2,2,2) under M3 over speeds 1..28
found=0; total=0
for S in combinations(range(1,29),5):
    g=0
    for v in S: g=gcd(g,v)
    if g!=1: continue
    total+=1
    if score(m3(S),5)==(2,2,2,2,2):
        found+=1
        if found<=3: print("REGULAR via M3:",S,flush=True)
print(f"M3 over {total} primitive 5-sets in 1..28: regular found = {found}",flush=True)
