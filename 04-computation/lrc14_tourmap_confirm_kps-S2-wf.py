"""
lrc14_tourmap_confirm_kps-S2-wf.py

Final confirmation pass. Two goals:

(1) CROSS-CHECK: the residue-exhaustive M4 forbidden set (4 classes incl. regular)
    equals the forbidden set observed over ACTUAL LRC speed sets (speeds 1..N).
    If they match, the forbidden-class claim for M4 is VERIFIED both ways.

(2) M3 (halfplane-majority over OPTIMAL LONELY TAUS) forbidden set at n=5, with a
    LARGER speed range, to see if the regular class is forbidden by genuine lonely-time
    data too (this is the most LRC-intrinsic map: it uses the actual gap optimum).
    Also report n=4 to confirm full coverage there.

Exact arithmetic throughout.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations, permutations, product

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def all_opt_taus(S):
    b = F(0)
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v
    return [t for t in cand(S) if gmin(S, t) == b], b

def canon_key(adj, m):
    best = None
    for perm in permutations(range(m)):
        bits = tuple(adj[perm[i]][perm[j]] for i in range(m) for j in range(m) if i != j)
        if best is None or bits < best: best = bits
    return best
def h_count(adj, m):
    return sum(1 for p in permutations(range(m)) if all(adj[p[k]][p[k+1]] for k in range(m-1)))
def score_seq(adj, m):
    return tuple(sorted(sum(adj[i][j] for j in range(m) if j != i) for i in range(m)))
def num_3cycles(adj, m):
    c=0
    for a,b,cc in combinations(range(m),3):
        if adj[a][b] and adj[b][cc] and adj[cc][a]: c+=1
        if adj[a][cc] and adj[cc][b] and adj[b][a]: c+=1
    return c
def primitive(S):
    g=0
    for v in S: g=gcd(g,v)
    return g==1
def all_iso_classes(m):
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

MOD=14
U=[a for a in range(1,MOD) if gcd(a,MOD)==1]
def method4(S):
    m=len(S); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            wi=wj=0
            for a in U:
                ri=(S[i]*a)%MOD; rj=(S[j]*a)%MOD
                di=min(ri,MOD-ri); dj=min(rj,MOD-rj)
                if di>dj: wi+=1
                elif dj>di: wj+=1
            if wi>wj: adj[i][j]=1
            elif wj>wi: adj[j][i]=1
            else: adj[i][j]=1 if S[i]<S[j] else 0; adj[j][i]=1-adj[i][j]
    return adj
def method3(S):
    taus,_=all_opt_taus(S); m=len(S); adj=[[0]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            sc=0
            for t in taus:
                fi=(S[i]*t)%1; fj=(S[j]*t)%1
                wi=0 if fi<F(1,2) else 1; wj=0 if fj<F(1,2) else 1
                if wi!=wj: sc+= (1 if wi<wj else -1)
                else: sc+= (1 if fi<fj else -1) if fi!=fj else 0
            if sc>0: adj[i][j]=1
            elif sc<0: adj[j][i]=1
            else: adj[i][j]=1 if S[i]<S[j] else 0; adj[j][i]=1-adj[i][j]
    return adj

import sys
def flush(): sys.stdout.flush()

free={m:all_iso_classes(m) for m in [4,5]}

print("="*70, flush=True)
print("(1) M4 cross-check: forbidden over ACTUAL LRC speed sets vs residue-exhaustive")
for maxsp in [20, 30, 40]:
    realized=set(); cnt=0
    for S in combinations(range(1,maxsp+1),5):
        if not primitive(S): continue
        cnt+=1
        realized.add(canon_key(method4(S),5))
    forb=set(free[5])-realized
    print(f"  speeds 1..{maxsp}: {cnt} prim 5-sets; realized={len(realized)}/12; forbidden={len(forb)}", flush=True)
    print(f"     forbidden: {sorted(free[5][k] for k in forb)}", flush=True)
print("  (residue-exhaustive M4 forbade exactly: H9(1,1,2,3,3), H13(1,2,2,2,3), H15(1,2,2,2,3), H15-REGULAR)")

print()
print("="*70, flush=True)
print("(2) M3 (lonely-tau halfplane majority) forbidden set, n=5, growing speed range")
for maxsp in [14, 18, 22]:
    realized=set(); cnt=0
    for S in combinations(range(1,maxsp+1),5):
        if not primitive(S): continue
        cnt+=1
        realized.add(canon_key(method3(S),5))
    forb=set(free[5])-realized
    reg_forbidden = any(free[5][k][2]==(2,2,2,2,2) for k in forb)
    print(f"  speeds 1..{maxsp}: {cnt} prim 5-sets; realized={len(realized)}/12; forbidden={len(forb)}; regular-forbidden={reg_forbidden}", flush=True)
    print(f"     forbidden: {sorted(free[5][k] for k in forb)}", flush=True)

print()
print("="*70, flush=True)
print("(3) n=4 full coverage check for M3 and M4 (should be all 4 classes => no n=4 constraint)")
for name,fn in [("M3",method3),("M4",method4)]:
    realized=set()
    for S in combinations(range(1,30),4):
        if not primitive(S): continue
        realized.add(canon_key(fn(S),4))
    print(f"  {name} n=4: realized {len(realized)}/4  forbidden={sorted(free[4][k] for k in (set(free[4])-realized))}", flush=True)
print("DONE", flush=True)
