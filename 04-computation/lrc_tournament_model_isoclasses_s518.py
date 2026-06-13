#!/usr/bin/env python3
"""oracle-2026-06-01-S518: realizable circular (half-turn) tournament iso-classes;
count = 2*Fib(m-2) (Fibonacci menu). See reflection lrc-as-a-tournament-and-its-fibonacci-isoclass-menu-s518.md."""
import random
from functools import lru_cache
from itertools import permutations, combinations

def half_turn(pts):
    n=len(pts); adj=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j: continue
            d=(pts[i]-pts[j])%1.0
            if 0<d<0.5: adj[i][j]=1
    return tuple(map(tuple,adj))

def canon(adj):
    n=len(adj); best=None
    for p in permutations(range(n)):
        f=tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f<best: best=f
    return best

def H(adj):
    n=len(adj); full=(1<<n)-1
    @lru_cache(None)
    def dp(m,l):
        if m==full: return 1
        return sum(dp(m|(1<<x),x) for x in range(n) if not (m>>x)&1 and adj[l][x])
    return sum(dp(1<<s,s) for s in range(n))

def n3(adj):
    n=len(adj);c=0
    for a,b,d in combinations(range(n),3):
        e=adj[a][b]+adj[b][d]+adj[d][a]
        if e in (0,3): c+=1
    return c
def reg(adj): return len(set(sum(r) for r in adj))==1
A000568={1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880}
rng=random.Random(1)
print("m | #realizable half-turn iso-classes | A000568(m) | regular? | H-values")
for m in range(3,8):
    raw=set()
    for _ in range(60000):
        pts=sorted(rng.random() for _ in range(m))
        raw.add(half_turn(pts))
    classes={canon(a) for a in raw}
    # per class info via one rep each
    info=[]
    seen=set()
    for a in raw:
        c=canon(a)
        if c in seen: continue
        seen.add(c); info.append((H(a), n3(a), reg(a)))
    Hs=sorted(set(h for h,_,_ in info))
    nreg=sum(1 for _,_,r in info if r)
    print(f"{m} | {len(classes):>3} | {A000568[m]:>4} | reg={nreg} | H={Hs}")

# --- m=8,9 confirmation ---
from itertools import permutations, combinations
from functools import lru_cache
def half_turn(pts):
    n=len(pts); adj=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j and 0<(pts[i]-pts[j])%1.0<0.5: adj[i][j]=1
    return tuple(map(tuple,adj))
def canon(adj):
    n=len(adj); best=None
    for p in permutations(range(n)):
        f=tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f<best: best=f
    return best
rng=random.Random(2)
for m in (8,9):
    raw=set()
    for _ in range(200000):
        raw.add(half_turn(sorted(rng.random() for _ in range(m))))
    classes={canon(a) for a in raw}
    import math
    # 2*Fib(m-2)
    def fib(k):
        a,b=1,1
        for _ in range(k-1): a,b=b,a+b
        return a
    pred=2*fib(m-2)
    print(f"m={m}: realizable iso-classes = {len(classes)}  | 2*Fib({m-2}) = {pred}  | match={len(classes)==pred}")
