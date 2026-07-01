#!/usr/bin/env python3
"""
flip_rank_settle_n6_klein.py  --  klein-2026-07-01-S71 (part 3)

Settle rho(6) and test whether the achiever is still BIPARTITE (unbalanced split).
(C) test every vertex-split (a,b) of n=6,7: does fixing that cut + freeing the sides realize all classes?
    -> the minimal f over realizing splits = the "bipartite flip-rank".
(D) EXHAUSTIVE: does ANY 6-edge subcube realize G_6? (settles rho(6)=6 vs 7 rigorously).
"""
import itertools, math

def edges(n): return [(i,j) for i in range(n) for j in range(i+1,n)]
def perm_maps(n,E):
    idx={e:k for k,e in enumerate(E)}; maps=[]
    for pi in itertools.permutations(range(n)):
        m=[(idx[(pi[i],pi[j])],0) if pi[i]<pi[j] else (idx[(pi[j],pi[i])],1) for (i,j) in E]
        maps.append(m)
    return maps
def canon(t,maps,ne):
    best=None
    for m in maps:
        v=0
        for k in range(ne):
            bit=(t>>k)&1; tk,fl=m[k]
            if fl: bit^=1
            v|=bit<<tk
        if best is None or v<best: best=v
    return best
def classes_of(n):
    E=edges(n); ne=len(E); maps=perm_maps(n,E); cid={}; cls=[0]*(1<<ne); nc=0
    for t in range(1<<ne):
        c=canon(t,maps,ne)
        if c not in cid: cid[c]=nc; nc+=1
        cls[t]=cid[c]
    return E,ne,cls,nc

def split_realizes(n,a,E,ne,cls,nc):
    idx={e:k for k,e in enumerate(E)}
    A=set(range(a))
    cross=[idx[e] for e in E if (e[0] in A) != (e[1] in A)]
    within=[idx[e] for e in E if (e[0] in A)==(e[1] in A)]
    f=len(within)
    for fa in range(1<<len(cross)):
        fb=0
        for bit,e in enumerate(cross):
            if (fa>>bit)&1: fb|=1<<e
        seen=set()
        for wa in range(1<<f):
            t=fb
            for bit,e in enumerate(within):
                if (wa>>bit)&1: t|=1<<e
            seen.add(cls[t])
            if len(seen)==nc: break
        if len(seen)==nc: return f,True
    return f,False

def exhaustive_k(n,k,E,ne,cls,nc):
    """does ANY k-edge free-set (with some fixed orientation) realize all nc classes?"""
    alledges=list(range(ne)); full=(1<<ne)-1
    for free in itertools.combinations(alledges,k):
        fmask=0
        for e in free: fmask|=1<<e
        keep=full & ~fmask
        buckets={}
        ok=False
        for t in range(1<<ne):
            fixed=t & keep
            s=buckets.get(fixed)
            if s is None: s=set(); buckets[fixed]=s
            s.add(cls[t])
            if len(s)==nc: ok=True; break
        if ok: return True, free
    return False, None

if __name__=="__main__":
    print("(C) bipartite over ALL splits (a,b): minimal f that realizes")
    for n in [6,7]:
        E,ne,cls,nc=classes_of(n)
        lb=math.ceil(math.log2(nc))
        print(f"  n={n} |G|={nc} LB={lb}:")
        best=None
        for a in range(1,n//2+1):
            f,ok=split_realizes(n,a,E,ne,cls,nc)
            print(f"    split {a}+{n-a}: f={f}  realizes={ok}")
            if ok and (best is None or f<best): best=f
        print(f"    -> minimal BIPARTITE flip-rank for n={n}: {best} (LB={lb}, naive C(n-1,2)={math.comb(n-1,2)})")
    print("\n(D) EXHAUSTIVE: does any 6-edge subcube realize G_6? (rho(6)=6 vs 7)")
    E,ne,cls,nc=classes_of(6)
    ok,free=exhaustive_k(6,6,E,ne,cls,nc)
    print(f"    n=6 k=6 realizing subcube exists: {ok}" + (f" (free={[E[e] for e in free]})" if ok else " => rho(6) >= 7 (so rho(6)=7, LB+1)"))
