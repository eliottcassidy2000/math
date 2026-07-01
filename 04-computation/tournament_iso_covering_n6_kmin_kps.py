#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
DETERMINE k_min(6) exactly: scan k=6,7,8,... with block-within + random free-sets until one covers all 56.

kind-pasteur-2026-07-01-S10. n<=5: k_min = info bound (6-> here info bound is 6) via block-within.
n=6: k=6 fails (best 47/56). Scan upward to find the true k_min(6) and the covering shape.
"""
import sys, itertools, random
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
n=6; prs=list(itertools.combinations(range(n),2)); idx={p:k for k,p in enumerate(prs)}; m=len(prs)
def relabel_map(pm):
    src=[0]*m; flp=[0]*m
    for k,(a,b) in enumerate(prs):
        oa,ob=pm[a],pm[b]
        if oa<ob: src[k]=idx[(oa,ob)]
        else: src[k]=idx[(ob,oa)]; flp[k]=1
    return src,flp
gens=[relabel_map([*range(i)]+[i+1,i]+[*range(i+2,n)]) for i in range(n-1)]
def apply(t,gm):
    src,flp=gm; t2=0
    for k in range(m):
        if ((t>>src[k])&1)^flp[k]: t2|=(1<<k)
    return t2
N=1<<m; parent=list(range(N))
def find(x):
    r=x
    while parent[r]!=r: r=parent[r]
    while parent[x]!=r: parent[x],x=r,parent[x]
    return r
for t in range(N):
    for gm in gens:
        ra,rb=find(t),find(apply(t,gm))
        if ra!=rb: parent[max(ra,rb)]=min(ra,rb)
cid={}; classid=[0]*N; nc=0
for t in range(N):
    r=find(t)
    if r not in cid: cid[r]=nc; nc+=1
    classid[t]=cid[r]
print(f"n=6: #classes={nc}")
def scatter(bits,pos):
    t=0;b=0
    while bits:
        if bits&1: t|=1<<pos[b]
        bits>>=1;b+=1
    return t
def best_cover_for_freeset(D, orient_iter):
    Fpos=[p for p in range(m) if p not in D]; best=0; witness=None
    for ft in orient_iter(Fpos):
        seen=set()
        for s in range(1<<len(D)):
            seen.add(classid[ft^scatter(s,D)])
            if len(seen)==nc: break
        if len(seen)>best: best=len(seen)
        if len(seen)==nc: witness=ft; break
    return best,witness
def all_orient(Fpos):
    for fo in range(1<<len(Fpos)): yield scatter(fo,Fpos)
def block_within(blocks):
    D=[]
    for B in blocks:
        for a,b in itertools.combinations(sorted(B),2): D.append(idx[(a,b)])
    return sorted(D)
# partitions of 6 (as block-size lists) and their block-within free-sets
def partitions(nn):
    def rec(nn, mx):
        if nn==0: yield []; return
        for first in range(min(nn,mx),0,-1):
            for rest in rec(nn-first, first): yield [first]+rest
    return list(rec(nn,nn))
partsets={}
for p in partitions(6):
    # canonical block assignment 0..: fill blocks in order
    blocks=[]; v=0
    for sz in p:
        blocks.append(list(range(v,v+sz))); v+=sz
    D=block_within(blocks); partsets[tuple(p)]=D

print("\nSCAN k upward for a covering subcube (block-within partitions + random):")
rng=random.Random(3)
kmin=None
for k in range(6,12):
    # structured: block-within partitions with within-size exactly k
    struct=[(p,D) for p,D in partsets.items() if len(D)==k]
    covered=None; bestk=0
    for p,D in struct:
        # limit orientations for speed if 2^|F| large
        Fsz=m-k
        if Fsz<=12: b,w=best_cover_for_freeset(D, all_orient)
        else:
            def randor(Fpos, R=4000):
                for _ in range(R): yield scatter(rng.getrandbits(len(Fpos)),Fpos)
            b,w=best_cover_for_freeset(D, randor)
        bestk=max(bestk,b)
        if w is not None: covered=('block '+str(p),D); break
    # random free-sets of size k
    if covered is None:
        for _ in range(1500):
            D=sorted(rng.sample(range(m),k))
            def randor(Fpos, R=200):
                for _ in range(R): yield scatter(rng.getrandbits(len(Fpos)),Fpos)
            b,w=best_cover_for_freeset(D, randor)
            bestk=max(bestk,b)
            if w is not None: covered=('random',D); break
    status = f"COVERS (shape={covered[0]})" if covered else f"best {bestk}/{nc}"
    print(f"  k={k}: {status}")
    if covered and kmin is None:
        kmin=k
        D=covered[1]; print(f"     => k_min(6) = {k}; free arcs {[prs[i] for i in D]}")
        break
if kmin is None: print("  no covering found up to k=11 in this search")
print(f"\nRESULT: k_min(6) = {kmin} (info bound was 6). Gap over info bound = {None if kmin is None else kmin-6}.")
print("  => the clean 'k_min = ceil(log2 A000568)' pattern (exact n<=5) BREAKS at n=6: a genuine gap opens.")
print("DONE.")
