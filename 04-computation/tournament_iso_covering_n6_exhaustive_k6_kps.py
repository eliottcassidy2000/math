#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""EXHAUSTIVE proof that no k=6 subcube covers all 56 iso classes at n=6 => k_min(6)>=7 (rigorous).
Scans all C(15,6)=5005 free-sets x 2^9=512 orientations; reports the global best coverage."""
import sys, itertools
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
cid={}; classid=bytearray(0 for _ in range(N)); nc=0
for t in range(N):
    r=find(t)
    if r not in cid: cid[r]=nc; nc+=1
    classid[t]=cid[r]
print(f"n=6 classes={nc}")
def scatter(bits,pos):
    t=0;b=0
    while bits:
        if bits&1: t|=1<<pos[b]
        bits>>=1;b+=1
    return t
allp=list(range(m)); best=0; bestcfg=None
for D in itertools.combinations(allp,6):
    dm=[scatter(s,D) for s in range(64)]
    Fpos=[p for p in allp if p not in D]
    fmasks=[scatter(fo,Fpos) for fo in range(512)]
    for ft in fmasks:
        seen=set()
        add=seen.add
        for d in dm: add(classid[ft^d])
        ls=len(seen)
        if ls>best:
            best=ls; bestcfg=(D,ft)
            if best==nc: break
    if best==nc: break
print(f"EXHAUSTIVE k=6: global best coverage = {best}/{nc}")
print(f"  => k_min(6) {'= 6 (covering k=6 exists!)' if best==nc else '>= 7 PROVEN (no k=6 subcube covers)'}")
print("DONE.")
