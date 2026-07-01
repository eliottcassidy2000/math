#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
SETTLE k_min(6): does the clean pattern (k_min = info bound, via block-within free-sets) survive to n=6?

kind-pasteur-2026-07-01-S10. Exhaustive n<=5 gave k_min=info-bound via block-within free-sets. But the
balanced (3,3) block-within free-set covers only 42/56 at n=6 (companion script). So EITHER some other
6-arc free-set covers (k_min(6)=6, but NOT block-within => the clean shape breaks at n=6), OR no 6-set
covers (k_min(6)>=7 => info bound not achieved at n=6).  Settle via fast union-find iso-classes.
"""
import sys, itertools, random
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass
n=6; prs=list(itertools.combinations(range(n),2)); idx={p:k for k,p in enumerate(prs)}; m=len(prs)  # m=15
# generators of S_n: adjacent transpositions; precompute bit relabel (src,flip) per generator
def relabel_map(pm):
    src=[0]*m; flp=[0]*m
    for k,(a,b) in enumerate(prs):
        oa,ob=pm[a],pm[b]
        if oa<ob: src[k]=idx[(oa,ob)]; flp[k]=0
        else:     src[k]=idx[(ob,oa)]; flp[k]=1
    return src,flp
gens=[]
for i in range(n-1):
    pm=list(range(n)); pm[i],pm[i+1]=pm[i+1],pm[i]; gens.append(relabel_map(pm))
def apply(t,gm):
    src,flp=gm; t2=0
    for k in range(m):
        if ((t>>src[k])&1)^flp[k]: t2|=(1<<k)
    return t2
# union-find over all 2^m tournaments; orbits = iso classes
N=1<<m; parent=list(range(N))
def find(x):
    r=x
    while parent[r]!=r: r=parent[r]
    while parent[x]!=r: parent[x],x=r,parent[x]
    return r
for t in range(N):
    for gm in gens:
        t2=apply(t,gm)
        ra,rb=find(t),find(t2)
        if ra!=rb: parent[max(ra,rb)]=min(ra,rb)
cid={}; classid=[0]*N; nc=0
for t in range(N):
    r=find(t)
    if r not in cid: cid[r]=nc; nc+=1
    classid[t]=cid[r]
print(f"n=6: #iso classes = {nc} (expect A000568(6)=56: {nc==56})")

def scatter(bits,pos):
    t=0; b=0
    while bits:
        if bits&1: t|=(1<<pos[b])
        bits>>=1; b+=1
    return t
def subcube_covers(Dpos,fixed_t):
    seen=set(); k=len(Dpos)
    for s in range(1<<k):
        t=fixed_t^scatter(s,Dpos)
        c=classid[t]
        if c not in seen:
            seen.add(c)
            if len(seen)==nc: return True
    return False

# --- k=6 search: try to find ANY free-set + orientation covering all 56 ---
print("\nSEARCH k=6 (lower bound: 2^5=32<56 => k>=6). Looking for a covering 6-arc subcube...")
allpos=list(range(m))
# structured free-sets: block-within for partitions with within-count 6: (3,3) and (4,1,1)
def block_within(blocks):
    D=[]
    for B in blocks:
        for a,b in itertools.combinations(sorted(B),2): D.append(idx[(a,b)])
    return sorted(D)
structured=[
  ("(3,3) blocks {0,1,2}{3,4,5}", block_within([[0,1,2],[3,4,5]])),
  ("(4,1,1) block {0,1,2,3}",     block_within([[0,1,2,3],[4],[5]])),
]
found=None; best=0
for label,D in structured:
    Fpos=[p for p in allpos if p not in D]
    for fo in range(1<<len(Fpos)):
        ft=scatter(fo,Fpos)
        seen=set()
        for s in range(1<<len(D)):
            seen.add(classid[ft^scatter(s,D)])
        if len(seen)>best: best=len(seen)
        if len(seen)==nc: found=(label,D,ft); break
    if found: break
    print(f"  structured {label}: best coverage over all orientations = {best}/{nc}")
# random 6-arc free-sets
if not found:
    rng=random.Random(0)
    for trial in range(3000):
        D=sorted(rng.sample(allpos,6)); Fpos=[p for p in allpos if p not in D]
        ft=scatter(rng.getrandbits(len(Fpos)),Fpos)
        seen=set()
        for s in range(1<<6):
            seen.add(classid[ft^scatter(s,D)])
        if len(seen)>best: best=len(seen)
        if len(seen)==nc: found=("random",D,ft); break
    print(f"  random 6-arc search (3000 tries): best coverage = {best}/{nc}")
if found:
    label,D,ft=found
    print(f"  ==> COVERING 6-arc subcube FOUND ({label}); k_min(6)=6. free arcs {[prs[i] for i in D]}")
else:
    print(f"  ==> NO covering 6-arc subcube found (best {best}/{nc}); strong evidence k_min(6)>=7.")
    # confirm k=7 achievable
    rng=random.Random(1); ok7=False
    for trial in range(2000):
        D=sorted(rng.sample(allpos,7)); Fpos=[p for p in allpos if p not in D]
        ft=scatter(rng.getrandbits(len(Fpos)),Fpos)
        seen=set()
        for s in range(1<<7):
            seen.add(classid[ft^scatter(s,D)])
            if len(seen)==nc: break
        if len(seen)==nc: ok7=True; break
    print(f"  k=7 achievable? {ok7}  => k_min(6) = {'7' if (not found and ok7) else '6 or 7 (see above)'}")
print("DONE.")
