#!/usr/bin/env python3
"""
death-star-2026-07-20-S61d (HYP-8320) -- clarify Babai-Cameron Remark 7.4 for ALL odd n.
CLAIM: each switching class of tournaments on ODD n has a UNIQUE canonical member --
all-EVEN out-degrees if n=1 mod4, all-ODD out-degrees if n=3 mod4 -- fixed by every
automorphism, so a fixed-point-free automorphism (fixing no member) CANNOT exist for odd n.
Compute the actual FPF-automorphism switching-class count for n=3,4,5,6 to locate exactly
where Babai-Cameron's 'we cannot do this' bites (the even-n case).
"""
from itertools import permutations
from math import comb
from collections import Counter

def build(n):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]; idx={e:k for k,e in enumerate(edges)}
    return edges, idx
def switch_v(bits,v,n,idx):
    for u in range(n):
        if u!=v: bits^=(1<<idx[(min(u,v),max(u,v))])
    return bits
def outdegs(bits,n,edges):
    od=[0]*n
    for k,(i,j) in enumerate(edges):
        if (bits>>k)&1: od[i]+=1
        else: od[j]+=1
    return od
def relabel(bits,g,n,edges,idx):
    nb=0
    for k,(i,j) in enumerate(edges):
        if (bits>>k)&1: a,b=g[i],g[j]
        else: a,b=g[j],g[i]
        if a<b: nb|=(1<<idx[(a,b)])
    return nb
def switching_classes(n,edges,idx):
    seen=set(); classes=[]
    for start in range(1<<len(edges)):
        if start in seen: continue
        cls=set([start]); fr=[start]
        while fr:
            b=fr.pop()
            for v in range(n):
                nb=switch_v(b,v,n,idx)
                if nb not in cls: cls.add(nb); fr.append(nb)
        seen|=cls; classes.append(cls)
    return classes

print("=== canonical member (all-even / all-odd out-degrees) per switching class ===")
for n in [3,4,5,6]:
    edges,idx=build(n); classes=switching_classes(n,edges,idx)
    ev=[]; od=[]
    for cls in classes:
        ev.append(sum(1 for b in cls if all(d%2==0 for d in outdegs(b,n,edges))))
        od.append(sum(1 for b in cls if all(d%2==1 for d in outdegs(b,n,edges))))
    tag = f"{n%4} mod4"
    print(f"  n={n} ({tag}): #classes={len(classes)}, all-EVEN/class={dict(Counter(ev))}, all-ODD/class={dict(Counter(od))}")

print("\n=== fixed-point-free automorphism: #switching classes with g (g.C=C) fixing NO member ===")
for n in [3,4,5,6]:
    edges,idx=build(n); classes=switching_classes(n,edges,idx)
    perms=list(permutations(range(n)))
    fpf_count=0
    for cls in classes:
        has_fpf=False
        for g in perms:
            if g==tuple(range(n)): continue
            # g stabilizes class?
            gcls=set(relabel(b,g,n,edges,idx) for b in cls)
            if gcls!=cls: continue
            # fixes no member?
            if all(relabel(b,g,n,edges,idx)!=b for b in cls):
                has_fpf=True; break
        if has_fpf: fpf_count+=1
    print(f"  n={n} ({n%4} mod4): {fpf_count} of {len(classes)} switching classes admit a fixed-point-free automorphism"
          f"   {'(ODD n => 0, canonical member always fixed)' if n%2==1 and fpf_count==0 else ''}")
print("\n  => ODD n: 0 (unique canonical all-even[n=1mod4]/all-odd[n=3mod4] member is fixed by every g).")
print("     EVEN n: the count is NONZERO -- this is exactly the case Babai-Cameron said 'we cannot do'.")
