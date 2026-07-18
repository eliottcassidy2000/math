from fractions import Fraction as F
from functools import reduce
from math import gcd
from itertools import combinations
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))
def tight(V): return uncovered(V)==0
def canon(V):
    g=reduce(gcd,V); return tuple(sorted(x//g for x in V))
BASE=tuple(range(1,14)); T2=canon([1,2,3,4,5,6,7,8,9,10,11,13,24])
base=list(range(1,14))
print("(7) TWO-SPEED MULTIPLICATION on {1,...,13}:  i->i*m, j->j*m'")
found=set()
for i,j in combinations(base,2):
    for m in [2,3,4,5]:
        for mp in [2,3,4,5]:
            V=sorted([x for x in base if x not in (i,j)]+[i*m,j*mp])
            if len(set(V))!=13: continue
            if tight(V): found.add((i,m,j,mp,canon(V)))
print(f"    tight results: {len(found)}")
for i,m,j,mp,c in sorted(found):
    tag = "{1..13}" if c==BASE else ("{1..11,13,24}" if c==T2 else "*** NEW ***")
    print(f"      {i}->{i*m}, {j}->{j*mp}   canon = {list(c)}  [{tag}]")

print()
print("(8) ONE MULTIPLE + ONE FREE REPLACEMENT (r <= 30)")
found2=set()
for i in base:
    for m in [2,3]:
        for j in base:
            if j==i: continue
            for r in range(2,31):
                V=sorted([x for x in base if x not in (i,j)]+[i*m,r])
                if len(set(V))!=13: continue
                if tight(V): found2.add(canon(V))
print(f"    distinct tight canon classes: {len(found2)}")
for c in sorted(found2):
    tag = "{1..13}" if c==BASE else ("{1..11,13,24}" if c==T2 else "*** NEW ***")
    print(f"      {list(c)}  [{tag}]")
print()
print("(9) SUMMARY: distinct tight families found across ALL S377 searches")
allc = {BASE, T2} | {c for *_,c in found} | found2
new = [c for c in allc if c not in (BASE,T2)]
print(f"    total distinct (up to dilation): {len(allc)}")
print(f"      {list(BASE)}   -- classical")
print(f"      {list(T2)}   -- S376")
print(f"    genuinely additional: {len(new)}  {[list(c) for c in new]}")
