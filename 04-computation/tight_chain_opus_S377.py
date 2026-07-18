# opus-2026-07-17-S377 -- CAN THE 12->24 OPERATION BE ITERATED OR VARIED?
from fractions import Fraction as F
from functools import reduce
from math import gcd
import random
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
BASE=tuple(range(1,14))

print("(4) ITERATE: single substitutions from the NEW tight family {1..11,13,24}")
T2=[1,2,3,4,5,6,7,8,9,10,11,13,24]
hits=[]
for i in T2:
    for r in range(2,150):
        V=sorted([x for x in T2 if x!=i]+[r])
        if len(set(V))!=13: continue
        if tight(V): hits.append((i,r,tuple(V)))
print(f"    substitutions keeping tightness: {len(hits)}")
seen=set()
for i,r,V in hits:
    c=canon(V)
    tag = "= {1..13} up to dilation" if c==BASE else ("= {1..11,13,24} (back)" if c==canon(T2) else "NEW")
    if (i,r) not in seen:
        print(f"      {i:3d} -> {r:3d}:  {list(V)}   [{tag}]")
        seen.add((i,r))

print()
print("(5) DOES THE OPERATION WORK ON DILATES?  2*{1..13}: replace 24 by 48?")
for k in [2,3]:
    B=[k*i for i in range(1,14)]
    print(f"    {k}*{{1..13}} tight: {tight(B)}")
    found=[]
    for i in B:
        for m in [2,3]:
            V=sorted([x for x in B if x!=i]+[i*m])
            if len(set(V))!=13: continue
            if tight(V): found.append((i,i*m,canon(V)==BASE))
    print(f"      tight single-multiplications: {[(a,b,'dilate' if d else 'NEW') for a,b,d in found]}")

print()
print("(6) BROAD HILL-CLIMB for tight families (minimise uncovered to exactly 0)")
random.seed(377)
tights=set()
for trial in range(30):
    V=sorted(random.sample(range(1,70),13))
    if reduce(gcd,V)!=1: continue
    cur=uncovered(V); stall=0
    for step in range(300):
        W=list(V); i=random.randrange(13)
        W[i]=max(1,W[i]+random.choice([-3,-2,-1,1,2,3]))
        W=sorted(set(W))
        if len(W)!=13: continue
        u=uncovered(W)
        if u<=cur:
            if u<cur: stall=0
            V,cur=W,u
        else: stall+=1
        if cur==0: break
        if stall>200: break
    if cur==0: tights.add(canon(V))
print(f"    distinct tight families (up to dilation) found by hill-climb: {len(tights)}")
for c in sorted(tights):
    print(f"      {list(c)}   {'= {1..13}' if c==BASE else ('= {1..11,13,24}' if c==canon(T2) else 'NEW!')}")
