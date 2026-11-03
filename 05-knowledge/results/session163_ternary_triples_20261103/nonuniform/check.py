#!/usr/bin/env python3
"""Finite endpoint-subset controls; run normally and with python -O."""
from itertools import product
from collections import Counter

def count(t, bits, split=False):
    n=3*t+1; x=n-1
    def edge(a,b):
        if a==x: return bool(bits[b])
        if b==x: return not bits[a]
        if a//3!=b//3: return a//3 < b//3
        return (a%3,b%3) in {(0,1),(1,2),(2,0)}
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v,c in enumerate(dp[mask]):
            if c:
                for w in range(n):
                    if not(mask>>w&1) and edge(v,w): dp[mask|1<<w][w]+=c
    return sum(dp[-1])

def split_counts(bits):
    import itertools
    def edge(a,b):
        if a//3!=b//3:return a//3<b//3
        return (a%3,b%3) in {(0,1),(1,2),(2,0)}
    result=[0]*7
    for p in itertools.permutations(range(7)):
        def e(a,b):
            if a==6:return bool(bits[b])
            if b==6:return not bits[a]
            return edge(a,b)
        if all(e(a,b) for a,b in zip(p,p[1:])):result[p.index(6)]+=1
    return result

def require(p,msg):
    if not p: raise RuntimeError(msg)

hist1=Counter(count(1,b) for b in product((0,1),repeat=3))
require(hist1=={3:2,5:6},str(hist1))
hist2=Counter(count(2,b) for b in product((0,1),repeat=6))
require(count(2,(1,0,0,1,0,0))==41,'41')
require(split_counts((1,0,0,1,0,0))==[3,5,7,8,7,5,6],'split')
hist3=Counter(count(3,b) for b in product((0,1),repeat=9))
require(all(h%3==0 for h in hist3),'three triangles')
print('one triangle:', sorted(hist1.items()))
print('two triangles:', sorted(hist2.items()))
print('three triangles:', len(hist3),'distinct values, all divisible by 3')

# Independently count the two-slot profile of the sharpness gadget and
# the actual nine-vertex strong wrapper, not merely the product formula.
def hp_subset(adj):
    n=len(adj); N=1<<n; dp=[[0]*n for _ in range(N)]
    for v in range(n):dp[1<<v][v]=1
    for m in range(1,N):
        for v,c in enumerate(dp[m]):
            if c:
                for w in range(n):
                    if not(m>>w&1) and adj[v][w]: dp[m|1<<w][w]+=c
    return [1]+[sum(dp[m]) for m in range(1,N)]

def gadget():
    n=7;x=6;adj=[[False]*n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            if a==b:continue
            if a==x:adj[a][b]=b in (0,3)
            elif b==x:adj[a][b]=a not in (0,3)
            elif a//3!=b//3:adj[a][b]=a//3<b//3
            else:adj[a][b]=(a%3,b%3) in {(0,1),(1,2),(2,0)}
    return adj

g=gadget(); h=hp_subset(g); F=sum(h[m]*h[127^m] for m in range(128))
require((h[-1],F)==(41,712),'profile')
N=9;adj=[[False]*N for _ in range(N)]
for a in range(N):
    for b in range(N):
        if a==b:continue
        if a<7 and b<7:adj[a][b]=g[a][b]
        elif a<7:adj[a][b]=(b==7)
        elif b<7:adj[a][b]=(a==8)
        else:adj[a][b]=(a,b)==(7,8)
require(hp_subset(adj)[-1]==753,'strong wrapper')
print('41 gadget profile (H,F2)=(41,712); wrapper H=753')
