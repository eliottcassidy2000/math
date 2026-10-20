#!/usr/bin/env python3
"""Independent ascending-path classifier; exhaustive vertex permutations."""
import itertools as it

def adjacency(n,x,free):
    a=[[0]*n for _ in range(n)]
    for i in range(n-1):a[i][i+1]=1
    for k,(i,j) in enumerate(free):
        if x&(1<<k):a[i][j]=1
        else:a[j][i]=1
    return a

def classify(n):
    free=[(i,j) for i in range(n) for j in range(i+2,n)]
    pairs=list(it.combinations(range(n),2))
    perms=list(it.permutations(range(n)))
    def code(a,p):
        return sum(a[p[i]][p[j]]*2**k for k,(i,j) in enumerate(pairs))
    mask=2**len(pairs)-1
    keys=[];selfconverse=[];As=[]
    for x in range(2**len(free)):
        a=adjacency(n,x,free); As.append(a)
        orbit=[code(a,p) for p in perms]
        ck=min(orbit);co=mask-max(orbit)
        keys.append(min(ck,co));selfconverse.append(ck==co)
    codes=sorted(set(keys));pi=[codes.index(k) for k in keys]
    return free,codes,pi,selfconverse,As

