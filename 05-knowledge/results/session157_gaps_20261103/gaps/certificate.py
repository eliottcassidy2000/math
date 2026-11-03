#!/usr/bin/env python3
"""Exact finite spectrum certificate; standard library only."""
from itertools import combinations, permutations
from collections import Counter

def graph(n, mask):
    # Directed path 0->1->...->n-1 is fixed. All other pairs vary.
    edges = [(i,j) for i in range(n) for j in range(i+2,n)]
    a = [[0]*n for _ in range(n)]
    for i in range(n-1): a[i][i+1] = 1
    for b,(i,j) in enumerate(edges):
        if mask>>b&1: a[i][j]=1
        else: a[j][i]=1
    return a

def strong(a):
    n=len(a)
    R=[sum(a[i][j]<<j for j in range(n)) | 1<<i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if R[i]>>k&1: R[i] |= R[k]
    return all(x==(1<<n)-1 for x in R)

def direct(a):
    n=len(a)
    return sum(all(a[p[i]][p[i+1]] for i in range(n-1)) for p in permutations(range(n)))

def odd_cycle_count(a):
    n=len(a)
    # Every odd cycle uses 3 or 5 vertices for n<=6. Fix minimum first,
    # then inspect all orders of the other vertices exactly once.
    c3=[]; c5=0
    for C in combinations(range(n),3):
        u,v,w=C
        if (a[u][v] and a[v][w] and a[w][u]) or (a[u][w] and a[w][v] and a[v][u]): c3.append(set(C))
    for C in combinations(range(n),5):
        u=C[0]
        for p in permutations(C[1:]):
            q=(u,)+p+(u,)
            if all(a[q[i]][q[i+1]] for i in range(5)): c5 += 1
    d=sum(x.isdisjoint(y) for x,y in combinations(c3,2))
    return 1+2*(len(c3)+c5)+4*d

for n in range(1,7):
    allh=Counter(); sh=Counter(); total=1<<((n-1)*(n-2)//2) if n>=2 else 1
    for mask in range(total):
        a=graph(n,mask); h=direct(a); h2=odd_cycle_count(a)
        if h!=h2: raise RuntimeError((n,mask,h,h2))
        allh[h]+=1
        if strong(a): sh[h]+=1
    print('n=',n,'codes=',total,'all=',sorted(allh.items()),'strong=',sorted(sh.items()),sep='')
    if any(h in allh for h in (7,21)): raise RuntimeError(('gap',n))
print('All checks passed.')
