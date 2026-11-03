#!/usr/bin/env python3
"""Exhaustive +18 edge-flip BFS from one specified labelled source."""
from itertools import combinations

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def source():
    # Z0 = transitive base 0,...,4; 0..4 -> 5 -> 6 -> 0..4.
    a=[0]*7
    for i in range(5):
        for j in range(i+1,5):a[i]|=1<<j
        a[i]|=1<<5;a[6]|=1<<i
    a[5]|=1<<6
    return tuple(a)

def hamilton(a):
    n=len(a);d=[[0]*n for _ in range(1<<n)]
    for i in range(n):d[1<<i][i]=1
    for mask in range(1,1<<n):
        for i,v in enumerate(d[mask]):
            if not v:continue
            avail=a[i]&~mask
            while avail:
                bit=avail&-avail; j=bit.bit_length()-1
                d[mask|bit][j]+=v;avail-=bit
    return sum(d[-1])

S={source()};need(hamilton(next(iter(S)))==33,'initial H')
for depth,target,expected in ((1,51,3),(2,69,14),(3,87,0)):
    N=set()
    for a in S:
        for i,j in combinations(range(7),2):
            b=list(a);b[i]^=1<<j;b[j]^=1<<i;b=tuple(b)
            if hamilton(b)==target:N.add(b)
    need(len(N)==expected,f'depth {depth}: expected {expected}, saw {len(N)}')
    print(f'depth {depth}: count {target}, distinct labelled states {len(N)}')
    S=N
print('PASS: all edge reversals; no strongness filter; seven vertices throughout')
