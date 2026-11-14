#!/usr/bin/env python3
"""Finite replay of the all-order SCC vertex/edge/converse identities."""
from collections import Counter, defaultdict
from functools import cache
from itertools import permutations, product
from math import comb
import sys

@cache
def canon(r):
    d = defaultdict(list)
    for i,v in enumerate(r):
        d[v.bit_count()].append(i)
    return min(relabel(r, sum(p, ())) for p in
               product(*(permutations(d[k]) for k in sorted(d))))

def relabel(r,p):
    return tuple(sum(((r[i]>>j)&1)<<k for k,j in enumerate(p)) for i in p)

def dual(r):
    F=(1<<len(r))-1
    return tuple(F^v^(1<<i) for i,v in enumerate(r))

def add(r,u):
    n=len(r)
    return tuple(x if u>>i&1 else x|(1<<n) for i,x in enumerate(r))+(u,)

@cache
def h(r):
    n=len(r)
    if n==0: return 1
    out=[[0]*n for _ in range(1<<n)]
    for i in range(n): out[1<<i][i]=1
    for s,v in enumerate(out):
        for i,a in enumerate(v):
            if a:
                for j in range(n):
                    if r[i]>>j&1 and not (s>>j&1):
                        out[s|(1<<j)][j]+=a
    return sum(out[-1])

def blocks(r):
    g=[a|1<<i for i,a in enumerate(r)]
    n=len(r)
    for i in range(n):
        for j in range(n):
            if g[j]>>i&1: g[j]|=g[i]
    d=defaultdict(list)
    for i,v in enumerate(g): d[v].append(i)
    return sorted(d.values(), key=lambda v:-g[v[0]].bit_count())

def join(w):
    r=()
    for b in w:
        new=((1<<len(b))-1)<<len(r)
        r=tuple(v|new for v in r)+tuple(a<<len(r) for a in b)
    return r

def pointing(r):
    aut=[]; anti=[]
    for p in permutations(range(len(r))):
        t=relabel(r,p)
        if t==r: aut.append(p)
        if t==dual(r): anti.append(p)
    for i in range(len(r)):
        if min(p[i] for p in aut)==i:
            fixed=any(p[i]==i for p in anti)
            yield i, fixed

def plus(c,d,fn):
    for a,x in c.items():
        for b,y in d.items(): yield fn(a,b),x*y

def main(N):
    V=[Counter({1:1})]+[Counter() for _ in range(N)]
    S=[Counter() for _ in range(N+1)]
    C=[Counter() for _ in range(N+1)]
    W=[Counter() for _ in range(N+1)]
    R=[Counter() for _ in range(N+1)]
    RC=[Counter() for _ in range(N+1)]
    layer={()}
    for n in range(1,N+1):
        nxt={canon(add(x,u)) for x in layer for u in range(1<<(n-1))}
        layer=nxt
        for x in layer:
            t=h(x); V[n][t]+=1; B=blocks(x); strong=len(B)==1
            if strong:S[n][t]+=1
            if canon(dual(x))==x:C[n][t]+=1
            for i,fixed in pointing(x):
                p=canon(relabel(x,tuple(j for j in range(n) if j!=i)))
                a=h(p); assert a<=t
                if a==t: assert any(b==[i] for b in B)
                else: assert not any(b==[i] for b in B)
                W[n][a,t,fixed]+=1
                if strong:
                    R[n][a,t]+=1
                    if fixed: RC[n][a,t]+=1
    P=[Counter({1:1})]+[Counter() for _ in range(N)]
    for n in range(1,N+1):
        for i in range(1,n+1):
            for k,v in plus(S[i],P[n-i],lambda x,y:x*y):P[n][k]+=v
        assert P[n]==V[n]
        e=Counter(); ec=Counter()
        for i in range(n):
            for j in range(n-i):
                for (a,b),m in R[n-i-j].items():
                    for x,l in P[i].items():
                        for y,r in P[j].items():
                            e[a*x*y,b*x*y]+=m*l*r
        for i in range((n-1)//2+1):
            for (a,b),m in RC[n-2*i].items():
                for x,l in P[i].items(): ec[a*x*x,b*x*x]+=m*l
        assert ec==Counter({(a,b):v for (a,b,f),v in W[n].items() if f})
        assert e==Counter({(a,b):W[n][a,b,True]+W[n][a,b,False]
                           for a,b,f in W[n]})
        assert all((v+ec[k])%2==0 for k,v in e.items())
        for target in (1,3,5,9,15):
            # Independent gap-position count against the class catalog.
            ans=0
            def words(t,b,r):
                nonlocal ans
                if t==target and b<=n: ans+=comb(n-b+r,r)
                for i in range(3,n-b+1):
                    for k,m in S[i].items():
                        if target%(t*k)==0:
                            old=ans
                            words(t*k,b+i,r+1)
                            ans=old+(ans-old)*m
            words(1,0,0)
            assert ans==V[n][target]
        print(n,sum(V[n].values()),len(V[n]),sum(e.values()),sum(ec.values()))
    if N>=6:
        for d, expected in [
            (2,[(3,(1,3),1),(4,(3,5),2)]),
            (4,[(4,(1,5),2),(5,(5,9),5)])]:
            got=[(n,k,v) for n in range(1,7) for k,v in R[n].items() if k[1]-k[0]==d]
            assert sorted(got)==expected
    for x,t in [((14,28,24,16,1),[1,5,5,5,1]),((10,28,25,16,1),[1,3,5,5,3])]:
        assert h(x)==9
        assert [h(relabel(x,tuple(j for j in range(5) if j!=i))) for i in range(5)]==t

if __name__=='__main__':
    print('order classes distinct-H pointed-edges self-converse-pointed-edges')
    main(int(sys.argv[1]) if len(sys.argv)>1 else 6)
