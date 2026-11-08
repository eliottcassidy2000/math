#!/usr/bin/env python3
"""Independent tournament/path enumerator; imports no producer code."""
from itertools import combinations, permutations, product
from functools import lru_cache
from fractions import Fraction as F
from collections import Counter
import json

def need(x,m):
    if not x: raise RuntimeError(m)
def det(M):
    n=len(M);out=0
    for p in permutations(range(n)):
        val=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        for i,j in enumerate(p): val*=M[i][j]
        out+=val
    return out

def graphs(n):
    pairs=list(combinations(range(n),2))
    for choices in product((0,1),repeat=len(pairs)):
        edges=set((i,j) if x else (j,i) for (i,j),x in zip(pairs,choices))
        yield edges

def data(n,edges,r):
    V=frozenset(range(n))
    @lru_cache(None)
    def tail(S,last):
        if not S:return 1
        return sum(tail(S-{v},v) for v in S if (last,v) in edges)
    @lru_cache(None)
    def hp(S):return sum(tail(S-{v},v) for v in S) if S else 1
    def sub(S):
        seq=sorted(S)
        for k in range(len(seq)+1):
            for a in combinations(seq,k):yield frozenset(a)
    def starts(S):return tail(S-{r},r) if r in S else 0
    def ends(S):
        return sum(all((p[i],p[i+1]) in edges for i in range(len(p)-1)) for p in permutations(S) if p[-1]==r) if r in S else 0
    H=hp(V);E=ends(V);U=sum(starts(S)*hp(V-S) for S in sub(V));q=hp(V-{r})
    Q=sum(hp(S)*hp(V-S) for S in sub(V));b=sum(hp(S)*hp(V-{r}-S) for S in sub(V-{r}))
    alpha=U;delta=sum(ends(S)*hp(V-S) for S in sub(V))
    T=sum(starts(S)*hp(W)*hp(V-S-W) for S in sub(V) for W in sub(V-S))
    return (H,E,U,q),(H,Q,starts(V),E,alpha,delta,T,q,b)

def matrix(d,n):
    H,Q,s,e,a,z,T,q,b=d
    return [[H,-H,Q-H,2*H-Q],[e,-e,z-e,e-z],[a,s-a,2*T-a,3*a-s-2*T], [q,-q,b-q,2*q-b] if n>1 else [1,0,0,0]]
def apply(M,x):return [sum(a*b for a,b in zip(row,x)) for row in M]
def reachable(n,edges,r):
    S={r}
    while True:
        T=S|{j for i,j in edges if i in S}
        if T==S:return len(S)==n
        S=T

def recover(st,x):
    H,Q,s,e,a,z,T,q,b=st;D1,E1,U1,d1=map(F,x);delta=q*Q-H*b
    t=(q*D1-H*d1)/delta;w=((Q-H)*d1-(b-q)*D1)/delta
    p=(E1-(z-e)*t)/e;d=w-p;E=(U1-a*p-(2*T-a)*t-(2*a-s)*d)/s
    return [p+E,E,t+d,d]

bank=[(1,set()),(2,{(0,1)}),(3,{(0,1),(1,2),(2,0)}),(3,{(0,1),(0,2),(1,2)})]
X=[list(data(n,e,0)[0]) for n,e in bank]
need(det(X)==-4,'physical span')
C=[[1,0,0,0],[1,-1,1,0],[1,0,2,0],[1,-1,3,-2]]
need(det(C)==4,'observer')
counts=Counter();minfull=100;signs=Counter()
for n in range(1,6):
    for edges in graphs(n):
        for r in range(n):
            x,st=data(n,edges,r);H,Q,s,e,a,z,T,q,b=st;M=matrix(st,n);delta=q*Q-H*b
            need(bool(s)==reachable(n,edges,r),'start reachability')
            need(bool(e)==reachable(n,{(j,i) for i,j in edges},r),'end reachability')
            need((x[2]==x[3]) == (not any(i==r for i,j in edges)), 'sink split')
            strong=all(reachable(n,edges,k) for k in range(n))
            need(bool(s*e)==strong,'strong equivalence')
            if n>1:
                need(det(M)==s*e*delta,'determinant')
                if strong and not delta:
                    need(any(det([[M[i][j] for j in cols] for i in rows]) for rows in combinations(range(4),3) for cols in combinations(range(4),3)),'rank three')
                if det(M):
                    minfull=min(minfull,n)
                    for j in range(4):
                        v=[int(i==j) for i in range(4)]
                        need(recover(st,apply(M,v))==v,'inverse all-space basis')
                        counts['inverse_basis_vectors']+=1
                signs[f'n={n},strong={strong},delta_sign={(delta>0)-(delta<0)}']+=1
            else:need(det(M)==0,'singleton singular')
            if n<=4:
                for (m,source),v in zip(bank,X):
                    total=source|{(i+m,j+m) for i,j in edges}|{(i+m,0) for i in range(n)}|{(i,j+m) for i in range(1,m) for j in range(n)}
                    out,_=data(m+n,total,m+r)
                    need(list(out)==apply(M,v),'independent literal graft')
                    counts['literal_grafts_all_roots']+=1
            counts['rooted_graphs']+=1
need(minfull==4,'minimal full-rank order')
print(json.dumps({'status':'PASS','counts':dict(counts),'minimum_full_rank_order':minfull,'physical_basis':X,'physical_determinant':det(X),'observer_determinant':det(C),'all_root_sign_census':dict(sorted(signs.items()))},indent=2,sort_keys=True))
