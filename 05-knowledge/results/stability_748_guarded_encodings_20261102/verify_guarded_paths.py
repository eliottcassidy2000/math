#!/usr/bin/env python3
"""Exhaustive graph/permutation versus linear-forest expansion controls."""
from itertools import combinations,permutations
from math import factorial,comb
import json

def need(ok,message):
    if not ok:raise RuntimeError(message)

def polynomial(n,E):
    counts=[0]*n
    for pi in permutations(range(n)):
        b=sum(tuple(sorted((u,v))) not in E for u,v in zip(pi,pi[1:]))
        counts[b]+=1
    return counts

def forests(n,E):
    E=sorted(E);C=[0]*(n+1)
    for mask in range(1<<len(E)):
        parent=list(range(n));deg=[0]*n;used=0;ok=True
        def root(v):
            while parent[v]!=v:v=parent[v]
            return v
        for i,(u,v) in enumerate(E):
            if not mask>>i&1:continue
            deg[u]+=1;deg[v]+=1;used+=1
            ru,rv=root(u),root(v)
            if max(deg[u],deg[v])>2 or ru==rv:ok=False;break
            parent[ru]=rv
        if ok:
            k=n-used;s=deg.count(0);C[k]+=2**(k-s)
    return C

def check(n,E):
    P=polynomial(n,E);C=forests(n,E);rhs=[0]*n
    for k in range(1,n+1):
        for j in range(n-k+1):
            rhs[k-1+j]+=factorial(k)*C[k]*comb(n-k,j)*(-1)**j
    need(P==rhs,'linear-forest expansion')
    pc=next(k for k in range(1,n+1) if C[k])
    need(next(i for i,x in enumerate(P) if x)==pc-1,'minimum-defect valuation')
    need(P[pc-1]==factorial(pc)*C[pc],'leading coefficient')
    complement=set(combinations(range(n),2))-E
    need(polynomial(n,complement)==P[::-1],'complement reciprocity')
    need(sum(P)==factorial(n),'fair-tournament collapse')
    return P,pc

def main():
    total=0
    for n in range(1,6):
        es=list(combinations(range(n),2))
        for mask in range(1<<len(es)):
            E={e for i,e in enumerate(es) if mask>>i&1};check(n,E);total+=1
    star={(0,1),(0,2),(0,3)};triangle={(1,2),(1,3),(2,3)}
    Ps,cs=check(4,star);Pt,ct=check(4,triangle)
    need(Ps==Pt==[0,12,12,0] and cs==ct==2,'complementary three-edge shapes')
    line={(0,1),(1,2),(2,3)}
    need(check(4,line)[0]==[2,10,10,2],'self-complementary path')
    attached=[check(5,E|{(0,4)}) for E in (star,triangle)]
    need([pc for P,pc in attached]==[3,2],'same port operation separates')
    path15=[8,1,15,10,6,3,13,12,4,5,11,14,2,7,9]
    need(sorted(path15)==list(range(1,16)),'inherited square-sum source labels')
    need(all(u+v in {9,16,25} for u,v in zip(path15,path15[1:])),'square-sum path15')
    print(json.dumps({'status':'passed','exhaustive_labelled_graphs_orders1_to5':total,
        'star_and_complement_polynomial':Ps,'path4_polynomial':[2,10,10,2],
        'marked_leaf_attachment_path_covers':[pc for P,pc in attached],
        'marked_leaf_attachment_polynomials':[P for P,pc in attached],
        'scope':'Exact finite controls of all-size identities proved in the note; no new infinite square-sum existence claim.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
