#!/usr/bin/env python3
"""Bounded independent controls for the one-slack theorem; no poset census."""
import json
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from math import factorial

checks=0
def require(test,message):
    global checks
    checks+=1
    if not test: raise RuntimeError(message)

def closure(n,edges):
    pred=[0]*n
    for a,b in edges: pred[b]|=1<<a
    for k in range(n):
        for b in range(n):
            if pred[b]>>k&1: pred[b]|=pred[k]
    require(all(not(pred[v]>>v&1) for v in range(n)), 'cycle')
    return pred

def structure(pred):
    n=len(pred)
    up=[sum(1<<j for j in range(n) if pred[j]>>i&1) for i in range(n)]
    inc=[((1<<n)-1) & ~(pred[i]|up[i]|1<<i) for i in range(n)]
    return up,inc

def chain(mask,pred):
    return all(pred[a]>>b&1 or pred[b]>>a&1 for a,b in combinations([i for i in range(len(pred)) if mask>>i&1],2))

def vg(a,b,pred):
    up,inc=structure(pred)
    return bool(inc[a]>>b&1) and ((pred[a]==pred[b] and chain(up[a]&~up[b],pred) and chain(up[b]&~up[a],pred)) or (up[a]==up[b] and chain(pred[a]&~pred[b],pred) and chain(pred[b]&~pred[a],pred)))

def count(pred):
    full=(1<<len(pred))-1
    @lru_cache(None)
    def visit(done):
        if done==full:return 1
        return sum(visit(done|1<<i) for i,p in enumerate(pred) if not(done>>i&1) and p&done==p)
    return visit(0)

def main():
    # A defect is absent, lower label 0..2, or upper label 0..2.
    # Sorting loses only a permutation of S, irrelevant to all checked claims.
    records=[]
    for m in range(3,9):
        representatives=labeled=0
        for cfg in combinations_with_replacement(range(7),m):
            D=[7 if z==0 or z>=4 else 7^(1<<(z-1)) for z in cfg]
            U=[7 if z<=3 else 7^(1<<(z-4)) for z in cfg]
            cliques=([i for i in range(m) if D[i]==7],[i for i in range(m) if U[i]==7])
            chosen=max(cliques,key=len)
            require(len(chosen)>=(m+1)//2,'clique size')
            for i,j in combinations(chosen,2):
                require((D[i]==D[j] and (U[i]&~U[j]).bit_count()<=1 and (U[j]&~U[i]).bit_count()<=1) or (U[i]==U[j] and (D[i]&~D[j]).bit_count()<=1 and (D[j]&~D[i]).bit_count()<=1),'very good defect pair')
            mult=factorial(m)
            for z in range(7):mult//=factorial(cfg.count(z))
            labeled+=mult;representatives+=1
        require(labeled==7**m,'permutation orbit completeness')
        records.append({'antichain_size':m,'orbit_representatives':representatives,'labeled_defect_assignments':labeled})
    # 14-vertex (6,6,0) root core: x=0,y=1,L=2..7,T=8..13.
    L=list(range(2,8));T=list(range(8,14))
    edges=[(a,b) for a in L for b in T]+[(a,1) for a in L]+[(0,b) for b in T]
    q=closure(14,edges)
    # Exterior a=14 lies below both roots and all but L[0].
    ext=edges+[(14,v) for v in [0,1]+T+L[1:]]
    p=closure(15,ext)
    up,inc=structure(p)
    require(max(x.bit_count() for x in inc)<=7,'ambient seven-thin')
    require((inc[0]|inc[1])==((1<<14)-1),'complete root core')
    require(pred_subset(q,p),'induced core unchanged')
    require(p[0]&~p[1]==0 and up[1]&~up[0]==0,'root critical')
    eq=count(q); ep=count(p)
    nq=count(closure(14,edges+[(L[0],L[1])]))
    np=count(closure(15,ext+[(L[0],L[1])]))
    require(2*nq==eq,'core twin half')
    require(2*np!=ep,'ambient weight differs')
    require(all(vg(a,b,p) for a,b in combinations(L[1:],2)),'ambient clique fixture')
    # An arbitrary-height extension: prepend a 5-chain and append a 4-chain.
    n=24
    high=ext+[(i,i+1) for i in range(15,19)]+[(19,v) for v in range(15)]+[(v,20) for v in range(15)]+[(i,i+1) for i in range(20,23)]
    h=closure(n,high)
    _,hi=structure(h)
    require(max(x.bit_count() for x in hi)<=7,'higher-height thinness')
    require(all(vg(a,b,h) for a,b in combinations(L[1:],2)),'higher-height ambient clique')
    print(json.dumps({'status':'ALL_CONTROLS_PASSED','guards':checks,'defect_orbits':records,'nonuniform_exterior_example':{'core_extensions':eq,'core_pair_count':nq,'core_probability':str(Fraction(nq,eq)),'ambient_extensions':ep,'ambient_pair_count':np,'ambient_probability':str(Fraction(np,ep))},'higher_height_vertices':n},indent=2))

def pred_subset(q,p):
    return all(p[i]&((1<<len(q))-1)==q[i] for i in range(len(q)))

if __name__=='__main__':main()
