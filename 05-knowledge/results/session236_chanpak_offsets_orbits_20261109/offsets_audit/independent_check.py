#!/usr/bin/env python3
"""Independent audit: no imports from the submitted verifier or repository.
Hamilton counts by first-vertex deletion recursion; rooted response by all
partial permutations; periodic drift by cyclic integer rotations.
"""
from functools import lru_cache
from itertools import combinations, permutations, product
from fractions import Fraction as F
import json
from pathlib import Path

def tour(n,z):
    arcs=set()
    for i,(a,b) in enumerate(combinations(range(n),2)):
        arcs.add((a,b) if z&(1<<i) else (b,a))
    return n,frozenset(arcs)

def counter(t):
    n,arcs=t
    @lru_cache(None)
    def start(v,S):
        if not S:return 1
        return sum(start(w,S-{w}) for w in S if (v,w) in arcs)
    def h(S=None):
        if S is None:S=frozenset(range(n))
        if not S:return 1
        return sum(start(v,S-{v}) for v in S)
    return h

def response(t,u):
    n,arcs=t;h=counter(t);ans=[0,0]
    for size in range(n):
        for tail in permutations([v for v in range(n) if v!=u],size):
            p=(u,)+tail
            if all((a,b) in arcs for a,b in zip(p,p[1:])):
                ans[len(p)%2]+=h(frozenset(range(n))-set(p))
    return h(),*ans

def padding(t,m):
    n,e=t
    return n+m,e|frozenset((i,j) for i in range(n+m) for j in range(max(n,i+1),n+m))

def joined(t,u,s,w):
    n,e=t;k,f=s
    arcs=set(e)|{(a+n,b+n) for a,b in f}|{(a,b+n) for a in range(n) for b in range(k)}
    arcs.remove((u,w+n));arcs.add((w+n,u))
    return n+k,frozenset(arcs)

def oneapex(bits):
    k=len(bits);arcs={(i,j) for i in range(k) for j in range(i+1,k)}
    arcs|={(k,i) if b else (i,k) for i,b in enumerate(bits)}
    return k+1,frozenset(arcs)

def wordcount(bits):
    # Count odd cycles through apex: choose an increasing odd-cardinality
    # chain subsequence, whose endpoints carry 1 and 0 respectively.
    # An odd cycle needs an even chain subset, thus an even size >=2.
    s=0
    for size in range(2,len(bits)+1,2):
        for ix in combinations(range(len(bits)),size):
            if bits[ix[0]]==1 and bits[ix[-1]]==0:s+=1
    return 1+2*s

def blockincrement(w,bits):
    change=0
    for b in bits:
        if b==0:change+=w+(w&1)
        w=2*w+b
    return change

def main():
    roots=[(tour(n,z),u) for n in range(1,4) for z in range(1<<(n*(n-1)//2)) for u in range(n)]
    padchecks=seamchecks=0
    for t,u in roots:
        h,e,o=response(t,u)
        for m in range(1,5):
            assert response(padding(t,m),u)==(h,2**(m-1)*(e+o),2**(m-1)*(e+o))
            padchecks+=1
        for s,w in roots:
            k,f=s
            hb,qe,qo=response((k,frozenset((b,a) for a,b in f)),w)
            for m in (1,2):
                assert counter(joined(padding(t,m),u,s,w))()==h*hb+2**m*(e+o)*(qe+qo)
                seamchecks+=1
    direct=0
    for k in range(0,9):
        for bits in product((0,1),repeat=k):
            assert counter(oneapex(bits))()==wordcount(bits)
            direct+=1
    # Determine affine coefficients from values at two admissible starting
    # integers, rather than the submitted coefficient-accumulation algorithm.
    checks=0; minima={}; rotationchecks=0
    for k in range(1,17):
        q=1<<k;den=q-1;least=None
        for d in range(q):
            bits=tuple(int(c) for c in f'{d:0{k}b}')
            w=bits[-1];base=blockincrement(w,bits)
            alpha=(blockincrement(w+2,bits)-base)//2
            beta=base-alpha*w
            numerator=0;rot=d
            for i,b in enumerate(bits):
                if b==0:numerator+=den*bits[i-1]-rot
                rot=((rot<<1)&den)|(rot>>(k-1))
            assert rot==d
            assert numerator==den*beta-alpha*d
            if d not in (0,den):
                assert numerator>0
                least=numerator if least is None else min(least,numerator)
            else:assert numerator==0
            checks+=1
            # Cyclic drift must be cut-independent, including wrapping runs.
            if k<=10:
                rotbits=bits[1:]+bits[:1];v=rotbits[-1]
                c=blockincrement(v,rotbits);a=(blockincrement(v+2,rotbits)-c)//2
                rot_d=((d<<1)&den)|(d>>(k-1))
                assert den*(c-a*v)-a*rot_d==numerator
                rotationchecks+=1
        if least is not None:minima[k]=str(F(least,den))
    assert counter(joined(tour(1,0),0,tour(1,0),0))()==1
    assert response(tour(2,1),0)==(1,1,1)
    assert response(tour(2,1),1)==(1,0,1)
    output=dict(status='PASS mathematical claims; citation corrections documented in AUDIT.md',padding_checks=padchecks,bridge_checks=seamchecks,independent_apex_counts=direct,periodic_words=checks,max_period_length=16,cyclic_rotation_checks=rotationchecks,minimal_mixed_drifts=minima)
    print(json.dumps(output,indent=2))
if __name__=='__main__':main()
