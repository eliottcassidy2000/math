#!/usr/bin/env python3
"""Exact structural controls for the arithmetic half-turn tournament lift."""
from fractions import Fraction as Q
from itertools import combinations, product, permutations
from math import gcd
import json

def sign(k,t):
    f=(k*t)%1
    if f==0 or f==Q(1,2): raise ValueError('wall has no tournament orientation')
    return 1 if f<Q(1,2) else -1

def walls(freq):
    return sorted({Q(j,2*d) for d in freq for j in range(2*d)}|{Q(1)})

def correlation(a,b):
    d=gcd(a,b); p,q=a//d,b//d
    return Q(1,p*q) if p%2 and q%2 else Q(0)

def tri_formula(a,b):
    return {(x,y,z):Q(1,8)*(1+x*y*correlation(a,b)+x*z*correlation(a,a+b)+y*z*correlation(b,a+b)) for x,y,z in product((-1,1),repeat=3)}

def tri_direct(a,b):
    ans={s:Q(0) for s in product((-1,1),repeat=3)}
    w=walls((a,b,a+b))
    for l,r in zip(w,w[1:]):
        t=(l+r)/2
        ans[tuple(sign(d,t) for d in (a,b,a+b))]+=r-l
    return ans

def word(V,K,t):
    return tuple(sign((v-u)*2**k,t) for k in range(K) for u,v in combinations(V,2))

def h_from_word(label,n):
    arcs={}
    for (i,j),sgn in zip(combinations(range(n),2),label):
        arcs[i,j]=sgn==1; arcs[j,i]=sgn==-1
    return sum(all(arcs[p[k],p[k+1]] for k in range(n-1)) for p in permutations(range(n)))

def check():
    triples=0; lo=[]; hi=[]
    for a in range(1,31):
        for b in range(1,31):
            exact=tri_direct(a,b); formula=tri_formula(a,b)
            if exact!=formula: raise RuntimeError((a,b,exact,formula))
            cyc=sum(v for (x,y,z),v in exact.items() if x==y and y==-z)
            if not Q(1,6)<=cyc<=Q(1,2): raise RuntimeError('cycle range')
            if cyc==Q(1,6): lo.append((a,b))
            if cyc==Q(1,2): hi.append((a,b))
            if (cyc==Q(1,6)) != (a==2*b or b==2*a): raise RuntimeError('lower equality')
            if (cyc==Q(1,2)) != (a==b): raise RuntimeError('upper equality')
            triples+=1
    languages=[]
    for V in ((0,1,2),(0,1,3,7),(0,2,5,11),(0,6,10)):
        ds=tuple(v-u for u,v in combinations(V,2))
        base=len({Q(j,d) for d in ds for j in range(d)})
        partition=walls(ds)
        cells=list(zip(partition,partition[1:]))
        labels=[word(V,1,(l+r)/2) for l,r in cells]
        lengths=[r-l for l,r in cells]
        adjacency=[]
        for l,r in cells:
            row=[]
            for c,e in cells:
                y=(c+e)/2
                row.append(any(l<x<r for x in (y/2,(y+1)/2)))
            adjacency.append(row)
        for i,row in enumerate(adjacency):
            if sum(lengths[j] for j,e in enumerate(row) if e)!=2*lengths[i]: raise RuntimeError('Markov right eigenvector')
        for j in range(len(cells)):
            if sum(row[j] for row in adjacency)!=2: raise RuntimeError('two preimages')
        state_words=[{label} for label in labels]
        for K in range(1,7):
            if K>1:
                state_words=[{prefix+labels[j] for i,row in enumerate(adjacency) if row[j] for prefix in state_words[i]} for j in range(len(cells))]
            w=walls(tuple(d*2**(K-1) for d in ds))
            observed={word(V,K,(l+r)/2) for l,r in zip(w,w[1:])}
            if set().union(*state_words)!=observed: raise RuntimeError('finite automaton language')
            if len(w)-1!=base*2**K: raise RuntimeError('cell scaling')
            if not 2**K<=len(observed)<=base*2**K: raise RuntimeError('language sandwich')
            languages.append({'speeds':V,'K':K,'cells':len(w)-1,'words':len(observed),'base_cells':base})
    margins=0
    for V in ((1,2,3),(1,3,7),(2,5,11)):
        for K in range(1,9):
            D=2**K
            for j in range(1,101):
                t=Q(j,211)
                m=[]
                for v in V:
                    bits=[int(sign(v*2**k,t)<0) for k in range(K)]
                    digit=sum(b*2**(K-1-k) for k,b in enumerate(bits))
                    if digit!=int(D*((v*t)%1)): raise RuntimeError('phase decoder')
                    m.append(digit)
                low=min(Q(min(x,D-x-1),D) for x in m)
                upp=min(Q(max(min(x,D-x),min(x+1,D-x-1)),D) for x in m)
                F=min(min((v*t)%1,1-(v*t)%1) for v in V)
                if not low<=F<=upp or upp-low>Q(1,D): raise RuntimeError('margin enclosure')
                margins+=1
    hostile=[]
    V=(0,1,3,7)
    for l,r,child in ((Q(1,14),Q(1,12),1),(Q(1,12),Q(1,8),5)):
        ds=[v-u for u,v in combinations(V,2)]
        cuts=sorted({l,r}|{x for x in walls([2*d for d in ds]) if l<x<r})
        for x,y in zip(cuts,cuts[1:]):
            t=(x+y)/2
            if h_from_word(word(V,1,t),4)!=5 or h_from_word(word(V,1,2*t),4)!=child: raise RuntimeError('nonlumpability hostile')
        hostile.append({'cell':[str(l),str(r)],'H':5,'next_H':child})
    return {'status':'PASS','scalar_hostile':hostile,'triple_gap_pairs':triples,'minimum_cycle_probability':'1/6','maximum_cycle_probability':'1/2','minimum_equality_pairs':len(lo),'maximum_equality_pairs':len(hi),'language_controls':languages,'margin_controls':margins,'scope':'Exact finite controls; all-order proofs in PROOF.md. No LRC or Chan-Pak closure.'}
if __name__=='__main__': print(json.dumps(check(),indent=2))
