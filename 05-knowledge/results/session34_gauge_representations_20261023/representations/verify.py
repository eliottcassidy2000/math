#!/usr/bin/env python3
"""Exact class-character computation; independent labelled-mask audit n<=4."""
from itertools import permutations,combinations
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from math import factorial
import json
@lru_cache(None)
def parts(n,cap=None):
    if not n:return ((),)
    if cap is None:cap=n
    return tuple((a,)+b for a in range(min(cap,n),1-1,-1) for b in parts(n-a,a))
def cells(p):return {(i,j) for i,a in enumerate(p) for j in range(a)}
@lru_cache(None)
def chi(lam,mu):
    if not mu:return int(not lam)
    k=mu[0]; out=0; L=cells(lam)
    for nu in parts(sum(lam)-k):
        N=cells(nu)
        if not N<=L:continue
        D=L-N
        if any({(i,j),(i+1,j),(i,j+1),(i+1,j+1)}<=D for i,j in D):continue
        seen={next(iter(D))}
        while True:
            new=seen|{(i+di,j+dj) for i,j in seen for di,dj in ((1,0),(-1,0),(0,1),(0,-1)) if (i+di,j+dj) in D}
            if new==seen:break
            seen=new
        if seen!=D:continue
        out+=(-1)**(len({i for i,j in D})-1)*chi(nu,mu[1:])
    return out
def perm(mu):
    p=[];start=0
    for a in mu:p+=list(range(start+1,start+a))+[start];start+=a
    return tuple(p)
def z(mu):
    o=1
    for a,m in Counter(mu).items():o*=a**m*factorial(m)
    return o
def fixed(n,p,converse):
    pairs=list(combinations(range(n),2)); remaining=set(pairs); cycles=0
    while remaining:
        start=next(iter(remaining));cur=start;flip=0
        while True:
            remaining.remove(cur);i,j=map(p.__getitem__,cur)
            flip^=(i>j)^converse;cur=tuple(sorted((i,j)))
            if cur==start:break
        if flip:return 0
        cycles+=1
    return 2**cycles
def brute(n,p,c):
    edges=list(combinations(range(n),2));ix={e:i for i,e in enumerate(edges)};total=0
    for mask in range(1<<len(edges)):
        good=True
        for k,(i,j) in enumerate(edges):
            a,b=p[i],p[j];l=ix[tuple(sorted((a,b)))];
            if (((mask>>k)^(mask>>l))&1)!=((a>b)^c):good=False;break
        total+=good
    return total
def transpose(lam):return tuple(sum(a>j for a in lam) for j in range(lam[0]))
def marked_brute(n):
    # Direct normal-form masks: C reverses order and reverses all arcs,
    # so it reflects the chord choices without complementing their bits.
    chords=[e for e in combinations(range(n),2) if e[1]-e[0]>1];ix={e:i for i,e in enumerate(chords)}
    return sum(all(((mask>>i)&1)==((mask>>ix[(n-1-b,n-1-a)])&1) for i,(a,b) in enumerate(chords)) for mask in range(1<<len(chords)))
rows=[]
for n in range(1,9):
    pp=parts(n);fixeds={mu:(fixed(n,perm(mu),False),fixed(n,perm(mu),True)) for mu in pp}
    for la in pp:
        for lb in pp:
            assert sum(Fraction(chi(la,mu)*chi(lb,mu),z(mu)) for mu in pp)==int(la==lb)
    if n<=4:
        for p in permutations(range(n)):
            for c in (False,True):assert fixed(n,p,c)==brute(n,p,c)
    mult={}
    for la in pp:
        a=sum(Fraction(chi(la,mu)*fixeds[mu][0],z(mu)) for mu in pp)
        b=sum(Fraction(chi(la,mu)*fixeds[mu][1],z(mu)) for mu in pp)
        assert a.denominator==b.denominator==1
        plus=(a+b)/2;minus=(a-b)/2
        assert plus.denominator==minus.denominator==1 and min(plus,minus)>=0
        mult[la]=(int(a),int(b),int(plus),int(minus))
    sign=(-1)**(n//2)
    for la in pp:
        other=mult[transpose(la)];a,b,_,_=mult[la]
        assert other[:2]==(a,sign*b)
    assert sum(chi(la,(1,)*n)*mult[la][0] for la in pp)==2**(n*(n-1)//2)
    m=(n-1)*(n-2)//2;f=(n-1)//2;q=(m+f)//2;w=(2,)*(n//2)+(1,)*(n%2)
    if n<=6:assert marked_brute(n)==2**q
    marked={}
    for la in pp:
        d=chi(la,(1,)*n);a=2**m*d;b=2**q*chi(la,w)
        assert (a+b)%2==0 and (a-b)%2==0
        marked[str(la)]={'dimension_multiplicity':a,'converse_trace':b,'plus':(a+b)//2,'minus':(a-b)//2}
    rows.append({'n':n,'unlabelled':mult[(n,)][0],'self_converse':mult[(n,)][1],'tournament_irreps':{str(k):{'multiplicity':v[0],'converse_trace':v[1],'plus':v[2],'minus':v[3]} for k,v in mult.items()},'marked_chords':m,'marked_fixed_chord_orbits':q,'marked_irreps':marked})
print(json.dumps(rows,indent=2))
