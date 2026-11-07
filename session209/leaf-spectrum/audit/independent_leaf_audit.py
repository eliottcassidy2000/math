#!/usr/bin/env python3
"""Independent exact controls for the leaf-spectrum merger. Stdlib only."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from math import comb, factorial
import json

def need(p, msg):
    if not p: raise RuntimeError(msg)

@lru_cache(None)
def words(n):
    if n == 0: return ((),)
    if n < 0: return ()
    return tuple((a,) + w for a in (1,3) for w in words(n-a))

def matrix(word):
    n=sum(word); out=[[False]*n for _ in range(n)]; start=0
    for a in word:
        for u in range(start,start+a):
            for v in range(start+a,n): out[u][v]=True
        if a==3:
            for u in range(3): out[start+u][start+(u+1)%3]=True
        if a==4:
            for u,v in ((0,1),(1,2),(2,3),(3,0),(2,0),(3,1)):
                out[start+u][start+v]=True
        start+=a
    return out

def raw_counts(out):
    n=len(out); hp=aut=0
    for p in permutations(range(n)):
        hp += all(out[p[i]][p[i+1]] for i in range(n-1))
        aut += all(out[u][v]==out[p[u]][p[v]] for u in range(n) for v in range(u+1,n))
    return hp,aut

def closed(n,k):
    s=n-3*k
    if s<0:return 0
    pal=0 if s%2 and k%2 else comb(s//2+k//2,k//2)
    return (comb(s+k,k)-pal)//2

counts=[]; tested=0; matrix_cases=0
for n in range(31):
    # Literal words and literal reversal comparison, independent of GF.
    pure=Counter(w.count(3) for w in words(n) if w<w[::-1])
    mixed=Counter()
    if n>=4 and n%2==0:
        mixed=Counter(2*w.count(3) for w in words((n-4)//2))
    for k in range(n//3+1):
        need(pure[k]==closed(n,k),('closed coefficient',n,k))
        pred=0
        if k%2==0 and n>=4+3*k and (n-4-3*k)%2==0:
            s=(n-4-3*k)//2;pred=comb(s+k//2,k//2)
        need(mixed[k]==pred,('mixed coefficient',n,k))
        tested+=1
    mass=sum((Fraction(2*v,3**k) for k,v in pure.items()),Fraction())
    mass+=sum((Fraction(v,3**k) for k,v in mixed.items()),Fraction())
    labelled=mass*factorial(n)
    need(labelled.denominator==1,('labelled integer',n))
    if n<=8:
        for w in words(n):
            if not w<w[::-1]:continue
            hp,aut=raw_counts(matrix(w));k=w.count(3)
            need((hp,aut)==(3**k,3**k),('literal pure',w,hp,aut))
            matrix_cases+=1
        if n>=4 and n%2==0:
            for w in words((n-4)//2):
                ww=w+(4,)+w[::-1];hp,aut=raw_counts(matrix(ww));k=2*w.count(3)
                need((hp,aut)==(5*3**k,3**k),('literal mixed',ww,hp,aut))
                matrix_cases+=1
    counts.append({'n':n,'pure':dict(sorted(pure.items())),
                   'mixed':dict(sorted(mixed.items())),
                   'leaves':sum(pure.values())+sum(mixed.values()),
                   'labelled_lifts':int(labelled)})

# Univariate rational denominator checked by independent literal coefficients.
ls=[row['leaves'] for row in counts]
for n in range(31):
    coeff=ls[n]
    for shift,sign in ((1,-1),(2,-1),(5,1),(6,-1),(7,1),(9,1)):
        if n>=shift:coeff+=sign*ls[n-shift]
    need(coeff=={4:2,5:-1,7:-1}.get(n,0),('rational numerator',n,coeff))

print(json.dumps({'status':'PASS','literal_word_max_order':30,
                  'closed_formula_cells':tested,
                  'literal_H_and_Aut_max_order':8,
                  'literal_H_and_Aut_cases':matrix_cases,
                  'rows':counts},indent=2))
