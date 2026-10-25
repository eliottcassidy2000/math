#!/usr/bin/env python3
"""Independent orbit audit: permutation Burnside versus onto-block inclusion-exclusion.
Does not import check.py. Reads only saved consequences, never their internals.
"""
from itertools import permutations
from math import comb, factorial, prod
from collections import Counter
from fractions import Fraction
from pathlib import Path
import json

checks=0
def need(p,m):
    global checks
    checks+=1
    if not p:raise RuntimeError(m)

def fixed_hist(r):
    return Counter(sum(p[i]==i for i in range(r)) for p in permutations(range(r)))

def burnside(q,k,w,h):
    num=sum(count*(k+f)**w for f,count in h.items())
    need(num%factorial(q-k)==0,'Burnside divisibility')
    return num//factorial(q-k)

def onto(q,k,w):
    total=0
    for b in range(min(q-k,w)+1):
        num=sum((-1)**(b-j)*comb(b,j)*(k+j)**w for j in range(b+1))
        need(num%factorial(b)==0,'onto-block divisibility')
        total+=num//factorial(b)
    return total

cases=0
for q in range(2,8):
    for k in range(q+1):
        h=fixed_hist(q-k)
        for w in range(9):
            need(burnside(q,k,w,h)==onto(q,k,w),'independent orbit formula')
            path=sum(c*(k+f)**w*(2+k+f) for f,c in h.items())//factorial(q-k)
            cycle=sum(c*(k+f)**(w+1)*(2+k+f) for f,c in h.items())//factorial(q-k)
            need(path==2*onto(q,k,w)+onto(q,k,w+1),'independent path budget')
            need(cycle==2*onto(q,k,w+1)+onto(q,k,w+2),'independent anchor budget')
            cases+=1

producer=json.loads(Path(__file__).with_name('normal.json').read_text())
h=fixed_hist(6)
need(producer['q6_w6_budgets']['path_orbits']==
     sum(c*f**6*(2+f) for f,c in h.items())//factorial(6),'published path number')
need(producer['q6_w6_budgets']['cycle_orbits']==
     sum(c*f**7*(2+f) for f,c in h.items())//factorial(6),'published cycle number')

# Leibniz determinants independently verify the sharp q+1 lower bound.
determinants=[]
for q in range(2,8):
    def entry(i,j):
        return 1 if i==0 or j==0 else int(i!=j)
    determinant=sum((-1)**sum(p[i]>p[j] for i in range(q+1) for j in range(i+1,q+1)) *
                    prod(entry(i,p[i]) for i in range(q+1))
                    for p in permutations(range(q+1)))
    need(determinant==(-1)**q,'independent Leibniz Hankel determinant')
    determinants.append([q,determinant])

# Explicit inverse of the matrix of singleton continuation observations.
for q in range(2,11):
    A=[[int(i!=j) for j in range(q)] for i in range(q)]
    B=[[Fraction(1,q-1)-int(i==j) for j in range(q)] for i in range(q)]
    for i in range(q):
        for j in range(q):
            need(sum(A[i][t]*B[t][j] for t in range(q))==int(i==j), 'explicit singleton response inverse')

print(json.dumps({'status':'PASS','checks':checks,'independent_orbit_cases':cases,
                 'methods':['permutation Burnside','onto-block inclusion-exclusion',
                            'Leibniz determinants','explicit rational inverse'],
                 'whole_word_determinants':determinants},indent=2))
