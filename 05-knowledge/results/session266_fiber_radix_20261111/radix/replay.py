#!/usr/bin/env python3
"""Small exact controls of the new lane's formulas, not a tournament census."""
import json
from itertools import product
from math import comb, factorial
from pathlib import Path

def require(ok, label):
    if not ok:
        raise RuntimeError(label)

def simulate(a,b,c,t):
    for _ in range(t):
        c += a+b
        b += 2*a
    return c

def compositions(m,r):
    if r==1:
        yield (m,)
    else:
        for j in range(m+1):
            for rest in compositions(m-j,r-1):
                yield (j,)+rest

def diffs(v):
    ans=[]
    while v:
        ans.append(v[0])
        v=[b-a for a,b in zip(v,v[1:])]
    return ans

def fixed_formula(m,k):
    if k%2:
        return comb(m//2+(k-1)//2,(k-1)//2) if m%2==0 else 0
    r=k//2
    if r==0:
        return 1
    return sum(comb(j+r-1,r-1) for j in range(m//2+1))

checks=0
for a,b,c in product(range(10),repeat=3):
    for t in range(16):
        require(simulate(a,b,c,t)==a*t*t+b*t+c,'timing polynomial')
        checks+=1
    cs=[simulate(a,b,c,t) for t in (3,4,5)]
    aa=(cs[2]-2*cs[1]+cs[0])//2
    bb=cs[1]-cs[0]-aa*7
    cc=cs[0]-aa*9-bb*3
    require((aa,bb,cc)==(a,b,c),'three-observation inverse')
for B in range(2,15):
    for t in range(1,18):
        images={a*t*t+b*t+c for a,b,c in product(range(B),repeat=3)}
        require((len(images)==B**3)==(t>=B),'injectivity threshold')
        if t==B:
            require(images==set(range(B**3)),'radix contiguous range')
        checks+=1
for k in range(7):
    for m in range(12):
        gaps=list(compositions(m,k+1))
        require(len(gaps)==comb(m+k,k),'gap census')
        fixed=sum(g==g[::-1] for g in gaps)
        require(fixed==fixed_formula(m,k),'converse gap formula')
        orbits={min(g,g[::-1]) for g in gaps}
        require(len(orbits)==(len(gaps)+fixed)//2,'Burnside')
        vals=[comb(t+k,k) for t in range(k+2)]
        require(diffs(vals)[-1]==0,'nilpotent recurrence')
        checks+=1
# Catalog tomography with deliberately overlapping core sizes and lengths.
catalog=[(3,1,1),(4,1,2),(6,2,1),(7,2,3),(9,3,1)]
def count(n):
    return sum(w*comb(n-s+k,k) for s,k,w in catalog if n>=s)
K=max(k for s,k,w in catalog)
n0=max(s for s,k,w in catalog)
jet=diffs([count(n0+t) for t in range(K+1)])
for t in range(100):
    require(sum(d*comb(t,j) for j,d in enumerate(jet) if j<=t)==count(n0+t),'tomography')
    checks+=1
for B in range(2,16):
    values={a*B*B-b*B+c for a,b,c in product(range(B),repeat=3)}
    require(values==set(range(-B*B+B,B**3-B*B+B)), 'negative radix interval')
    checks+=1
for i,j,t in product(range(7),range(7),range(15)):
    rhs=sum(factorial(k)//(factorial(k-i)*factorial(k-j)*factorial(i+j-k))*comb(t,k)
            for k in range(max(i,j),i+j+1) if k<=t)
    require(comb(t,i)*comb(t,j)==rhs, 'positive binomial product')
    checks+=1
base_examples={str(B):[427//(B*B),(427//B)%B,427%B] for B in [8,9,10,11,12,13,16]}
receipt={'status':'PASS','scope':'plant recurrence and radix injectivity; finite gap enumeration; converse Burnside; synthetic finite-catalog tomography','checks':checks,'decimal_solution':[4,2,7],'integer_427_base_examples':base_examples,'catalog_jet':jet,'not_claimed':'No exhaustive tournament census, Collatz convergence, or new attribution of the original puzzle.'}
print(json.dumps(receipt,indent=2))
if '--save' in __import__('sys').argv:
    Path(__file__).with_name('receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
