#!/usr/bin/env python3
"""Independent regular-translation monomial boundaries for r,k<=4."""
from itertools import product
from math import comb
import json
checks=0

def need(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)

def compositions(n,r):
 if r==1:yield (n,);return
 for j in range(n+1):
  for tail in compositions(n-j,r-1):yield (j,)+tail

def valuation(x):
 if not x:return 10**6
 v=0
 while x%3==0:x//=3;v+=1
 return v

def axis_power(j):
 c=[1,0,0]
 for _ in range(j):c=[c[(a+1)%3]-c[a] for a in range(3)]
 return c

def kron_vector(v,w):return [a*b for a in v for b in w]
reports=[]
for r,k in product(range(1,5),repeat=2):
 n=3**r;m=3**k;cutoff=2*(k+r)-1;points=list(product(range(3),repeat=r));index={a:i for i,a in enumerate(points)}
 shifts=[]
 for axis in range(r):
  shifts.append([index[tuple((x+1)%3 if i==axis else x for i,x in enumerate(a))] for a in points])
 shell=[];max_j=cutoff+1;smallpowers=[axis_power(j) for j in range(max_j+1)]
 for degree in (cutoff-1,cutoff,cutoff+1):
  monomials=0;survivors=[]
  for js in compositions(degree,r):
   monomials+=1
   tensor=[1]
   for j in js:tensor=kron_vector(tensor,smallpowers[j])
   direct=[int(i==0) for i in range(n)]
   for axis,j in enumerate(js):
    for _ in range(j):direct=[direct[z]-direct[a] for a,z in enumerate(shifts[axis])]
   need(tensor==direct,'regular translation and tensor delta columns agree')
   v=min(valuation(c) for c in direct);pred=sum((j-1)//2 for j in js if j)
   need(v==pred,'exact integral monomial content valuation')
   reduced=[a%m for a in direct]
   survives=any(reduced)
   need(survives==(pred<k),'exact nonzero modulus criterion')
   if survives:
    survivors.append(list(js))
    need(degree==cutoff-1 and all(j>0 and j%2==0 for j in js),'only even full-support top monomials survive')
    need(reduced==[((-1)**(k-1)*3**(k-1))%m]*n,'top-layer norm coefficient vector')
  expected=comb(k+r-2,r-1) if degree==cutoff-1 else 0
  need(len(survivors)==expected,'complete sharp boundary survivor count')
  shell.append(dict(degree=degree,monomials=monomials,surviving=len(survivors),surviving_exponents=survivors))
 # Norm action on each delta is the same nonzero constant; translations cannot enlarge this line.
 norm=[3**(k-1)]*n
 need(all([norm[z] for z in shift]==norm for shift in shifts),'top norm translation-invariant')
 need(any(x%m for x in norm) and all(3*x%m==0 for x in norm),'top ideal additive order exactly three')
 reports.append(dict(r=r,k=k,modulus=m,regular_dimension=n,nilpotency_index=cutoff,boundaries=shell))
print(json.dumps(dict(status='PASS',checks=checks,scope='regular representation only; arbitrary modules may have smaller index',reports=reports),sort_keys=True,indent=2))
