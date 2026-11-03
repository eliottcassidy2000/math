#!/usr/bin/env python3
"""Independent exact audit: rational rank, prime guards, exponent-state simulator."""
from fractions import Fraction as F
from itertools import product
from math import gcd,isqrt,prod
import json
N=0
def need(t,msg):
 global N
 N+=1
 if not t: raise RuntimeError(msg)
def prime(n): return n>1 and all(n%d for d in range(2,isqrt(n)+1))
def rank(a):
 a=[[F(v) for v in r] for r in a]; k=0
 for j in range(len(a[0])):
  pivot=next((i for i in range(k,len(a)) if a[i][j]),None)
  if pivot is None: continue
  a[k],a[pivot]=a[pivot],a[k]; v=a[k][j]; a[k]=[x/v for x in a[k]]
  for i in range(k+1,len(a)):
   if a[i][j]:
    v=a[i][j]; a[i]=[x-v*y for x,y in zip(a[i],a[k])]
  k+=1
  if k==len(a): break
 return k
ranks=[]
for m in (3,5,7,11,15,21,35,145):
 s=[tuple(int(gcd(h*h-d*d,m)==1) for d in range(m)) for h in range(m)]
 unique=list(dict.fromkeys(s)); factors=[p for p in range(3,m+1,2) if prime(p) and m%p==0]
 expected=prod((p+1)//2 for p in factors)
 need(len(unique)==expected,'distinct row count')
 need(rank(unique)==expected,'rational row rank')
 for h in range(m):
  need(sum(s[h])==prod(p-1 if h%p==0 else p-2 for p in factors),'weighted row size')
  for d in range(m): need(s[h][d]==int(gcd(h-d,m)==gcd(h+d,m)==1),'port mixing identity')
 ranks.append([m,expected])
# Weighted local spectral decomposition is checked on explicit full-residue vectors.
for p in (3,5,7,11,29):
 s=[[int((h*h-d*d)%p!=0) for d in range(p)] for h in range(p)]
 def act(v): return [sum(a*b for a,b in zip(r,v)) for r in s]
 for j in range(1,(p+1)//2):
  v=[0]*p; v[j]=1; v[-j]=-1
  need(act(v)==[0]*p,'odd zero eigenvector')
 for j in range(2,(p+1)//2):
  v=[0]*p; v[j]=v[-j]=1; v[1]=v[-1]=-1
  need(act(v)==[-2*x for x in v],'even zero-sum eigenvector')
 e0=[1]+[0]*(p-1); e1=[0]+[1]*(p-1)
 need(act(e0)==e1,'remaining block column0')
 need(act(e1)==[p-1]+[p-3]*(p-1),'remaining block column1')
 need((p-1)//2+(p-3)//2+2==p,'basis dimensions')
# Complete odd-pair triangle through sum600, with the height/self-hit sieve.
row=[]; unitrow=[]
for h in range(3,301):
 for d in range(h-2):
  if (h+d)%2!=1: continue
  p,q=h-d,h+d
  need(p>=3 and p<=q and (p+q)//2==h and (q-p)//2==d,'triangle inverse')
  tested=[ell for ell in range(3,isqrt(2*h-3)+1,2) if prime(ell)]
  accepted=all((p%ell!=0 or p==ell) and (q%ell!=0 or q==ell) for ell in tested)
  need(accepted==(prime(p) and prime(q)),'full row sieve')
  if h==145 and accepted: row.append([p,q])
  if prime(p) and prime(q) and p<q:
   a,b,c=p*q,(q*q-p*p)//2,(q*q+p*p)//2
   need(a*a+b*b==c*c and gcd(a,b)==1,'primitive PPT')
for p in range(1,146,2):
 q=290-p
 if gcd(p*q,145)==1: unitrow.append([p,q])
need(len(row)==10 and len(unitrow)==56 and unitrow[0]==[1,289],'row290 cardinalities')
# Conway list, simulated independently as guarded exponent translations.
primes=(2,3,5,7,11,13,17,19,23,29)
program=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
def exponents(n):
 v=[]
 for p in primes:
  k=0
  while n%p==0: n//=p; k+=1
  v.append(k)
 need(n==1,'typed source support')
 return tuple(v)
rules=[(exponents(a),exponents(b)) for a,b in program]
def number(v): return prod(p**k for p,k in zip(primes,v))
def advance(v):
 for i,(a,b) in enumerate(rules):
  if all(x>=y for x,y in zip(v,b)):
   return i,tuple(x+y-z for x,y,z in zip(v,a,b))
 raise RuntimeError('missing last rule')
# Exhaust all zero/one exponent vectors, including unreachable control combinations.
for v in product((0,1),repeat=10):
 n=number(v); i,w=advance(v)
 literal=next(j for j,(a,b) in enumerate(program) if n%b==0)
 need(i==literal and number(w)==n//program[i][1]*program[i][0],'priority and integer realization')
 a,u,ctrl=v[0],prod(p**k for p,k in zip(primes[1:4],v[1:4])),prod(p**k for p,k in zip(primes[4:],v[4:]))
 h,d=(u+ctrl)//2,(ctrl-u)//2
 av,bv=rules[i]; rho=F(prod(p**k for p,k in zip(primes[1:4],av[1:4])),prod(p**k for p,k in zip(primes[1:4],bv[1:4])))
 sigma=F(prod(p**k for p,k in zip(primes[4:],av[4:])),prod(p**k for p,k in zip(primes[4:],bv[4:])))
 hh=((rho+sigma)*h+(sigma-rho)*d)/2; dd=((sigma-rho)*h+(rho+sigma)*d)/2
 need(hh.denominator==dd.denominator==1,'chart integral on legal cell')
 need(2**w[0]*(hh*hh-dd*dd)==number(w),'chart literal target')
trialcount=0
for n in range(2,17):
 for d in range(1,n):
  p,q=2*d+1,2*(n-d)+1
  need((n%d==0)==((q-1)%(p-1)==0),'ordered divisor ray')
  v=exponents(5**n*7**d*13); divisible=n%d==0
  target=2**n*7**(d-1) if divisible else 5**n*7**(d-1)*13
  expected=4*n+2*(n//d)+2 if divisible else 6*n+2*(n//d)+2
  for step in range(expected):
   _,v=advance(v)
   if step+1<expected: need(number(v)!=target,'first macro arrival')
  need(number(v)==target,'macro endpoint and cost'); trialcount+=1
need(advance(exponents(39))[0]==9 and number(advance(exponents(33))[1])==29,'priority defeats reciprocal word')
need(number(advance(exponents(13))[1])==11 and number(advance(exponents(11))[1])==13,'legal reciprocal control')
print(json.dumps({'status':'PASS','gates':N,'ranks':ranks,'row290':row,'unit_pair_count':len(unitrow),
 'universes':{'odd_pair_sums_through':600,'exponent_cube_states':1024,'complete_division_trials':trialcount,
              'local_weighted_spectra':[3,5,7,11,29]}},indent=2,sort_keys=True))
