#!/usr/bin/env python3
"""Independent character engine: Jacobi-Trudi + permutation character cycle allocations.
Producer uses rim hooks. This checker imports no producer source.
"""
import itertools as it,math,json
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
@lru_cache(None)
def parts(n,m=None):
 if not n:return ((),)
 if m is None:m=n
 return tuple((i,)+p for i in range(min(n,m),0,-1) for p in parts(n-i,i))
def transpose(l):return tuple(sum(i>k for i in l) for k in range(l[0]))
@lru_cache(None)
def allocations(mu,a):
 if min(a,default=0)<0:return 0
 if not mu:return int(not any(a))
 return sum(allocations(mu[1:],a[:i]+(v-mu[0],)+a[i+1:]) for i,v in enumerate(a) if v>=mu[0])
def chi(l,mu):
 if len(l)>l[0]:return (-1)**(sum(mu)-len(mu))*chi(transpose(l),mu)
 L=len(l);s=0
 for p in it.permutations(range(L)):
  sign=(-1)**sum(p[i]>p[j] for i in range(L) for j in range(i+1,L))
  s+=sign*allocations(mu,tuple(l[i]-i+p[i] for i in range(L)))
 return s
def z(mu):
 from collections import Counter
 return math.prod(k**v*math.factorial(v) for k,v in Counter(mu).items())
def graded(mu):
 p=[]
 for v in mu:
  a=len(p);p+=list(range(a+1,a+v))+[a]
 E=set(it.combinations(range(len(p)),2));coeff=[1]
 while E:
  e=next(iter(E));d=e;L=0;s=1
  while True:
   E.remove(d);a,b=p[d[0]],p[d[1]];s*=1 if a<b else -1;d=tuple(sorted((a,b)));L+=1
   if d==e:break
  new=coeff+[0]*L
  for k,c in enumerate(coeff):new[k+L]+=s*c
  coeff=new
 return coeff
root=Path(__file__).resolve().parents[1];expected=json.loads((root/'representations/results.json').read_text());spect=json.loads((root/'representations/spectrum.json').read_text())
summary=[]
for n in range(1,9):
 pp=parts(n);G={mu:graded(mu) for mu in pp};characters={(la,mu):chi(la,mu) for la in pp for mu in pp};N=n*(n-1)//2
 for la in pp:
  v=[sum(Fraction(characters[la,mu]*G[mu][k],z(mu)) for mu in pp) for k in range(N+1)]
  assert all(x.denominator==1 and x>=0 for x in v)
  assert v==spect[n-1]['walsh_degree_multiplicity'][str(la)]
  a=sum(v);b=sum((-1)**k*x for k,x in enumerate(v));r=expected[n-1]['tournament_irreps'][str(la)]
  assert (a,b,(a+b)/2,(a-b)/2)==tuple(r[k] for k in ['multiplicity','converse_trace','plus','minus'])
  m=(n-1)*(n-2)//2;q=(m+(n-1)//2)//2;w=(2,)*(n//2)+(1,)*(n%2)
  a=2**m*characters[la,(1,)*n];b=2**q*characters[la,w];r=expected[n-1]['marked_irreps'][str(la)]
  assert (a,b,(a+b)//2,(a-b)//2)==tuple(r[k] for k in ['dimension_multiplicity','converse_trace','plus','minus'])
 summary.append(dict(n=n,partitions=len(pp),characters_checked=len(pp)**2,all_multiplicities_match=True))
Path(__file__).with_name('character_certificate.json').write_text(json.dumps(summary,indent=2)+'\n')
print('PASS independent Jacobi-Trudi character engine n<=8, every graded/converse/marked multiplicity')
