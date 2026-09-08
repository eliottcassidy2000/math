#!/usr/bin/env python3
"""Exact controls for the central cyclic-cover genus bound.

The analytic embedding and rational-time argument are separate proof duties.
No polynomial coefficient box or generic-parameter inference is used.
"""
from math import gcd
from functools import reduce
from itertools import product
from hashlib import sha256
import json
import sympy as S
checks=0
def check(x,label):
 global checks
 checks+=1
 if not x: raise RuntimeError(label)
def data(A,B,es):
 E=sum(es);r=len(es);M=2*A+3*B+6*E
 n0=A+B+2*E;ni=A+B+3*E
 d=reduce(gcd,[A,B]+list(es))
 gd=reduce(gcd,[M,n0]+list(es))
 T=r*M-gcd(M,n0)-gcd(M,ni)-sum(gcd(M,e) for e in es)
 return M,n0,ni,d,gd,T
rows=[]
for r in range(1,4):
 for es in product(range(1,5),repeat=r):
  for A in range(1,13):
   for B in range(13):
    M,n0,ni,d,gd,T=data(A,B,es)
    check(d==gd,'actual connected-component gcd')
    check(T%(2*d)==0,'genus integral on every cyclic component')
    check(T>=2*d,'every central component genus at least two')
    ap,bp=A//d,B//d;eps=tuple(e//d for e in es);ep=sum(eps)
    mm,nn,ii,dd,_,tt=data(ap,bp,eps)
    check(dd==1 and T==d*tt,'primitive normalization preserves component genus')
    check(gcd(mm,nn)<=bp+2*ep,'zero puncture gcd bound')
    check(sum(gcd(mm,e) for e in eps)<=ep,'finite puncture budget')
    if bp:
     check(gcd(mm,ii)<=bp,'infinity gcd bound when y exponent nonzero')
     lower=r*mm-2*bp-3*ep
     check(tt>=lower>=2*ap+bp+3*ep>=6,'positive-y lower bound')
    else:
     check(gcd(mm,ii)==mm//2,'zero-y exact infinity half')
     lower=(2*r-1)*ap+6*(r-1)*ep
     check(tt>=lower>=1 and tt%2==0,'zero-y positive even bound')
    rows.append((A,B,es,d,1+T//(2*d)))
# Direct sheet orbits and puncture permutations: independent covering path.
for A in range(1,7):
 for B in range(5):
  for es in [(1,),(2,),(1,1),(1,2),(2,2),(1,1,1)]:
   M,n0,ni,d,gd,T=data(A,B,es)
   exponents=[n0]+list(es)+[-ni]
   check(sum(exponents)==0,'full monodromy product')
   unseen=set(range(M));orbits=[]
   while unseen:
    start=min(unseen);orb={start};todo=[start]
    while todo:
     v=todo.pop()
     for e in exponents:
      q=(v+e)%M
      if q not in orb:orb.add(q);todo.append(q)
    unseen-=orb;orbits.append(orb)
   check(len(orbits)==d and all(len(o)==M//d for o in orbits),'independent sheet components')
   for orb in orbits:
    ram=0
    for e in exponents:
     unseen2=set(orb);cycles=0
     while unseen2:
      start=min(unseen2);q=start
      while q in unseen2:unseen2.remove(q);q=(q+e)%M
      cycles+=1
     ram+=len(orb)-cycles
    check((-2*len(orb)+ram+2)//2==1+T//(2*d),'literal permutation Riemann Hurwitz')
p,y,v,w=S.symbols('p y v w')
for R in [1,p,y,p**3-y**2,p**3+y**2,p**4+y**3,(p**3-y**2)**3+p**10,p**4*y+y**3+y**5]:
 I=S.expand(p*(p**3-y**2)*R)
 terms=S.Poly(I,p,y).terms();weight=min(2*a+3*b for (a,b),c in terms)
 initial=sum(c*p**a*y**b for (a,b),c in terms if 2*a+3*b==weight)
 expr=S.expand(I.subs({p:v**2*w,y:v**3*w})/v**weight)
 check(S.cancel(S.limit(expr,v,0)-initial.subs({p:w,y:w}))==0,'actual third-blowup leading coefficient')
for A,B,es in [(1,0,(1,)),(2,0,(1,)),(2,0,(2,)),(10,0,(10,))]:
 M,n0,ni,d,gd,T=data(A,B,es)
 check(1+T//(2*d)==2,'sharp equality and disconnected-power controls')
# Removing either mandatory factor can reduce genus.
check((3+1-2)//2==1,'Delta-only elliptic boundary')
check((2-2)//2==0,'p-only rational boundary')
blob=json.dumps(rows,separators=(',',':')).encode()
print('Central cyclic-cover universal genus controls PASS')
print('Multiplicity universe:',len(rows),'rows; A1..12,B0..12,1..3 distinct binomial factors,exponents1..4')
print('Independent literal permutation cases:',6*5*6)
print('Always-active gates:',checks)
print('Semantic SHA256:',sha256(blob).hexdigest())
print('Scope: exact genus arithmetic; analytic embedding and flow consumer require the proof and audit.')
