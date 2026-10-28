#!/usr/bin/env python3
"""Exact C2-equivariant tournament/prime indices; no supplied prime census input."""
from fractions import Fraction as F
from collections import Counter
from math import factorial,gcd
import json,sys
N=int(sys.argv[1]) if len(sys.argv)>1 else 15

def parts(n,m=1):
 if not n: yield ()
 for d in range(m,n+1):
  for q in parts(n-d,d): yield (d,)+q

def add(a,b,s=1):
 c=a.copy()
 for p,v in b.items():
  c[p]=c.get(p,0)+s*v
  if not c[p]: del c[p]
 return c

def mul(a,b):
 c={}
 for p,v in a.items():
  for q,w in b.items():
   if sum(p)+sum(q)<=N:
    k=tuple(sorted(p+q));c[k]=c.get(k,0)+v*w
 return {p:v for p,v in c.items() if v}

def psi(a,d):return {tuple(d*j for j in p):v for p,v in a.items() if d*sum(p)<=N}
def comp(p,ae,at=None):
 cache={d:psi(at if at is not None and d%2 else ae,d) for d in range(1,N+1)}
 c={}
 for q,v in p.items():
  term={():v}
  for d in q:term=mul(term,cache[d])
  c=add(c,term)
 return c

def inv(a):
 f=add(a,{():1},-1);r={():F(1)};t=r.copy()
 for k in range(1,N+1):t=mul(t,f);r=add(r,t,(-1)**k)
 return r

def fixed(p,anti=False):
 if anti:
  if p.count(1)>1 or any(d!=1 and d%4!=2 for d in p):return 0
 else:
  if any(d%2==0 for d in p):return 0
 e=sum(d//2 for d in p)+sum(gcd(d,e) for i,d in enumerate(p) for e in p[i+1:])
 return 2**e

def z(p):
 r=1
 for d,m in Counter(p).items():r*=d**m*factorial(m)
 return r

def recover():
 be={};bt={}
 for n in range(N+1):
  for p in parts(n):
   for a,anti in [(be,False),(bt,True)]:
    v=fixed(p,anti)
    if v:a[p]=F(v,z(p))
 ae=add(be,{():1},-1);at=add(bt,{():1},-1)
 se=add({():1},inv(be),-1)
 st=add(mul(bt,add({():1},psi(se,2),-1)),{():1},-1)
 re=add(se,{(1,):1},-1);rt=add(st,{(1,):1},-1)
 pe={};pt={}
 for n in range(3,N+1):
  for p,r,anti in [(pe,re,False),(pt,rt,True)]:
   residual=add(r,comp(p,ae,at if anti else None),-1)
   p.update({q:v for q,v in residual.items() if sum(q)==n})
 if comp(pe,ae)!=re or comp(pt,ae,at)!=rt:raise RuntimeError('composition')
 if any(v<0 for p in (pe,pt) for v in p.values()):raise RuntimeError('negative prime coefficient')
 if any(not fixed(q,True) for q in pt):raise RuntimeError('anti support')
 for a in (be,bt,se,st,pe,pt):
  if any(sum(v for p,v in a.items() if sum(p)==n).denominator!=1 for n in range(N+1)):raise RuntimeError('noninteger total')
 return be,bt,se,st,pe,pt

def serial(a):return {','.join(map(str,p)):str(v) for p,v in sorted(a.items(),key=lambda x:(sum(x[0]),x[0]))}
def counts(a):return [int(sum(v for p,v in a.items() if sum(p)==n)) for n in range(N+1)]
if __name__=='__main__':
 be,bt,se,st,pe,pt=recover()
 print(json.dumps({'N':N,'all':counts(be),'self_converse':counts(bt),'strong_self_converse':counts(st),'prime':counts(pe),'prime_self_converse':counts(pt),'prime_anti_index':serial(pt),'checks':'PASS coefficientwise mixed and ordinary plethysm; nonnegativity; integral totals; anti-support'},indent=2))
