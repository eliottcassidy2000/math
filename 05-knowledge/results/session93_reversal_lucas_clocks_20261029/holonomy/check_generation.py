#!/usr/bin/env python3
"""Independent explicit-permutation check of recursive triangle generators."""
import json
from collections import deque
from sympy.combinatorics import Permutation, PermutationGroup
checks=0
def check(b,why):
 global checks
 checks+=1
 if not b: raise RuntimeError(why)
def ident(n):return tuple(range(n))
def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inv(a):
 b=[0]*len(a)
 for i,j in enumerate(a):b[j]=i
 return tuple(b)
def power(a,e):
 r=ident(len(a))
 while e:
  if e&1:r=mul(r,a)
  a=mul(a,a);e//=2
 return r
def root(a):
 b=[None]*len(a);seen=set()
 for i in range(len(a)):
  if i in seen:continue
  cyc=[];j=i
  while j not in seen:cyc.append(j);seen.add(j);j=a[j]
  l=len(cyc);check(l%2==1,'odd cycle')
  for q,j in enumerate(cyc):b[j]=cyc[(q+(l+1)//2)%l]
 return tuple(b)
def digits(i,k):
 a=[0]*k
 for j in range(k-1,-1,-1):a[j]=i%3;i//=3
 return a
def number(a):
 r=0
 for x in a:r=3*r+x
 return r
def alpha(k):return tuple(number([(-x)%3 for x in digits(i,k)]) for i in range(3**k))
def digit_rotation(k,j):
 out=[]
 for i in range(3**k):
  a=digits(i,k)
  if not any(a[:j]):a[j]=(a[j]+1)%3
  out.append(number(a))
 return tuple(out)
def lift(a,k):
 n=3**k;s=n//3
 return tuple(a[i] if i<s else i for i in range(n))
def block(a,b,c):
 s=len(a);return tuple(a)+tuple(s+x for x in b)+tuple(2*s+x for x in c)
def edge(a,b):return root(mul(a,b))
def triangle(s,t,u):return mul(mul(inv(edge(u,s)),edge(u,t)),edge(t,s))
def triples(k):
 if k<2:return []
 s=alpha(k);r=digit_rotation(k,0);out=[]
 for j in range(k-1):
  v=digit_rotation(k-1,j);y=lift(v,k)
  t=mul(power(r,2),s);u=mul(power(y,2),s)
  out.append((s,t,u,block(ident(len(v)),v,inv(v))))
 for a,b,c,h in triples(k-1):
  # Marks outside block0 retain alpha, not identity.
  beta=alpha(k-1)
  def mark(p):return mul(lift(mul(p,beta),k),s)
  out.append((mark(a),mark(b),mark(c),lift(h,k)))
 return out
def arc(a,b,k):
 da=digits(a,k);db=digits(b,k)
 for x,y in zip(da,db):
  if x!=y:return (y-x)%3==1
 raise ValueError('diagonal')
def closure(gens,n):
 e=ident(n);seen={e};todo=deque([e])
 while todo:
  a=todo.popleft()
  for g in gens:
   b=mul(a,g)
   if b not in seen:seen.add(b);todo.append(b)
 return seen
rows=[]
for k in range(1,6):
 n=3**k;s=alpha(k);gens=[]
 for a,b,c,expected in triples(k):
  for mark in (a,b,c):
   check(mul(mark,mark)==ident(n),'involution')
   for i in range(n):
    for j in range(i):check(arc(mark[i],mark[j],k)!=arc(i,j,k),'anti-isomorphism')
  h=triangle(a,b,c)
  check(h==expected,'exact triangle coordinates')
  check(mul(h,s)==mul(s,h),'centralizer')
  gens.append(h)
 expected_order=3**((n-2*k-1)//4)
 group=PermutationGroup([Permutation(g) for g in gens] or [Permutation(list(range(n)))])
 check(group.order()==expected_order,'Schreier-Sims order')
 check(len(gens)==k*(k-1)//2,'minimal generator count')
 if k<=3:check(len(closure(gens,n))==expected_order,'independent BFS closure')
 rows.append(dict(k=k,n=n,generators=len(gens),centralizer_order=str(expected_order)))
print(json.dumps(dict(status='PASS',checks=checks,rows=rows),indent=2))
