#!/usr/bin/env python3
"""Independent original gap-coordinate matrices; no producer imports or transformed basis."""
from itertools import product,permutations
import json,sys
from fractions import Fraction as Q

def need(x,m):
 if not x:raise RuntimeError(m)
def identity(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def multiply(A,B):
 return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def minus(A,c):return [[x-c*(i==j) for j,x in enumerate(r)] for i,r in enumerate(A)]
def letter(k,p,b):
 states=[s+(f,) for s in product(range(3),repeat=k) for f in (0,2)];ix={s:i for i,s in enumerate(states)}
 A=[[0]*len(states) for _ in states]
 for s in states:
  for g in range(k+1):
   if s[g]==0 and g and not (b>>p[g-1]&1):continue
   t=list(s);t[g]=2 if g==k or not(b>>p[g]&1) else 1
   A[ix[tuple(t)]][ix[s]]+=1
 return states,A
if '--fail' in sys.argv:need(False,'intentional sentinel')
tests=0
for k in range(3):
 for p in permutations(range(k)):
  for ell in range(1,4):
   for w in product(range(1<<k),repeat=ell):
    states,A=letter(k,p,w[0]);M=identity(len(A))
    for b in w:M=multiply(letter(k,p,b)[1],M)
    R=identity(len(M))
    for j in range(k+2):
     power=k+1 if j==0 else k-j+2
     for _ in range(power):R=multiply(minus(M,j**ell),R)
    need(not any(any(row) for row in R),'original basis annihilator')
    tests+=1
def scalar(k,w,m):
 total=0
 for p in permutations(range(k)):
  states,A=letter(k,p,w[0]);v=[int(s==(0,)*(k+1)) for s in states]
  for b in w*m:
   A=letter(k,p,b)[1];v=[sum(x*y for x,y in zip(row,v)) for row in A]
  for st,z in zip(states,v):
   if all(st[g]==2 or (st[g]==0 and (g==0 or p[g-1]<p[g])) for g in range(k)):total+=z
 return total
examples=[]
for k,w in ((1,(1,0)),(2,(3,2,1))):
 vals=[]
 for m in range(14):
  got=scalar(k,w,m)
  if k==1:expected=Q(4*4**m+6*m+5,9)
  else:expected=Q(236833,811538)+Q(2650,4459)*m+Q(13,49)*m*m+(Q(426044,866761)+Q(1684,6517)*m)*8**m+Q(26433,122018)*27**m
  need(got==expected,'physical scalar recurrence certificate')
  vals.append(got)
 examples.append(dict(k=k,mask_block=w,initial_counts=vals))
print(json.dumps({'examples':examples,'status' :'PASS','original_basis_annihilator_blocks':tests,'scope':'k=0,1,2; all apex permutations; all nonempty blocks through length three'},indent=2,sort_keys=True))
