#!/usr/bin/env python3
"""Independent exception-safe Session56 audit. No producer imports."""
from itertools import product,permutations
from collections import Counter
from fractions import Fraction
from math import factorial,comb,inf
from pathlib import Path
import json,hashlib,subprocess
ROOT=Path('/tmp/session56-publication/clean-repo')
OUT=Path(__file__).parent

def check(v,msg):
 if not v: raise RuntimeError(msg)
def val(x,p):
 if x==0:return inf
 k=0
 while x%p==0:x//=p;k+=1
 return k
def parent_tree(a):
 # Stack construction independent of producer recursive splitting.
 parent=[None]*len(a);stack=[]
 for i in range(len(a)):
  last=None
  while stack and (a[i],i)<(a[stack[-1]],stack[-1]):last=stack.pop()
  if stack:parent[i]=stack[-1]
  if last is not None:parent[last]=i
  stack.append(i)
 return tuple(parent)
def lca(p,i,j):
 ancestors=set()
 while i is not None:ancestors.add(i);i=p[i]
 while j not in ancestors:j=p[j]
 return j
def hook(p):
 sizes=[1]*len(p)
 for i in range(len(p)):
  j=p[i]
  while j is not None:sizes[j]+=1;j=p[j]
 z=1
 for s in sizes:z*=s
 return z
checks=0
for n in range(8):
 signatures={};trees={}
 for a in product(range(3),repeat=n):
  t=parent_tree(a);s=tuple(min(range(i,j+1),key=lambda k:(a[k],k)) for i in range(n) for j in range(i,n))
  check(s==tuple(lca(t,i,j) for i in range(n) for j in range(i,n)),'stack LCA');checks+=1
  check(t not in trees or trees[t]==s,'tree implies signature');trees[t]=s
  check(s not in signatures or signatures[s]==t,'signature implies tree');signatures[s]=t
fibres=[]
for n in range(9):
 c=Counter(parent_tree(a) for a in permutations(range(n)))
 check(len(c)==comb(2*n,n)//(n+1),'Catalan count')
 check(all(k*hook(t)==factorial(n) for t,k in c.items()),'hook fibres')
 check(sum((Fraction(1,hook(t)) for t in c),Fraction())==1,'weighted mass')
 fibres.append([n,len(c),sum(c.values())])
vp_cases=0
for p in (2,3,5,7):
 for n in range(5):
  for a in product((-5,-1,0,1,5),repeat=n):
   values=set(a);H=0
   while values and p**H<=max(values)-min(values):H+=1
   order=sorted(values,key=lambda x:tuple((x//p**k)%p for k in range(H)))
   adj=[val(y-x,p) for x,y in zip(order,order[1:])]
   for i,x in enumerate(order):
    for j in range(i+1,len(order)):check(val(order[j]-x,p)==min(adj[i:j]),'LCP min')
   for x,y in product(a,repeat=2):
    for h in range(H+2):check(((x-y)%p**h==0)==(val(x-y,p)>=h),'residue hierarchy')
   if len(a)!=len(values):check(any(val(x-y,p)==inf for i,x in enumerate(a) for y in a[i+1:]),'duplicate infinite')
   vp_cases+=1
check(val(2,2)!=min(val(1,2),val(1,2)),'natural order hostile')
W=(0,9,1,16,15,14,29,21,30)
for p,q in ((2,3),(2,5),(3,7),(5,7)):
 for a,b in product(range(5),repeat=2):
  for x,y in product(W,repeat=2):check(((x-y)%(p**a*q**b)==0)==(val(x-y,p)>=a and val(x-y,q)>=b),'CRT exact')
A={0,4,8};B={2,4,6};check(A&B and not A<=B and not B<=A,'crossed clusters')
# Native parity first descent and affine map via rational pairs.
x=27;alpha=Fraction(1);beta=Fraction(0);word=[]
while True:
 e=x%2;word.append(e)
 if e:x=(3*x+1)//2;alpha=3*alpha/2;beta=(3*beta+1)/2
 else:x//=2;alpha/=2;beta/=2
 if x<27:break
q=2**len(word);a=int(alpha*q);b=int(beta*q);d=q-a
check((len(word),sum(word),x,b,d)==(59,37,23,1100931843921811423,126176846412426125),'Collatz source')
check(d%2 and b%d,'noninteger odd-denominator fixed point')
fuel=0
for h in range(1,9):
 m=2**(59*h);r=b*pow(d,-1,m)%m
 for lift in range(6):
  n=r+lift*m;v=val(d*n-b,2);K=v//59;z=n
  for _ in range(K):
   for bit in word:
    check(z%2==bit,'native itinerary');z=(3*z+1)//2 if bit else z//2
  check(val(d*z-b,2)<59,'actual exhaustion')
  for cap in (0,1,58,59,v,v+1):
   residue=b*pow(d,-1,2**cap)%(2**cap)
   check(min(cap,val(n-residue,2))==min(cap,v),'finite cap lower bound')
  fuel+=1
check(val(d*27-b,2)//59==1 and val(d*23-b,2)//59==0,'coverage guard')
# Simple atlas algebra controls.
def mm(A,B):return tuple(min(A[2*i+k]+B[2*k+j] for k in (0,1)) for i in (0,1) for j in (0,1))
mat=list(product((0,1),repeat=4))
for A,B,C in product(mat,repeat=3):check(mm(mm(A,B),C)==mm(A,mm(B,C)),'tropical associativity')
check(mm((1,)*4,(1,)*4)==(2,)*4,'nonidempotence')
for r in range(3,9):
 for s in range(r+1):check(min(abs(k)*(r-1)+s*abs(1-k)+(r-s)*abs(k) for k in range(-20,21))==min(s,r-1),'higher arity')
# Exact integer capacity rounding, no floating-point logarithms.
for n in range(101):
 N=n*(n-1)//2;c=comb(2*n,n)//(n+1);den=c*factorial(n);B=max(0,N-(den.bit_length()-1))
 check(den*2**B>=2**N and (B==0 or den*2**(B-1)<2**N),'capacity bound')
check(comb(14,7)//8==429 and 456>429,'OEIS calibration')
print(json.dumps({'status':'PASS','tied_arrays':checks,'fibres':fibres,'valuation_record_cases':vp_cases,'Collatz_replays':fuel,'finite_cap_cases':fuel*6,'tropical_triples':4096,'capacity_n':[0,100],'hostiles':['tree loses min value','natural order fails','duplicate depth infinite','finite cap only lower bound','crossed CRT no one hierarchy','fuel exhaustion not coverage'],'limits':'Finite controls supplement proofs; no global convergence, tournament runtime bound, or atlas implementation promotion.'},indent=2,sort_keys=True))
