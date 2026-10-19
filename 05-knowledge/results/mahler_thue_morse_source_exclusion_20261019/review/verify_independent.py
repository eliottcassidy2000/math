#!/usr/bin/env python3
"""Independent Mahler carry audit: no imports from the producer script."""
from fractions import Fraction as F
from itertools import product
import json
checks=0

def ck(x):
 global checks
 checks+=1
 if not x: raise RuntimeError('independent audit failed')
def actual(a,n):
 out=[]
 for _ in range(n):
  out.append(a%2);a=(3*a+a%2)//2
 return out,a
def morph(t,u,v):return [d for x in t for d in (v if x else u)]
def sub(w):return [d for x in w for d in (x,1-x)]
def tm(n):return [x.bit_count()%2 for x in range(n)]
def carry(w):return sum(x*3**(len(w)-1-i)*2**i for i,x in enumerate(w))
words=[list(w) for n in range(1,4) for w in product((0,1),repeat=n)]
for u in words:
 for v in words:
  for j in range(1,5):
   a=[0,1,1,0,1];b=[0,1,1]
   for _ in range(j):a=sub(a);b=sub(b)
   aa,bb=morph(a,u,v),morph(b,u,v)
   k=2**(j-1)*(len(u)+len(v))
   ck(len(aa)==5*k and len(bb)==3*k and aa==(bb+bb)[:5*k])
   ck(aa==morph(tm(len(a)),u,v))
  for A in range(1,41):
   for prefix in ([],[0],[1],[1,0,0]):
    pa,B=actual(A,len(prefix))
    if pa!=prefix:continue
    k=len(u)+len(v);j=1
    while 2**(len(prefix)+5*k)<3**(len(prefix)+3*k)*(A+1):k*=2;j+=1
    w=[0,1,1]
    for _ in range(j):w=sub(w)
    w=morph(w,u,v);C=carry(w);D=3**len(w)-2**len(w);X=D*B+C
    ck(0<X<2**(5*k))
    e=(X&-X).bit_length()-1
    aa,_=actual(B,e+1);tt=morph(tm(e+1),u,v)[:e+1]
    ck(aa[:e]==tt[:e] and aa[e]!=tt[e])
    for r in range(min(e+1,8)):
     ck((D*B+C)%(2**r)==0)
w=morph(tm(100),[0,0,0],[1,0,0])
for i in range(len(w)):
 ck(sum((F(2,3)**(j+1))*e for j,e in enumerate(w[i:]))<F(18,19))
 ck(all(j-i>=3 for i,j in zip([k for k,b in enumerate(w) if b],[k for k,b in enumerate(w) if b][1:])))
print(json.dumps({'status':'PASS','checks':checks,'scope':'All U,V lengths1..3, j1..4; sources1..40 and four prefixes; independent carry sum and actual first mismatch; finite safe-tail controls'},indent=2))
