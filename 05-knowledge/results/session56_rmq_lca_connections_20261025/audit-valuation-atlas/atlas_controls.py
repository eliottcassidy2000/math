from itertools import product
from math import isqrt
from fractions import Fraction
import json

def ck(v,s):
 if not v:raise RuntimeError(s)
# Charge via exact radical floor, independent of greedy atom implementation.
def B(n):return (3*(n+1)-isqrt(5*(n+1)**2)-1)//2
def Q(n):return (n-2*B(n),B(n))
for a,b in product(range(100),repeat=2):
 qa,qb,qc=Q(a),Q(b),Q(a+b);d=qa[1]+qb[1]-qc[1]
 ck(d in (-1,0,1) and qa[0]+qb[0]-qc[0]==-2*d,'carry exact')
for a,b,c in product(range(20),repeat=3):
 delta=lambda x,y:B(x)+B(y)-B(x+y)
 ck(delta(a,b)+delta(a+b,c)==delta(b,c)+delta(a,b+c),'carry cocycle')
na=0
for n in range(1,6):
 for masks in product(range(4),repeat=n):
  for l in range(n):
   for r in range(l,n):
    w=2**((r-l+1).bit_length()-1)
    def meet(inds):
     bits=3
     for i in inds:bits&=masks[i]
     return bits
    ck(meet(range(l,r+1))==meet(range(l,l+w))&meet(range(r-w+1,r+1)),'overlap AND');na+=1
for k in range(2,6):
 one=sum(c[0]!=c[1] for c in product(range(k),repeat=3))
 two=sum(c[0]!=c[1] and c[0]!=c[2] for c in product(range(k),repeat=3))
 ck((one,two)==(k*k*(k-1),k*(k-1)**2),'named ports')
for r in range(12):
 h=(33*9**r-9)//8;hn=(33*9**(r+1)-9)//8
 ck(hn==9*h+9 and 3*(r+1)+4-(3*r+4)==3,'address scales')
for n in range(1,101):
 x=n;a=Fraction(1);b=Fraction(0)
 for t in range(12):
  odd=x%2
  aa,bb=(3,1) if odd else (Fraction(1,2),0)
  a,b=aa*a,aa*b+bb;x=3*x+1 if odd else x//2
  ck(a*n+b==x,'standard itinerary')
print(json.dumps({'status':'PASS','carry_pairs':10000,'carry_cocycles':8000,'mask_interval_queries':na,'named_port_k':[2,5],'paired_address_checks':12,'standard_Collatz_steps':1200,'scope':'Atlas elementary controls only; inherited theorems checked as pinned scoped claims, proposals remain OPEN.'},indent=2,sort_keys=True))
