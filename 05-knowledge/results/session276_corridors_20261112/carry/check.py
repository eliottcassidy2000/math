from functools import lru_cache
from math import gcd
import json
@lru_cache(None)
def g(x):
 if x==0:return 0
 a=x//2
 return g(a)+(0 if x%2 else a+a%2)
def F(A,B):return 15381+3357*(g(A)+g(B))+50004*A*B+25915*A+20732*B+1865*(A%2)+1492*(B%2)
def D(k,d):return (2**(k+1)-3)*d+2**k-d%2
checks=0
for k in range(1,14):
 for d in range(256):
  z=2**k*(2*d+1)
  assert g(z)-g(z-1)==D(k,d);checks+=1
for a in range(256):
 for b in range(32):
  assert abs(F(a+1,b)-F(a,b))>=6
  assert abs(F(a,b+1)-F(a,b))>=6
  checks+=2
# Entire nearest-neighbor overlap classification, bounded replay of the all-input proof.
for a in range(128):
 for b in range(1,128):
  assert F(a+1,b+1)-F(a,b)>4
  d=F(a+1,b-1)-F(a,b)
  if abs(d)<=4:assert (a%2,b%2,d)==(1,0,4)
  checks+=2
examples=[]
for j in [0,1,2,3,17,100,1000,10**10,10**30]:
 a=29527+1362656*j;b=29140+1344752*j
 assert (a+1)%16==8 and b%8==4
 assert F(a+1,b-1)-F(a,b)==4
 examples.append(dict(j=j,A=a,B=b,N=F(a,b),counts=[F(a,b)+i for i in [0,2,4,6]]))
 checks+=1
# Formula4, all parity/carry strata in an independent finite matrix of inputs.
for k in range(1,6):
 for l in range(1,6):
  ak=(3357-50004)*2**(k+1)-3*3357;al=(3357-50004)*2**(l+1)-3*3357
  for d in range(12):
   for e in range(12):
    a=2**k*(2*d+1)-1;b=2**l*(2*e+1)
    r,alpha=divmod(d,2);t,beta=divmod(e,2)
    c=(50004-3357)*(2**l-2**k)+4810+ak*alpha-al*beta+3357*(beta-alpha)
    assert F(a+1,b-1)-F(a,b)==2*ak*r-2*al*t+c
    checks+=1
print(json.dumps(dict(status='PASS',checks=checks,examples=examples),indent=2))
