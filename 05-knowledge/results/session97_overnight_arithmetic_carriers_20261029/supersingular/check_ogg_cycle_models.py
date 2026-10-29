#!/usr/bin/env python3
"""Explicit F5 hyperelliptic models for the odd-order Ogg cycle lengths."""
from math import isqrt
CHECKS=0
def check(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(msg)
def rem(a,b,p=5):
 a=list(a)
 while len(a)>=len(b):
  t=a[-1]*pow(b[-1],-1,p)%p;d=len(a)-len(b)
  for j,x in enumerate(b):
   if x:a[j+d]=(a[j+d]-t*x)%p
  while a and not a[-1]:a.pop()
 return a
def mul(a,b,f,p=5):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:c[i+j]=(c[i+j]+x*y)%p
 return rem(c,f,p)
def power(a,n,f,p=5):
 z=[1]
 while n:
  if n&1:z=mul(z,a,f,p)
  a=mul(a,a,f,p);n//=2
 return z
def gcd_poly(a,b,p=5):
 while b:a,b=b,rem(a,b,p)
 return a
def order2(n):
 d=1;x=2%n
 while x!=1:x=2*x%n;d+=1
 return d
OGG=[2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]
print('Ogg value, binary multiplicative order, kernel count for cycle(1,value)')
for n in OGG:
 if n==2:print(2,'not applicable: odd-cycle theorem');continue
 d=order2(n);print(n,d,0 if d%2==0 else 2**((n-1)//(2*d)))
check([n for n in OGG if n!=2 and order2(n)%2]==[7,23,31,47,71],'Ogg odd-order subset')
# A finite, reproducible low-degree coefficient search, stopping at first certificate.
for n in [7,23,31,47,71]:
 check(all(n%a for a in range(2,isqrt(n)+1)),'prime degree')
 found=None
 for mask in range(1,3125):
  a=[];m=mask
  for _ in range(5):a.append(m%5);m//=5
  if not a[0]:continue
  f=a+[0]*(n-5)+[1]
  if len(gcd_poly(f,[0,4,0,0,0,1]))!=1:continue
  if power([0,1],5**n,f)==[0,1]:found=(mask,a,f);break
 check(found is not None,('model found',n))
 mask,a,f=found
 check(len(gcd_poly(f,[0,4,0,0,0,1]))==1,'Rabin proper-divisor check')
 check(power([0,1],5**n,f)==[0,1],'Rabin full-degree check')
 print('MODEL',n,'genus',(n-1)//2,'low coefficients',a,'candidate mask',mask,'kernels',2**((n-1)//(2*order2(n))))
print('CHECKS',CHECKS)
