#!/usr/bin/env python3
"""Upper bound in the two-cluster cyclic row template."""
from math import prod,factorial
from fractions import Fraction as F
s=56;k=s//2;n=s*s
# In each near-regular template exactly half the rows have each degree.
deg=[k//2]* (k//2)+[k//2-1]*(k//2)
# In cyclic order with ties sent forward, row 0 is high, row 27 low.
t=7;a=3
e=F(27182818,10**7)
assert sum(F(1,factorial(i)) for i in range(13))>e
def test(mate):
 rows=[]
 for side in (0,1):
  for i,d in enumerate(deg):
   for j in range(s):
    x=(s-2)*d+k-1+s-1-j + (F(k*(s-1),t) if not side else k*t)
    x+=F(1 if not side else t,a) + 5*(1 if not side else t)
    if mate:
     if (side,i,j)==(0,0,0):x-=1
     if (side,i,j)==(0,k-1,s-1):x+=1-F(1,t)
     if (side,i,j)==(1,0,0):x+=t-1
     if (side,i,j)==(1,k-1,s-1):x+=1
    rows.append(x)
 rows.append(s*k*(a+F(a,t))+5*a)
 B=prod(rows)/e**(n+1)
 bs=B.numerator.bit_length()-B.denominator.bit_length()
 bs+=B.numerator>=B.denominator<<bs
 print(mate,bs)
 assert B.numerator < B.denominator<<26795 and bs==26795
 return bs
print([test(m) for m in (False,True)])
