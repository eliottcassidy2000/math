"""Exact standard-library controls; the all-order theorem is proved in PROOF.md."""
from fractions import Fraction as F
from math import comb,factorial

def falling(x,r):
 out=1
 for i in range(r):out*=x-i
 return out

def elementary(values):
 e=[1]+[0]*len(values)
 for a in values:
  for r in range(len(values),0,-1):e[r]+=a*e[r-1]
 return e

cases=0
for ell in range(0,11):
 for size in range(-3,15):
  values=list(range(size-1,size-ell-1,-1));direct=elementary(values);newton=[F(1)]
  for r in range(1,min(ell,5)+1):
   newton.append(sum((-1)**(a-1)*newton[r-a]*sum(v**a for v in values) for a in range(1,r+1))/r)
   assert newton[r]==direct[r];cases+=1
print('PASS exact Newton versus direct subset coefficient controls:',cases)
for b in range(2,31):
 q=elementary(list(range(3*b-1,2*b-1,-1)))[2]
 q-=2*b*(b-1)*sum(range(3*b-2,2*b-1,-1))
 q+=5*falling(b,3)+2*falling(b,4)
 expected=F(b*(b-1)*(3*b*b+17*b+2),24)
 assert q==expected,(b,q,expected)
assert F(2*1*(12+34+2),24)/factorial(2)==2
assert F(3*2*(27+51+2),24)/factorial(3)==F(10,3)
print('PASS real-inventory pure-three Q2 controls b=2..30 and constant-boundary examples.')
for j in range(7):
 for k in range(max(1,j),16):
  for d in range(j+1):
   ell=k-d;r=j-d
   assert 0<=r<=ell
   for b in range(k+1):
    c=k-b
    for U in range(2*d+1):
     for V in range(2*d-U+1):
      if U>b or V>c:assert falling(b,U)*falling(c,V)==0
print('PASS coefficient-domain and infeasible-consumption zero-factor controls.')
