from fractions import Fraction as F
from math import gcd
import json
INF=float('inf')
def v5(n):
 if n==0:return INF
 if isinstance(n,F):return v5(n.numerator)-v5(n.denominator)
 a=0
 while n%5==0:a+=1;n//=5
 return a
def mod5(x):
 x=F(x);return x.numerator*pow(x.denominator,-1,5)%5
def parameters(m,b):
 t=next(t for t in (1,2,4) if pow(mod5(m),t,5)==1)
 return t,(v5(m**t-1) if t>1 else min(v5(m-1),v5(b)))
def depth(m,b,n,c):
 t,r=parameters(m,b); C=v5((m-1)*c-b*(n-1))
 if t!=1:
  if pow(mod5(n),t,5)!=1:return 0
  beta=v5(n**t-1)
  return min(C,beta if beta<r else INF)
 s,h=v5(m-1),v5(b)
 alpha=v5(n-1) if s<=h else v5(c)
 return max(0,min(C-r,alpha if alpha<r else INF))
def group(m,b,q):
 s={(1,0)}; a,c=1,0
 while True:
  a,c=a*m%q,(c*m+b)%q
  if (a,c)==(1,0):return s
  if (a,c) in s:raise RuntimeError('bad orbit')
  s.add((a,c))
def require(p,msg):
 if not p:raise RuntimeError(msg)
# Lift all affine pairs mod25 to positive nonidentity integer multipliers.
units=[n for n in range(1,25) if n%5]
checks=0; orderchecks=0
for mr in units:
 m=mr if mr!=1 else 26
 for br in range(25):
  b=br if br else 25
  for q,e in ((5,1),(25,2)):
   G=group(m,b,q);t,r=parameters(m,b)
   require(len(G)==t*5**max(0,e-r),('order',m,b,q))
   orderchecks+=1
   for nr in units:
    n=nr if nr!=1 else 26
    for cr in range(25):
     c=cr if cr else 25
     require(((n%q,c%q) in G)==(e<=depth(m,b,n,c)),('membership',m,b,n,c,q))
     checks+=1
# Higher-depth adversaries: unit-pivot and translation-pivot stalls.
for m,b in [(26,1),(126,5),(626,25),(6,25),(126,625),(2,5),(7,25),(57,125)]:
 for q,e in [(125,3),(625,4)]:
  G=group(m,b,q);t,r=parameters(m,b)
  require(len(G)==t*5**max(0,e-r),('order high',m,b,q));orderchecks+=1
  for n in range(1,150):
   if n%5==0:continue
   if n==1:n=626
   for c in range(0,75,3):
    require(((n%q,c%q) in G)==(e<=depth(m,b,n,c)),('high',m,b,n,c,q));checks+=1
print(json.dumps({'status':'PASS','arbitrary_affine_membership_checks':checks,'order_checks':orderchecks,'scope':'All affine residue pairs mod25 at depths1/2, plus stalled-pivot families at depths3/4; positive integer lifts.'},indent=2))
