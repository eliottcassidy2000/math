from fractions import Fraction as F
from math import gcd
import json
checks=0
def gate(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def exact_intersection(V,r):
 parts=[(F(0),F(1))]
 for v in V:
  nextparts=[]
  for a,b in parts:
   # Only integer bands whose safe interval could meet [a,b].
   for k in range(max(0,(v*a-r).__floor__()),min(v-1,(v*b-r).__floor__())+1):
    lo=max(a,(k+r)/v);hi=min(b,(k+1-r)/v)
    if lo<=hi:nextparts.append((lo,hi))
  parts=sorted(set(nextparts))
 return parts
def dist(v,t):
 z=(v*t)%1;return min(z,1-z)
rows=[]
for q in range(3,81):
 C=tuple(range(1,q))+(q+1,)
 units=[a for a in range(1,q) if gcd(a,q)==1]
 s=next(a for a in units if a>=2)
 owners=[]
 for a in units:
  p=pow(a,-1,q);n=q-p
  if p==1:p=q+1
  if n==1:n=q+1
  owners.append((a,p,n))
 gate(min(min(p,n) for a,p,n in owners)==s,('least owner',q))
 for weight in (F(1,7),F(1,2),F(6,7),F(1)):
  r=F(1,q+2)+(F(1,q)-F(1,q+2))*weight
  delta=F(1,q)-r
  predicted=[(F(a,q)-delta/p,F(a,q)+delta/n) for a,p,n in owners]
  gate(exact_intersection(C,r)==predicted,('complete window',q,r))
 for m in range(1,s+4):
  M=max(F(1,q+2),F(m,q*m+s))
  if 2*m<=s:
   gate(q*m%(q+2)!=0,('floor tail nonzero',q,m))
   gate(min(dist(v,F(1,q+2)) for v in C+(q*m,))>=M,('floor witness',q,m))
  else:
   endpoints=[F(a,q)-(F(1,q)-M)/p for a,p,n in owners if p==s]+[F(a,q)+(F(1,q)-M)/n for a,p,n in owners if n==s]
   for t in endpoints:gate(min(dist(v,t) for v in C+(q*m,))==M,('upper witness',q,m,t))
 rows.append({'q':q,'s':s,'primitive_centers':len(units)})
print(json.dumps(dict(status='PASS',checks=checks,cores='q=3..80',radii_per_core=4,method='independent exact interval intersection versus inverse-residue owner windows',rows=rows),indent=2))
