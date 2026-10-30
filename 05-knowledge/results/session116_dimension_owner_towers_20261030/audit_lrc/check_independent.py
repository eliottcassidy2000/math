#!/usr/bin/env python3
"""Independent literal safe-interval and affine-breakpoint audit; no producer imports."""
from fractions import Fraction as F
from math import gcd, lcm
from collections import Counter
import json
checks=0

def check(a,b,label):
 global checks
 checks+=1
 if a!=b: raise RuntimeError((label,a,b))

def intersect(X,Y):
 out=[]; i=j=0
 while i<len(X) and j<len(Y):
  a,b=X[i];c,d=Y[j]
  if max(a,c)<=min(b,d): out.append((max(a,c),min(b,d)))
  if b<d:i+=1
  else:j+=1
 return out

def good(V,r):
 X=[(F(0),F(1))]
 for v in sorted(set(V)):
  X=intersect(X,[(F(i)+r,F(i+1)-r) for i in range(v)]) if v==1 else intersect(X,[((F(i)+r)/v,(F(i+1)-r)/v) for i in range(v)])
 return X

def owners(q,A):return [u+q*(u in A) for u in range(1,q) if gcd(u,q)==1]
def core(q,A):return list(range(1,q))+[q+b for b in A]
def windows(q,A,r):
 ep=F(1,q)-r; E={e%q:e for e in owners(q,A)}
 return [(F(a,q)-ep/E[pow(a,-1,q)],F(a,q)+ep/E[q-pow(a,-1,q)]) for a in range(1,q) if gcd(a,q)==1]

def dist(x):return min(x%1,1-x%1)
def literal_max(V):
 # Every positive peak is intersection of a rising and a falling affine branch.
 V=sorted(set(V)); pts={F(0),F(1)}
 for i,v in enumerate(V):
  for w in V[i:]:pts.update(F(j,v+w) for j in range(1,v+w))
 vals={t:min(dist(v*t) for v in V) for t in pts}
 r=max(vals.values());return r, sorted(t for t,z in vals.items() if z==r)

def predictions(q,A,m):
 E=owners(q,A); n=q+1
 pts=[];B=F(0)
 for e in E:
  a=pow(e,-1,q)
  for j in range(1,m//e+1):
   if gcd(j,n)!=1:continue
   t=F(a,q)-F(j,q*n*m);pts.extend([t,1-t])
   ell=pow(j,-1,n);B+=F(2,q*m)*(F(1,ell)+F(1,n-ell))
 return sorted(pts),B

def literal_B(V,pts,r):
 ans=F(0)
 for t in pts:
  # Active increasing constraints bound extension to the left; decreasing to right.
  left=[];right=[]
  for v in set(V):
   x=(v*t)%1
   if x==r:left.append(v)
   if x==1-r:right.append(v)
  if not left or not right:raise RuntimeError('non-isolated maximum')
  ans+=F(1,max(left))+F(1,max(right))
 return ans

atlas_total=0
for q in range(3,10):
 U=[u for u in range(1,q) if gcd(u,q)==1]; atlas=Counter();mins=Counter()
 for bits in range(1<<(q-2)):
  A={1}|{b for b in range(2,q) if bits>>(b-2)&1}; E=owners(q,A)
  atlas[tuple(E)]+=1;mins[min(E)]+=1;atlas_total+=1
  for r in [F(1,q), (F(1,q)+F(1,q+2))/2, F(1,q+2)+F(1,100*q*(q+2))]:check(good(core(q,A),r),windows(q,A,r),('windows',q,A,r))
 check(len(atlas),2**(len(U)-1),('atlas count',q))
 check(set(atlas.values()),{2**(q-1-len(U))},('fiber',q))
 expected={u:2**sum(v>u for v in U)*2**(q-1-len(U)) for u in U if u>1};expected[q+1]=2**(q-1-len(U))
 check(dict(mins),expected,('min census',q))

# Literal complete maxima and terminal slope at every atlas in q<=5.
harmonic_cases=0
for q in range(3,6):
 for bits in range(1<<(q-2)):
  A={1}|{b for b in range(2,q) if bits>>(b-2)&1};s=min(owners(q,A))
  for m in sorted({1,s,max(1,s-1),s+1,2*s}):
   V=core(q,A)+[q*m*j for j in range(1,q+1)]
   # A literal safe interval at 1/(q+1) already identifies all preserving maxima.
   intervals=good(V,F(1,q+1)); predicted,B=predictions(q,A,m)
   check(intervals,[(t,t) for t in predicted],('complete maximizing set',q,A,m))
   check(bool(intervals),m>=s,('preservation iff',q,A,m))
   if m>=s:check(literal_B(V,predicted,F(1,q+1)),B,('mass',q,A,m))
   harmonic_cases+=1
  for k in [1,q]:
   m=s//2+1; V=core(q,A)+[q*m*j for j in range(1,k+1)]
   r,pts=literal_max(V);want=min(F(m,q*m+s),F(1,k+1));check(r,want,('harmonic max',q,A,m,k))
   if F(m,q*m+s)<F(1,k+1):
    a=pow(s,-1,q);h=(a*s-1)//q;t=F(a*m+h,q*m+s)
    check(pts,sorted([t,1-t]),('two endpoints',q,A,m,k))
    check(literal_B(V,pts,r),2*(F(1,s)+F(1,q*m)),('strict slope',q,A,m,k))

check(literal_max([1,2,3,4,5,6,9])[0],F(1,7),'hostile lower chamber')
# Surgery compares literal survivor sets, not merely floor counts.
for q in range(3,9):
 for u in range(2,q):
  if gcd(u,q)!=1:continue
  A={1}; Ap={1,u};m=3*(q+1);n=q+1
  old,Bo=predictions(q,A,m);new,Bn=predictions(q,Ap,m)
  lost=[j for j in range(m//(q+u)+1,m//u+1) if gcd(j,n)==1]
  check(len(old)-len(new),2*len(lost),('surgery count',q,u))
  check(set(new)<=set(old),True,('surgery inclusion',q,u))
  check(Bo-Bn,sum((F(2,q*m)*(F(1,pow(j,-1,n))+F(1,n-pow(j,-1,n))) for j in lost),F(0)),('surgery mass',q,u))
# q7 phase census and character resolution, including composite min and all promoted.
for A in [{1},{1,2,3},set(range(1,7))]:
 E=owners(7,A);s=min(E)
 for m in range(s,3*s+2):
  pts,B=predictions(7,A,m);phases={(7*m*t)%1 for t in pts}
  check(phases=={F(1,8),F(3,8),F(5,8),F(7,8)},m>=3*s,('four eighths',A,m))
  character=sum(1 if j%8 in (1,7) else -1 for e in E for j in range(1,m//e+1,2))
  check(B,F(88*len(pts)+64*character,735*m),('character',A,m))
# Exact quasiperiod increments in small moduli; direct phase enumeration.
for q in range(3,6):
 for A in [{1},set(range(1,q))]:
  E=owners(q,A); n=q+1;L=n*lcm(*E);m=max(E)
  pts,B=predictions(q,A,m);pts2,B2=predictions(q,A,m+L)
  HE=sum((F(1,e) for e in E),F(0));HU=sum((F(1,u) for u in range(1,n) if gcd(u,n)==1),F(0));phi=sum(gcd(u,n)==1 for u in range(1,n))
  check(len(pts2)-len(pts),2*L*phi*HE/n,('period count',q,A))
  check((m+L)*B2-m*B,4*L*HE*HU/(q*n),('period mass',q,A))
print(json.dumps({'status' :'PASS','checks':checks,'literal_core_subsets':atlas_total,'harmonic_cases':harmonic_cases},sort_keys=True))
