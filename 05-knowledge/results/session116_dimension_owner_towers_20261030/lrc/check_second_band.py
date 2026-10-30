#!/usr/bin/env python3
"""Exact independent interval intersections for second-band owner promotion."""
from fractions import Fraction as F
from math import gcd,lcm
from itertools import combinations_with_replacement
import json
checks=0

def gate(ok,label):
 global checks
 checks+=1
 if not ok: raise RuntimeError(label)

def survivor(V,r):
 parts=[(F(0),F(1))]
 for v in sorted(V):
  nxt=[]
  for a,b in parts:
   lo=max(0,(v*a-r).__ceil__()-1); hi=min(v-1,(v*b-r).__floor__())
   for j in range(lo,hi+1):
    c=max(a,(j+r)/v); d=min(b,(j+1-r)/v)
    if c<=d:nxt.append((c,d))
  parts=sorted(set(nxt))
 return parts

def maximum(V):
 best=F(0)
 for den in {v+w for v,w in combinations_with_replacement(V,2)}:
  for a in range(1,den//2+1):
   r=min(min(a*v%den,(-a*v)%den) for v in V)
   if r*best.denominator>best.numerator*den:best=F(r,den)
 return best

def owners(q,A):return [u+q*(u in A) for u in range(1,q) if gcd(u,q)==1]
def core(q,A):return set(range(1,q))|{q+b for b in A}
def atlas(q,A,r):
 ep=F(1,q)-r; out=[]
 for a in range(1,q):
  if gcd(a,q)>1:continue
  p=pow(a,-1,q);n=q-p
  out.append((F(a,q)-ep/(p+q*(p in A)),F(a,q)+ep/(n+q*(n in A))))
 return out

def prediction(q,A,m):
 points=set();B=F(0);n=q+1
 for e in owners(q,A):
  for j in range(1,m//e+1):
   if gcd(j,n)>1:continue
   t=F(pow(e,-1,q),q)-F(j,q*n*m);points|={t,1-t}
   ell=pow(j,-1,n);B+=F(2,q*m)*(F(1,ell)+F(1,n-ell))
 return sorted(points),B

def slope(V,r,pts):
 B=F(0)
 for t in pts:
  P=[v for v in V if v*t%1==r];N=[v for v in V if v*t%1==1-r]
  gate(bool(P and N),('missing binding direction',V,t))
  B+=F(1,max(P))+F(1,max(N))
 return B

def periodic_stats(q,A,m):
 n=q+1;N=0;W=F(0)
 weights={j:F(1,pow(j,-1,n))+F(1,n-pow(j,-1,n)) for j in range(1,n) if gcd(j,n)==1}
 for e in owners(q,A):
  cycles,rem=divmod(m//e,n)
  N+=cycles*len(weights)+sum(j<=rem for j in weights)
  W+=cycles*sum(weights.values(),F(0))+sum((w for j,w in weights.items() if j<=rem),F(0))
 return 2*N,2*W/(q*m)

atlas_rows=harmonic_rows=strict_rows=0
for q in range(3,10):
 seen=set();counts={}
 for mask in range(1<<(q-2)):
  A={1}|{u for u in range(2,q) if mask>>(u-2)&1};E=owners(q,A);s=min(E);C=core(q,A)
  seen.add(tuple(E));counts[s]=counts.get(s,0)+1
  for alpha in (F(1,20),F(1,3),F(2,3),F(19,20),F(1)):
   r=F(1,q+2)+alpha*(F(1,q)-F(1,q+2))
   gate(survivor(C,r)==atlas(q,A,r),('atlas',q,sorted(A),r));atlas_rows+=1
  for m in sorted({1,max(1,s-1),s,s+1,2*s+1}):
   V=C|{q*m*j for j in range(1,q+1)};r=F(1,q+1)
   got=survivor(V,r);pts,B=prediction(q,A,m)
   gate(bool(got)==(m>=s),('preservation iff',q,A,m))
   gate(got==[(t,t) for t in pts],('census',q,A,m))
   if pts:
    gate(slope(V,r,pts)==B,('mass owners',q,A,m))
    eps=F(1,4*(q+1)*max(V));near=survivor(V,r-eps)
    gate(sum(b-a for a,b in near)==eps*B,('literal mass',q,A,m))
   harmonic_rows+=1
  if q<=6:
   for m in sorted({s//2+1,s,max(1,s-1)}):
    if 2*m<=s:continue
    for k in sorted({1,q-1,q}):
     V=C|{q*m*j for j in range(1,k+1)};want=min(F(m,q*m+s),F(1,k+1))
     gate(maximum(V)==want,('global max',q,A,m,k));strict_rows+=1
     if F(m,q*m+s)<F(1,k+1):
      a=pow(s,-1,q);h=(a*s-1)//q;t=F(a*m+h,q*m+s);pts=sorted([t,1-t])
      gate(survivor(V,want)==[(x,x) for x in pts],('strict census',q,A,m,k))
      gate(slope(V,want,pts)==2*(F(1,s)+F(1,q*m)),('strict slope',q,A,m,k))
 phi=sum(gcd(u,q)==1 for u in range(1,q))
 gate(len(seen)==2**(phi-1),('atlas classification',q))
 for s,c in counts.items():
  want=2**(q-1-phi)*(1 if s==q+1 else 2**sum(gcd(u,q)==1 for u in range(s+1,q)))
  gate(c==want,('minimum classification',q,s))
# Surgery, period and asymptotic coefficient identities for all atlas patterns q<=8.
for q in range(3,9):
 U=[u for u in range(1,q) if gcd(u,q)==1];n=q+1
 for mask in range(1<<(len(U)-1)):
  A={1}|{u for i,u in enumerate(U[1:]) if mask>>i&1};E=owners(q,A);s=min(E)
  HE=sum((F(1,e) for e in E),F(0));HU=sum((F(1,u) for u in range(1,n) if gcd(u,n)==1),F(0));phin=sum(gcd(u,n)==1 for u in range(1,n))
  L=n*lcm(*E)
  # Period sums calculated by complete-period arithmetic; avoids enormous physical speed enumeration.
  for m in (s,s+1,2*s+3):
   pts,B=prediction(q,A,m)
   N0,B0=periodic_stats(q,A,m);N1,B1=periodic_stats(q,A,m+L)
   gate((N0,B0)==(len(pts),B),('periodic independent census',q,A,m))
   deltaN=N1-N0
   deltaMB=(m+L)*B1-m*B0
   gate(deltaN==2*L*phin*HE/n,('N period',q,A,m))
   gate(deltaMB==4*L*HE*HU/(q*n),('mB period',q,A,m))
   for u in U:
    if u in A:continue
    A2=A|{u};m2=max(m,min(owners(q,A2)))
    p1,b1=prediction(q,A,m2);p2,b2=prediction(q,A2,m2)
    js=[j for j in range(m2//(q+u)+1,m2//u+1) if gcd(j,n)==1]
    gate(len(p1)-len(p2)==2*len(js),('surgery count',q,A,u,m2))
    lost=F(2,q*m2)*sum((F(1,pow(j,-1,n))+F(1,n-pow(j,-1,n)) for j in js),F(0))
    gate(b1-b2==lost,('surgery weight',q,A,u,m2))
# All q=7 promotion patterns retain the eighth-character signature.
for mask in range(32):
 A={1}|{u for u in range(2,7) if mask>>(u-2)&1};s=min(owners(7,A))
 for m in range(s,3*s+2):
  pts,B=prediction(7,A,m);phases={7*m*t%1 for t in pts}
  gate((len(phases)==4)==(m>=3*s),('four phases threshold',A,m))
  charsum=sum((1 if j%8 in (1,7) else -1) for e in owners(7,A) for j in range(1,m//e+1,2))
  gate(B==F(88*len(pts)+64*charsum,735*m),('eighth character',A,m))
# Named hostile gates: composite minimum, full promotion, and no inherited lower floor.
gate(min(owners(7,{1,2,3}))==4,'composite minimum')
gate(survivor(core(3,{1}),F(1,5))!=atlas(3,{1},F(1,5)),'excluded chamber endpoint')
gate(min(owners(7,set(range(1,7))))==8,'full promotion')
V=core(3,{1,2})|{3,6,9}
gate(maximum(V)==F(1,7),'lower-floor counterexample')
receipt={'status':'PASS','checks':checks,'atlas_rows':atlas_rows,'harmonic_rows':harmonic_rows,'global_maximum_rows':strict_rows,'arithmetic':'exact Fraction; explicit gates; complete subsets q=3..9','scope':'explicit second-band family, open top chamber; no general LRC closure'}
print(json.dumps(receipt,sort_keys=True,indent=2))
