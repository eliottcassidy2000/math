#!/usr/bin/env python3
"""Exact harmonic-completion controls; explicit gates survive python -O."""
from fractions import Fraction as F
from math import gcd,lcm
from itertools import combinations_with_replacement
import json
checks=0

def gate(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)

def E(q):return [e for e in range(2,q) if gcd(e,q)==1]+[q+1]
def speeds(q,m,k):return sorted(set(range(1,q))|{q+1}|{q*m*j for j in range(1,k+1)})
def dist(v,t):
 z=v*t%1;return min(z,1-z)
def survivor(V,r):
 parts=[(F(0),F(1))]
 for v in sorted(V):
  nxt=[]
  for a,b in parts:
   low=max(0,(v*a-r).__ceil__()-1);high=min(v-1,(v*b-r).__floor__())
   for j in range(low,high+1):
    c=max(a,(j+r)/v);d=min(b,(j+1-r)/v)
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

def prediction(q,m):
 n=q+1;points=set();B=F(0)
 for e in E(q):
  a=pow(e,-1,q)
  for j in range(1,m//e+1):
   if gcd(j,n)>1:continue
   t=F(a,q)-F(j,q*n*m);points|={t,1-t}
   ell=pow(j,-1,n)
   B+=F(2,q*m)*(F(1,ell)+F(1,n-ell))
 return sorted(points),B

def count(q,m):return 2*sum(sum(gcd(j,q+1)==1 for j in range(1,m//e+1)) for e in E(q))
def slope(V,r,points):
 B=F(0)
 for t in points:
  P=[v for v in V if v*t%1==r];N=[v for v in V if v*t%1==1-r]
  gate(bool(P and N),('owners',q,m,t))
  B+=F(1,max(P))+F(1,max(N))
 return B

rows=[]
# Every harmonic-top preserving row in this complete parameter rectangle.
for q in range(3,18):
 s=min(E(q));n=q+1
 for m in range(1,2*q+1):
  V=speeds(q,m,q);r=F(1,n);parts=survivor(V,r)
  gate(bool(parts)==(m>=s),('threshold',q,m))
  if m<s:continue
  pts,B=prediction(q,m)
  gate(parts==[(x,x) for x in pts],('maximizer census',q,m))
  gate(len(pts)==count(q,m),('count',q,m))
  gate(slope(V,r,pts)==B,('physical slope',q,m))
  # Small guarded relaxation: all event radii have denominator <=2max(V).
  eps=F(1,4*n*max(V));relaxed=survivor(V,r-eps)
  gate(sum(b-a for a,b in relaxed)==eps*B,('literal mass germ',q,m))
  rows.append({'q':q,'m':m,'points':len(pts),'B':str(B)})
# Independent full maximum engine on a smaller universe, including the unclassified lower chamber.
maxrows=0
for q in range(3,9):
 s=min(E(q))
 for m in range(1,7):
  for k in sorted({1,2,q-1,q}):
   V=speeds(q,m,k);got=maximum(V)
   if 2*m>s:
    want=min(F(m,q*m+s),F(1,k+1));gate(got==want,('global maximum',q,m,k,got,want))
    if F(m,q*m+s)<F(1,k+1):
     a=pow(s,-1,q);h=(a*s-1)//q;t=F(a*m+h,q*m+s)
     parts=survivor(V,got)
     gate(parts==[(x,x) for x in sorted([t,1-t])],('strict maximizers',q,m,k))
     gate(slope(V,got,[t,1-t])==2*(F(1,s)+F(1,q*m)),('strict slope',q,m,k))
   if k==q:gate((got==F(1,q+1))==(m>=s),('global threshold',q,m))
   maxrows+=1
# Exact affine-period identities, not merely fitted asymptotics.
for q in range(3,13):
 n=q+1;L=n*lcm(*E(q));HE=sum(F(1,e) for e in E(q));HU=sum(F(1,a) for a in range(1,n) if gcd(a,n)==1)
 phi=sum(gcd(a,n)==1 for a in range(1,n))
 for m in [min(E(q)),min(E(q))+1,2*q+1]:
  _,B0=prediction(q,m);_,B1=prediction(q,m+L)
  gate(count(q,m+L)-count(q,m)==2*L*F(phi,n)*HE,('count shift',q,m))
  gate((m+L)*B1-m*B0==4*L*HE*HU/(q*n),('weight shift',q,m))
# The actual eighth phase carrier and derivative hostile.
odd8=[]
for m in range(2,17):
 pts,B=prediction(7,m);phases=sorted({(7*m*t)%1 for t in pts})
 gate(phases==([F(1,8),F(7,8)] if m<6 else [F(a,8) for a in (1,3,5,7)]),('eighth phases',m))
 odd8.append({'m':m,'points':len(pts),'B':str(B),'phases':[str(x) for x in phases]})
for m in range(2,122):
 A=sum(sum((-1)**((j*j-1)//8) for j in range(1,m//e+1,2)) for e in E(7))
 _,B=prediction(7,m)
 gate(abs(A)<=6,('character bound',m))
 gate(B==F(88*count(7,m)+64*A,735*m),('character correction',m))
for j in range(1,65,2):
 ell=pow(j,-1,8);chi=(-1)**((j*j-1)//8)
 gate(F(1,ell)+F(1,8-ell)==F(88+32*chi,105),('character weight',j))
gate(prediction(7,2)[1]==F(8,49),'critical owner hostile')
gate(4*sum(F(1,e) for e in E(7))*sum(F(1,a) for a in (1,3,5,7))/56==F(33,175),'eighth limiting B')
print(json.dumps({'status':'PASS','checks':checks,'top_rows':len(rows),'top_rectangle':'q3..17,m1..2q','full_maximum_rows':maxrows,'full_maximum_rectangle':'q3..8,m1..6,k in {1,2,q-1,q}','period_checks':'q3..12, three m each','odd_eighth_rows':odd8},indent=2))
