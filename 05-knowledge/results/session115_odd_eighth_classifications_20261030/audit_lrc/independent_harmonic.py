"""Independent exact finite audit; imports no producer implementation."""
from fractions import Fraction as F
from math import gcd,lcm
import json
checks=0

def gate(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
def owners(q):return [e for e in range(2,q) if gcd(e,q)==1]+[q+1]
def formula(q,m):
 E=owners(q);n=q+1
 pts=set();B=F(0)
 for e in E:
  a=pow(e,-1,q)
  for j in range(1,m//e+1):
   if gcd(j,n)!=1:continue
   t=F(a,q)-F(j,q*n*m)
   pts.update((t,1-t))
   b=pow(j,-1,n)
   B+=F(2,q*m)*(F(1,b)+F(1,n-b))
 return pts,B

def literal_top(q,m):
 n=q+1;D=q*m*n;core=list(range(1,q))+[q+1]
 V=core+[q*m*k for k in range(1,q+1)]
 pts=[];B=F(0)
 for a in range(1,D):
  if gcd(a,n)!=1:continue
  if min(min(a*v%D,(-a*v)%D) for v in core)<q*m:continue
  pts.append(F(a,D))
  pp=[v for v in V if a*v%D==q*m]
  nn=[v for v in V if a*v%D==D-q*m]
  gate(bool(pp and nn),('owners',q,m,a))
  B+=F(1,max(pp))+F(1,max(nn))
 return set(pts),B

def max_clearance(V):
 best=F(0);points=set()
 for d in {v+w for v in V for w in V}:
  for a in range(1,d):
   num=min(min(a*v%d,(-a*v)%d) for v in V)
   value=F(num,d)
   if value>best:best=value;points={F(a,d)}
   elif value==best:points.add(F(a,d))
 return best,points

rows=0;boundary=[]
for q in range(3,16):
 s=min(owners(q))
 for m in range(1,31):
  pts,B=literal_top(q,m);want,wB=formula(q,m)
  gate(bool(pts)==(m>=s),('iff',q,m))
  gate(pts==want,('points',q,m))
  gate(B==wB,('weights',q,m,B,wB))
  if m==s:
   gate(len(pts)==2 and B==F(2*(q+1),q*q*s),('collision',q,m))
   boundary.append({'q':q,'m':m,'count':len(pts),'B':str(B)})
  rows+=1
maximum_rows=0
for q in range(3,8):
 s=min(owners(q))
 for m in range(1,6):
  if 2*m<=s:continue
  for k in range(1,q+1):
   V=sorted(set(range(1,q))|{q+1}|{q*m*j for j in range(1,k+1)})
   got,P=max_clearance(V);want=min(F(m,q*m+s),F(1,k+1))
   gate(got==want,('maximum',q,m,k,got,want))
   if F(m,q*m+s)<F(1,k+1):gate(len(P)==2,('strict points',q,m,k))
   maximum_rows+=1
# Formula-level long-period identities, computed from complete residue blocks.
period_rows=0
for q in range(3,10):
 E=owners(q);s=min(E);n=q+1;L=n*lcm(*E)
 H_E=sum((F(1,e) for e in E),F(0));U=[a for a in range(1,n) if gcd(a,n)==1]
 H_U=sum((F(1,a) for a in U),F(0))
 for m in (s,s+1,2*s+3):
  P,B=formula(q,m);Q,C=formula(q,m+L)
  gate(len(Q)-len(P)==2*L*len(U)*H_E/n,('count period',q,m))
  gate((m+L)*C-m*B==4*L*H_E*H_U/(q*n),('weight period',q,m))
  period_rows+=1
E=owners(7)
def chi(j):return 1 if j%8 in (1,7) else -1
for m in range(2,101):
 P,B=formula(7,m)
 A=sum(chi(j) for e in E for j in range(1,m//e+1,2))
 A2=sum(chi(j) for e in E for j in range(1,(m+960)//e+1,2))
 gate(abs(A)<=6 and A==A2,('character period',m))
 gate(B==F(88*len(P)+64*A,735*m),('character weight',m))
print(json.dumps({'status':'PASS','checks':checks,'literal_top_rows':rows,'complete_maximum_rows':maximum_rows,'period_rows':period_rows,'critical_rows':boundary},indent=2))
