#!/usr/bin/env python3
from itertools import product,combinations
from functools import lru_cache
from math import prod,comb
import json
from pathlib import Path

def require(test,msg):
 if not test: raise RuntimeError(msg)
def factor(n):
 out=[]; p=2
 while p*p<=n:
  a=0
  while n%p==0:n//=p;a+=1
  if a:out.append((p,a))
  p+=1
 if n>1:out.append((n,1))
 return out
@lru_cache(None)
def enumerate_box(a):
 points=[]
 for v in product(*(range(-x,x+1) for x in a)):
  first=next((x for x in v if x),0)
  if first>0:points.append(v)
 count=primitive=allprim=0
 for triple in combinations(points,3):
  maxima=[max(abs(v[i]) for v in triple) for i in range(len(a))]
  prim=tuple(maxima)==a
  allprim+=prim
  balanced=all(sum(abs(v[i])==maxima[i] for v in triple)>=2 for i in range(len(a)))
  count+=balanced;primitive+=balanced and prim
 return count,primitive,comb(len(points),3),allprim
Q=[[1,18,5],[1,4,-1],[1,3],[1,1]]; R=[3,3,2,2]; C=[10,3,2,1];weights=[1,-6,-3,14]
def series_coeff(q,r,a):
 # Multiply polynomial by r independent geometric series; no binomial formula.
 coeff=q[:a+1]+[0]*max(0,a+1-len(q))
 for _ in range(r):
  for j in range(1,a+1): coeff[j]+=coeff[j-1]
 return coeff[a]
def component(n,q,r,prim):
 result=1
 for p,a in factor(n):
  if p%4!=1:
   if prim:return 0
  else: result*=series_coeff(q,r-int(prim),a)
 return result
B=5000;N=[0]*(B+1);P=N.copy();T=N.copy();PT=N.copy()
for n in range(1,B+1):
 fac=factor(n); a=tuple(sorted(a for p,a in fac if p%4==1)); en,ep,et,ept=enumerate_box(a)
 N[n]=en;T[n]=et
 if all(p%4==1 for p,a in fac): P[n]=ep;PT[n]=ept
 got=sum(w*component(n,q,r,False) for w,q,r in zip(weights,Q,R))-6
 gotp=sum(w*component(n,q,r,True) for w,q,r in zip(weights,Q,R))-6*(n==1)
 require(got==48*N[n],('total',n,got,N[n]))
 require(gotp==48*P[n],('primitive',n,gotp,P[n]))
for data,pr in [(N,P),(T,PT)]:
 accum=[0]*(B+1)
 for d in range(1,B+1):
  for n in range(d,B+1,d):accum[n]+=pr[d]
 require(accum==data,'ancestry')
# Independent local linear checks and all-triple cubic local coefficients.
for q,r,c in zip(Q,R,C):
 require(series_coeff(q,r,0)==1,'constant')
 require(series_coeff(q,r,1)==2*c+1,'total slope')
 require(series_coeff(q,r-1,1)==2*c,'primitive slope')
for a in range(101):
 require(series_coeff([1,23,23,1],4,a)==(2*a+1)**3,('cube',a))
# Primitive n=1 correction and stripping inert/ramified primes are load-bearing.
require(sum(weights)==6,'unit correction')
require(N[65]==P[65]==4,'first split pair')
require(N[130]==4 and P[130]==0,'ramified scaling')
require(N[195]==4 and P[195]==0,'inert scaling')
report={'status':'PASS','bound':B,'direct_box_count':enumerate_box.cache_info().currsize,'direct_triple_totals':{'N':sum(N),'P':sum(P),'T':sum(T),'PT':sum(PT)},'checks':['direct signed-orbit triple enumeration vs Euler coefficients','total and all-triple divisor ancestry','ramified/inert primitive stripping','unit correction','all-triple cubic local series through exponent 100'],'samples':[{ 'n':n,'N':N[n],'P':P[n],'T':T[n],'PT':PT[n]} for n in [1,5,25,65,130,195,325,845,1105,4225]]}
print(json.dumps(report,indent=2));Path(__file__).with_name('certificate.json').write_text(json.dumps(report,indent=2)+'\n')
