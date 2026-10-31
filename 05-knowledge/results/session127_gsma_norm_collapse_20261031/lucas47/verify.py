#!/usr/bin/env python3
"""Session127 Lucas47 candidate junction. Exact arithmetic; no disabled asserts.
Only the Lucas-index interval 0..100 is searched. All positive q are covered
per index by the monotonicity proof in README.md, not by a denominator cap.
"""
import hashlib, json
from pathlib import Path

checks=0
def require(ok,label):
 global checks
 checks+=1
 if not ok: raise RuntimeError(label)
def mul(A,B):
 return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def power(A,n):
 B=[[int(i==j) for j in range(len(A))] for i in range(len(A))]
 while n:
  if n&1:B=mul(B,A)
  A=mul(A,A);n//=2
 return B
def det3(A):
 a,b,c=A[0];d,e,f=A[1];g,h,i=A[2]
 return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
def norm3(p,q):return p**3-p*p*q-p*q*q-q**3
def shear(r,q):return r**3+2*r*r*q-2*q**3
Q=[[0,1],[1,1]];Q8=power(Q,8)
require(Q8==[[13,21],[21,34]],'Q8')
require(sum(Q8[i][i] for i in range(2))==47,'Lucas trace')
require(Q8[0][0]*Q8[1][1]-Q8[0][1]*Q8[1][0]==1,'unit norm')
# Deliberately explicit guard: trace47 is not the unit-minus-one norm45.
require((Q8[0][0]-1)*(Q8[1][1]-1)-Q8[0][1]*Q8[1][0]==-45,'unit-minus-one norm')
require(47**2-5*21**2==4,'Pell')
require(2**7-3**4==47,'Pillai gap')
require(103-56==47,'convergent difference')
require(norm3(103,56)==-1,'Tribonacci norm')
require(56**4-9*103**3==-47,'collapsed residual')
M=[[0,0,1],[1,0,1],[0,1,1]]
for p in range(-8,9):
 for q in range(-8,9):
  require(norm3(p,q)==det3([[p*int(i==j)-q*M[i][j] for j in range(3)] for i in range(3)]),'determinant norm')
  require(norm3(p,q)==shear(p-q,q),'unimodular shear')
  require((p-q)*norm3(p,q)==q**4-(2*q-p)*p**3,'quartic collapse')
a,b=2,1;rows=[];solutions=[]
for n in range(101):
 r=a;a,b=b,a+b
 lo,hi=r,2*r
 require(shear(r,lo)==r**3 and shear(r,hi)==-11*r**3,'bracket')
 while hi-lo>1:
  mid=(lo+hi)//2
  if shear(r,mid)>0:lo=mid
  else:hi=mid
 require(hi==lo+1 and shear(r,lo)>0 and shear(r,hi)<0,'strict root bracket')
 row=[n,r,lo,hi,shear(r,lo),shear(r,hi)];rows.append(row)
 for q in (lo,hi):
  value=norm3(r+q,q)
  if abs(value)==1:solutions.append([n,r+q,q,value])
require(solutions==[[1,2,1,1],[8,103,56,-1]],'bounded Lucas/Thue intersection')
# Classical complete-solution list: the exact witnesses are checked here;
# completeness is the CITED Delone-Nagell theorem, not a finite-search claim.
plus=[(0,-1),(1,0),(2,1),(-103,-56)]
require(all(norm3(p,q)==1 for p,q in plus),'four global plus-one witnesses')
minus=[(-p,-q) for p,q in plus]
require(all(norm3(p,q)==-1 for p,q in minus),'four global minus-one witnesses')
filtered=sorted((p,q,norm3(p,q)) for p,q in plus+minus if q>0 and p-q>0)
require(filtered==[(2,1,1),(103,56,-1)],'global classification filter')
# Hostiles: adjacent Lucas indices and the next doubling trace fail.
require(all(n not in [7,9,16] for n,_,_,_ in solutions),'neighbor and doubling hostiles')
result={'status':'FINITE-EXACT','checks':checks,'lucas_index_range':[0,100],
 'positive_denominators':'all, by monotone cubic root isolation per index',
 'solutions_n_p_q_norm':solutions,'cited_global_filter_p_q_norm':filtered,'Q8':Q8,'trace':47,'norm_phi8_minus_one':-45,
 'row_digest':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest(),
 'rows_n_r_floor_ceil_gfloor_gceil':rows}
print(json.dumps(result,sort_keys=True,indent=2))
