#!/usr/bin/env python3
"""Independent Lucas/Thue audit. Decimal proposes; integers certify.
Does not import producer and obtains Lucas numbers by a binomial closed form.
"""
from decimal import Decimal, localcontext
from itertools import permutations
from math import comb
import hashlib,json
checks=0
def check(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
def lucas(n):
 if n==0:return 2
 # 2^-n[(1+sqrt5)^n+(1-sqrt5)^n]
 numerator=2*sum(comb(n,2*j)*5**j for j in range(n//2+1))
 check(numerator%(2**n)==0,'Lucas integrality')
 return numerator//(2**n)
def determinant(A):
 total=0
 for s in permutations(range(len(A))):
  term=(-1)**sum(s[i]>s[j] for i in range(len(s)) for j in range(i+1,len(s)))
  for i,j in enumerate(s):term*=A[i][j]
  total+=term
 return total
def N(p,q):
 return determinant([[p,0,-q],[-q,p,-q],[0,-q,p-q]])
with localcontext() as ctx:
 ctx.prec=120;t=Decimal('1.84')
 for _ in range(30):t-=(t*t*t-t*t-t-1)/(3*t*t-2*t-1)
 ratio=1/(t-1)
 rows=[];sol=[]
 for n in range(101):
  r=lucas(n);q=int(Decimal(r)*ratio)
  left,right=N(q+r,q),N(q+1+r,q+1)
  check(left>0 and right<0,'exact adjacent sign certificate')
  rows.append([n,r,q,q+1,left,right])
  for denominator,value in ((q,left),(q+1,right)):
   if abs(value)==1:sol.append([n,r+denominator,denominator,value])
check(sol==[[1,2,1,1],[8,103,56,-1]],'finite intersection')
check(lucas(8)==47 and lucas(4)**2-2==47,'trace duplication')
check(N(103,56)==-1,'unit witness')
# Native element 56tau-103 has norm +1 and inverse in Z[tau].
# The multiplication table proves this without numerical approximation.
def product(a,b):
 c=[0]*5
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 for d in (4,3):
  v=c[d]
  for j in range(3):c[d-3+j]+=v
 return c[:3]
u=[1,0,0]
for _ in range(16):u=product(u,[0,1,0])
check(product(u,[-103,56,0])==[1,0,0],'tau^-16 affine unit')
check(N(-103,-56)==1,'inverse power norm')
all_solutions=[(0,-1),(1,0),(2,1),(-103,-56),(0,1),(-1,0),(-2,-1),(103,56)]
check(all(abs(N(p,q))==1 for p,q in all_solutions),'eight cited witnesses')
check(sorted((p,q) for p,q in all_solutions if q>0 and p>q)==[(2,1),(103,56)],'positive difference filter')
result={'status':'FINITE-EXACT independent audit','checks':checks,
 'solutions_n_p_q_norm':sol,'tau16_coordinates_1_tau_tau2':u,
 'row_digest':hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest()}
print(json.dumps(result,sort_keys=True,indent=2))
