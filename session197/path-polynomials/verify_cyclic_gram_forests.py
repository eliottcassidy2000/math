#!/usr/bin/env python3
from verify import *
from math import factorial,prod
from random import Random

def ordered(a,L):return [(cover_polynomial(a)+[0]*L)[r]*factorial(r) for r in range(1,L+1)]
def jacobi(f):
 m=len(f)-1;L=m+1;f=f+[0,0];J=[[0]*L for _ in range(L)]
 for r in range(1,L+1):
  J[r-1][r-1]=f[r-1]+3*f[r]+f[r+1]
  if r<L:J[r-1][r]=J[r][r-1]=f[r]+f[r+1]
 return J
def form(A,x,y):return sum(x[i]*A[i][j]*y[j] for i in range(len(x))for j in range(len(y)))
def forest(f):
 m=len(f)-1;w=[f[r]+(f[r+1]if r<m else 0)for r in range(1,m+1)];lam=[2*f[1]]+f[2:]+[0];out=0
 for mask in range(1<<m):
  weight=1;root_sum=lam[0]
  for i in range(m):
   if mask>>i&1:weight*=w[i];root_sum+=lam[i+1]
   else:weight*=root_sum;root_sum=lam[i+1]
  out+=weight*root_sum
 return out

def run():
 kernels=forests=minorcount=0;rng=Random(1971349)
 for m in range(1,7):
  for _ in range(8):
   f=[0]+[rng.randrange(1,12)for r in range(m)];J=jacobi(f);L=m+1
   literal=[[0]*L for _ in range(L)]
   for i,j,k in product(range(1,L+1),range(1,L+1),range(1,m+1)):
    N=sum(tuple(tuple((s+d)%3 for d in range(i+j+k)).count(c)for c in range(3))==(i,j,k) for s in range(3))
    literal[i-1][j-1]+=N*f[k]
   require(literal==J,'literal cyclic run Jacobi matrix');kernels+=1
   x=[rng.randrange(-4,5)for _ in range(L)];sq=2*f[1]*x[0]**2+sum(f[r]*x[r-1]**2 for r in range(2,m+1))+sum((f[r]+(f[r+1] if r<m else 0))*(x[r-1]+x[r])**2 for r in range(1,m+1))
   require(form(J,x,x)==sq,'exact weighted square decomposition')
   require(det(J)==forest(f)>0,'rooted path forest determinant');forests+=1
   for k in range(1,min(3,L)+1):
    for rs in combinations(range(L),k):
     for cs in combinations(range(L),k):require(det([[J[i][j]for j in cs]for i in rs])>=0,'Jacobi total nonnegativity');minorcount+=1
 inputs=[tournament(n,c)for n in range(1,4)for c in range(1<<(n*(n-1)//2))];physical=cschecks=0
 for V in inputs:
  m=len(V);f=[0]+ordered(V,m);J=jacobi(f);L=m+1
  require(det(J)==forest(f),'physical cover-weight forest identity')
  for T,U in product(inputs,repeat=2):
   x=ordered(T,L);y=ordered(U,L);h=h_subsets(cyclic(T,U,V))[-1]
   require(h==form(J,x,y),'actual three-input graph/kernel')
   lhs=h*h;rhs=form(J,x,x)*form(J,y,y);require(lhs<=rhs,'physical Cauchy-Schwarz')
   proportional=all(x[i]*y[0]==y[i]*x[0]for i in range(L));require((lhs==rhs)==proportional,'exact Gram equality scope')
   physical+=1;cschecks+=1
  bank=[ordered(transitive(n),L)for n in range(1,L+1)];gram=[[form(J,x,y)for y in bank]for x in bank]
  require(det(gram)==det(J)*prod(factorial(r)**2 for r in range(1,L+1)),'full physical transitive Gram determinant')
 J1=jacobi([0,1]);J2=jacobi([0,1,2]);require(J1==[[3,1],[1,1]]and det(J1)==2,'singleton Jacobi kernel')
 K1=[[J1[i][j]*factorial(i+1)*factorial(j+1)for j in range(2)]for i in range(2)]
 K2=[[J2[i][j]*factorial(i+1)*factorial(j+1)for j in range(3)]for i in range(3)]
 require(K1==[[3,2],[2,4]] and K2==[[5,6,0],[6,28,24],[0,24,72]]and det(J2)==32 and det(K2)==4608,'exact first two unordered kernels')
 need_norm=0
 for h in range(1,101,2):
  for u in range(101):
   H=2*h*h+(h+2*u)**2
   if H%7==0:require(h%7==0 and(h+2*u)%7==0 and H%49==0 and H>49,'norm obstruction including exceptional49')
   need_norm+=1
 print(json.dumps({'status':'PASS','literal_color_kernels':kernels,'positive_weight_forest_controls':forests,'exact_TN_minors':minorcount,'actual_triple_graphs':physical,'physical_Gram_inequalities':cschecks,'norm_controls':need_norm,'first_Jacobi_kernels':[J1,J2],'first_unordered_kernels':[K1,K2],'determinants':[det(J1),det(J2),det(K1),det(K2)]},sort_keys=True,indent=2))
if __name__=='__main__':run()
