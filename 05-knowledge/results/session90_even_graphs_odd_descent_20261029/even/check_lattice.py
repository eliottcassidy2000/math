#!/usr/bin/env python3
"""Integral partial-Walsh audit for decorated tournaments. Stdlib only."""
import json
from collections import Counter
from fractions import Fraction
from check_decorated_walsh import setup
checks=0
def check(test,what):
 global checks
 checks+=1
 if not test:raise RuntimeError(what)
def valuation(n,p=2):
 k=0
 while n%p==0:n//=p;k+=1
 return k
def factor(n):
 ans=Counter();p=2
 while p*p<=n:
  while n%p==0:ans[p]+=1;n//=p
  p+=1
 if n>1:ans[n]+=1
 return ans
def bareiss(a):
 a=[r[:] for r in a];n=len(a);sgn=1;last=1
 for k in range(n-1):
  if not a[k][k]:
   j=next((j for j in range(k+1,n) if a[j][k]),None)
   if j is None:return 0
   a[k],a[j]=a[j],a[k];sgn=-sgn
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    num=a[i][j]*pivot-a[i][k]*a[k][j]
    check(num%last==0,'nonexact Bareiss division')
    a[i][j]=num//last
  for i in range(k+1,n):a[i][k]=0
  last=pivot
 return sgn*a[-1][-1] if n else 1

def data(n,matrix=True):
 edges,perms,maps=setup(n);m=len(edges);L=1<<m;V=1<<n
 def action(z,k,tour=False):
  x,r,d=z;gm,vm,rev,c=maps[k]
  return (gm[x]^(c if tour else 0),gm[r],vm[d])
 def reps(tour):
  done=set();out=[]
  for x in range(L):
   for r in range(L):
    for d in range(V):
     z=(x,r,d)
     if z in done:continue
     orbit={action(z,k,tour) for k in range(len(perms))};done.update(orbit)
     st=[k for k in range(len(perms)) if action(z,k,tour)==z]
     if tour or all((x&maps[k][2]).bit_count()%2==0 for k in st):out.append((z,orbit,st))
  return out
 src=reps(True);dst=reps(False);r=len(src);check(r==len(dst),'dimension')
 sf=Counter();df=Counter()
 for z,orb,st in src:sf.update(factor(len(st)));check(len(st)%2==1,'source odd')
 for z,orb,st in dst:df.update(factor(len(st)))
 check({p:k for p,k in sf.items() if p!=2}=={p:k for p,k in df.items() if p!=2},'odd stabilizer product')
 twice_k=m*r-df[2];check(twice_k>=0 and twice_k%2==0,'integral determinant exponent')
 out={'n':n,'rank':r,'bits':m,'predicted_det_exponent':twice_k//2,'source_stabilizer_product_factors':dict(sf),'target_stabilizer_product_factors':dict(df)}
 if not matrix:return out
 columns=[]
 for z,orb,st in dst:
  signs={action(z,k):(-1)**((z[0]&maps[k][2]).bit_count()%2) for k in range(len(perms))}
  columns.append([sum(sign*((-1)**((x&t).bit_count()%2)) for (x,r1,d1),sign in signs.items() if (r1,d1)==(rr,dd)) for (t,rr,dd),_,_ in src])
 M=[list(row) for row in zip(*columns)];det=bareiss(M)
 check(abs(det)==1<<(twice_k//2),'determinant')
 # N=2^m M^-1 via weighted transpose. Confirm all divisions exact.
 N=[]
 for j in range(r):
  row=[]
  for i in range(r):
   num=len(src[i][1])*M[i][j];den=len(dst[j][1])
   check(num%den==0,'integer inverse numerator')
   row.append(num//den)
  N.append(row)
 for i in range(r):
  for j in range(r):
   expected=L if i==j else 0
   check(sum(M[i][k]*N[k][j] for k in range(r))==expected,'MN')
   check(sum(N[i][k]*M[k][j] for k in range(r))==expected,'NM')
 # Exact cokernel exponent from least denominator of entries of inverse.
 denexp=max((Fraction(v,L).denominator.bit_length()-1 for row in N for v in row),default=0)
 check(denexp<=m,'cokernel exponent')
 out.update({'determinant':det,'cokernel_exponent_power_of_two':denexp})
 return out
if __name__=='__main__':
 results=[data(n,n<=3) for n in range(5)]
 print(json.dumps({'results':results,'checks':checks,'status':'PASS'},indent=2))
