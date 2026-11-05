from itertools import combinations,product
from fractions import Fraction
from random import Random
from math import prod,factorial
import json

def need(v,s):
 if not v:raise RuntimeError(s)
def det(a):
 a=[list(row)for row in a];n=len(a);prev=1;sgn=1
 if not n:return 1
 for k in range(n-1):
  z=next((i for i in range(k,n)if a[i][k]),None)
  if z is None:return 0
  if z!=k:a[k],a[z]=a[z],a[k];sgn=-sgn
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    v=a[i][j]*pivot-a[i][k]*a[k][j];need(v%prev==0,'exact Bareiss');a[i][j]=v//prev
  for i in range(k+1,n):a[i][k]=0
  prev=pivot
 return sgn*a[-1][-1]
def jac(f):
 m=len(f)-1;L=m+1;g=f+[0,0];a=[[0]*L for _ in range(L)]
 for i in range(L):
  r=i+1;a[i][i]=g[r-1]+3*g[r]+g[r+1]
  if i+1<L:a[i][i+1]=a[i+1][i]=g[r]+g[r+1]
 return a
rng=Random(1971401);controls=minors=selections=0;first=[]
for m in range(1,7):
 for case in range(4):
  f=[0]+[rng.randrange(1,16)for _ in range(m)];J=jac(f);L=m+1;w=[f[r]+(f[r+1]if r<m else 0)for r in range(1,m+1)];lam=[2*f[1]]+f[2:]+[0]
  # Build the incidence factor independently and sum every square minor.
  columns=[];weights=[]
  for i in range(m):
   col=[0]*L;col[i]=1;col[i+1]=-1;columns.append(col);weights.append(w[i])
  for i in range(L):
   if lam[i]:
    col=[0]*L;col[i]=1;columns.append(col);weights.append(lam[i])
  conjugated=[[(-1)**(i+j)*J[i][j]for j in range(L)]for i in range(L)]
  need(conjugated==[[sum(weights[k]*columns[k][i]*columns[k][j]for k in range(len(columns)))for j in range(L)]for i in range(L)],'incidence factor')
  CB=0;nonzero=0
  for I in combinations(range(len(columns)),L):
   v=det([[columns[j][i]for j in I]for i in range(L)]);need(v in(-1,0,1),'incidence unimodularity');CB+=v*v*prod(weights[j]for j in I);selections+=1;nonzero+=bool(v)
  need(CB==det(J)>0,'weighted Cauchy Binet sum')
  forest=0
  for mask in range(1<<m):
   comps=[[0]];weight=1
   for i in range(m):
    if mask>>i&1:comps[-1].append(i+1);weight*=w[i]
    else:comps.append([i+1])
   for roots in product(*comps):forest+=weight*prod(lam[i]for i in roots)
  need(forest==CB,'individual rooted forests, zero last root weight')
  literal=[[0]*L for _ in range(L)]
  for length in range(3,3*m+3):
   for start in range(3):
    counts=[0,0,0]
    for j in range(length):counts[(start+j)%3]+=1
    a,b,c=counts
    if min(counts)>0 and a<=L and b<=L and c<=m:literal[a-1][b-1]+=f[c]
  need(literal==J,'exact maximal run-color words')
  for k in range(1,L+1):
   for rs in combinations(range(L),k):
    for cs in combinations(range(L),k):need(det([[J[i][j]for j in cs]for i in rs])>=0,'all TN minors');minors+=1
  piv=Fraction(J[0][0]);need(piv>0,'first pivot')
  for i in range(1,L):piv=J[i][i]-Fraction(J[i][i-1]**2,piv);need(piv>0,'exact LDL positive pivot')
  controls+=1
for f in ([0,1],[0,1,2]):
 J=jac(f);K=[[factorial(i+1)*J[i][j]*factorial(j+1)for j in range(len(J))]for i in range(len(J))];first.append({'J':J,'K':K,'detJ':det(J),'detK':det(K)})
need([x['detJ']for x in first]==[2,32]and[x['detK']for x in first]==[8,4608],'factorial determinant normalization')
for h in range(7):
 for u in range(7):
  if (2*h*h+(h+2*u)**2)%7==0:need(h==0 and(h+2*u)%7==0,'complete norm residue test')
print(json.dumps({'status':'PASS','positive_coefficient_kernels':controls,'all_finite_TN_minors':minors,'incidence_square_minor_selections':selections,'complete_norm_residue_pairs':49,'first_kernels':first},sort_keys=True,indent=2))
