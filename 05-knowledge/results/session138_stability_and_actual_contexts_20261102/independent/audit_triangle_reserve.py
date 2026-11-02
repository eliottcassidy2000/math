"""Exact rational independent checks for triangle + inherited-variance bounds."""
from fractions import Fraction as F
from itertools import combinations,product
import random,json,hashlib

def require(p,message):
 if not p:raise RuntimeError(message)
def matrix(n,p):
 A=[[F(0) for _ in range(n)] for _ in range(n)]
 for (i,j),x in zip(combinations(range(n),2),p):A[i][j]=x;A[j][i]=1-x
 return A
def det(A):
 B=[list(map(F,r)) for r in A];out=F(1)
 for k in range(len(B)):
  r=next((r for r in range(k,len(B)) if B[r][k]),None)
  if r is None:return F(0)
  if r!=k:B[k],B[r]=B[r],B[k];out=-out
  x=B[k][k];out*=x
  for j in range(k,len(B)):B[k][j]/=x
  for i in range(k+1,len(B)):
   x=B[i][k]
   for j in range(k,len(B)):B[i][j]-=x*B[k][j]
 return out
def uv(A):
 n=len(A);return det([[A[i][j]+(i==j) for j in range(n)] for i in range(n)]),(-1)**n*det(A)
def all_h(A):
 n=len(A);d=[[F(0)]*n for _ in range(2**n)];h=[F(0)]*2**n;h[0]=1
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,2**n):
  h[mask]=sum(d[mask])
  for i in range(n):
   if d[mask][i]:
    for j in range(n):
     if not mask>>j&1:d[mask|(1<<j)][j]+=d[mask][i]*A[i][j]
 return h
def triangles(A):
 return {(i,j,k):A[i][j]*A[j][k]*A[k][i]+A[i][k]*A[k][j]*A[j][i] for i,j,k in combinations(range(len(A)),3)}
def pow2(e):return F(2**e) if e>=0 else F(1,2**(-e))
def c(n):return (24-pow2(3-n)*(n*n+2*n+3))/(n*(n-1))
def beta(n):return (24-pow2(4-n)*(n*n-n+3))/(n*(n-1))
counts={'matrices':0,'weak_roots':0,'strong_roots':0,'strong_matrices':0};digest=hashlib.sha256()
def audit(A):
 n=len(A);h=all_h(A);u,v=uv(A);delta=h[-1]-u;gap=delta-abs(v);tris=triangles(A);tau=sum(tris.values(),F(0));V=sum(A[i][j]*A[j][i] for i,j in combinations(range(n),2))
 if n>=2:
  require(delta>=tau+c(n)*V,'weak global triangle reserve')
  for x in range(n):
   vs=[i for i in range(n) if i!=x];B=[[A[i][j] for j in vs] for i in vs];up,vp=uv(B);dparent=h[(2**n-1)^(1<<x)]-up
   incident=sum(A[x][i]*A[i][x] for i in vs);through=sum(w for triple,w in tris.items() if x in triple)
   require(delta>=dparent+(n-1)*pow2(2-n)*incident+through,'weak root triangle reserve');counts['weak_roots']+=1
   if n>=4:
    require(gap>=dparent+pow2(3-n)*incident+through,'strong root triangle reserve');counts['strong_roots']+=1
 if n>=4:
  require(gap>=tau+beta(n)*V,'strong global triangle reserve');counts['strong_matrices']+=1
 counts['matrices']+=1;digest.update(str((n,delta,gap,tau,V)).encode())
for n in range(1,5):
 for p in product([F(0),F(1,2),F(1)],repeat=n*(n-1)//2):audit(matrix(n,p))
rng=random.Random(138)
for n in range(5,9):
 for _ in range(32):
  p=[]
  for e in combinations(range(n),2):
   den=rng.choice([2,3,5,7,11]);p.append(F(rng.randrange(den+1),den))
  audit(matrix(n,p))
A=matrix(3,[F(1),F(0),F(1)]);u,v=uv(A);require(all_h(A)[-1]-u-abs(v)==0 and sum(triangles(A).values())==1,'C3 excludes n3 strong extension')
# n5 stronger sharp 7/8 joint reserve is false: fair pair ordinal fair triple.
A=matrix(5,[F(1,2) if (i<2 and j<2) or (i>=2 and j>=2) else F(1) for i,j in combinations(range(5),2)])
h=all_h(A)[-1];u,v=uv(A);tau=sum(triangles(A).values());V=sum(A[i][j]*A[j][i] for i,j in combinations(range(5),2));gap=h-u-abs(v)
require((h,u,v,tau,V,gap-tau)==(F(3,2),F(3,8),F(1,16),F(1,4),F(1),F(13,16)),'sharp+triangle n5 hostile coordinates')
print(json.dumps({'scope':'all ternary matrices n1..4 and32 fixed-seed rational matrices each n5..8','counts':counts,'semantic_sha256':digest.hexdigest(),'n3_strong_hostile':'C3:Delta0,tau1','n5_7_over_8_joint_hostile':{'H':str(h),'u':str(u),'v':str(v),'tau':str(tau),'V':str(V),'Delta_minus_tau':str(gap-tau)},'verdict':'PASS'},indent=2,sort_keys=True))
