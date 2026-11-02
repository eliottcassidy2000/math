"""Exact rational independent checks for triangle + inherited-variance bounds."""
from fractions import Fraction as F
from itertools import combinations,product
import random,json,hashlib,functools

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

def gamma(n,k):return 4*(n-k+2-F(n+k+2,2**k))/F(n*(n-1))
def gap(A):
 u,v=uv(A);return all_h(A)[-1]-u-abs(v)
@functools.lru_cache(None)
def binary_gap(n,code):return gap(matrix(n,[F(code>>i&1) for i in range(n*(n-1)//2)]))
def expectation(A):
 n=len(A);edges=list(combinations(range(n),2));fixed=0;uncertain=[]
 for e,(i,j) in enumerate(edges):
  p=A[i][j]
  if p==1:fixed|=1<<e
  elif p:uncertain.append((e,p))
 value=F(0)
 for bits in range(2**len(uncertain)):
  code=fixed;prob=F(1)
  for k,(e,p) in enumerate(uncertain):
   if bits>>k&1:code|=1<<e;prob*=p
   else:prob*=1-p
  value+=prob*binary_gap(n,code)
 return value
def uncertainty_cover(A):
 n=len(A);unc=[(i,j) for i,j in combinations(range(n),2) if A[i][j] not in (0,1)]
 return min(C.bit_count() for C in range(1<<n) if all((C>>i&1) or (C>>j&1) for i,j in unc))
counts={'matrices':0,'binary_core_bounds':0,'multistar_bounds':0};digest=hashlib.sha256()
def audit(A):
 n=len(A);u,v=uv(A);c=uncertainty_cover(A)
 require(u-abs(v)>=F(1,2**c),'binary-core margin lemma');counts['binary_core_bounds']+=1
 if n>=3:
  D=all_h(A)[-1]-u-abs(v);E=expectation(A);V=sum(A[i][j]*A[j][i] for i,j in combinations(range(n),2))
  for k in range(1,n+1):require(D>=E+gamma(n,k)*V,'multistar bound');counts['multistar_bounds']+=1
  digest.update(str((n,D,E,V,c)).encode())
 counts['matrices']+=1
for n in range(1,5):
 for p in product([F(0),F(1,2),F(1)],repeat=n*(n-1)//2):audit(matrix(n,p))
rng=random.Random(138534)
for n in range(5,9):
 for trial in range(12):
  edges=list(combinations(range(n),2));p=[F(rng.randrange(2)) for e in edges]
  for e in rng.sample(range(len(edges)),min(6,len(edges))):
   den=rng.choice([3,5,7,11]);p[e]=F(rng.randrange(1,den),den)
  audit(matrix(n,p))
for n in range(3,201):
 for k in range(1,n+1):
  direct=F(2,n*(n-1))*sum(F(n-k+2*s,2**s) for s in range(k))
  require(gamma(n,k)==direct,'closed-form coefficient sum')
  if k<n:require((gamma(n,k+1)-gamma(n,k))*F(n*(n-1),4)==-1+F(n+k+1,2**(k+1)),'successive numerator difference')
 require(gamma(n,1)==F(2,n),'one-star base')
# Independent direct values of all 8*27 triangle Bernstein coefficients.
triangle_min=None;triangle_zeros=0
for signs in product([-1,1],repeat=3):
 for inds in product([0,1,2],repeat=3):
  A=matrix(3,[(1+F(s*k,2))/2 for s,k in zip(signs,inds)])
  tau=sum(triangles(A).values());active=sum(k<2 for k in inds)
  coeff=tau-F(active*(active-1),32)
  require(coeff>=0,'triangle Bernstein coefficient')
  triangle_min=coeff if triangle_min is None else min(triangle_min,coeff);triangle_zeros+=coeff==0
# Order-two is an explicit prohibited extension.
A=matrix(2,[F(1,2)]);require(gap(A)==0 and gamma(2,1)*F(1,4)>0,'n2 boundary hostile')
print(json.dumps({'status':'PASS','scope':'all ternary n1..4, twelve six-fractional-edge rational matrices each n5..8; coefficient identities n3..200','counts':counts,'binary_gap_cache':binary_gap.cache_info().currsize,'triangle_coefficients':216,'triangle_minimum':str(triangle_min),'triangle_zero_coefficients':triangle_zeros,'semantic_sha256':digest.hexdigest()},sort_keys=True,indent=2))
