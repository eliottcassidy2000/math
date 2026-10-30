"""Independent exact matrix audit of full sedenion spectrum and bad primes."""
from check_tower import mul, basis, norm, dot, conj, add, scale
from itertools import product
from pathlib import Path
import json, random
from fractions import Fraction
GATES=0
def need(v,s):
 global GATES
 GATES+=1
 if not v:raise RuntimeError(s)
def mm(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def determinant(a):
 a=[r[:] for r in a];n=len(a);old=1;sign=1
 for k in range(n-1):
  pivot=next((i for i in range(k,n) if a[i][k]),None)
  if pivot is None:return 0
  if pivot!=k:a[k],a[pivot]=a[pivot],a[k];sign=-sign
  q=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    z=a[i][j]*q-a[i][k]*a[k][j]
    need(z%old==0,'Bareiss exact division');a[i][j]=z//old
  for i in range(k+1,n):a[i][k]=0
  old=q
 return sign*a[-1][-1]
def rank_mod(a,p):
 a=[[v%p for v in r] for r in a];rank=0
 for j in range(len(a[0])):
  z=next((i for i in range(rank,len(a)) if a[i][j]),None)
  if z is None:continue
  a[rank],a[z]=a[z],a[rank];r=a[rank];iv=pow(r[j],-1,p);a[rank]=[v*iv%p for v in r]
  for i in range(len(a)):
   if i!=rank:
    k=a[i][j];a[i]=[(u-k*v)%p for u,v in zip(a[i],a[rank])]
  rank+=1
 return rank
rng=random.Random(116)
samples=[tuple(rng.randrange(-3,4) for _ in range(16)) for _ in range(32)]
samples += [basis(8,i)+basis(8,j) for i in range(1,8) for j in range(1,8) if i!=j]
# Non-coordinate zero divisors from rational Householder orthogonal frames.
for _ in range(12):
 v=tuple(rng.randrange(-3,4) for _ in range(7));r=norm(v)
 A=tuple(r*int(i==0)-2*v[0]*v[i] for i in range(7))
 B=tuple(r*int(i==1)-2*v[1]*v[i] for i in range(7))
 samples.append((0,)+A+(0,)+B)
samples += [(0,)*16, basis(16,0), basis(8,1)+basis(8,1), basis(8,0)+basis(8,1)]
for x in samples:
 cols=[mul(x,basis(16,i)) for i in range(16)];L=[list(r) for r in zip(*cols)]
 G=mm(cols,L);N=norm(x);a,b=x[:8],x[8:];D=norm(a[1:])*norm(b[1:])-dot(a[1:],b[1:])**2
 B=[[G[i][j]-N*int(i==j) for j in range(16)] for i in range(16)];B2=mm(B,B);B3=mm(B2,B)
 need(B3==[[4*D*v for v in r] for r in B],'Gram cubic polynomial')
 need(sum(B[i][i] for i in range(16))==0,'Gram centered trace')
 need(sum(B2[i][i] for i in range(16))==32*D,'Gram squared trace')
 det=determinant(L);need(det==N**4*(N*N-4*D)**2,'determinant formula')
 for p in (2,3,5,7,11,13):need((rank_mod(L,p)<16)==(N*(N*N-4*D)%p==0),'bad-prime criterion')
 if N and N*N==4*D:
  need(a[0]==b[0]==0 and norm(a)==norm(b) and dot(a,b)==0,'real zero-divisor cone')
  need(rank_mod(L,1000003)==12,'four-dimensional kernel sampled')
  r=norm(a);ab=mul(a,b)
  for j in range(1,8):
   z=basis(8,j)
   c=add(add(scale(r*r,z),scale(-r*dot(z,a),a)),add(scale(-r*dot(z,b),b),scale(-dot(z,ab),ab)))
   need(dot(c,a)==dot(c,b)==dot(c,ab)==c[0]==0,'quaternion perpendicular construction')
   d=scale(Fraction(1,r),mul(b,mul(c,a)))
   need(mul(x,c+d)==(0,)*16,'full kernel parameter formula')
out={'status':'PASS','gates':GATES,'matrix_cases':len(samples),'primes':[2,3,5,7,11,13],'method':'recursive16D multiplication -> Bareiss determinant, Gram cubic and traces, independent modular row reduction'}
Path(__file__).with_name('spectrum_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
