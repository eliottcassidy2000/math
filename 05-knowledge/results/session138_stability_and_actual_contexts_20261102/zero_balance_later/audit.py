"""Independent fixed-module profiles and finite identities for an all-a,b proof.

This implementation imports no producer/repository counter. It computes all
induced H values by backward endpoint DP, then unordered path covers by an
anchored set-partition DP, multiplying by r! for ordered covers.
"""
import json,math,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent

def covers(a):
 n=len(a);N=1<<n;E=[[0]*n for _ in range(N)];hs=[0]*N;hs[0]=1
 for S in range(1,N):
  for i in range(n):
   if not(S>>i&1):continue
   R=S^(1<<i)
   E[S][i]=sum(E[R][j] for j in range(n) if R>>j&1 and a[j]>>i&1) if R else 1
  hs[S]=sum(E[S])
 G=[[0]*(n+1) for _ in range(N)];G[0][0]=1
 for S in range(1,N):
  z=S&-S;R=S^z;U=R
  while True:
   A=U|z;B=S^A;weight=hs[A]
   for r in range(1,S.bit_count()+1):G[S][r]+=weight*G[B][r-1]
   if U==0:break
   U=(U-1)&R
 return [G[-1][r]*math.factorial(r) for r in range(1,n+1)]

def matrix(d):
 dd=[0]+d+[0,0];N=len(d)+1;K=[[0]*(N+1) for _ in range(N+1)]
 for s in range(1,N+1):
  K[s][s]=dd[s-1]+3*dd[s]+dd[s+1]
  if s<N:K[s][s+1]=K[s+1][s]=dd[s]+dd[s+1]
 return K

def f(a,r):return sum((-1)**(r-j)*math.comb(r,j)*(j*j*j+2*j)**a for j in range(r+1))
def bilinear(K,A,B):return sum(K[i][j]*A[i]*B[j] for i in range(1,len(K)) for j in range(1,len(K)))
data=json.loads((P/'mixed_profile_inputs.json').read_text())
for a in [data['original_masks']]+[row['child_masks'] for row in data['rows']]:
 if len(a)!=12 or any(a[i]>>i&1 for i in range(12)) or any(((a[i]>>j)&1)+((a[j]>>i)&1)!=1 for i in range(12) for j in range(i+1,12)):raise RuntimeError('literal tournament typing')
orig=covers(data['original_masks'])
for row in data['rows']:
 if row['original_H']!=orig[0] or row['child_masks'][-1]!=row['cut']:raise RuntimeError('marked module identity')
 for i in range(11):
  for j in range(11):
   if ((row['child_masks'][i]>>j)&1)!=((data['original_masks'][i+1]>>(j+1))&1):raise RuntimeError('same literal deletion card')
base=[0,9,126,684,1656,1800,720]+[0]*7
if [f(2,r) for r in range(14)]!=base:raise RuntimeError('base profile')
log_slacks=[base[r]**2-base[r-1]*base[r+1] for r in range(2,6)]
if min(log_slacks)<0:raise RuntimeError('base log concavity')
rows=[]
for source in data['rows']:
 cut=source['cut'];got=covers(source['child_masks']);d=[b-a for a,b in zip(orig,got)]
 if d!=source['difference']:raise RuntimeError(('independent profile',cut,d))
 K=matrix(d);N=len(K)-1
 if cut==996:
  if any(K[i][j]<0 for i in range(1,N+1) for j in range(1,N+1)) or K[1][1]!=135618:raise RuntimeError('positive kernel')
  boundary=None
 else:
  boundary={867:4,995:5}[cut]
  if any(K[i][j]<0 for i in range(1,boundary) for j in range(1,N+1)):raise RuntimeError('positive prefix rows')
  if any(K[i][j]>0 for i in range(boundary+1,N+1) for j in range(1,N+1)):raise RuntimeError('negative suffix rows')
  if not (K[boundary][boundary-1]>0 and K[boundary][boundary]<0 and K[boundary][boundary+1]<0):raise RuntimeError('single mixed row')
  if min(d[:4])<=0:raise RuntimeError('thin-context positivity')
  if bilinear(K,base,base)>=0:raise RuntimeError('negative base')
 # Bounded algebra/control checks, not the all-order argument.
 for a in range(2,13):
  F=[f(a,r) for r in range(14)]
  if any(F[i]*base[j]>F[j]*base[i] for i in range(1,7) for j in range(i+1,7)):raise RuntimeError('MLR control')
 row={'cut':cut,'independent_profile_difference':d,'tridiagonal_diagonal':[K[s][s] for s in range(1,N+1)],'tridiagonal_upper':[K[s][s+1] for s in range(1,N)],'one_crossing_boundary':boundary,'gap_2_2':bilinear(K,base,base),'gap_1_1':bilinear(K,[f(1,r) for r in range(14)],[f(1,r) for r in range(14)])};rows.append(row);print(json.dumps(row,sort_keys=True))
report={'original_profile':orig,'base_profile':base[1:7],'log_concavity_slacks':log_slacks,'rows':rows,'all_fixed_module_profiles_independently_match':True}
print('PASS '+json.dumps({'source_sha256':data['source_sha256'],'semantic_sha256':hashlib.sha256(json.dumps(report,separators=(',',':')).encode()).hexdigest()},sort_keys=True))
