#!/usr/bin/env python3
"""Exact controls for compatible fibers; run from any working directory."""
import json,sys,math
from fractions import Fraction as F
from itertools import permutations,combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from explore import coefficients,packet,roots,allrows

def check(p,msg):
 if not p:raise RuntimeError(msg)

def det(a):
 a=[list(map(F,r)) for r in a];d=F(1)
 for k in range(len(a)):
  p=next((p for p in range(k,len(a)) if a[p][k]),None)
  if p is None:return 0
  if p!=k:a[k],a[p]=a[p],a[k];d=-d
  pivot=a[k][k];d*=pivot
  for j in range(k,len(a)):a[k][j]/=pivot
  for i in range(k+1,len(a)):
   c=a[i][k]
   for j in range(k,len(a)):a[i][j]-=c*a[k][j]
 return int(d) if d.denominator==1 else str(d)

def rank(a):
 a=[list(map(F,row)) for row in a];i=0
 for j in range(len(a[0]) if a else 0):
  p=next((p for p in range(i,len(a)) if a[p][j]),None)
  if p is None:continue
  a[i],a[p]=a[p],a[i];v=a[i][j];a[i]=[x/v for x in a[i]]
  for p in range(len(a)):
   if p!=i:
    v=a[p][j];a[p]=[x-v*y for x,y in zip(a[p],a[i])]
  i+=1
  if i==len(a):break
 return i

def attach(r,mask):
 return [mask<<1]+[(v<<1)|(0 if mask>>i&1 else 1) for i,v in enumerate(r)]

def direct_packet(r,brute=False):
 """No use of attachment coefficients. Literal subset paths then rooted covers."""
 n=len(r);N=1<<n;ends=[[0]*n for _ in range(N)];starts=[[0]*n for _ in range(N)]
 if brute:
  for length in range(1,n+1):
   for p in permutations(range(n),length):
    if all(r[a]>>b&1 for a,b in zip(p,p[1:])):
     S=sum(1<<i for i in p);ends[S][p[-1]]+=1;starts[S][p[0]]+=1
 else:
  for i in range(n):ends[1<<i][i]=starts[1<<i][i]=1
  for S in range(1,N):
   for i in range(n):
    if not S>>i&1:continue
    for j in range(n):
     if S>>j&1:continue
     if r[i]>>j&1:ends[S|1<<j][j]+=ends[S][i]
     if r[j]>>i&1:starts[S|1<<j][j]+=starts[S][i]
 H=[sum(row) for row in ends];H[0]=1;full=N-1;D=full^1
 U=V=d=e=0
 for S in range(1,N,2):
  c=H[full^S];sg=1 if S.bit_count()%2 else -1
  U+=starts[S][0]*c;V+=ends[S][0]*c
  d+=sg*starts[S][0]*c;e+=sg*ends[S][0]*c
 Q=sum(H[S]*H[D^S] for S in range(N) if S&D==S)
 return H[-1],H[D],Q,U,V,d,e

def transitive(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def cyclic(n,S):return [sum(1<<((i+j)%n) for j in S) for i in range(n)]

out={};count=0
for n in range(1,5):
 for code,r in allrows(n):
  seen={};points=[]
  for mask,response,H in roots(r):
   P=packet(r,[mask>>i&1 for i in range(n)])
   check(P==direct_packet(attach(r,mask),True),'literal small packet')
   if response in seen:check(seen[response]==H,'unexpected lower central switch')
   seen[response]=H;count+=1;points.append([P[j] for j in [0,3,4,5,6]])
  h,Q,u,v,a,b,B,s,t=coefficients(r)
  dim=rank([[x-y for x,y in zip(p,points[0])] for p in points[1:]])
  check(dim==(rank([u,v,a,b])+1 if n>=2 else 1),'exact affine dimension')
out['literal_all_labelled_core_tests_through4']=count
r=transitive(4);ids=[1,2,3,4,8];base=packet(r,[0]*4)
mat=[[packet(r,[mask>>i&1 for i in range(4)])[j]-base[j] for j in [0,3,4,5,6]] for mask in ids]
check(det(mat)==24,'affine hull determinant');out['affine_minor']=24
for n in range(4,13):
 r=transitive(n)
 cc=coefficients(r);h,Q,u,v,a,b,B,s,t=cc
 check(h==1 and Q==2**n,'transitive totals')
 check(u==[2**(n-1-i) for i in range(n)] and v==[2**i for i in range(n)],'dyadic linear rows')
 check(a==[0]*(n-1)+[-1] and b==[-1]+[0]*(n-1),'boundary rows')
 for i in range(n):
  for j in range(n):
   w=0 if i==j else 2 if abs(i-j)==1 else 2**(abs(i-j)-1)
   check(B[i][j]+B[j][i]==w,'weight formula')
 for l in range(n-4):
  y=[0]*n;y[l+1:l+4]=[2,-5,2]
  check(sum(u[i]*y[i] for i in range(n))==sum(v[i]*y[i] for i in range(n))==0,'kernel filter')
 K=[[([32,-25,11,-2][abs(i-j)] if abs(i-j)<4 else 0) for j in range(n-4)] for i in range(n-4)]
 for i in range(n-4):
  for j in range(n-4):
   val=sum([2,-5,2][p]*[2,-5,2][q]*(B[i+p+1][j+q+1]+B[j+q+1][i+p+1]) for p in range(3) for q in range(3))
   check(F(-val,2)==K[i][j],'quadratic restriction')
out['curvature_leading_determinants']=[det([[([32,-25,11,-2][abs(i-j)] if abs(i-j)<4 else 0) for j in range(m)] for i in range(m)]) for m in range(1,6)]
check(out['curvature_leading_determinants']==[32,399,2646,8988,-286],'leading minors')
z=[10,25,32,25,10];q=sum(z[i]*z[j]*([32,-25,11,-2][abs(i-j)] if abs(i-j)<4 else 0) for i in range(5) for j in range(5))
check(q==-2,'negative witness');out['negative_witness_value']=q
out['symmetry_controls']=[]
for n,S in [(3,[1]),(5,[1,2]),(7,[1,2,4]),(7,[1,2,3]),(9,[1,7,6,5]),(11,[1,3,4,5,9])]:
 r=cyclic(n,S);h,Q,u,v,a,b,B,s,t=coefficients(r)
 check(len(set(u))==len(set(v))==len(set(a))==len(set(b))==1,'constant symmetry response')
 check(sum(u)==Q-h and sum(a)==-h,'response sums')
 banks={};w=sorted(set(B[0][j]+B[j][0] for j in range(1,n)))
 for mask,response,H in roots(r):
  k=mask.bit_count();C=[i for i in range(n) if mask>>i&1]
  check(response==(h+k*(Q-h)//n,h+(n-k)*(Q-h)//n,h-k*h//n,k*h//n),'symmetry packet')
  check(H==h+(Q-2*h)*k//n-sum(B[i][j]+B[j][i] for i,j in combinations(C,2)),'energy formula')
  banks.setdefault(k,set()).add(H)
  if n<=7:check(packet(r,[mask>>i&1 for i in range(n)])==direct_packet(attach(r,mask),False),'literal symmetry packet')
 if len(w)==1:check(all(len(vals)==1 for vals in banks.values()),'isotropic scalar rigidity')
 out['symmetry_controls'].append(dict(n=n,S=S,h=h,Q=Q,weights=w,bank_sizes={k:len(v) for k,v in banks.items()}))
r=cyclic(9,[1,7,6,5]);pair=[]
for mask in [21,11]:
 P=direct_packet(attach(r,mask),False);check(P==packet(r,[mask>>i&1 for i in range(9)]),'independent9packet');pair.append(P)
check(pair==[(14207,3357,50004,18906,34455,2238,1119),(14209,3357,50004,18906,34455,2238,1119)],'new unit pair')
out['cyclic9_unit_pair']={'deletion_rows':r,'cuts':[21,11],'packets':pair,'gcd_with_175':math.gcd(3357,175)}
print(json.dumps(out,sort_keys=True,indent=2))
