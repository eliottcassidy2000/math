"""Independent n=5 Bernstein audit using LOCAL FACTORS, not power coefficients.
Every tensor factor is scaled by 4. Entire numerator tensor is denominator 4^m
for Q_s=8*(H-det(I+A)+s*(-1)^n*det(A))-7*sum p(1-p).
A directed arc has degree-two Bernstein vector (2,2+eps,2+2eps)/4;
a reciprocal pair has vector (1,1,0)/4; absent edges have vector (4,4,4)/4.
"""
import numpy as np
from itertools import combinations,permutations,product
from fractions import Fraction as F
import importlib.util,json,hashlib,os
ROOT=os.path.dirname(__file__)
def require(p,message="audit condition failed"):
 if not p:raise RuntimeError(message)

n=5;E=list(combinations(range(n),2));m=len(E);EI={e:i for i,e in enumerate(E)};P=list(permutations(range(n)));D=4**m


def encode(M): return sum(int(M[i][j])<<k for k,(i,j) in enumerate(E))
def decode(o):
 M=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(E):M[i][j]=(o>>k)&1;M[j][i]=1-M[i][j]
 return M

def orbits():
 seen=set();out=[]
 for o in range(1<<m):
  if o in seen:continue
  M=decode(o);orbit={encode([[M[p[i]][p[j]] for j in range(n)] for i in range(n)]) for p in P}
  require(not (orbit&seen),'orientation orbits overlap')
  seen|=orbit;out.append((min(orbit),len(orbit)))
 require(seen==set(range(1<<m)),'orientation coverage')
 return out

def term(arcs,o,dtype):
 factors=[[4,4,4] for _ in E]
 byedge={}
 for i,j in arcs:byedge.setdefault(EI[tuple(sorted((i,j)))],[]).append(1 if i<j else -1)
 for k,dirs in byedge.items():
  if len(dirs)==2:
   require(sorted(dirs)==[-1,1],'reciprocal arc factor');factors[k]=[1,1,0]
  else:
   require(len(dirs)==1,'single arc factor')
   eps=dirs[0]*(1 if (o>>k)&1 else -1);factors[k]=[2,2+eps,2+2*eps]
 v=np.array([1],dtype=dtype)
 for f in factors:v=np.kron(np.array(f,dtype=dtype),v)
 return v

def direct(o,dtype=np.int64):
 h=np.zeros(3**m,dtype=dtype);u=h.copy();v=h.copy();q=h.copy()
 for p in P:
  h+=term(list(zip(p,p[1:])),o,dtype)
  sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
  arcs=[(i,p[i]) for i in range(n) if p[i]!=i]
  vec=term(arcs,o,dtype);u+=sign*vec
  if len(arcs)==n:v+=(-1)**n*sign*vec
 for i,j in E:q+=term([(i,j),(j,i)],o,dtype)
 return [8*(h-u+s*v)-7*q for s in [-1,1]]

def digits(i):
 out=[]
 for _ in E:out.append(int(i%3));i//=3
 return tuple(out)

def equality_faces(b):
 # An open face coordinate is 0,1,or '*' where each Bernstein factor is positive
 # at respectively index0,index2,or all three indices. If support stays all-zero,
 # nonnegative-coefficient polynomial vanishes identically on the face.
 zeros={digits(i) for i in np.flatnonzero(b==0)};faces=[]
 for face in product([0,1,2],repeat=m): # 0=z0, 1=z1, 2=free interior
  support=1
  for f in face:
   if f==2:support*=3
  if support>len(zeros):continue
  if all(tuple(t) in zeros for t in product(*[([0] if f==0 else [2] if f==1 else [0,1,2]) for f in face])):faces.append(face)
 maximal=[f for f in faces if not any(all(gk==2 or gk==fk for fk,gk in zip(f,g)) and f!=g for g in faces)]
 return zeros,maximal

def determinant(A):
 if not A:return F(1)
 B=[list(map(F,r)) for r in A];det=F(1)
 for k in range(len(B)):
  r=next((r for r in range(k,len(B)) if B[r][k]),None)
  if r is None:return F(0)
  if r!=k:B[k],B[r]=B[r],B[k];det=-det
  z=B[k][k];det*=z
  for j in range(k,len(B)):B[k][j]/=z
  for i in range(k+1,len(B)):
   z=B[i][k]
   for j in range(k,len(B)):B[i][j]-=z*B[k][j]
 return det

def evaluate_direct(z,o,s):
 A=[[F(0)]*n for _ in range(n)]
 for k,(i,j) in enumerate(E):A[i][j]=(1+(1 if o>>k&1 else -1)*z[k])/2;A[j][i]=1-A[i][j]
 h=sum((__import__('functools').reduce(lambda a,b:a*b,(A[p[i]][p[i+1]] for i in range(n-1)),F(1)) for p in P),F(0))
 u=determinant([[A[i][j]+(i==j) for j in range(n)] for i in range(n)]);v=(-1)**n*determinant(A)
 V=sum(A[i][j]*A[j][i] for i,j in E)
 return 8*(h-u+s*v)-7*V

def eval_bern(b,z):
 vals=[F(int(x),D) for x in b]
 for t in z:
  vals=[vals[i]*(1-t)**2+vals[i+1]*2*t*(1-t)+vals[i+2]*t*t for i in range(0,len(vals),3)]
 return vals[0]


O=orbits();ref={(r['orientation'],r['sign']):r for r in map(json.loads,open(os.path.join(ROOT,'..','stability','triangle5_13over16.jsonl')))};rows=[]
for o,size in O:
 bs=direct(o);q=np.zeros(3**m,dtype=np.int64);t=q.copy()
 for i,j in E:q+=term([(i,j),(j,i)],o,np.int64)
 for i,j,k in combinations(range(n),3):
  t+=term([(i,j),(j,k),(k,i)],o,np.int64)+term([(i,k),(k,j),(j,i)],o,np.int64)
 for sign,bold in zip([-1,1],bs):
  b=2*bold+q-16*t;r=ref[o,sign]
  require(int(b.min())>=0,'negative Bernstein coefficient')
  require(int(b.min())==32*r['min'] and int((b==0).sum())==r['zero'] and int(b.argmin())==r['argmin'],'producer joint receipt mismatch')
  row={'orientation':o,'orbit_size':size,'sign':sign,'minimum':int(b.min()),'zeros':int((b==0).sum()),'maximum':int(b.max()),'sha256':hashlib.sha256(b.astype('<i8').tobytes()).hexdigest()}
  if row['zeros']:
   zeros,faces=equality_faces(b);row['maximal_zero_faces']=[list(f) for f in faces]
  rows.append(row);print(json.dumps(row),flush=True)
summary={'n':n,'k':'13/16','subtracted_triangle':True,'direct_scale_for_margin':16*D,'certificates':len(rows),'orientation_coverage':sum(c for _,c in O),'rows':rows}
open(ROOT+'/audit_joint_n5_summary.json','w').write(json.dumps(summary,indent=2)+'\n')
print('PASS',json.dumps({k:v for k,v in summary.items() if k!='rows'}),flush=True)
