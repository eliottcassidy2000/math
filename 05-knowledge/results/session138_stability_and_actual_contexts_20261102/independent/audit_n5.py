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
spec=importlib.util.spec_from_file_location('producer',os.path.join(ROOT,'..','stability','bernstein.py')); mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

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

O=orbits();require(len(O)==12,'twelve orientation classes')
h,u,v,V=mod.polys(n);rows=[]
for o,size in O:
 bs=direct(o)
 for s,b in zip([-1,1],bs):
  bp=mod.bernstein(8*(h-u+s*v)-7*V,m,o)
  require(np.array_equal(b,32*bp),'producer coefficient array mismatch')
  require(int(b.min())>=0,'negative Bernstein coefficient')
  rows.append(dict(orientation=o,orbit_size=size,sign=s,minimum=int(b.min()),zeros=int((b==0).sum()),maximum=int(b.max()),sha256=hashlib.sha256(b.astype('<i8').tobytes()).hexdigest()))
  print(json.dumps(rows[-1]),flush=True)
 if o==0:
  for s,b in zip([-1,1],bs):
   zeros,faces=equality_faces(b)
   print('EQUALITY',s,json.dumps({'zero_indices':[list(x) for x in sorted(zeros)],'maximal_faces':[list(x) for x in faces]}),flush=True)
  require(all(np.array_equal(x,y) for x,y in zip(bs,direct(o,object))),'arbitrary-precision replay mismatch')
 # Independent Gaussian elimination Fraction evaluations at nonbinary points.
 z=[F(((i*5+o)%11),11) for i in range(m)]
 for s,b in zip([-1,1],bs):require(eval_bern(b,z)==evaluate_direct(z,o,s),'Fraction point check')
# Absolute bound accounts for every unsimplified local tensor term and final sum.
B=8*3*len(P)*D+7*m*D
require(B<2**63,'int64 bound')
summary={'n':n,'m':m,'orbits':O,'orientation_coverage':sum(c for _,c in O),'certificates':len(rows),'scale_for_Q':D,'scale_for_margin':8*D,'conservative_int64_bound':B,'all_rows_match_producer':True,'transitive_object_dtype_match':True,'fraction_evaluations':2*len(O),'rows':rows}
open(ROOT+'/audit_n5_summary.json','w').write(json.dumps(summary,indent=2)+'\n')
print('PASS',json.dumps({k:v for k,v in summary.items() if k!='rows'}),flush=True)
