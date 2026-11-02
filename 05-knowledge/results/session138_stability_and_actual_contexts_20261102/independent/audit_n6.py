"""Independent n6 sharp-margin certificate: subset-DP polynomial assembly,
blockwise Bernstein transform, independently generated orientation orbits.
No producer code is imported. Only the producer's JSON receipt is compared.
"""
from itertools import combinations,permutations
import numpy as np,json,math,time,hashlib,os
N=6;E=list(combinations(range(N),2));M=len(E);EI={p:k for k,p in enumerate(E)};PW=[3**k for k in range(M)];SIZE=3**M
ROOT=os.path.dirname(__file__)
def require(p,message="audit condition failed"):
 if not p:raise RuntimeError(message)

def add(dst,src,mul=1):
 for e,c in src.items():
  nc=dst.get(e,0)+mul*c
  if nc:dst[e]=nc
  elif e in dst:del dst[e]
def arc(p,i,j):
 k=EI[tuple(sorted((i,j)))];s=1 if i<j else -1;out=p.copy()
 for e,c in p.items():
  ee=e+PW[k];nc=out.get(ee,0)+s*c
  if nc:out[ee]=nc
  elif ee in out:del out[ee]
 return out

def hp():
 states={(1<<i,i):{0:1} for i in range(N)}
 for mask in range(1,1<<N):
  for j in range(N):
   cur=states.get((mask,j))
   if cur is None:continue
   for k in range(N):
    if mask>>k&1:continue
    dest=states.setdefault((mask|(1<<k),k),{})
    add(dest,arc(cur,j,k))
 out={}
 for j in range(N):add(out,states[(2**N-1,j)],2)
 return out

def detpoly(identity):
 states={0:{0:1}}
 for mask in range(1<<N):
  if mask not in states:continue
  i=mask.bit_count()
  if i==N:continue
  cur=states[mask]
  for j in range(N):
   if mask>>j&1:continue
   if i==j:
    if not identity:continue
    pp={e:2*c for e,c in cur.items()}
   else:pp=arc(cur,i,j)
   s=(-1)**((mask>>(j+1)).bit_count())
   add(states.setdefault(mask|(1<<j),{}),pp,s)
 return states[2**N-1]

def orbits():
 P=list(permutations(range(N)));unseen=set(range(1<<M));reps=[]
 while unseen:
  o=min(unseen);a=[[0]*N for _ in range(N)]
  for k,(i,j) in enumerate(E):a[i][j]=(o>>k)&1;a[j][i]=1-a[i][j]
  orb=set()
  for p in P:orb.add(sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(E)))
  require(orb<=unseen,'orientation orbit overlap')
  unseen-=orb;reps.append((o,len(orb)))
 require(len(reps)==56 and sum(x[1] for x in reps)==32768,'56 classes cover32768 orientations')
 return reps

def bern(poly,o):
 b=np.zeros(SIZE,dtype=np.int64)
 for e,c in poly.items():
  q=e;par=0
  for k in range(M):
   if q%3==1 and not(o>>k&1):par^=1
   q//=3
  b[e]=(-1 if par else 1)*c
 stride=1
 for k in range(M):
  v=b.reshape((-1,3,stride));x=v[:,0,:].copy();y=v[:,1,:].copy()
  v[:,2,:]=2*(x+y+v[:,2,:]);v[:,1,:]=2*x+y;v[:,0,:]=2*x
  stride*=3
 return b

def digits(i):
 out=[]
 for _ in E:out.append(int(i%3));i//=3
 return out

start=time.time();h=hp();u=detpoly(True);v=detpoly(False);V={0:M*2**(N-2)}
for p in PW:V[2*p]=-2**(N-2)
print(json.dumps({'stage':'assembly','seconds':time.time()-start,'support_sizes':[len(h),len(u),len(v),len(V)],'sums_abs':[sum(map(abs,p.values())) for p in [h,u,v,V]]}),flush=True)
O=orbits();bound=2**M*(4*3*2**N*math.factorial(N)+3*M*2**(N-1));require(bound<2**63,'int64 bound')
reference={(r['orientation'],r['sign']):r for r in map(json.loads,open(os.path.join(ROOT,'..','stability','bernstein6.jsonl')))};require(len(reference)==112,'112 reference rows')
rows=[]
for o,orbsize in O:
 for s in [-1,1]:
  q={};add(q,h,4);add(q,u,-4);add(q,v,4*s*(-1)**N);add(q,V,-3)
  b=bern(q,o);row={'orientation':o,'orbit_size':orbsize,'sign':s,'min':int(b.min()),'zero':int((b==0).sum()),'max':int(b.max()),'argmin':int(b.argmin()),'sha256':hashlib.sha256(b.astype('<i8').tobytes()).hexdigest()}
  require(row['min']>=0,'negative Bernstein coefficient')
  ref=reference[o,s]
  for key in ['min','zero','argmin']:require(row[key]==ref[key],(key,row,ref))
  if o==0:row['zero_indices']=[digits(i) for i in np.flatnonzero(b==0)]
  rows.append(row);print(json.dumps(row),flush=True)
summary={'n':N,'m':M,'scale_for_margin':4*2**(N+M),'conservative_all_intermediate_bound':bound,'orbits':O,'certificates':len(rows),'orientation_coverage':sum(x[1] for x in O),'seconds':time.time()-start,'rows':rows}
open(ROOT+'/audit_n6_summary.json','w').write(json.dumps(summary,indent=2)+'\n');print('PASS',json.dumps({k:v for k,v in summary.items() if k!='rows'}),flush=True)
