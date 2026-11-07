import itertools, collections, json, importlib.util, random
from pathlib import Path
P=Path(__file__).parent
spec=importlib.util.spec_from_file_location('producer',P/'frozen/probe.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
# Independent coordinate order: permutation followed by base-three interior
# state and a binary final gap. Transition reconstructed by physical boundaries.
def compiler(k,a):
 perms=list(itertools.permutations(range(k)))
 states=[s+(z,) for s in itertools.product(range(3),repeat=k) for z in (0,2)]
 keys=[(p,s) for p in perms for s in states]; ix={q:i for i,q in enumerate(keys)}
 initial=tuple(int(not any(s)) for p,s in keys)
 accept=[];mat=[[] for _ in range(1<<k)]
 for i,(p,s) in enumerate(keys):
  accept.append(all(x!=1 for x in s) and all(s[g] or a[p[g-1]][p[g]] for g in range(1,k)))
  for m in range(1<<k):
   row=[]
   for g in range(k+1):
    if not s[g] and g and not ((m>>p[g-1])&1):continue
    q=list(s);q[g]=2 if g==k or not ((m>>p[g])&1) else 1
    row.append(ix[p,tuple(q)])
   mat[m].append(row)
 def step(v,m,cap=None):
  u=[0]*len(keys)
  for i,c in enumerate(v):
   for j in mat[m][i]:u[j]+=c
  return tuple(min(x,cap) for x in u) if cap is not None else tuple(u)
 def value(v):return sum(c for c,b in zip(v,accept) if b)
 return initial,step,value,keys

def adjacency(k,a,w):
 n=k+len(w);out=[0]*n
 for i in range(n):
  for j in range(n):
   if i==j:continue
   b=a[i][j] if i<k and j<k else ((w[j-k]>>i)&1) if i<k else not((w[i-k]>>j)&1) if j<k else i<j
   if b:out[i]|=1<<j
 return out

def permutations_count(out):
 return sum(all(out[x]>>y&1 for x,y in zip(p,p[1:])) for p in itertools.permutations(range(len(out))))

def inverse(h,prune):
 v,step,val,_=compiler(2,[[0,1],[0,0]])
 todo=collections.deque([v]);seen={v};found=False;expanded=0
 while todo:
  x=todo.popleft();expanded+=1
  if val(x)==h:found=True;break
  if prune and val(x)>h:continue
  for m in range(4):
   y=step(x,m,h+1)
   if y not in seen:seen.add(y);todo.append(y)
 return {'found':found,'states':len(seen),'expanded':expanded}

rng=random.Random(216);controls=0
for k in range(5):
 arcs=[]
 if k<=3:
  for bits in itertools.product((0,1),repeat=k*(k-1)//2):
   a=[[0]*k for _ in range(k)];b=iter(bits)
   for i in range(k):
    for j in range(i+1,k):a[i][j]=next(b);a[j][i]=1-a[i][j]
   arcs.append(a)
 else:
  for _ in range(3):
   a=[[0]*k for _ in range(k)]
   for i in range(k):
    for j in range(i+1,k):a[i][j]=rng.randrange(2);a[j][i]=1-a[i][j]
   arcs.append(a)
 for a in arcs:
  v,step,val,keys=compiler(k,a); pv,pstep,pval=mod.make(k,a)
  words=[()]
  for length in range(1,min(4,8-k)):
   words += list(itertools.product(range(1<<k),repeat=length)) if k<=2 else [tuple(rng.randrange(1<<k) for _ in range(length)) for z in range(12)]
  for w in words:
   x=v;y=pv
   for m in w:x=step(x,m);y=pstep(y,m)
   expected=permutations_count(adjacency(k,a,w))
   if not val(x)==pval(y)==expected:raise RuntimeError((k,a,w,val(x),pval(y),expected))
   for h in (0,1,3,7,21):
    z=v
    for m in w:z=step(z,m,h+1)
    if z!=tuple(min(c,h+1) for c in x):raise RuntimeError('saturation')
    if (val(z)==h)!=(expected==h):raise RuntimeError('readout')
   controls+=1
results={str(h):inverse(h,True) for h in (0,1,3,5,7,9,21,23)}
if results['7']['states']!=708 or results['21']['states']!=4626:raise RuntimeError(results)
# No-pruning replay ensures hidden large-state deletion is harmless.
results['7_unpruned']=inverse(7,False)
if results['7_unpruned']['found']:raise RuntimeError('unexpected unpruned witness')
# Validate domain and consistency of both saved witness formats.
data=json.loads((P/'frozen/witnesses.json').read_text()); text={}
for line in (P/'frozen/witnesses.txt').read_text().splitlines():
 h,w=line.split();text[h]=[] if w=='-' else list(map(int,w))
if data!=text:raise RuntimeError('witness format mismatch')
if set(map(int,data))!=set(range(1,10001,2))-{7,21}:raise RuntimeError('witness coverage')
if any(len(w)>10 or any(m not in range(4) for m in w) for w in data.values()):raise RuntimeError('witness domain')
# Re-run producer resource and boundary checks; ordinary negative inputs outside domain.
spec=importlib.util.spec_from_file_location('inverse',P/'frozen/inverse.py');inv=importlib.util.module_from_spec(spec)
import sys;sys.path.insert(0,str(P/'frozen'));spec.loader.exec_module(inv)
resource={str(h):inv.compile_target(h,1) for h in (7,21,23)}
if any(x['status']!='INCOMPLETE' for x in resource.values()):raise RuntimeError('resource false absence')
print(json.dumps({'permutation_controls':controls,'saturation_caps':[0,1,3,7,21],'inverse':results,'resource':resource,'witnesses':len(data),'max_vertices':2+max(map(len,data.values()))},indent=2))
