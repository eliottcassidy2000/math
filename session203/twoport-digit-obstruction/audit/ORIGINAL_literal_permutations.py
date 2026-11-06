import itertools,sys,pathlib,json
from collections import Counter
PTH=pathlib.Path(__file__).resolve().parents[1]/'session203-twoport-digit-obstruction';sys.path.insert(0,str(PTH));import verify as v

def graphs(n):
 es=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(es)):
  a=[0]*n
  for k,(i,j)in enumerate(es):
   u,w=(i,j)if bits>>k&1 else(j,i);a[u]|=1<<w
  yield a

def forests(a,ports):
 n=len(a);found=set();lab={u:i for i,u in enumerate(ports)}
 for perm in itertools.permutations(range(n)):
  for cuts in range(1<<(n-1)):
   paths=[];start=0
   for i in range(n):
    if i==n-1 or cuts>>i&1:
     p=perm[start:i+1];start=i+1
     if not all(a[x]>>y&1 for x,y in zip(p,p[1:])):break
     paths.append(p)
   else:
    if sum(p[0]not in lab for p in paths)<=2:found.add(tuple(sorted(paths)))
 out=Counter()
 for f in found:out[tuple(sorted((lab.get(p[0],-1),lab.get(p[-1],-1))for p in f))]+=1
 return out

def direct_read(a,root):
 # Independent Hamilton/cover count from vertex words and one split.
 h=s=L=0;n=len(a)
 for p in itertools.permutations(range(n)):
  if all(a[x]>>y&1 for x,y in zip(p,p[1:])):h+=1;s+=p[0]==root
  cut=p.index(root)+1
  if all(a[x]>>y&1 for x,y in zip(p[:cut],p[1:cut]))and all(a[x]>>y&1 for x,y in zip(p[cut:],p[cut+1:])):L+=1
 return h+L-s,s
ct=cc=inputs=0
for n in range(2,5):
 for a in graphs(n):
  inputs+=1;st=forests(a,[0,1]);x=v.vector(st)
  if n>2:assert x[v.BULK]>=1 and not any(x[i]for i in v.BAD)
  for oi,(mask,drop,swap)in enumerate(v.OPS):
   b=a[:]+[0]
   for i in range(n):
    if i<2 and mask>>i&1:b[n]|=1<<i
    else:b[i]|=1<<n
   ports=[0,1,n];ports.pop(drop)
   if swap:ports.reverse()
   assert v.col_times(v.M[oi],x)==v.vector(forests(b,ports));ct+=1
  for ci in range(6):
   if ci<2:b=a;root=ci
   else:
    b=a[:]+[0];root=n
    for i in range(n):
     if i<2 and (ci-2)>>i&1:b[n]|=1<<i
     else:b[i]|=1<<n
   got=tuple(sum(t*y for t,y in zip(r,x))for r in v.C[ci]);assert got==direct_read(b,root);cc+=1
print(json.dumps(dict(status='PASS',method='Independent vertex permutations, all cut sets, unordered actual path forests; direct Hamilton/two-bin closing readouts',inputs=inputs,literal_transitions=ct,literal_closures=cc,scope='All labeled n2..4 seeds, all updates/closures, including empty bulk'),indent=2))
