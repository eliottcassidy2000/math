"""Small complete structural control, standard library only; theorem is all-order."""
from itertools import permutations
from math import comb
from collections import defaultdict
import json

def readings(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)];hs=[0]*(1<<n);hs[0]=1
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v)
    if rest:dp[mask][v]=sum(dp[rest][u] for u in range(n) if rest>>u&1 and rows[u]>>v&1)
  hs[mask]=sum(dp[mask])
 allmask=(1<<n)-1
 return hs[-1],sum(hs[m]*hs[allmask^m] for m in range(1<<n))

def components(rows):
 n=len(rows);reach=[row|(1<<i) for i,row in enumerate(rows)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 unseen=set(range(n));cs=[]
 while unseen:
  i=min(unseen);c={j for j in unseen if reach[i]>>j&1 and reach[j]>>i&1}
  cs.append(sorted(c));unseen-=c
 return cs

def induced(rows,vs):return tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(vs)) for v in vs)
def canonical(rows):
 return min(induced(rows,p) for p in permutations(range(len(rows))))

classes={};checks=0
for n in range(1,6):
 reps={}
 for bits in range(1<<(n*(n-1)//2)):
  rows=[0]*n;i=0
  for u in range(n):
   for v in range(u+1,n):
    if bits>>i&1:rows[u]|=1<<v
    else:rows[v]|=1<<u
    i+=1
  cs=components(rows);vs=[v for c in cs if len(c)>1 for v in c];core=induced(rows,vs)
  h,p=readings(rows);hc,pc=readings(core)
  if h!=hc or p!=pc*2**(n-len(core)) or (h>1 and len(core)>h):raise ArithmeticError('core')
  reps.setdefault(canonical(rows),(h,p,cs,tuple(rows)))
  checks+=1
 classes[n]=reps
cores={}
for n,reps in classes.items():
 for key,(h,p,cs,rows) in reps.items():
  if all(len(c)>1 for c in cs):cores[key]=(h,p,n,len(cs))
controls=0
for n,reps in classes.items():
 observed=defaultdict(int);predicted=defaultdict(int)
 for h,p,cs,rows in reps.values():observed[(h,p)]+=1
 predicted[(1,2**n)]=1
 for h,p,R,r in cores.values():
  if R<=n:predicted[(h,p*2**(n-R))]+=comb(n-R+r,r)
 if observed!=predicted:raise ArithmeticError(('multiplicity',n,dict(observed),dict(predicted)))
 controls+=len(observed)
print(json.dumps({'labelled_tournaments_checked':checks,'isomorphism_class_totals':{n:len(r) for n,r in classes.items()},'exact_joint_multiplicity_controls':controls,'nontrivial_core_classes_through_5':len(cores)},indent=2))
