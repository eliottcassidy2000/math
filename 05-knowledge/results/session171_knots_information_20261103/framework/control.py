import importlib.util,itertools,json,random
from pathlib import Path
spec=importlib.util.spec_from_file_location('c',Path(__file__).with_name('check.py'));c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
rng=random.Random(171);cases=0;orders=0
for n in range(1,5):
 for trial in range(12):
  a=[[rng.choice([-2,-1,0,1,2]) for _ in range(n)] for _ in range(n)]
  want={}
  for k in range(1,n+1):
   for p in itertools.permutations(range(n),k):
    w=1
    for i,j in zip(p,p[1:]):w*=a[i][j]
    s=sum(1<<i for i in p)
    want[s]=want.get(s,0)+w
  want={s:w for s,w in want.items() if w}
  m=c.augmented(a);full=c.star(m,n+1)
  if full[0][-1]!=want:raise Exception('weighted path mismatch')
  for imask in range(1<<n):
   inds=[i+1 for i in range(n) if imask>>i&1];red,keep=c.eliminate(m,inds,n+1)
   if c.star(red,n+1)!=[[full[i][j] for j in keep] for i in keep]:raise Exception('full boundary identity')
   for order in itertools.permutations(inds):
    curr=m;labels=list(range(n+2))
    for old in order:
     index=labels.index(old);curr,k=c.eliminate(curr,[index],n+1);labels=[labels[i] for i in k]
    if curr!=red:raise Exception('all orders')
    orders+=1
   cases+=1
# Independent exhaustive response-radical finite-horizon control over F2^2.
def act(A,x):return sum(((row&x).bit_count()%2)<<i for i,row in enumerate(A))
matrices=list(itertools.product(range(4),repeat=2));obschecks=0
for A,B in itertools.product(matrices,repeat=2):
 for obs in range(4):
  blind={x for x in range(4) if (obs&x).bit_count()%2==0};history=[blind]
  while True:
   nxt={x for x in blind if act(A,x) in blind and act(B,x) in blind}
   if nxt==blind:break
   history.append(nxt);blind=nxt
  horizon=len(history)-1
  if horizon>2-bool(obs):raise Exception('horizon')
  if any(act(M,x) not in blind for M in [A,B] for x in blind):raise Exception('invariance')
  # Every linear subspace inside ker C invariant under both lies in radical.
  for S in [{0},{0,1},{0,2},{0,3},{0,1,2,3}]:
   if S<=history[0] and all(act(M,x) in S for M in [A,B] for x in S) and not S<=blind:raise Exception('universality')
  obschecks+=1
print(json.dumps({'weighted_boundary_cases':cases,'all_single_vertex_orders':orders,'exhaustive_F2_d2_observer_systems':obschecks,'status':'PASS'},indent=2))
