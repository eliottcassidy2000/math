import itertools,random,json
from pathlib import Path
# Direct permutation enumeration, independently computing response and signed contrast.
def core(k,mask):
 q=[[0]*k for _ in range(k)]
 for i,(u,v) in enumerate(itertools.combinations(range(k),2)):q[u][v]=(mask>>i)&1;q[v][u]=1-q[u][v]
 return q
def experiment(k,mask,u,v,fixedmask,perm=None):
 q=core(k,mask);n=2*k-1;remaining=[w for w in range(k) if w not in (u,v)]
 if perm is not None:remaining=list(perm)
 R=[v]
 for i in range(k-1):
  R.append(k+i)
  if i<k-2:R.append(remaining[i])
 edges=list(zip(R,R[1:]));extra=(u,R[-1]);randomized=edges+[extra]
 sets={frozenset(e) for e in randomized};fixed=[(a,b) for a in range(k) for b in range(k,n) if frozenset((a,b)) not in sets]
 if len(edges)!=2*k-3 or len(sets)!=2*k-2:raise RuntimeError('incorrect edge set')
 t=[[0]*n for _ in range(n)]
 for a in range(k):
  for b in range(k):t[a][b]=q[a][b]
 for a in range(k,n):
  for b in range(a+1,n):t[a][b]=1
 for i,(a,b) in enumerate(fixed):t[a][b]=(fixedmask>>i)&1;t[b][a]=1-t[a][b]
 paths=list(itertools.permutations(range(n)));total=0;survivors={}
 for assignment in range(1<<len(randomized)):
  for i,(a,b) in enumerate(randomized):t[a][b]=(assignment>>i)&1;t[b][a]=1-t[a][b]
  sign=(-1)**(len(edges)-sum((assignment>>i)&1 for i in range(len(edges))))
  for p in paths:
   if all(t[p[i]][p[i+1]] for i in range(n-1)):
    total+=sign;survivors[p]=survivors.get(p,0)+sign
 nonzero={p:z for p,z in survivors.items() if z}
 allowed={tuple([u]+R),tuple(R+[u]),tuple([u]+R[::-1]),tuple(R[::-1]+[u])}
 if set(nonzero)-allowed:raise RuntimeError('unexpected survivor')
 expected=2*(2*q[u][v]-1)
 if total!=expected:raise RuntimeError((k,mask,u,v,fixedmask,total,expected))
 return len(fixed)
counts={};counts[2]=0;counts[3]=0
for k in [2,3]:
 for mask in range(1<<(k*(k-1)//2)):
  for u in range(k):
   for v in range(k):
    if u==v:continue
    for fixedmask in range(1<<((k-1)*(k-2))):
     experiment(k,mask,u,v,fixedmask);counts[k]+=1
rng=random.Random(2703);counts[4]=0
for _ in range(20):
 u,v=rng.sample(range(4),2);p=[w for w in range(4) if w not in (u,v)];rng.shuffle(p)
 experiment(4,rng.randrange(64),u,v,rng.randrange(64),p);counts[4]+=1
result={'verdict':'PASS','exhaustive_orders':[2,3],'contrast_checks_by_core_order':counts,'method':'full vertex permutations for every randomized edge assignment; per-word cancellation also checked','queries_per_contrast':'4^(k-1)','expected_contrast':'2(2Q_uv-1)'}
print(json.dumps(result,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
