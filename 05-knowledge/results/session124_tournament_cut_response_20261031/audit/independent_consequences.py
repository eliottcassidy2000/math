"""Independent structured and seeded packing/small-gain checks."""
from pathlib import Path
from itertools import combinations
from functools import lru_cache
import json,hashlib,random,runpy
D=runpy.run_path(str(Path(__file__).with_name('independent_check.py')))
graph,hamiltonian,compile_response,child=(D[x] for x in ('graph','hamiltonian','compile_response','child'))
def triangles(A):return [sum(1<<u for u in tri) for tri in combinations(range(len(A)),3) if all(sum(A[u][v] for v in tri)==1 for u in tri)]
def packing(A):
 ts=triangles(A)
 @lru_cache(None)
 def go(mask):
  return max([0]+[1+go(mask^t) for t in ts if mask&t==t])
 return go((1<<len(A))-1)
def ordinal(t):
 return [[(i//3<j//3) or (i//3==j//3 and (j-i)%3==1) for j in range(3*t)] for i in range(3*t)]
rng=random.Random(124034);parents=[graph(n,rng.randrange(1<<(n*(n-1)//2))) for n in range(6,10) for _ in range(8)]+[ordinal(3),ordinal(4)]
checks=dict(parents=0,pair_packing=0,positive_gain=0,crossing_deletions=0,small_gain_transitive=0,selected_child_dp=0);controls=[];digest=hashlib.sha256()
for A in parents:
 n=len(A);nu=packing(A);out=tuple(sum(1<<v for v in range(n) if A[u][v]) for u in range(n));r=compile_response(out)
 if r.h!=hamiltonian(A) or r.h<3**nu:raise RuntimeError('parent packing')
 checks['parents']+=1;pair={}
 for u,v in combinations(range(n),2):
  left=[x for x in range(n) if x not in (u,v)];B=[[A[x][y] for y in left] for x in left];h=hamiltonian(B)
  if h<3**max(nu-2,0):raise RuntimeError('pair packing')
  pair[u,v]=pair[v,u]=(h,not triangles(B));checks['pair_packing']+=1
 positive=[]
 for S,F in enumerate(r.all_values()):
  gain=F-r.h
  if not gain:continue
  if gain<2*3**max(nu-2,0):raise RuntimeError('gain packing')
  positive.append(gain);checks['positive_gain']+=1
  for u in range(n):
   for v in range(n):
    if S>>u&1 and not S>>v&1 and A[u][v]:
     h,transitive=pair[u,v]
     if 2*h>gain:raise RuntimeError('crossing deletion')
     checks['crossing_deletions']+=1
     if gain<=4:
      if h!=1 or not transitive:raise RuntimeError('small gain transitive')
      checks['small_gain_transitive']+=1
 for S in (1,((1<<n)-1)//3):
  if r.value(S)!=hamiltonian(child(A,S)):raise RuntimeError('selected actual child')
  checks['selected_child_dp']+=1
 if n in (9,12) and A==ordinal(n//3):controls.append(dict(triangles=nu,H=r.h,minimum_positive_gain=min(positive)))
 digest.update(json.dumps([out,r.h,nu,min(positive)]).encode()+b'\n')
print(json.dumps(dict(status='PASS',checks=checks,seed=124034,structured_controls=controls,data_sha256=digest.hexdigest()),sort_keys=True))
