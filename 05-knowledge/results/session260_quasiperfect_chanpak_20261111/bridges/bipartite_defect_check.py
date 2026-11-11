from itertools import combinations,permutations
from collections import Counter
from pathlib import Path
import json,random

def require(q,m):
 if not q:raise RuntimeError(m)
def transitive(out,vs):
 return not any(out[a]>>b&1 and out[b]>>c&1 and out[c]>>a&1 for a in vs for b in vs for c in vs if len({a,b,c})==3)
def probe(out):
 n=len(out);ps=[p for p in permutations(range(n)) if all(out[p[i]]>>p[i+1]&1 for i in range(n-1))];H=len(ps);cases=[]
 for z in range(n):
  V=[i for i in range(n) if i!=z];A={i for i in V if out[z]>>i&1};B=set(V)-A
  adj={i:set() for i in V}
  for a in A:
   for b in B:
    if out[a]>>b&1:adj[a].add(b);adj[b].add(a)
  fibres=Counter(frozenset(p[:p.index(z)]) for p in ps);D=0;nonnegative=True
  for mask in range(1<<len(V)):
   L={v for j,v in enumerate(V) if mask>>j&1};S=(A-L)|(B&L)
   eps=all(v in S or bool(adj[v]&S) for v in V);D+=eps
   require(fibres[frozenset(L)]>=eps,('fibre',out,z,L))
  acyclic=transitive(out,V)
  require(H>=D and (H==D)==acyclic,('global',out,z,H,D,acyclic))
  if not acyclic:require(H>=D+2,('gap',out,z))
  cases.append((H,D))
 return cases
rootcases=0;graphs=0;equal=0
for n in range(1,6):
 pairs=list(combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  out=[0]*n
  for bit,(i,j) in enumerate(pairs):
   if mask>>bit&1:out[i]|=1<<j
   else:out[j]|=1<<i
  cs=probe(out);rootcases+=len(cs);graphs+=1;equal+=sum(h==d for h,d in cs)
rng=random.Random(260);sample=0
for _ in range(100):
 n=6;out=[0]*n
 for i,j in combinations(range(n),2):
  if rng.randrange(2):out[i]|=1<<j
  else:out[j]|=1<<i
 probe(out);sample+=n
# Two independent graph edges: A={0,1},B={2,3}; z=4.
match=[]
for x,y in [(0,0),(0,1),(1,0),(1,1)]:
 out=[0]*5
 for i,j in combinations(range(5),2):
  if j==4:win,lose=(j,i) if i<2 else (i,j)
  elif i<2<=j:win,lose=(i,j) if j-i==2 else (j,i)
  elif (i,j)==(0,1):win,lose=(i,j) if x else(j,i)
  else:win,lose=(i,j) if y else(j,i)
  out[win]|=1<<lose
 cs=probe(out);match.append({'shore_bits':[x,y],'H':cs[4][0],'d':cs[4][1]})
res={'status':'FINITE-EXACT','all_labelled_tournaments_orders_1_to_5':graphs,'all_root_cases':rootcases,'equality_cases':equal,'additional_seeded_order6_root_cases':sample,'two_K2_all_shore_orientations':match}
Path(__file__).with_name('bipartite_defect_receipt.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
