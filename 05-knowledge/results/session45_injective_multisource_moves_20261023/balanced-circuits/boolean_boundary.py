from itertools import combinations
from pathlib import Path
import json
H=Path(__file__).resolve().parent;pool=json.loads((H/'pool.json').read_text());L={int(s):tuple(v) for s,v in pool['labels'].items()};cases=json.loads((H/'candidates.json').read_text())['found'];reports=[]
for key,c in sorted(cases.items(),key=lambda z:int(z[0])):
 ids=c['plus']+c['minus'];k=len(c['plus']);n=len(ids);full=(1<<n)-1;counts=[0]*(1<<n);p=[0]*(1<<n);x=[0]*(1<<n)
 for bits in range(1,1<<n):
  low=bits&-bits;j=low.bit_length()-1;prev=bits^low;sg=1 if j<k else -1;counts[bits]=counts[prev]+sg;p[bits]=p[prev]^L[ids[j]][0];x[bits]=x[prev]+sg*L[ids[j]][1]
 for modulus in (4,128):
  S={b for b in range(1<<n) if counts[b]==0 and p[b]==0 and x[b]%modulus==0}
  if {b^full for b in S}!=S or not {0,full}<=S:raise RuntimeError('complement/constant hypothesis')
  groups=[];unused=set(range(n))
  while unused:
   i=min(unused);C={j for j in unused if all(((b>>i)&1)==((b>>j)&1) for b in S)};unused-=C;groups.append(sorted(C))
  pair_closure={b for b in range(1<<n) if all(len({(b>>j)&1 for j in C})==1 for C in groups)}
  # Literal binary projections, not only the implied equality shortcut.
  pairproj={(i,j):{((b>>i)&1,(b>>j)&1) for b in S} for i,j in combinations(range(n),2)}
  if any(A not in ({(0,0),(1,1)},{(0,0),(0,1),(1,0),(1,1)}) for A in pairproj.values()):raise RuntimeError('projection hypothesis')
  witness=next(([0,a,b,a&b] for a,b in combinations(sorted(S),2) if (a&b) not in S),None)
  if (pair_closure==S)!=(witness is None):raise RuntimeError('exact Boolean algebra criterion')
  strong=None
  if n==16 and modulus==4:
   strong=[0,1806,6500,260]
   if not all(v in S for v in strong[:3]) or (1806&6500)!=260 or 260 in S or counts[260]!=0:raise RuntimeError('current-only majority hostile')
  reports.append({'current_only_majority_hostile':strong,'slots':n,'paired_X_modulus':modulus,'accepted':len(S),'equality_components':groups,'pairwise_closure_size':len(pair_closure),'exact_pairwise_representation':pair_closure==S,'majority_hostile':witness})
print(json.dumps({'status':'PASS','reports':reports},indent=2))
