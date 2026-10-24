import itertools as it
import json
from functools import lru_cache
from decoder import decode
ROWS=[[102,108,408,657,771,220,412,787,555,103],[102,108,664,401,771,220,412,787,551,107]]
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def forward(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if not dp[mask][v]:continue
   for w in range(n):
    if not(mask>>w&1) and rows[v]>>w&1:dp[mask|1<<w][w]+=dp[mask][v]
 return sum(dp[-1])
def backwards(rows):
 n=len(rows)
 @lru_cache(None)
 def visit(mask,v):
  if mask==(1<<v):return 1
  sub=mask^(1<<v)
  return sum(visit(sub,w) for w in range(n) if sub>>w&1 and rows[w]>>v&1)
 return sum(visit((1<<n)-1,v) for v in range(n))
def literal_aut(rows):
 groups={}
 for v,r in enumerate(rows):groups.setdefault(r.bit_count(),[]).append(v)
 cells=list(groups.values());domain=sum(cells,[]);out=[]
 for pieces in it.product(*(it.permutations(C) for C in cells)):
  g=dict(zip(domain,sum(pieces,())))
  if all(((rows[v]>>w)&1)==((rows[g[v]]>>g[w])&1) for v in range(10) for w in range(v)):out.append([g[v] for v in range(10)])
 return out
ans=[]
for rows in ROWS:
 d=decode(rows);a=literal_aut(rows);h=forward(rows)
 require(h==backwards(rows),'two Hamiltonian algorithms')
 require(len(a)==1==d['group_order'],'independent full ordinary automorphisms')
 require(d['error']==40,'root error')
 require(sorted(map(int.bit_count,rows))==[4]*5+[5]*5,'same scores')
 require(d['cells']==[list(range(5)),list(range(5,10))],'literal stable cells')
 require(d['edges']==[(1,0,'rigid')],'same rigid crossing graph')
 ans.append({'rows':rows,'hamiltonian_paths':h,'automorphism_order':len(a),'refinement_error':d['error']})
require([r['hamiltonian_paths'] for r in ans]==[15465,15549],'explicit counts')
print(json.dumps({'status':'PASS','witnesses':ans},sort_keys=True))
