from itertools import permutations,combinations
from fractions import Fraction as F
from collections import Counter
import json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def graphkey(n,E):
 return min(tuple(sorted(tuple(sorted((p[a],p[b]))) for a,b in E)) for p in permutations(range(n)))
def tournamentkey(n,B):
 arcs={(b,a) if (a,b) in B else (a,b) for a,b in combinations(range(n),2)}
 return min(tuple(int((p[a],p[b]) not in arcs) for a,b in combinations(range(n),2)) for p in permutations(range(n)))
def components(K):
 unseen=set(range(len(K)));out=[]
 while unseen:
  a=min(unseen);todo=[a];unseen.remove(a);block=[]
  while todo:
   i=todo.pop();block.append(i)
   for j in list(unseen):
    if K[i][j]:unseen.remove(j);todo.append(j)
  out.append(sorted(block))
 return out
witnesses=[{(0,1),(0,2)},{(0,1),(1,2)},{(0,1),(2,3)}]
compressed=[graphkey(4,E) for E in witnesses]
full=[graphkey(5,{(a,b+1) for a,b in E}) for E in witnesses]
need(compressed[0]==compressed[1] and full[0]!=full[1],'compressed does not determine full')
need(full[0]==full[2] and compressed[0]!=compressed[2],'full does not determine compressed')
reports=[]
n=4;CE=list(combinations(range(n-1),2));FE=set(combinations(range(n),2))
for mode in ['full','compressed']:
 for paired in [False,True]:
  cells=Counter()
  for x in range(1<<len(CE)):
   C={e for k,e in enumerate(CE) if x>>k&1};B={(a,b+1) for a,b in C}
   gs=graphkey(n,B) if mode=='full' else graphkey(n-1,C);ts=tournamentkey(n,B)
   if paired:
    gs=min(gs,graphkey(n,FE-B) if mode=='full' else graphkey(n-1,set(CE)-C));ts=min(ts,tournamentkey(n,FE-B))
   cells[gs,ts]+=1
  G=sorted({g for g,t in cells});T=sorted({t for g,t in cells});C=[[cells[g,t] for t in T] for g in G];r=list(map(sum,C));c=[sum(row[j] for row in C) for j in range(len(T))]
  K=[[sum((F(C[i][j]*C[k][j],r[i]*c[j]) for j in range(len(T))),F()) for k in range(len(G))] for i in range(len(G))]
  comp=components(K);need(len(comp)==(1 if mode=='compressed' and paired else 2),'pairing components')
  if mode=='compressed' and paired:need(len(K)==2 and K[0][0]+K[1][1]-1==F(7,15),'compressed paired nonconstant eigenvalue')
  reports.append({'mode':mode,'paired':paired,'row_masses':r,'components':comp,'incidence':C,'transition':[[str(v) for v in row] for row in K]})
print(json.dumps({'status':'PASS','incomparable_partition_witnesses':[sorted(x) for x in witnesses],'order_four_comparisons':reports},indent=2,sort_keys=True))
