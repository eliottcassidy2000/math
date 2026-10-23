#!/usr/bin/env python3
"""Independent edge-subset audit: connected acyclic max-degree-two n-1-edge subsets."""
import json, math
from pathlib import Path

def solve(n):
 E=[(a,b) for b in range(2,n+1) for a in range(1,b) if math.isqrt(a+b)**2==a+b]
 # Force scarce vertices early; correctness independent of ordering.
 degrees=[0]*(n+1)
 for a,b in E: degrees[a]+=1;degrees[b]+=1
 E.sort(key=lambda e:(min(degrees[e[0]],degrees[e[1]]),sum(degrees[x] for x in e),e))
 m=len(E); deg=[0]*(n+1); rem=degrees[:]; parent=list(range(n+1)); chosen=[]; counts=[0]*m; endpoint={}; total=0; calls=0
 def root(a):
  while parent[a]!=a:a=parent[a]
  return a
 def dfs(i):
  nonlocal total,calls
  calls+=1
  k=len(chosen)
  if k>n-1 or k+m-i<n-1:return
  if any(deg[v]==0 and rem[v]==0 for v in range(1,n+1)) and n>1:return
  if i==m:
   if k!=n-1:return
   # Acyclic maintained: n-1 edges implies connected.
   total+=1
   for j in chosen:counts[j]+=1
   ends=tuple(v for v in range(1,n+1) if deg[v]==1)
   endpoint[ends]=endpoint.get(ends,0)+1
   return
  a,b=E[i];rem[a]-=1;rem[b]-=1
  ra,rb=root(a),root(b)
  if deg[a]<2 and deg[b]<2 and ra!=rb:
   parent[ra]=rb;deg[a]+=1;deg[b]+=1;chosen.append(i)
   dfs(i+1)
   chosen.pop();deg[a]-=1;deg[b]-=1;parent[ra]=ra
  dfs(i+1)
  rem[a]+=1;rem[b]+=1
 if n==1:total=1;endpoint[(1,)]=1
 else:dfs(0)
 if sum(counts)!=(n-1)*total: raise RuntimeError("edge participation mass failure")
 return {'n':n,'edges':m,'paths_up_to_reversal':total,'edge_counts':{f'{a}:{b}':c for (a,b),c in sorted(zip(E,counts))},'endpoint_counts':{':'.join(map(str,p)):c for p,c in sorted(endpoint.items())},'search_nodes':calls}
if __name__=='__main__':
 rows=[solve(n) for n in range(1,26)]
 out={'method':'edge-subset acyclic degree-constrained exhaustive search, independent of path DFS','rows':rows}
 Path(__file__).with_name('square-edge-subsets-receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print([(r['n'],r['paths_up_to_reversal']) for r in rows]);print('search nodes',sum(r['search_nodes'] for r in rows))
