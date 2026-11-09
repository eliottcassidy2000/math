#!/usr/bin/env python3
"""Independent permutation scatter, not endpoint DP. Run from repository root."""
import itertools,json
from pathlib import Path
OUT=Path(__file__).parent

def paths(a):
 return [p for p in itertools.permutations(range(len(a))) if all(a[u][v] for u,v in zip(p,p[1:]))]

def decode(raw,zero='empty'):
 n=len(raw); a=[[0]*n for _ in range(n)]
 for i,j in itertools.combinations(range(n),2):
  x,y=raw[i][j],raw[j][i]
  if x!=y: a[i][j],a[j][i]=x,y
  elif x==0:
   if zero=='reciprocal': a[i][j]=a[j][i]=1
   elif zero=='invalid': raise ValueError('00 disallowed')
 return a

def scan(n):
 pairs=list(itertools.combinations(range(n),2)); idx={p:k for k,p in enumerate(pairs)}; m=len(pairs)
 totals=[0]*(1<<m); uses=[[0]*m for _ in totals]
 for p in itertools.permutations(range(n)):
  fixed=bits=0; used=[]
  for u,v in zip(p,p[1:]):
   k=idx[min(u,v),max(u,v)]; fixed|=1<<k; bits|=(u<v)<<k; used.append(k)
  free=((1<<m)-1)^fixed; s=free
  while True:
   t=bits|s; totals[t]+=1
   for k in used: uses[t][k]+=1
   if s==0: break
   s=(s-1)&free
 result={'n':n,'tournaments':len(totals),'max_tournament_H':max(totals),'targets':{}}
 for h in (7,21):
  hits=[(t,k) for t in range(1<<m) for k in range(m) if totals[t]-uses[t][k]==h]
  item={'tournament_arc_deletions':len(hits),'distinct_oriented_graphs':len(hits)//2}
  # Each resulting missing pair has exactly two tournament completions.
  if hits:
   t,k=hits[0]; a=[[0]*n for _ in range(n)]
   for q,(i,j) in enumerate(pairs):
    if q!=k: a[i][j],a[j][i]=(1,0) if t>>q&1 else (0,1)
   pp=paths(a); assert len(pp)==h
   i,j=pairs[k]; completion=[]
   for u,v in ((i,j),(j,i)):
    a[u][v]=1; completion.append(len(paths(a))); a[u][v]=0
   raw=[r[:] for r in a]; raw[i][j]=raw[j][i]=1
   assert decode(raw,'empty')==a and decode(raw,'invalid')==a
   assert decode(raw,'reciprocal')==a # all other pairs unequal
   item.update(rows=[''.join(map(str,r)) for r in a],raw_11_empty_rows=[''.join(map(str,r)) for r in raw],missing_pair=pairs[k],completion_counts=completion,paths=[''.join(map(str,p)) for p in pp])
  result['targets'][h]=item
 return result
def oriented_scan(n):
 pairs=list(itertools.combinations(range(n),2)); idx={p:k for k,p in enumerate(pairs)}; m=len(pairs)
 counts=[0]*(3**m)
 for p in itertools.permutations(range(n)):
  used=set(); base=0
  for u,v in zip(p,p[1:]):
   k=idx[min(u,v),max(u,v)]; used.add(k); base+=(1 if u<v else 2)*3**k
  free=[k for k in range(m) if k not in used]
  for choices in itertools.product(range(3),repeat=len(free)):
   counts[base+sum(x*3**k for x,k in zip(choices,free))]+=1
 return {'n':n,'universe':3**m,'max_H':max(counts),'H7':counts.count(7),'H21':counts.count(21)}

if __name__=='__main__':
 result={'method':'Scatter every permutation over compatible tournaments; remove all paths using a specified arc. Independently enumerate witness paths.','scans':[scan(n) for n in range(2,7)],'all_oriented_through_five':[oriented_scan(n) for n in range(1,6)]}
 (OUT/'independent.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))
