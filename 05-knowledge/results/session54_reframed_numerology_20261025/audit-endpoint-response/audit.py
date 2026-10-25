#!/usr/bin/env python3
"""Fresh backward endpoint counts, independently constructed deletion sets and pool census."""
from pathlib import Path
from collections import Counter,deque
from itertools import combinations,permutations
import json,hashlib
D=Path(__file__).resolve().parent
P=D/'inputs/producer'
def require(b,m):
 if not b: raise ValueError(m)
M=json.loads((P/'MANIFEST.json').read_text())
require(set(M)=={str(f.relative_to(P)) for f in P.rglob('*') if f.is_file() and f.name!='MANIFEST.json'},'producer closure')
for f,h in M.items(): require(hashlib.sha256((P/f).read_bytes()).hexdigest()==h,'producer pin '+f)
for s in json.loads((P/'SOURCES.json').read_text()): require(hashlib.sha256((P/s['bundled']).read_bytes()).hexdigest()==s['sha256'],'source pin')
records={x['n']:x for x in json.loads((P/'inputs/result.json').read_text())['classification']}
reports=[]
for n in (7,11):
 V=set(range(n)); sq={i*i%n for i in range(1,n)}
 E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in sq]
 incoming=[sum(1<<u for u,v in E if v==t) for t in range(n)]
 # Backward recurrence for each initial vertex; no symmetry assumed in DP.
 counts=[]
 for start in range(n):
  dp=[[0]*n for _ in range(1<<n)];dp[1<<start][start]=1
  for mask in range(1,1<<n):
   if not(mask>>start&1) or mask==1<<start:continue
   for t in range(n):
    if t==start or not(mask>>t&1):continue
    prev=mask^(1<<t);pred=prev&incoming[t]
    while pred:
     b=pred&-pred;dp[mask][t]+=dp[prev][b.bit_length()-1];pred-=b
  counts.append(dp[-1])
 ab=(1,8) if n==7 else (628,1101)
 require(all(counts[s][t]==(0 if s==t else ab[int((t-s)%n not in sq)]) for s in range(n) for t in range(n)),'endpoint counts')
 rows=[frozenset(e for e in E if e[0]==v) for v in range(n)]
 cols=[frozenset(e for e in E if e[1]==v) for v in range(n)]
 basic={rows[u]|rows[v]:'outgoing' for u,v in combinations(range(n),2)}
 basic.update({cols[u]|cols[v]:'incoming' for u,v in combinations(range(n),2)})
 basic.update({rows[v]|cols[v]:'star' for v in range(n)})
 families=dict(basic)
 if n==7:
  for c in range(n):
   predecessors={u for u,v in E if v==c}
   for S,T in ((predecessors,{c}),(V-{c},predecessors|{c})):
    B=frozenset((u,v) for u,v in E if u in S and v not in T)
    require(B not in families,'Hall distinct');families[B]='Hall'
 masks={sum(1<<j for j,e in enumerate(E) if e in B) for B in families}
 require(masks==set(records[n]['all_minimum_blockers']),'independent family reconstruction')
 paths=[]
 if n==7:
  # Full permutation filter, not producer's adjacency DFS.
  for p in permutations(range(n)):
   es=tuple(zip(p,p[1:]))
   if all((v-u)%n in sq for u,v in es):paths.append((p,es))
  require(len(paths)==189,'P7 total')
 histograms=Counter();pathchecks=0
 for B,fam in families.items():
  # Independent signed graph solve anchored at head zero.
  graph=[[] for _ in range(2*n)]
  for u,v in E:
   z=int((u,v) in B);graph[u].append((n+v,-z));graph[n+v].append((u,z))
  potentials={n:0};q=deque([n])
  while q:
   u=q.popleft()
   for v,z in graph[u]:
    wanted=potentials[u]+z
    if v in potentials:require(potentials[v]==wanted,'potential consistency')
    else:potentials[v]=wanted;q.append(v)
  require(len(potentials)==2*n,'connected')
  deficiency=sum(potentials[u]-potentials[n+u] for u in range(n))
  require(deficiency==2,'deficiency')
  h=Counter()
  for s in range(n):
   for t in range(n):
    if s!=t:h[deficiency-potentials[t]+potentials[n+s]]+=counts[s][t]
  if paths:
   direct=Counter()
   for p,es in paths:
    value=sum(e in B for e in es);direct[value]+=1
    require(value==2-potentials[p[-1]]+potentials[n+p[0]],'P7 pathwise identity');pathchecks+=1
   require(direct==h,'P7 direct histogram')
  wanted=({1:57,2:129,3:3} if fam=='Hall' else {1:54,2:135}) if n==7 else {1:17290,2:77805}
  require(h==wanted,'histogram')
  histograms[(fam,tuple(sorted(h.items())))]+=1
 reports.append({'n':n,'endpoint_counts':counts,'total_paths':sum(map(sum,counts)),'blockers':len(families),'direct_path_checks':pathchecks,'histograms':[[f,list(h),k] for (f,h),k in sorted(histograms.items())]})

def labels(n):return [('O',p) for p in combinations(range(n),2)]+[('I',p) for p in combinations(range(n),2)]+[('S',(v,)) for v in range(n)]
def response(label,s,t):
 f,U=label
 return 2-int(t in U) if f=='O' else 2-int(s in U) if f=='I' else 2-int(s in U)-int(t in U)
def word(label,pool):return tuple(response(label,s,t) for s,t in pool)
finite=[]
for n in range(5,32):
 labs=labels(n);pool=[(0,v) for v in range(1,n)]+[(v,0) for v in range(1,n)]
 fullwords=[sum((2-r)<<i for i,r in enumerate(word(b,pool))) for b in labs]
 fullmask=(1<<len(pool))-1
 for u in range(1,n):
  for v in range(1,n):
   selected=[q for q in pool if q not in ((0,u),(v,0))]
   mask=fullmask^(1<<(u-1))^(1<<(n-2+v))
   sig=[w&mask for w in fullwords]
   require((len(set(sig))==n*n)==(u==v),'all double omissions')
   if u!=v:require(word(('O',tuple(sorted((u,v)))),selected)==word(('S',(v,)),selected),'unequal omissions witness')
  for v in range(u+1,n):
   for typ in ('O','I'):
    missing=((0,u),(0,v)) if typ=='O' else ((u,0),(v,0))
    selected=[q for q in pool if q not in missing]
    require(word((typ,(0,u)),selected)==word((typ,(0,v)),selected),'lower bound collision')
 finite.append({'n':n,'matching_designs':n-1,'queries':2*n-4})
# Independent full subset census: compare every pair of response words by bitmask.
exhaustive=[]
for n in range(5,8):
 pool=[(0,v) for v in range(1,n)]+[(v,0) for v in range(1,n)]
 words=[word(b,pool) for b in labels(n)]
 differences={sum((1<<i) for i,(x,y) in enumerate(zip(a,b)) if x!=y) for a,b in combinations(words,2)}
 census=Counter()
 for selected in range(1<<len(pool)):
  if all(selected&d for d in differences):census[selected.bit_count()]+=1
 require(min(census)==2*n-4 and census[2*n-4]==n-1,'exhaustive optimum and designs')
 exhaustive.append({'n':n,'separating_subsets_by_size':dict(sorted(census.items()))})
require(word(('O',(1,2)),[(0,1),(0,2),(1,0),(2,0)])==word(('I',(0,3)),[(0,1),(0,2),(1,0),(2,0)]),'n4 formal collision')
print(json.dumps({'status':'PASS','paley':reports,'formal_signature_checks':finite,'exhaustive_pool_census':exhaustive,'n4_formal_countermodel':True},indent=2,sort_keys=True))
