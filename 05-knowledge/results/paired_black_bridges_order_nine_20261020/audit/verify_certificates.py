"""Independent edge-identity Tarjan checks of sparse Q certificate and blue graph."""
import json, mmap, struct
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
DATA=json.loads((ROOT/'census_n9.json').read_text())
N=DATA['merged'];n=9;m=28;full=(1<<m)-1
fh=(ROOT/'census_n9.map').open('rb');mp=mmap.mmap(fh.fileno(),0,access=mmap.ACCESS_READ)
def node(t):return struct.unpack_from('<i',mp,4*t)[0]
coords=[(i,j) for j in range(n-1,1,-1) for i in range(j-1)]
sig=[coords.index((n-1-j,n-1-i)) for i,j in coords]
def refl(t):return sum(((t>>i)&1)<<sig[i] for i in range(m))
def tarjan(edges):
 adj=[[] for _ in range(N)]
 for eid,(t,a,b) in enumerate(edges):adj[a].append((b,eid));adj[b].append((a,eid))
 tin=[-1]*N;low=[-1]*N;parent=[-1]*N;parentedge=[-1]*N;idx=[0]*N;tm=0;cc=0;out=[]
 for root in range(N):
  if tin[root]>=0:continue
  cc+=1;tin[root]=low[root]=tm;tm+=1;stack=[root]
  while stack:
   v=stack[-1]
   if idx[v]<len(adj[v]):
    w,eid=adj[v][idx[v]];idx[v]+=1
    if eid==parentedge[v]:continue
    if tin[w]>=0:low[v]=min(low[v],tin[w]);continue
    parent[w]=v;parentedge[w]=eid;tin[w]=low[w]=tm;tm+=1;stack.append(w)
   else:
    stack.pop();p=parent[v]
    if p>=0:
     if low[v]>tin[p]:out.append(parentedge[v])
     low[p]=min(low[p],low[v])
 return cc,[edges[e] for e in out]
class DSU:
 def __init__(self):self.p=list(range(N));self.sz=[1]*N
 def find(self,x):
  while x!=self.p[x]:self.p[x]=self.p[self.p[x]];x=self.p[x]
  return x
 def add(self,x,y):
  x=self.find(x);y=self.find(y)
  if x==y:return False
  if self.sz[x]<self.sz[y]:x,y=y,x
  self.p[y]=x;self.sz[x]+=self.sz[y];return True
edges=[];forests=[DSU(),DSU()];counts=Counter();seen=set()
for line in (ROOT/'census_n9_sparse_certificate.tsv').read_text().splitlines()[1:]:
 t,a,b,k=map(int,line.split());assert t not in seen;seen.add(t)
 rt=refl(t);assert t!=rt and t==min(t,rt,full^t,full^rt)
 assert (node(t),node(full^t))==(a,b)
 assert k in (1,2) and forests[k-1].add(a,b)
 edges.append((t,a,b));counts[k]+=1
cc,bridges=tarjan(edges)
expected=sorted((e['tiling'],*e['ends']) for e in DATA['paired_black']['bridges'])
assert sorted(bridges)==expected
assert cc==DATA['paired_black']['components']
sparse={'certificate_edges':len(edges),'forest_sizes':dict(counts),'components':cc,'bridges':sorted(bridges),'all_checks':'passed'}
# Generate fixed masks via the independent coordinate involution.
orbits=[];visited=set()
for i in range(m):
 if i not in visited:
  block={i,sig[i]};visited.update(block);orbits.append(sum(1<<j for j in block))
fixed=[0]
for mask in orbits:fixed += [x|mask for x in fixed]
blue_edges=[(t,node(t),node(full^t)) for t in fixed if t<full^t]
blue_degree=Counter()
for _,a,b in blue_edges:blue_degree[a]+=1;blue_degree[b]+=1
assert all(blue_degree[u]==r['blue'] for u,r in enumerate(DATA['nodes']))
cc,bridges=tarjan(blue_edges)
other=json.loads((ROOT/'blue_n9_bridge_audit.json').read_text())
assert cc==other['components']
assert sorted((t,*sorted((a,b))) for t,a,b in bridges)==sorted((r['tiling'],*r['ends']) for r in other['bridges'])
assert len(bridges)==24 and all(min(blue_degree[a],blue_degree[b])==1 for _,a,b in bridges)
blue={'edges':len(blue_edges),'support_pairs':len({tuple(sorted((a,b))) for _,a,b in blue_edges}),'components':cc,'bridge_count':len(bridges),'leaf_bridges':len(bridges),'all_checks':'passed'}
(ROOT/'audit/certificates_audit.json').write_text(json.dumps({'sparse_paired_black':sparse,'blue':blue},indent=2)+'\n')
print(sparse);print(blue)
# Check every Q leaf's literal strong-component word, not just its saved type.
leaf_rows=[]
for rec in DATA['nodes']:
 if rec['degree']!=1:continue
 t=rec['representative'];adj=[0]*n
 for i in range(n-1):adj[i]|=1<<(i+1)
 for k,(i,j) in enumerate(coords):
  if t>>k&1:adj[j]|=1<<i
  else:adj[i]|=1<<j
 reach=[a|(1<<i) for i,a in enumerate(adj)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 parts=[];left=set(range(n))
 while left:
  i=min(left);part={j for j in left if reach[i]>>j&1 and reach[j]>>i&1};parts.append(part);left-=part
 frozen=list(parts)
 parts=sorted(parts,key=lambda s:-sum(bool(adj[next(iter(s))]>>next(iter(t))&1) for t in frozen if s!=t))
 assert all(adj[u]>>v&1 for i,s in enumerate(parts) for t in parts[i+1:] for u in s for v in t)
 word=list(map(len,parts))
 assert all(a in [1,3] for a in word) and word!=word[::-1]
 assert not rec['SC'] and rec['blue']==0 and rec['mass']==2
 leaf_rows.append({'node':rec['node'],'representative':t,'SCC_word':word})
assert len(leaf_rows)==7
isolates=[r for r in DATA['nodes'] if r['degree']==0]
assert len(isolates)==6 and all(r['mass']==r['blue'] for r in isolates)
(ROOT/'audit/leaf_family_audit.json').write_text(json.dumps({'all_seven_leaves_nonpalindromic_singleton_C3_words':leaf_rows,'pure_blue_Q_isolates':6},indent=2)+'\n')
