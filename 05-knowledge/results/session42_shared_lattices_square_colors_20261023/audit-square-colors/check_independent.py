from pathlib import Path
import json,math,itertools,collections,hashlib,subprocess
from fractions import Fraction
ROOT=Path(__file__).resolve().parent.parent
P=ROOT/'square-colors' if (ROOT/'square-colors').is_dir() else ROOT/'session42-square-colors'
REPO=ROOT/'session35-publication/clean-repo'
if not (REPO/'.git').exists():
 REPO=next((a for a in [ROOT,*ROOT.parents] if (a/'.git').exists()),None)
 if REPO is None:raise RuntimeError('Pinned Git checkout required for source provenance audit')
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def parents(n):return [a for a in range(1,n) if math.isqrt(a+n)**2==a+n]
A={}
for n in range(1,1001):
 pp=parents(n)
 A[n]=set().union(*(A[a] for a in pp)) if pp else {n}
need([n for n in A if not parents(n)]==[1,2,4],'roots')
need(max(n for n in A if len(A[n])<3)==70,'saturation minimal')
need(A[70]=={1,4} and all(A[n]=={1,2,4} for n in range(71,100)),'induction certificate')
need(min(n for n in A if len(A[n])==3)==22,'first fully mixed')
pure={c:[n for n in A if A[n]=={c}] for c in (4,2,1)}
need(pure=={4:[4,5,11,12],2:[2,7,9,16,18],1:[1,3,6,8,10,15,17,19,30]},'complete pure founder islands')
# Enumerate distinct colors directly, weighting by compatible parents rather
# than enumerating choices and identifying colors afterward.
states={(1,):1}
for n in range(2,26):
 nxt=collections.defaultdict(int);pp=parents(n)
 for colors,w in states.items():
  multiplicities=collections.Counter(colors[p-1] for p in pp) if pp else {n:1}
  for c,k in multiplicities.items():nxt[colors+(c,)]+=w*k
 states=dict(nxt)
need(len(states)==192 and sum(states.values())==1024,'selector polynomial totals')
need(not any(c[12]==1 and c[22]==4 for c in states),'joint marginal hostile')
# Independent equality-constraint rank via connected components, no union-find.
def rank_graph(p,q):
 graph=collections.defaultdict(set)
 for a,b in zip(p,q):graph[a].add(b);graph[b].add(a)
 seen=set();components=0
 for a in graph:
  if a in seen:continue
  components+=1;stack=[a];seen.add(a)
  while stack:
   for b in graph[stack.pop()]:
    if b not in seen:seen.add(b);stack.append(b)
 return len(graph)-components
def ew(p):return [tuple(sorted((a,b))) for a,b in zip(p,p[1:])]
raw=json.loads((P/'inherited_paths.json').read_text());receipt=json.loads((P/'receipt.json').read_text())
source=subprocess.check_output(['git','-C',str(REPO),'show',raw['source_commit']+':'+raw['source_path']])
need(hashlib.sha256(source).hexdigest()==raw['source_sha256'],'source receipt hash')
source_json=json.loads(source)
source_rows={r['n']:r for r in source_json['rows']}
for row in raw['rows']:
 need(row['paths']==source_rows[row['n']]['paths'],'exact inherited path extraction')
# Retain canonical endpoint orientation while leaving source completeness inherited.
rows={r['n']:[p if p[0]<p[-1] else p[::-1] for p in r['paths']] for r in raw['rows']}
for n,ps in rows.items():
 need(len(set(map(tuple,ps)))==len(ps),'unique paths')
 for p in ps:
  need(sorted(p)==list(range(1,n+1)) and all(math.isqrt(a+b)**2==a+b for a,b in zip(p,p[1:])),'actual path')
 for kind in ('vertex','edge'):
  ranks=collections.Counter()
  for p,q in itertools.combinations(ps,2):
   if (p[0],p[-1])==(q[0],q[-1]):ranks[rank_graph(p,q) if kind=='vertex' else rank_graph(ew(p),ew(q))]+=1
  recorded=receipt['collision_probabilities'][f'{n}-{kind}']
  need({str(k):v for k,v in sorted(ranks.items())}==recorded['ranks'],'independent ranks')
  for q in (2,3):
   bound=sum((Fraction(c,q**r) for r,c in ranks.items()),Fraction())
   need(str(bound)==recorded['union_bounds'][str(q)],'union bound exact')
# Tiny full coloring enumeration verifies collision probability and joint-event rank.
paths=[(0,1,2,3),(0,2,1,3),(0,1,3,2)]
for q in (2,3):
 colorings=list(itertools.product(range(q),repeat=4));pairs=list(itertools.combinations(paths,2))
 actual=sum(len({tuple(c[i] for i in p) for p in paths})==3 for c in colorings)
 ie=0
 for mask in range(1<<len(pairs)):
  left=[];right=[]
  for i,(p,r) in enumerate(pairs):
   if mask>>i&1:left+=p;right+=r
  rank=rank_graph(left,right)
  ie+=(-1)**mask.bit_count()*q**(4-rank)
 need(actual==ie,'inclusion exclusion independent joint ranks')
# Independently certify first odd cycle / first triangle and 3-4-5 edge counts.
for n in range(1,16):
 graph={i:set() for i in range(1,n+1)}
 for b in graph:
  for a in parents(b):graph[a].add(b);graph[b].add(a)
 colors={};okay=True
 for a in graph:
  if a in colors:continue
  colors[a]=0;stack=[a]
  while stack:
   u=stack.pop()
   for v in graph[u]:
    if v not in colors:colors[v]=1-colors[u];stack.append(v)
    elif colors[v]==colors[u]:okay=False
 need(okay==(n<15),'first odd cycle')
tri=next((c,a,b) for c in range(1,31) for a,b in itertools.combinations(parents(c),2) if math.isqrt(a+b)**2==a+b)
need(tri==(30,6,19),'first triangle')
islands=[{12,4,5,11},{2,18,7,9,16},{10,6,3,1,8,17}]
need([sum(math.isqrt(a+b)**2==a+b for a,b in itertools.combinations(S,2)) for S in islands]==[3,4,5],'induced edge counts')
print(json.dumps({'status':'PASS','checks':checks,'distinct_colorings25':len(states),'selector_forests25':sum(states.values()),'tail_start':71,'independent_rank_controls':'all supplied paths25..32','scope':'source completeness inherited; no Q36/Q80 complete census rerun'},indent=2))
