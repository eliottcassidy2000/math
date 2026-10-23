import itertools,json,hashlib
from pathlib import Path

def need(x,msg):
 if not x: raise RuntimeError(msg)
base=Path(__file__).resolve().parent.parent
p=next((base/name for name in ('session42-correspondence','correspondence') if (base/name/'CORRESPONDENCE_AND_OBSTRUCTIONS.md').exists()),None)
need(p is not None,'source correspondence manifest not found beside audit')
m=json.loads((Path(__file__).resolve().parent/'SOURCE_PINS.json').read_text())
for f in m['files']:
 b=(p/f['path']).read_bytes();need(len(b)==f['bytes'] and hashlib.sha256(b).hexdigest()==f['sha256'],'pin mismatch '+f['path'])
summary=[]
for n in range(2,6):
 E=list(itertools.combinations(range(n),2));perms=list(itertools.permutations(range(n)));full=(1<<len(E))-1
 def bits(seq):return sum(v<<k for k,v in enumerate(seq))
 gs={};ts={}
 for w in range(full+1):
  G=[[False]*n for _ in range(n)];T=[[False]*n for _ in range(n)]
  for k,(a,b) in enumerate(E):
   G[a][b]=G[b][a]=bool(w>>k&1)
   T[a][b]=not G[a][b];T[b][a]=G[a][b]
  gs[w]=min(bits(G[q[a]][q[b]] for a,b in E) for q in perms)
  ts[w]=min(bits(T[q[b]][q[a]] for a,b in E) for q in perms)
 overlap={(min(gs[w],gs[w^full]),min(ts[w],ts[w^full])) for w in range(full+1)}
 row=[n,len(set(gs.values())),len({min(gs[w],gs[w^full]) for w in gs}),len(set(ts.values())),len({min(ts[w],ts[w^full]) for w in ts}),len(overlap)]
 summary.append(row)
need(summary==[[2,2,1,1,1,1],[3,4,2,2,2,3],[4,11,6,4,3,14],[5,34,18,12,10,120]],'independent census differs')
# Four-order triple contribution using direct directed edges; outside contribution cancels.
signs=[]
for ab,ac,bc in itertools.product((0,1),repeat=3):
 if ac==bc:continue
 g={frozenset((0,1)):ab,frozenset((0,2)):ac,frozenset((1,2)):bc}
 def cycles(order):
  rank={a:k for k,a in enumerate(order)}
  arrows={(a,b) for a in range(3) for b in range(3) if a!=b and ((rank[a]<rank[b])!=bool(g[frozenset((a,b))]))}
  return int(all(sum((a,b) in arrows for b in range(3))==1 for a in range(3)))
 c=[cycles(q) for q in ((0,1,2),(1,0,2),(2,0,1),(2,1,0))]
 signs.append((c[1]-c[0])-(c[3]-c[2]))
need(set(signs)=={-2,2},'local four-order identity differs')
print(json.dumps({'status':'PASS','source_files_pinned':len(m['files']),'independent_direct_matrix_counts':summary,'four_order_local_signs':signs,'scope':'independent finite controls plus separate proof audit'},indent=2))
