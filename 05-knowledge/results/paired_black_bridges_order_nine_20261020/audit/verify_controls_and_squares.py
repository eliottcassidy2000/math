"""Audit earlier-order controls and independently enumerate hostile n=5 squares."""
from pathlib import Path
import itertools,json,struct
ROOT=Path(__file__).resolve().parent.parent
OLD=ROOT.parent/'session14-metagraph-census'
OUT=ROOT/'audit'
controls={};cycles={}
for n in range(3,9):
 p=OLD/('n8' if n==8 else '')/f'census_n{n}.json'
 d=json.loads(p.read_text());N=len(d['nodes'])
 if n>=7:
  b=json.loads((ROOT/f'census_n{n}.json').read_text())
  mp=struct.unpack('<'+'i'*len(d['pi']),(ROOT/f'census_n{n}.map').read_bytes())
  assert list(mp)==d['pi']
  assert b['classes']==len(d['classes']) and b['merged']==N
  for u,v in zip(d['nodes'],b['nodes']):
   for k in ['node','representative','mass','blue','SC','c3']:assert u[k]==v[k]
   assert u['reflection_pairs']==v['degree']
  for name in ['blue','paired_black']:
   c=d['graphs'][name];e=b[name]
   for k in ['edges','loops']:assert c[k]==e[k]
   assert c['component_count']==e['components']
   if name=='blue':assert len(c['bridges'])==e['bridge_count']
   else:assert sorted(c['bridges'])==sorted(sorted(x['ends']) for x in e['bridges'])
  controls[n]={'map_entries_identical':len(mp),'nodes_identical':N,'ordinary_classes':b['classes'],'all_summary_checks':'pass'}
 adj=[set() for _ in range(N)];degrees=[0]*N
 counts={'loops':0,'digon_edges':0,'triangle_edges':0,'quadrilateral_edges':0,'leaf_edges':0,'uncovered_edges':0}
 for a,b,c in d['paired_black_edges']:
  degrees[a]+=c;degrees[b]+=c
  if a!=b:adj[a].add(b);adj[b].add(a)
 for a,b,c in d['paired_black_edges']:
  if a==b:counts['loops']+=c
  elif c>1:counts['digon_edges']+=c
  elif min(degrees[a],degrees[b])==1:counts['leaf_edges']+=1
  elif adj[a]&adj[b]:counts['triangle_edges']+=1
  elif any((adj[x]&adj[b])-{a,b} for x in adj[a]-{b}):counts['quadrilateral_edges']+=1
  else:counts['uncovered_edges']+=1
 assert counts['uncovered_edges']==0
 assert sum(counts.values())==d['graphs']['paired_black']['edges']
 cycles[n]=counts
(OUT/'control_comparison.json').write_text(json.dumps(controls,indent=2)+'\n')
(OUT/'lower_order_cycles.json').write_text(json.dumps(cycles,indent=2)+'\n')
# Entirely direct five-vertex permutation enumeration, without importing producer code.
n=5;coords=[(i,j) for j in range(n-1,1,-1) for i in range(j-1)];m=len(coords);full=(1<<m)-1
sig=[coords.index((n-1-j,n-1-i)) for i,j in coords]
def refl(t):return sum(((t>>i)&1)<<j for i,j in enumerate(sig))
def line(t):return min(t,full^t,refl(t),full^refl(t))
def black(t):return t!=refl(t)
def encode(adj,p):return sum((1-((adj[p[i]]>>p[j])&1))<<k for k,(i,j) in enumerate(coords))
def reverse(adj,p):
 a=adj[:]
 for x,y in zip(p,p[1:]):a[x]^=1<<y;a[y]^=1<<x
 return a
p=tuple(range(n));edgep={frozenset((a,b)) for a,b in zip(p,p[1:])};rows=[]
d=json.loads((OLD/'census_n5.json').read_text())
for t in range(full+1):
 if line(t)!=t or not black(t):continue
 a,b=d['pi'][t],d['pi'][full^t]
 if a==b or min(d['nodes'][a]['reflection_pairs'],d['nodes'][b]['reflection_pairs'])==1:continue
 adj=[0]*n
 for i in range(n-1):adj[i]|=1<<(i+1)
 for k,(i,j) in enumerate(coords):
  if t>>k&1:adj[j]|=1<<i
  else:adj[i]|=1<<j
 for q in itertools.permutations(range(n)):
  if any(not(adj[x]>>y&1) for x,y in zip(q,q[1:])):continue
  if any(frozenset((x,y)) in edgep for x,y in zip(q,q[1:])):continue
  masks=[t,encode(adj,q),encode(reverse(adj,p),q),encode(reverse(adj,q),p)]
  lines=list(map(line,masks))
  rows.append({'target':t,'path':q,'masks':masks,'all_four_black':all(map(black,masks)),'black_flags':list(map(black,masks)),'lines':lines,'target_occurrences':lines.count(t)})
assert not any(r['all_four_black'] and r['target_occurrences']==1 for r in rows)
assert len(rows)==5 and len({r['target'] for r in rows})==4
(OUT/'commuting_n5.json').write_text(json.dumps({'order':5,'total_edge_disjoint_path_witnesses':len(rows),'distinct_target_edges':4,'witnesses':rows},indent=2)+'\n')
print('Control, lower-order cycle, and n=5 hostile square checks passed.')
