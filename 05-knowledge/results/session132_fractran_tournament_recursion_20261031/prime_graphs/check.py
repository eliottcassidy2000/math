#!/usr/bin/env python3
"""Independent exact verifier: pairwise edge construction and unpruned DP.
No imports from the producer; no optimization-sensitive assert statements.
"""
from pathlib import Path
from math import isqrt
from collections import Counter,deque
from functools import lru_cache
import json,sys
ROOT=Path(__file__).resolve().parent
GATES=0

def need(value,message):
 global GATES
 GATES+=1
 if not value:raise RuntimeError(message)

def raw_graph(n):
 a=[set() for _ in range(n)]
 for u in range(n):
  for v in range(u+1,n):
   s=u+v+2
   if isqrt(s)**2==s:a[u].add(v);a[v].add(u)
 return a

def components(a,deleted=frozenset(),removed_edge=None):
 unseen=set(range(len(a)))-set(deleted);out=[]
 while unseen:
  s=min(unseen);unseen.remove(s);seen={s};todo=[s]
  while todo:
   v=todo.pop()
   for w in a[v]:
    if removed_edge is not None and {v,w}==removed_edge:continue
    if w in unseen:unseen.remove(w);seen.add(w);todo.append(w)
  out.append(sorted(seen))
 return out

def girth(a):
 best=None
 for u in range(len(a)):
  for v in a[u]:
   if u>v:continue
   d={u:0};todo=deque([u])
   while todo and v not in d:
    x=todo.popleft()
    for y in a[x]:
     if {x,y}=={u,v}:continue
     if y not in d:d[y]=d[x]+1;todo.append(y)
   if v in d:
    length=d[v]+1
    if best is None or length<best:best=length
 return best

def ham_counts(a):
 n=len(a);full=(1<<n)-1;bits=[sum(1<<j for j in vs) for vs in a]
 leaves=[i for i,vs in enumerate(a) if len(vs)==1]
 if len(components(a))>1 or len(leaves)>2:pc=0
 else:
  @lru_cache(None)
  def walk(v,mask):
   if mask==full:return 1
   rem=bits[v]&~mask;ans=0
   while rem:
    b=rem&-rem;rem-=b;ans+=walk(b.bit_length()-1,mask|b)
   return ans
  # A fixed leaf is an endpoint of every Hamiltonian path; otherwise
  # sum all starts and quotient reversal (except the singleton).
  starts=leaves[:1] if leaves else range(n)
  pc=sum(walk(s,1<<s) for s in starts)
  if not leaves and n>1:pc//=2
  walk.cache_clear()
 if n<3 or min(map(len,a))<2:cc=0
 else:
  @lru_cache(None)
  def closed(v,mask):
   if mask==full:return int(0 in a[v])
   rem=bits[v]&~mask;ans=0
   while rem:
    b=rem&-rem;rem-=b;ans+=closed(b.bit_length()-1,mask|b)
   return ans
  cc=closed(0,1)//2
  closed.cache_clear()
 return pc,cc

def witness(vertices,n,cycle=False):
 need(len(vertices)==n and set(vertices)==set(range(1,n+1)),'witness vertices')
 edges=list(zip(vertices,vertices[1:]))
 if cycle:edges.append((vertices[-1],vertices[0]))
 for u,v in edges:need(u!=v and isqrt(u+v)**2==u+v,'witness square edge')

def main():
 if '--test-failing-gate' in sys.argv:need(False,'intentional failure control')
 data=json.loads((ROOT/'graphs.json').read_text());cert=json.loads((ROOT/'witnesses.json').read_text())
 primes=[n for n in range(2,98) if all(n%d for d in range(2,isqrt(n)+1))]
 need(data['first_25_primes']==primes and len(primes)==25,'first 25 primes')
 groups={'A':[2,3,5,7,11,13],'B':[17,19,23,31,47],'C':[29,41,59,71],'D':[37,43,53,61,71],'E':[73,79,83,89,97]}
 need(data['groups_as_supplied']==groups,'unaltered groups')
 flat=[p for g in groups.values() for p in g];counts=Counter(flat)
 need(len(flat)==25 and len(counts)==24 and counts[71]==2 and set(primes)-set(flat)=={67},'duplicate and omitted labels')
 for g,ps in groups.items():
  for p,typed in zip(ps,data['groups_typed'][g]):
   need(typed=={'prime_label':p,'prime_rank':primes.index(p)+1,'prime_order_graph':p,'rank_order_graph':primes.index(p)+1},'two typed interpretations')
 rows=data['all_graphs_1_through_97'];need([r['n'] for r in rows]==list(range(1,98)),'exact data universe')
 previous_edges=0
 for r in rows:
  n=r['n'];a=raw_graph(n);edges=[[u+1,v+1] for u in range(n) for v in sorted(a[u]) if u<v]
  need(edges==r['edge_list'],'pairwise edges n='+str(n))
  ds=Counter(map(len,a));cs=components(a);e=len(edges)
  need(r['edges']==e and r['degree_histogram']=={str(k):v for k,v in sorted(ds.items())},'degree histogram')
  need(r['min_degree']==min(ds) and r['max_degree']==max(ds) and r['leaves']==ds[1],'degree extremes')
  need(r['components']==len(cs) and r['component_sizes']==sorted(map(len,cs),reverse=True),'components')
  need(e-previous_edges==isqrt(2*n-1)-isqrt(n),'arrival edge formula')
  need(r['arrival_roots']==list(range(isqrt(n)+1,isqrt(2*n-1)+1)),'arrival roots')
  previous_edges=e
  need(e==sum(max(0,min((k*k-1)//2,n-(k*k)//2)) for k in range(2,isqrt(2*n-1)+1)),'root matching edge count')
  for v in range(1,n+1):
   expected=isqrt(n+v)-isqrt(v)-int(isqrt(2*v)**2==2*v)
   need(len(a[v-1])==expected,'vertex degree formula')
  ap=[v+1 for v in range(n) if len(components(a,{v}))>len(cs)]
  bridges=[[u+1,v+1] for u in range(n) for v in sorted(a[u]) if u<v and len(components(a,removed_edge={u,v}))>len(cs)]
  need(ap==r['articulation_points'] and bridges==r['bridges'],'deletion-based articulation/bridges')
  gr=girth(a);need(gr==r['girth'],'edge-deletion girth')
  need(r['odd_girth']==(None if n<15 else 5 if n<30 else 3),'odd-girth transition')
  if n<15:need(gr is None,'forest through 14')
  else:
   cyc=[1,3,6,10,15] if n<30 else [6,19,30]
   need(all(v-1 in a[u-1] for u,v in zip(cyc,cyc[1:]+cyc[:1])),'odd-girth explicit witness')
  odd=[{v for v in a[u] if isqrt(u+v+2)%2} for u in range(n)]
  ocs=components(odd);need(len(ocs)==r['odd_root_components'],'odd-root components')
  if n>=44:
   block={0:0,1:0,2:1,7:1,3:2,6:2,4:3,5:3}
   want=sorted(sorted(v-1 for v in range(1,n+1) if block[v%8]==b) for b in range(4))
   need(sorted(ocs)==want,'exact four residue components')
  need((len(cs)==1)==(n==1 or n>=14),'connected threshold')
  need((n>=3 and len(cs)==1 and not ap)==(n>=31),'2-connected threshold')
  if n<=32:
   pc,cc=ham_counts(a);c=cert['finite_hamilton'][n-1]
   need(c['n']==n and c['path']['count']==pc==r['hamilton_path_count'],'independent path count')
   need(c['cycle']['count']==cc==r['hamilton_cycle_count'],'independent cycle count')
   if pc:witness(c['path']['witness'],n)
   if cc:witness(c['cycle']['witness'],n,True)
 prime_cycles=cert['prime_cycles'];need([c['n'] for c in prime_cycles]==[p for p in primes if p>=37],'prime cycle universe')
 for c in prime_cycles:witness(c['cycle'],c['n'],True)
 # Open ears certify a spanning 2-connected subgraph, independently of cuts.
 ears=cert['ears31'];base=ears[0];need(base[0]==base[-1] and len(set(base[:-1]))==len(base)-1,'initial simple ear cycle')
 present=set(base[:-1])
 for j,ear in enumerate(ears):
  need(all(isqrt(u+v)**2==u+v for u,v in zip(ear,ear[1:])),'ear edges')
  if j:
   need(ear[0]!=ear[-1] and ear[0] in present and ear[-1] in present,'old distinct ear endpoints')
   need(len(set(ear))==len(ear) and not(set(ear[1:-1])&present),'new ear interior')
   present.update(ear)
 need(present==set(range(1,32)),'spanning open ears G31')
 hp_ranks=[n for n in range(2,26) if rows[n-1]['hamilton_path_count']]
 need(hp_ranks==[15,16,17,23,25],'first 25 nontrivial path orders')
 need([primes[i-1] for i in hp_ranks]==[47,53,59,83,97],'rank interpretation path prime labels')
 need(not any(rows[n-1]['hamilton_cycle_count'] for n in range(1,26)),'first 25 no cycles')
 # Minimal certificates and cross-group collision controls.
 need(len(components(raw_graph(19),{6}))==3,'G19 path cut obstruction at vertex7')
 need(rows[40]['arrival_roots']==rows[42]['arrival_roots']==rows[46]['arrival_roots']==[7,8,9],'cross-group shell collision')
 need(rows[66]['hamilton_cycle_exists'] and rows[70]['hamilton_cycle_exists'],'prime duplicate correction invisible to cycle flag')
 need(rows[18]['hamilton_path_count']==rows[19]['hamilton_path_count']==0,'rank duplicate correction invisible to path flag')
 # Verify table interpretation, graph identities, and vector panel edge sets.
 import csv,xml.etree.ElementTree as ET,subprocess,tempfile
 atlas_dir=tempfile.TemporaryDirectory(prefix='s132-prime-atlas-')
 subprocess.run([sys.executable,str(ROOT/'render_atlas.py'),'--output-dir',atlas_dir.name],check=True,capture_output=True)
 for typ,orders in [('prime_order',primes),('rank_order',list(range(1,26)))]:
  with (ROOT/(typ+'_first25.csv')).open() as f:table=list(csv.DictReader(f))
  need(len(table)==25,'25 table rows')
  for i,(n,row) in enumerate(zip(orders,table),1):
   p=primes[i-1];groups_at=','.join(g for g,ps in groups.items() if p in ps) or 'MISSING'
   need(int(row['i'])==i and int(row['prime_label'])==p and int(row['graph_n'])==n and row['supplied_groups']==groups_at,'typed table identity')
   need(int(row['edges'])==rows[n-1]['edges'] and int(row['components'])==rows[n-1]['components'],'table invariants')
  ns={'s':'http://www.w3.org/2000/svg'};root=ET.parse(Path(atlas_dir.name)/(typ+'_first25.svg')).getroot();panels=root.findall('s:g',ns)
  need(len(panels)==25,'25 SVG panels')
  for i,(n,panel) in enumerate(zip(orders,panels),1):
   need(int(panel.attrib['data-graph-n'])==n and int(panel.attrib['data-prime-label'])==primes[i-1] and int(panel.attrib['data-rank'])==i,'SVG panel identity')
   edges=[[int(e.attrib['data-u']),int(e.attrib['data-v'])] for e in panel.findall('s:line',ns)]
   vs=[int(v.attrib['data-vertex']) for v in panel.findall('s:circle',ns)]
   need(edges==rows[n-1]['edge_list'] and vs==list(range(1,n+1)),'SVG exact graph')
 atlas_dir.cleanup()
 result={'status':'PASS','gates':GATES,'graph_orders':'1..97','table_interpretations':['p_i -> G_p_i','p_i -> G_i'],'independent_hamilton_counts':'all orders 1..32, paths modulo reversal and cycles modulo dihedral symmetry','prime_cycle_witnesses':len(prime_cycles),'nontrivial_path_ranks_1_to_25':hp_ranks,'rank_path_prime_labels':[primes[i-1] for i in hp_ranks]}
 out=json.dumps(result,indent=2,sort_keys=True)+'\n'
 if '--output' in sys.argv:Path(sys.argv[sys.argv.index('--output')+1]).write_text(out)
 else:print(out,end='')
if __name__=='__main__':main()
