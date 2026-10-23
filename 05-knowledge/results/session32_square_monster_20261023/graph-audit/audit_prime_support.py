#!/usr/bin/env python3
"""Audit exactly the 105 pairs of the supplied 15-prime support, no group claims."""
from pathlib import Path
from itertools import combinations
from collections import Counter
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
parts=[[17,19,23,31,47],[29,41,59,71],[2,3,5,7,11,13]]
P=sorted(sum(parts,[]));squares={q*q for q in range(2,13)}
pairs=list(combinations(P,2));E=[e for e in pairs if sum(e) in squares]
checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
def component_sets(vertices,edges):
 p={v:v for v in vertices}
 def root(v):
  while p[v]!=v:p[v]=p[p[v]];v=p[v]
  return v
 for a,b in edges:p[root(a)]=root(b)
 blocks={}
 for v in vertices:blocks.setdefault(root(v),[]).append(v)
 return sorted(sorted(c) for c in blocks.values())
expected=[(2,7),(2,23),(2,47),(3,13),(5,11),(5,31),(5,59),(7,29),
          (13,23),(17,19),(17,47),(23,41),(29,71),(41,59)]
need(len(pairs)==105 and E==expected,'complete 105-pair edge census')
C=component_sets(P,E);need(C==[P] and len(E)==len(P)-1,'connected tree')
d=Counter(v for e in E for v in e);leaves=[v for v in P if d[v]==1]
need(leaves==[3,11,19,31,71],'five leaves')
rows=[{'vertices':block,'edges':[e for e in E if set(e)<=set(block)]} for block in parts]
need([row['edges'] for row in rows]==[[(17,19),(17,47)],[(29,71),(41,59)],[(2,7),(3,13),(5,11)]],
     'all three induced blocks')
B=sorted(set(P)-set(parts[1]));BE=[e for e in E if set(e)<=set(B)];BC=component_sets(B,BE)
need(BC==[[2,3,7,13,17,19,23,47],[5,11,31]],'Baby components')
need({(23,41),(41,59),(5,59),(7,29),(29,71)}<=set(E),'restoration paths')
producer=json.loads((ROOT/'check_prime_block_graphs.json').read_text())
need(E==[tuple(e) for e in producer['all_edges']] and C==producer['components'],'producer graph')
need(leaves==producer['leaves'] and producer['cycle_rank']==0,'producer invariants')
for row,old in zip(rows,producer['parts']):
 need(row['vertices']==old['vertices'] and row['edges']==[tuple(e) for e in old['edges']],
      'producer induced block')
sources=['PRIME_BLOCK_CARDINALITY_HOSTILE.md',
         'check_prime_block_graphs.py','check_prime_block_graphs.json']
out={'status':'PASS bounded prime-support graph audit','explicit_checks':checks,
     'scope':'105 unordered pairs on the supplied prime support; no re-verification of group-order prime lists',
     'pair_count':105,'vertices':P,'all_edges':E,'components':C,'leaves':leaves,'parts':rows,
     'Baby_vertices':B,'Baby_edges':BE,'Baby_components':BC,
     'source_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources}}
print(json.dumps(out,indent=2))
