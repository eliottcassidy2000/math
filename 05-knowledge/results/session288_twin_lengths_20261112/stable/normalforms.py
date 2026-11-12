from pathlib import Path
from itertools import permutations,product,combinations
from collections import defaultdict,Counter
import json,hashlib
HERE=Path(__file__).resolve().parent
SRC=HERE.parents[1]/'five_stress_connection_atlas_20261112/cyclic_port_receipt.json'
def need(x,m):
 if not x: raise RuntimeError(m)
def canon(rows,O):
 n=len(rows);cells=defaultdict(list)
 for v,r in enumerate(rows):cells[(bool(O>>v&1),r.bit_count(),(r&O).bit_count())].append(v)
 best=None
 for parts in product(*(permutations(cells[k]) for k in sorted(cells))):
  p=sum(parts,());rr=tuple(sum(1<<j for j,w in enumerate(p) if rows[v]>>w&1) for v in p)
  key=(sum(1<<j for j,v in enumerate(p) if O>>v&1),rr)
  if best is None or key<best:best=key
 return best

def contract(rows,O,a):
 p=[v for v in range(len(rows)) if v!=a]
 return [sum(1<<j for j,w in enumerate(p) if rows[v]>>w&1) for v in p],sum(1<<j for j,v in enumerate(p) if O>>v&1)
def run():
 source=json.loads(SRC.read_text());states={canon(r['rows'],r['selected']):r for r in source['positive_screen_certificates'] if r['individually_stable']}
 need(len(states)==161,'bank size')
 children={};arcs=0;diamonds=0
 for key,r in states.items():
  child=set();stable=set(r['individually_stable']);rows=r['rows'];O=r['selected'];n=len(rows)
  pairs=[]
  for a,b in combinations(range(n),2):
   if O>>a&1 or O>>b&1:continue
   if a not in stable and b not in stable:continue
   if all(bool(rows[a]>>v&1)==bool(rows[b]>>v&1) for v in range(n) if v not in (a,b)):
    cr,co=contract(rows,O,a);ck=canon(cr,co);need(ck in states,'child closure');child.add(ck);pairs.append((a,b))
  children[key]=child;arcs+=len(child)
  need(set(pairs)==set(map(tuple,r['stable_twin_contractions'])),'contraction list')
 terminal={};heights={}
 for key in sorted(states,key=lambda k:len(k[1])):
  ch=children[key];terminal[key]=set().union(*(terminal[c] for c in ch)) if ch else {key}
  need(len(terminal[key])==1,'nonunique terminal')
  heights[key]=len(key[1])-len(next(iter(terminal[key]))[1])
  for a,b in combinations(ch,2):
   need(bool(({a}|children[a])&({b}|children[b])),'one-step diamond')
   diamonds+=1
 byn={}
 for n in range(3,9):
  subset=[k for k in states if len(k[1])==n]
  byn[str(n)]={'states':len(subset),'terminal_births':sum(not children[k] for k in subset),'normal_form_order_histogram':dict(sorted(Counter(len(next(iter(terminal[k]))[1]) for k in subset).items())),'max_contraction_depth':max(heights[k] for k in subset)}
 out={'scope':'All 161 strong marked unit seeds through order8 with at least one exact stable unselected port, inherited census. Canonical quotient/DAG recomputed; no census completeness re-enumeration.','source_sha256':hashlib.sha256(SRC.read_bytes()).hexdigest(),'states':len(states),'distinct_reduction_arcs':arcs,'nontrivial_local_diamonds':diamonds,'unique_terminal_all':True,'levels':byn}
 print(json.dumps(out,sort_keys=True,indent=2));return out
if __name__=='__main__':run()
