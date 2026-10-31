#!/usr/bin/env python3
"""Compile rooted port orders through disjoint primitive path permutations."""
import itertools,json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'square_sum_spectrum_overnight_20261031f'))
from edge_atlas import Atlas,edges,edge,need
from native_edge_atlas import build

class Orders:
 def __init__(self,cube):
  self.cube=cube
  p=list(cube.reference);i=p.index(1);p=p[i:]+p[:i]
  if p[1]!=cube.anchor:p=[p[0]]+p[:0:-1]
  need(p[1]==cube.anchor,'canonical rooted orientation')
  self.reference=p;self.pos={v:j for j,v in enumerate(p)};self.regions=[];self.owner={}
  for bit,(removed,added) in cube.moves.items():
   adj={}
   for u,v in removed|added:adj.setdefault(u,set()).add(v);adj.setdefault(v,set()).add(u)
   todo=set(adj)
   while todo:
    seed=min(todo);stack=[seed];vs=set()
    while stack:
     v=stack.pop()
     if v in vs:continue
     vs.add(v);stack.extend(adj[v]-vs)
    todo-=vs
    need(len(vs) in (4,6) and all(len(adj[v])==2 for v in vs),'primitive variable circuit')
    lo,hi=min(self.pos[v] for v in vs),max(self.pos[v] for v in vs)
    old=self.reference[lo:hi+1];S=set(old);es=edges(old)
    local_removed={e for e in removed if set(e)<=vs};local_added={e for e in added if set(e)<=vs}
    need(len(local_removed)==len(vs)//2 and len(local_added)==len(vs)//2,'primitive matching sizes')
    es=(es-local_removed)|local_added;aa={v:[] for v in old}
    for u,v in es:aa[u].append(v);aa[v].append(u)
    need(set(aa)==S and len(aa[old[0]])==len(aa[old[-1]])==1,'same ordered ends')
    new=[old[0]];previous=None
    while new[-1]!=old[-1]:
     choices=[v for v in aa[new[-1]] if v!=previous]
     need(len(choices)==1,'local path degree');previous,v=new[-1],choices[0]
     need(v not in new,'local simple path');new.append(v)
    need(len(new)==len(old) and set(new)==S and new!=old,'same local set, changed order')
    region=len(self.regions)
    for v in old:need(v not in self.owner,'disjoint reorder intervals');self.owner[v]=region
    self.regions.append({'bit':bit,'lo':lo,'hi':hi,'old':old,'new':new})
  need(all(1 not in r['old'] and cube.anchor not in r['old'] for r in self.regions),'orientation anchor is outside moving regions')
 def observed(self,bits,labels):
  wanted=set(labels);return [v for v in self.cube.walk(bits) if v in wanted]
 def compile(self,wanted):
  need(len(set(wanted))==len(wanted) and all(v in self.pos for v in wanted),'distinct native port labels')
  S=set(wanted);reference=[v for v in self.reference if v in S]
  def key(v):return ('r',self.owner[v]) if v in self.owner else ('v',v)
  def compress(seq):return [k for k,_ in itertools.groupby(map(key,seq))]
  if compress(wanted)!=compress(reference):return None
  pins={}
  for ri in {self.owner[v] for v in wanted if v in self.owner}:
   r=self.regions[ri];local=[v for v in wanted if self.owner.get(v)==ri]
   allowed=[z for z,name in [(0,'old'),(1,'new')] if [v for v in r[name] if v in S]==local]
   if not allowed:return None
   if len(allowed)==1:
    bit=r['bit'];z=allowed[0]
    if bit in pins and pins[bit]!=z:return None
    pins[bit]=z
  return {'pins':pins,'free_bits':self.cube.bits-len(pins),'one_code':''.join(str(pins.get(b,0)) for b in range(self.cube.bits))}

def main():
 a=Atlas();records=[];tests=0
 for n,word in [(2812,'A'),(2861,'A'),(70312,'A'),(70312,'AA'),(137813,'AC')]:
  c=a.cycle(n,word);o=Orders(c)
  # Pick a pair whose relative order changes in each of two different bits.
  picks=[];used=set()
  for r in o.regions:
   if r['bit'] in used:continue
   pp={v:j for j,v in enumerate(r['new'])}
   pair=next((list(p) for p in itertools.combinations(r['old'],2) if pp[p[0]]>pp[p[1]]),None)
   need(pair is not None,'local ordering change');picks+=pair;used.add(r['bit'])
   if len(used)==2:break
  reference=[v for v in o.reference if v in picks];compiled=o.compile(reference)
  need(compiled is not None and len(compiled['pins'])==2,'four reference ports pay exactly two bits')
  # Every permutation of these four ports; compare all four relevant codes.
  bitlist=sorted(used);actual={}
  for values in itertools.product('01',repeat=2):
   bits=['0']*c.bits
   for b,z in zip(bitlist,values):bits[b]=z
   order=tuple(o.observed(''.join(bits),picks));actual[order]=actual.get(order,0)+1
  for order in itertools.permutations(picks):
   out=o.compile(order);need((out is not None)==(order in actual),'complete four-port order relation')
   if out is not None:
    need(2**(2-len(out['pins']))==actual[order],'exact local order fibre')
    need(o.observed(out['one_code'],picks)==list(order),'literal returned order witness')
   tests+=1
  # Marking one label in each region cannot reveal any switch.
  singles=[r['old'][1] for r in o.regions]
  observed=[v for v in o.reference if v in set(singles)]
  need(o.compile(observed)['free_bits']==c.bits,'one port per region is order-blind')
  records.append({'N':n,'schedule':word,'bits':c.bits,'primitive_regions':len(o.regions),'four_ports':reference,'free_bits':compiled['free_bits']})
 # Native port-preserving two-layer lift.
 _,c,*_=build(270112,'');o=Orders(c)
 r=o.regions[0];pp={v:j for j,v in enumerate(r['new'])};pair=next(list(p) for p in itertools.combinations(r['old'],2) if pp[p[0]]>pp[p[1]])
 P=[1,8,17,19,6,10,15,21,4,12,24];wanted=[v for v in o.reference if v in set(P+pair)]
 out=o.compile(wanted);need(out is not None and len(out['pins'])==1,'native prefix free, one inversion pair paid')
 need(o.observed(out['one_code'],wanted)==wanted,'native order witness')
 print(json.dumps({'status':'passed','complete_four_port_permutations':tests,'standard_cases':records,'native_case':{'N':270112,'bits':c.bits,'free_bits':out['free_bits'],'prefix':P,'additional_pair':pair},'scope':'Exact rooted order relation in the authenticated family, graph-complete only after primitive refinement. Any k reference ports cost at most floor(k/2) order bits.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
