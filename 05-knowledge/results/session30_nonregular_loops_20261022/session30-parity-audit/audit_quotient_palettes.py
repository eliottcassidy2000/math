"""Independent signed vertex-orbit and per-edge palette audit."""
import json,itertools,collections
import numpy as np
from pathlib import Path
BASE=Path(__file__).resolve().parent.parent
cores=json.loads((BASE/'session30-parity-basis/parity-kernels.json').read_text())
quotient=json.loads((BASE/'session30-parity-basis/quotient-kernels.json').read_text())
def canonical(core,idx):
 edges=[tuple(e) for e in core['edges']];target=collections.Counter(edges);candidates=[]
 for p in itertools.permutations(range(core['V'])):
  transformed=collections.defaultdict(list);sign=1
  for (a,b),k in zip(edges,idx):
   u,v=p[a],p[b]
   if u>v:
    u,v=v,u
    if k>=2:sign=-sign
   transformed[(u,v)].append(k)
  if collections.Counter({e:len(ks) for e,ks in transformed.items()})!=target:continue
  for ks in transformed.values():ks.sort(reverse=True)
  label=tuple(transformed[e].pop() for e in edges)
  candidates.append((label,sign))
 minimum=min(k for k,s in candidates);signs={s for k,s in candidates if k==minimum}
 return minimum,0 if len(signs)>1 else next(iter(signs))
def reduce(core):
 out=collections.defaultdict(int)
 for term in core['terms']:
  idx,sign=canonical(core,term['indices']);out[idx]+=sign*term['numerator']
 return {idx:v for idx,v in out.items() if v}
def is_bridge(core,j):
 if core['edges'][j][0]==core['edges'][j][1]:return False
 seen={0}
 while True:
  old=seen.copy()
  for i,(a,b) in enumerate(core['edges']):
   if i==j:continue
   if a in seen:seen.add(b)
   if b in seen:seen.add(a)
  if seen==old:return len(seen)<core['V']
summary=[];bridge_union=set();w1bridge_union=set();loop_union=set()
for core,q in zip(cores,quotient):
 ours=reduce(core);candidate=reduce(q)
 assert ours==candidate
 assert len(ours)==len(q['terms'])
 row={'weight':core['weight'],'edges':core['edges'],'quotient_terms':len(ours),'palettes':[]}
 for j,(a,b) in enumerate(core['edges']):
  palette={t['indices'][j] for t in core['terms']};bridge=is_bridge(core,j)
  row['palettes'].append({'edge':[a,b],'bridge':bridge,'basis':[core['basis'][k] for k in sorted(palette)]})
  if a==b:loop_union|=palette;assert palette<={0,1}
  if bridge:
   bridge_union|=palette;assert palette<={0,3}
   if core['weight']==1:w1bridge_union|=palette;assert palette<={0}
 summary.append(row)
assert sum(r['quotient_terms'] for r in summary)==263
# Exact integer contraction tests on arbitrary symmetric/skew matrices, not
# just rational functions of a common M. This tests the quotient's actual scope.
rng=np.random.default_rng(20261022);integer_checks=0
for n in (2,3,4):
 random=[rng.integers(-2,3,size=(n,n),dtype=np.int64) for _ in range(4)]
 B=[random[0]+random[0].T,random[1]+random[1].T,random[2]-random[2].T,random[3]-random[3].T]
 def numerator(core):
  spec=','.join('abcd'[a]+'abcd'[b] for a,b in core['edges'])+'->'
  return sum(t['numerator']*int(np.einsum(spec,*[B[k] for k in t['indices']],optimize=False)) for t in core['terms'])
 for core,q in zip(cores,quotient):
  assert numerator(core)==numerator(q)
  integer_checks+=1
print(json.dumps({'exact_integer_contraction_checks':integer_checks,'status' :'PASS','cores':31,'quotient_terms':263,'loop_palette':sorted(loop_union),'bridge_palette':sorted(bridge_union),'weight1_bridge_palette':sorted(w1bridge_union),'basis':['G','H','O','P'],'summary':summary},indent=2))
