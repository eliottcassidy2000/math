#!/usr/bin/env python3
"""Independent operator-block and kernel cycle checks; no producer imports."""
from collections import Counter
from itertools import product
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
HERE=Path(__file__).resolve().parent
GATES=0
def need(ok,msg):
 global GATES
 GATES+=1
 if not ok:raise RuntimeError(str(msg))
def main():
 source=HERE.parent/'pairings/normal.json';raw=source.read_bytes();audit=json.loads(raw);rows=[]
 for row in audit['rows']:
  s=row['source_data'];ordinary=s['ordinary_nodes'];merged=s['merged_nodes'];oc=[x['code'] for x in ordinary];mc=[x['code'] for x in merged]
  oi={c:i for i,c in enumerate(oc)};mi={c:i for i,c in enumerate(mc)};R={x['code']:x['converse'] for x in ordinary};d={x['code']:x['degree'] for x in ordinary}
  C={(a,b):v for a,b,v in s['C']};A={(a,b):v for a,b,v in s['A']};blue={(a,b):v for a,b,v in s['blueC']}
  fibres={c:[u for u in oc if min(u,R[u])==c] for c in mc};ch=[c for c in mc if len(fibres[c])==2]
  D={(a,b):C.get((a,b),0)-C.get((a,R[b]),0) for a in ch for b in ch}
  for x in oc:
   u=min(x,R[x]);need(d[x]==d[R[x]],'converse source masses')
   for v in mc:
    lhs=sum(C.get((x,y),0) for y in fibres[v]);need(F(lhs)==F(A.get((u,v),0),len(fibres[u])),'natural plus-block intertwining')
    need(F(lhs,d[x])==F(A.get((u,v),0),len(fibres[u])*d[x]),'stochastic natural plus block')
   for v in ch:
    lhs=C.get((x,v),0)-C.get((x,R[v]),0)
    rhs=0 if x==R[x] else (1 if x==u else -1)*D[u,v]
    need(lhs==rhs,'minus-block intertwining')
    need(F(lhs,d[x])==(0 if x==R[x] else (1 if x==u else -1)*F(D[u,v],d[u])),'stochastic minus block')
   for y in oc:
    need(C.get((x,y),0)==C.get((R[x],R[y]),0),'native converse commutation')
    if x!=R[x] or y!=R[y]:need(blue.get((x,y),0)==0,'blue confined to plus singleton sector')
  # Automorphism-invariant weighted vertex refinement, rebuilt from sparse A.
  colors=[(sum(A.get((u,v),0) for v in mc),A.get((u,u),0)) for u in mc]
  traces=[len(set(colors))]
  while True:
   signatures=[(colors[i],tuple(sorted(Counter((colors[j],A.get((u,v),0)) for j,v in enumerate(mc)).items()))) for i,u in enumerate(mc)]
   unique={x:k for k,x in enumerate(sorted(set(signatures)))}
   colors=[unique[x] for x in signatures];traces.append(len(set(colors)))
   if traces[-1]==traces[-2]:break
  need(traces==row['regimes']['A']['refinement_cell_counts'],'independent weighted vertex refinement')
  need(traces[-1]==len(mc) or row['n']==3,'vertex-rigidity certificate')
  # Independent component construction by repeated set union, not producer BFS.
  components=[{c} for c in ch]
  for a in ch:
   for b in ch:
    if a==b or not D[a,b]:continue
    ca=next(z for z in components if a in z);cb=next(z for z in components if b in z)
    if ca is not cb:ca.update(cb);components.remove(cb)
  need(sorted(map(sorted,components))==sorted(row['signed_support_components']),'independent signed components')
  fixed=Counter();burnside={2:0,3:0}
  for flips in product([False,True],repeat=len(components)):
   chosen=set().union(*(z for z,b in zip(components,flips) if b)) if any(flips) else set()
   p={x:(R[x] if min(x,R[x]) in chosen else x) for x in oc}
   for x in oc:
    for y in oc:need(C.get((x,y),0)==C.get((p[x],p[y]),0),'component action on ordinary matrix')
   f=sum(p[x]==x for x in oc);fixed[f]+=1
   seen=set();cycles=[]
   for x in oc:
    if x in seen:continue
    length=0;y=x
    while y not in seen:seen.add(y);length+=1;y=p[y]
    cycles.append(length)
   need(all(k in (1,2) for k in cycles),'kernel is involutive')
   for q in burnside:burnside[q]+=q**len(cycles)
  need(dict(fixed)=={int(k):v for k,v in row['physical_kernel_fixed_polynomial'].items()},'kernel fixed polynomial')
  single=sum(x==R[x] for x in oc)
  for q,total in burnside.items():
   expected=q**single
   for comp in components:expected*=q**(2*len(comp))+q**len(comp)
   need(total==expected,'cycle-index coloring numerator')
   need(total%(1<<len(components))==0,'Burnside integer')
  rows.append({'n':row['n'],'plus_dimension':len(mc),'minus_dimension':len(ch),
               'weighted_refinement':traces,'signed_component_sizes':sorted(map(len,components)), 'kernel_fixed_histogram':dict(sorted(fixed.items()))})
 return {'status':'PASS','always_active_checks':GATES,'scope':'operator intertwining and component-kernel cycle consequences on audited full matrices; not a new n7 census',
         'rows':rows,'input_sha256':hashlib.sha256(raw).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
