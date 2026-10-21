#!/usr/bin/env python3
import json, math
from pathlib import Path
from collections import Counter
from sympy import Matrix, Rational
ROOT=next(p for p in [*Path(__file__).resolve().parents, Path(__file__).resolve().parent.parent/'math'] if (p/'AGENTS.md').is_file() and (p/'01-canon').is_dir())
SRC=ROOT/'05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020'
records=[]
for n in range(3,7):
 d=json.loads((SRC/f'census_n{n}.json').read_text());N=len(d['nodes']);L=[[0]*N for _ in range(N)]
 for name in ['blue_edges','black_edges']:
  for u,v,w in d[name]:
   if u!=v:L[u][u]+=w;L[v][v]+=w;L[u][v]-=w;L[v][u]-=w
 A=Matrix([row[1:] for row in L[1:]]);inv=A.inv();cols=[[0]*N]+[[0]+list(inv[:,u-1]) for u in range(1,N)]
 checked=0;hostiles=[];minimum=None
 for u in range(N):
  if d['nodes'][u]['SC']:continue
  for v in range(N):
   if u==v:continue
   phi=[cols[u][i]-cols[v][i] for i in range(N)];den=math.lcm(*(int(x.q) if hasattr(x,'q') else 1 for x in phi));R=phi[u]-phi[v]
   if den%2:raise RuntimeError('NS divisor unexpectedly odd')
   if Matrix(L)*Matrix(phi)!=Matrix([int(i==u)-int(i==v) for i in range(N)]):raise RuntimeError('Voltage equation failed')
   checked+=1
   if int(R.q)%2==1 and len(hostiles)<3:hostiles.append({'source':u,'sink':v,'common_voltage_denominator':den,'effective_resistance':str(R),'potential_with_even_denominator':next({'node':i,'phi':str(x)} for i,x in enumerate(phi) if hasattr(x,'q') and int(x.q)%2==0)})
 r={'n':n,'ordered_NS_source_pairs_checked':checked,'odd_effective_resistance_denominator_hostiles':hostiles}
 if n==5:
  full=Matrix(L);boundary=[0,2];interior=[i for i in range(N) if i not in boundary]
  ext=-full.extract(interior,interior).inv()*full.extract(interior,boundary)
  schur=full.extract(boundary,boundary)+full.extract(boundary,interior)*ext
  if schur!=Rational(5137,6092)*Matrix([[1,-1],[-1,1]]):raise RuntimeError('Exact Kron mismatch')
  if ext[interior.index(1),1]!=Rational(11327,12184):raise RuntimeError('Harmonic extension mismatch')
  r['kron_boundary']=boundary;r['kron_conductance']='5137/6092';r['harmonic_extension_node1_source_coefficient']='11327/12184'
 records.append(r);print(json.dumps(r),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(records,indent=2)+'\n')
