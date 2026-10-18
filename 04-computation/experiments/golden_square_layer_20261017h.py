#!/usr/bin/env python3
"""Complete induced graph (A+uC+vC)+{0,v^2}, preserving eight ports.
Finite field only rejects; every retained difference gets an exact norm.
"""
import argparse,json,time
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import numpy as np
from pysat.solvers import Solver
from escaped_ports_20261017d import make,need,P
from golden_boundary_20261017c import PORT_EDGES,patterns


def run(mode):
 start=time.monotonic()
 old=json.loads(Path('04-computation/experiments/root_positive_transport_20261017g/golden_shared.json').read_text())
 f,z,j,alpha,ports,units,eta,conj,escape,reduce=make('sqrt3')
 bar=lambda a:f.apply(a,conj)
 Z=tuple(map(Q,old['contact_sum']));D=tuple(map(Q,old['quadratic_relation']))
 A=[tuple(map(Q,a)) for a in old['first_axis_coordinates']]
 C=[tuple(map(Q,a)) for a in old['spindle_coordinates']]
 v=(f.scale(Z,Q(-1,2)),f.scale(Z,Q(1,2)))
 def mul(x,y):
  a,b=x;c,e=y
  return f.add(f.mul(a,c),f.mul(D,f.mul(b,e))),f.add(f.mul(a,e),f.mul(b,c))
 def norm(x):
  a,b=x;aa,bb=bar(a),bar(b)
  return f.sub(f.mul(a,aa),f.mul(D,f.mul(b,bb))),f.sub(f.mul(b,aa),f.mul(a,bb))
 v2=mul(v,v);need(norm(v2)==(f.one,f.zero),'v squared is a unit')
 u=units[0];uc=[f.mul(u,c) for c in C];vc=[(f.mul(v[0],c),f.mul(v[1],c)) for c in C]
 base=[(f.add(f.add(A[a],uc[b]),vc[c][0]),vc[c][1]) for a,b,c in old['source_addresses']]
 coords=list(base);index={x:i for i,x in enumerate(coords)};second=[]
 for a,b in base:
  x=(f.add(a,v2[0]),f.add(b,v2[1])) if mode=='rung' else ((a,f.scale(b,-1)) if mode=='conjugate' else mul(v,(a,b)))
  if x not in index:index[x]=len(coords);coords.append(x)
  second.append(index[x])
 need(len(set(second))==len(base),'second copy is injective')
 if mode=='rung':need(len(coords)==2*len(base),'both translated layers distinct')
 dn,im=reduce(D);need(im==0,'real discriminant at conjugation-compatible place')
 rr=np.array([[*reduce(a),*reduce(b)] for a,b in coords],dtype=np.int64)
 print('coordinates',len(coords),'seconds',round(time.monotonic()-start,2),flush=True)
 edges=[];cache={};candidates=0
 for a in range(len(coords)):
  ds=(rr[a+1:]-rr[a])%P
  ok=((ds[:,0]**2+ds[:,1]**2-dn*(ds[:,2]**2+ds[:,3]**2))%P==1)&((ds[:,3]*ds[:,0]-ds[:,2]*ds[:,1])%P==0)
  for off in np.flatnonzero(ok):
   b=a+1+int(off);candidates+=1
   delta=(f.sub(coords[b][0],coords[a][0]),f.sub(coords[b][1],coords[a][1]))
   if delta not in cache:cache[delta]=norm(delta)==(f.one,f.zero)
   if cache[delta]:edges.append((a,b))
  if a%2000==0:print('sieve',a,'edges',len(edges),'norms',len(cache),'seconds',round(time.monotonic()-start,2),flush=True)
 n=len(base);e0=set(map(tuple,old['edges']));full=set(edges)
 need({e for e in full if e[1]<n}==e0,'first layer complete inherited graph')
 second_edges={tuple(sorted((second[a],second[b]))) for a,b in e0}
 need(second_edges<=full,'all rotated or translated copy edges')
 if mode=='rung':
  need({(a-n,b-n) for a,b in full if a>=n}==e0,'translated layer complete inherited graph')
  need(all((a,a+n) in full for a in range(n)),'all v squared unit rungs')
 need([e for e in edges if e[1]<8]==sorted(tuple(sorted(e)) for e in PORT_EDGES),'literal port edges unchanged')
 cross=[e for e in edges if e[0]<n<=e[1]];extra=[e for e in cross if e[1]-e[0]!=n]
 print('graph',len(coords),len(edges),'cross',len(cross),'nonrung',len(extra),flush=True)
 pats=patterns(8,PORT_EDGES);witnesses=[];recipe=None
 if mode=='rung':
  need(not extra,'Cartesian matching layers')
  recipe=dict(kind='cyclic_colour_shift',modulus=6,shift=1)
 elif mode=='conjugate':
  fold=list(range(n))+[None]*(len(coords)-n)
  for a,b in enumerate(second):
   if fold[b] is not None:need(fold[b]==a,'overlap fixed pointwise by sigma')
   fold[b]=a
  need(all(tuple(sorted((fold[a],fold[b]))) in e0 for a,b in edges),'graph fold preserves every edge')
  recipe=dict(kind='fold_to_inherited_graph',fold=fold)
 if recipe is not None:
  for pat,col0 in zip(pats,old['witnesses']):
   col=list(map(int,col0))
   if mode=='rung':col=col+[(c+1)%6 for c in col]
   else:col=[col[x] for x in fold]
   need(col[:8]==list(pat) and all(col[a]!=col[b] for a,b in edges),'recipe reconstructs boundary witness')
 else:
  with Solver(name='g3') as sol:
   for a in range(len(coords)):sol.add_clause([6*a+c+1 for c in range(6)])
   for a,b in edges:
    for c in range(6):sol.add_clause([-6*a-c-1,-6*b-c-1])
   for i,pat in enumerate(pats):
    need(sol.solve(assumptions=[6*a+c+1 for a,c in enumerate(pat)]),'blocked pattern requires separately audited UNSAT certificate')
    model={x for x in sol.get_model() if x>0}
    col=[next(c for c in range(6) if 6*a+c+1 in model) for a in range(len(coords))];col[:8]=list(pat)
    need(col[:8]==list(pat) and all(col[a]!=col[b] for a,b in edges),'full witness verified')
    witnesses.append(''.join(map(str,col)))
    if i%100==0:print('boundary',i,'seconds',round(time.monotonic()-start,2),flush=True)
 with Solver(name='g3') as sol:
  for a in range(len(coords)):sol.add_clause([4*a+c+1 for c in range(4)])
  for a,b in edges:
   for c in range(4):sol.add_clause([-4*a-c-1,-4*b-c-1])
  four=None
  if sol.solve():
   model={x for x in sol.get_model() if x>0};four=[next(c for c in range(4) if 4*a+c+1 in model) for a in range(len(coords))]
   need(all(four[a]!=four[b] for a,b in edges),'four colour witness')
 return dict(status='FINITE-EXACT, complete induced graph and positive boundary witnesses',vertices=len(coords),edges=edges,
  cross_edges=cross,nonrung_cross_edges=extra,sieve_prime=P,sieve_candidates=candidates,distinct_exact_norms=len(cache),
  patterns=pats,witnesses=witnesses,witness_recipe=recipe,four_colour=four,source=('G+{0,v^2}' if mode=='rung' else ('G union sigma(G)' if mode=='conjugate' else 'G union vG'))+'; G=A+uC+vC; same eight literal ports',mode=mode,second_copy=second,
  address_encoding='first inherited source_addresses; second_copy gives the index of each transformed address',
  v_squared=[[str(c) for c in a] for a in v2])

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);p.add_argument('--mode',choices=['rung','rotate','conjugate'],default='rung');a=p.parse_args()
 r=run(a.mode);Path(a.output).write_text(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
 print(json.dumps({k:v for k,v in r.items() if k not in ('edges','cross_edges','nonrung_cross_edges','patterns','witnesses','four_colour','v_squared','second_copy','witness_recipe')}))
