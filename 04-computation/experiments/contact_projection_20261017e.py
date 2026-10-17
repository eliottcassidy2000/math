#!/usr/bin/env python3
"""Port-fixing graph projections and independently lifted colour witnesses."""
import json,argparse
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
from pysat.solvers import Solver
from escaped_ports_20261017d import make,need,P
from golden_boundary_20261017c import PORT_EDGES,patterns


def run(path,kind):
 raw=json.loads(Path(path).read_text())
 f,z,j,alpha,ports,units,eta,conj,escape,reduce=make('sqrt3')
 A=[tuple(Q(c) for c in a) for a in raw['first_axis_coordinates']]
 M=[tuple(Q(c) for c in a) for a in raw['spindle_coordinates']]
 if kind=='straight':r=s=f.scale(eta,-1);shift=f.one
 else:r=f.apply(eta,conj);s=f.scalar(-1);shift=f.scale(eta,-1)
 norm=lambda a:f.mul(a,f.apply(a,conj))
 d=f.add(f.one,eta)
 need(norm(r)==norm(s)==f.one and norm(f.add(d,f.add(f.mul(r,shift),s)))==f.one,'all four direction families map to unit directions')
 points=list(ports);index={p:i for i,p in enumerate(points)};projection=[]
 for a,b,c in raw['source_addresses']:
  p=f.add(A[a],f.add(f.mul(r,M[b]),f.mul(s,M[c])))
  if p not in index:index[p]=len(points);points.append(p)
  projection.append(index[p])
 need(projection[:8]==list(range(8)),'literal ports fixed')
 residues=np.array([reduce(p) for p in points],dtype=np.int64);edges=[];cache={}
 for a in range(len(points)):
  delta=residues[a+1:]-residues[a]
  for off in np.flatnonzero((delta[:,0]**2+delta[:,1]**2)%P==1):
   b=a+1+int(off);v=f.sub(points[b],points[a])
   if v not in cache:cache[v]=norm(v)==f.one
   if cache[v]:edges.append((a,b))
 eset=set(edges)
 need(all(tuple(sorted((projection[a],projection[b]))) in eset for a,b in raw['edges']),'every source unit edge survives projection, none collapses')
 need([e for e in edges if e[1]<8]==sorted(tuple(sorted(e)) for e in PORT_EDGES),'induced original boundary unchanged')
 pats=patterns(8,PORT_EDGES);ws=[]
 with Solver(name='g3') as sol:
  for a in range(len(points)):
   sol.add_clause([6*a+c+1 for c in range(6)])
   for c in range(6):
    for b in range(c):sol.add_clause([-6*a-c-1,-6*a-b-1])
  for a,b in edges:
   for c in range(6):sol.add_clause([-6*a-c-1,-6*b-c-1])
  for pat in pats:
   need(sol.solve(assumptions=[6*a+c+1 for a,c in enumerate(pat)]),'all boundary patterns extend on projection image')
   model=set(v for v in sol.get_model() if v>0)
   col=[next(c for c in range(6) if 6*a+c+1 in model) for a in range(len(points))]
   need(col[:8]==list(pat) and all(col[a]!=col[b] for a,b in edges),'image colour witness')
   need(all(col[projection[a]]!=col[projection[b]] for a,b in raw['edges']),'independent lifted source witness')
   ws.append(''.join(map(str,col)))
 # Generic replay expects a valid relative-peeling receipt; compute it.
 adj=[set() for _ in points]
 for a,b in edges:adj[a].add(b);adj[b].add(a)
 alive=set(range(len(points)));peel=[]
 while True:
  batch=[a for a in sorted(alive) if a>=8 and len(adj[a]&alive)<=5]
  if not batch:break
  for a in batch:alive.remove(a);peel.append(a)
 return dict(status='FINITE-EXACT port-fixing projection; all570 image and lifted boundary witnesses',kind=kind,
  vertices=len(points),edges=edges,patterns=pats,witnesses=ws,source_projection=projection,
  source_vertices=raw['vertices'],source_edges=len(raw['edges']),
  coordinates=[[str(c) for c in p] for p in points],peeling_order=peel,relative_six_core=len(alive)-8)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--kind',choices=['straight','tilted'],required=True);p.add_argument('--output',required=True);a=p.parse_args()
 d=run(a.input,a.kind);Path(a.output).write_text(json.dumps(d,sort_keys=True,separators=(',',':'))+'\n')
 print(json.dumps(dict(kind=a.kind,vertices=d['vertices'],edges=len(d['edges']),relative_six_core=d['relative_six_core'])))
