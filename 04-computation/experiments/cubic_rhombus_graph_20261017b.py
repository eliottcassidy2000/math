#!/usr/bin/env python3
"""Full induced graph of the first rhombus closure, in degree-24 arithmetic.
The five old ports outside B(u) attach only to original ports, by the
proved cyclotomic-conjugate-circle classification in the companion note.
"""
import argparse,hashlib,json
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import numpy as np
import sympy as sp
from pysat.solvers import Solver
import cubic_hub_coupling_20261017a as H
import cubic_rhombus_contacts_20261017b as R
checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)

def graph(points,den):
 ids=sorted(points);P=[points[i] for i in ids];filters=[]
 for p in sp.primerange(503,10000):
  if den%p==0 or p%4!=1 or any(pow(d,(p-1)//2,p)!=1 for d in (7,11)):continue
  roots=[w for w in range(p) if (w**3+2*w*w-4*w-4)%p==0]
  if not roots:continue
  r=next(t for t in range(p) if t*t%p==7);s=next(t for t in range(p) if t*t%p==11);ii=next(t for t in range(p) if t*t%p==p-1);w=roots[0]
  if (w+3+ii*s)%p==0 or (w+3-ii*s)%p==0:continue
  ev=[]
  for sign in (1,-1):
   jj=sign*ii*s%p;u=(w+3-jj)*pow((w+3+jj)%p,-1,p)%p
   basis=[pow(u,k,p)*(r if j&1 else 1)*(s if j&2 else 1)*(sign*ii if j&4 else 1)%p for k in range(3) for j in range(8)]
   ev.append([sum(c*v for c,v in zip(x,basis))%p for x in P])
  filters.append((p,np.array(ev,dtype=np.int64).T))
  if len(filters)==3:break
 need(len(filters)==3,'three exact field homomorphisms')
 E=[];cache={};survivors=0
 for i in range(len(P)):
  js=np.arange(i+1,len(P),dtype=np.int64)
  for p,ev in filters:js=js[((ev[js,0]-ev[i,0])*(ev[js,1]-ev[i,1])-den*den%p)%p==0]
  survivors+=len(js)
  for j in js:
   delta=tuple(a-b for a,b in zip(P[i],P[int(j)]));key=min(delta,tuple(-a for a in delta))
   if key not in cache:
    c=tuple(F(a,den) for a in key);cache[key]=H.unit_direct(tuple(c[k*8:k*8+8] for k in range(3)))
   if cache[key]:E.append((ids[i],ids[int(j)]))
 return E,dict(primes=[p for p,ev in filters],survivors=survivors,exact_displacement_checks=len(cache))

def colouring(N,E,q):
 clauses=[[q*i+c+1 for c in range(q)] for i in range(N)]
 clauses.extend([[-q*i-c-1,-q*j-c-1] for i,j in E for c in range(q)])
 with Solver(name='g3',bootstrap_with=clauses) as solver:
  if not solver.solve():return None
  model=set(solver.get_model());col=[next(c for c in range(q) if q*i+c+1 in model) for i in range(N)]
 need(all(col[i]!=col[j] for i,j in E),'proper global colouring replay')
 return col

def boundary(N,E,col):
 pe=[(i,j) for i,j in E if i<8 and j<8];port_neighbours={i:[] for i in range(8,N)}
 for i,j in E:
  if i<8<=j:port_neighbours[j].append(i)
 need(all(set(v)<={0,1,6} for v in port_neighbours.values()),'only three old ports can touch the new carrier')
 double=[i for i,v in port_neighbours.items() if len(v)==2]
 need(len(double)<=2 and all(len(v)<=2 for v in port_neighbours.values()),'at most two double-port contacts')
 patterns=[]
 def visit(t):
  n=len(t)
  if n==8:patterns.append(tuple(t));return
  for c in range(min(6,max(t,default=-1)+2)):
   if all(t[i]!=c for i,j in pe if j==n):visit(t+[c])
 visit([])
 for beta in patterns:
  unused=[c for c in range(6) if c not in {beta[i] for i in (0,1,6)}];single=unused[:2];rest=[c for c in range(6) if c not in single]
  found=False
  for pair in combinations(rest,2):
   palettes=[pair,tuple(c for c in rest if c not in pair),(single[0],),(single[1],)]
   assignment=list(beta)
   for i in range(8,N):
    options=[c for c in palettes[col[i]] if all(c!=beta[j] for j in port_neighbours[i])]
    if not options:break
    assignment.append(options[0])
   if len(assignment)==N:
    need(all(assignment[i]!=assignment[j] for i,j in E),'every six-colour boundary extension explicitly replayed');found=True;break
  need(found,'four-colour source-palette theorem supplies extension')
 return dict(patterns=len(patterns),removed=[],double_port_points=double)

def run():
 old,V,VI,den,P,addr,E0,candidates,source_hash=R.setup();P[6]=(3*den//4,)+(0,)*23
 index={p:i for i,p in P.items()};attachments=[];nextid=old['vertices']
 for z,item in sorted(candidates.items()):
  if z not in index:index[z]=nextid;P[nextid]=z;nextid+=1;attachments.append(dict(vertex=index[z],witness=item['witnesses'][0]))
 E,meta=graph(P,den);restricted={tuple(e) for e in E if e[0] in addr and e[1] in addr}
 need(restricted==set(E0),'independent degree-24 graph reproduces every inherited grid edge')
 # Complete old-port graph is inherited byte-for-byte. The external-port
 # classification proves these are all additional possible incidences.
 E=sorted(set(E)|{tuple(e) for e in old['edges'] if e[0]<8 and e[1]<8})
 setE=set(E)
 for z,item in candidates.items():
  for p in item['neighbours']:need(tuple(sorted((index[z],p))) in setE,'every generated rhombus contact present in the full induced graph')
 col=colouring(nextid,E,3);bound=3
 if col is None:col=colouring(nextid,E,4);bound=4
 if col is None:raise RuntimeError('neither declared colour upper bound certified')
 # Test the inherited additive colour on old grid, with all attachment
 # colours free and all actual new-new edges present.
 q=3;clauses=[[q*i+c+1 for c in range(q)] for i in range(nextid)]
 clauses.extend([[-q*i-c-1,-q*j-c-1] for i,j in E for c in range(q)])
 f=old['product_base_colours'];fixed={i:sum(f[t] for t in a)%3 for i,a in addr.items()}
 with Solver(name='g3',bootstrap_with=clauses) as solver:
  saved_extends=solver.solve(assumptions=[q*i+c+1 for i,c in fixed.items()])
  if saved_extends:
   model=set(solver.get_model());values=[fixed.get(i,next(c for c in range(q) if q*i+c+1 in model)) for i in range(nextid)]
   need(all(values[i]!=values[j] for i,j in E),'fixed additive colouring extends with new-new edges')
 b=boundary(nextid,E,col)
 return dict(status='FINITE-EXACT induced graph plus proved external-port classification',vertices=nextid,edges=E,colouring=col,colour_bound=bound,
             saved_additive_colouring_extends=saved_extends,saved_additive_extension=values if saved_extends else None,boundary=b,geometry=meta,source_sha256=source_hash,
             coordinates={str(i):list(p) for i,p in P.items()},denominator=den,attachments=attachments,checks=checks)

if __name__=='__main__':print(json.dumps(run(),sort_keys=True,separators=(',',':')))
