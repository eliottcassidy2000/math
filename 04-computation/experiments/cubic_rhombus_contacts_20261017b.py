#!/usr/bin/env python3
"""Literal rhombus completion and exact additive-colouring CSP.
Candidate attachments inherit unit edges by z=p+q-r; source ports fixed.
"""
import argparse,hashlib,json
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations,product
from math import lcm
from pathlib import Path
import numpy as np
from pysat.solvers import Solver
import cubic_hub_coupling_20261017a as H

checks=0
def need(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)

def setup():
 raw=Path('04-computation/experiments/anchor_hub_20261017a/hub_valuation_grid_four.json').read_bytes();old=json.loads(raw)
 V=[tuple(F(c) for c in x) for x in old['base_vertices']];den=lcm(*(v.denominator for x in V for v in x));VI=[tuple(int(v*den) for v in x) for x in V]
 zero=(0,)*24;one=VI[1]+(0,)*16
 points={0:zero,1:one};addr={0:(0,0,0),1:(1,0,0)};index={zero:0,one:1};nextid=8
 for a in product(range(len(V)),repeat=3):
  p=sum((VI[i] for i in a),())
  if p in index:continue
  points[nextid]=p;addr[nextid]=a;index[p]=nextid;nextid+=1
 need(nextid==old['vertices'],'same literal ordered grid with the six non-grid ports reserved')
 E=[tuple(e) for e in old['edges'] if all(i in points for i in e)]
 adjacency=defaultdict(set)
 for i,j in E:adjacency[i].add(j);adjacency[j].add(i)
 candidates={}
 for r,ns in adjacency.items():
  for p,q in combinations(sorted(ns),2):
   z=tuple(a+b-c for a,b,c in zip(points[p],points[q],points[r]))
   if z in index:continue
   item=candidates.setdefault(z,dict(neighbours=set(),witnesses=[]))
   item['neighbours'].update((p,q));item['witnesses'].append((r,p,q))
 return old,V,VI,den,points,addr,E,candidates,hashlib.sha256(raw).hexdigest()

class ColourCSP:
 def __init__(self,V,E,addr):
  self.n=len(V);self.nextvar=3*3*self.n;self.clauses=[];self.addr=addr;self.attachments=0
  for position in range(3):
   for i in range(self.n):self.onehot(position*self.n+i)
  for position in range(3):self.clauses.append([self.var(position*self.n,0)])
  self.clauses.append([self.var(1,1)])
  unique=set()
  for i,j in E:
   terms=self.difference(addr[i],addr[j]);unique.add(tuple(sorted(terms.items())))
  for terms in sorted(unique):self.nonzero(dict(terms))
  self.solver=Solver(name='g3',bootstrap_with=self.clauses)
  self.constraint_count=len(unique)
 def var(self,index,c):return 3*index+c+1
 def onehot(self,index):
  self.clauses.append([self.var(index,c) for c in range(3)])
  self.clauses.extend([[-self.var(index,c),-self.var(index,d)] for c,d in combinations(range(3),2)])
 def difference(self,a,b):
  d=defaultdict(int)
  for p in range(3):d[p*self.n+a[p]]+=1;d[p*self.n+b[p]]-=1
  return {i:c%3 for i,c in d.items() if c%3}
 def nonzero(self,terms):
  keys=sorted(terms)
  for values in product(range(3),repeat=len(keys)):
   if sum(terms[i]*c for i,c in zip(keys,values))%3==0:
    self.clauses.append([-self.var(i,c) for i,c in zip(keys,values)])
 def solve(self):
  if not self.solver.solve():return None
  m=set(self.solver.get_model());base=[next(c for c in range(3) if self.var(i,c) in m) for i in range(3*self.n)]
  colours={i:sum(base[p*self.n+a[p]] for p in range(3))%3 for i,a in self.addr.items()}
  return base,colours
 def attach(self,neighbours):
  index=3*self.n+self.attachments;self.attachments+=1;start=len(self.clauses);self.onehot(index)
  for i in sorted(neighbours):
   a=self.addr[i];terms={p*self.n+a[p]:1 for p in range(3)};terms[index]=2;self.nonzero(terms)
  for clause in self.clauses[start:]:self.solver.add_clause(clause)

def run(limit):
 old,V,VI,den,P,addr,E,candidates,source_hash=setup()
 f=old['product_base_colours'];fixed={i:sum(f[t] for t in a)%3 for i,a in addr.items()}
 fixed_bad=[z for z,c in candidates.items() if len({fixed[i] for i in c['neighbours']})==3]
 csp=ColourCSP(V,E,addr);need(csp.solve() is not None,'inherited additive family satisfiable')
 selected=[];remaining=dict(candidates);status='LIMIT';models=[]
 for iteration in range(limit+1):
  model=csp.solve()
  if model is None:status='ADDITIVE_UNSAT';break
  base,colours=model;need(all(colours[i]!=colours[j] for i,j in E),'all grid edges replay current additive colouring')
  bad=[z for z,item in remaining.items() if len({colours[i] for i in item['neighbours']})==3]
  if not bad:status='ADDITIVE_SURVIVES_ALL_COMPLETIONS';models=base;break
  if iteration==limit:break
  z=max(bad,key=lambda z:(len(remaining[z]['neighbours']),z));item=remaining.pop(z)
  # Each reported contact is proved via the original two-edge witness and
  # independently checked in the cubic quotient for the selected attachment.
  for i in item['neighbours']:
   diff=tuple(F(a-b,den) for a,b in zip(z,P[i]));need(H.unit_direct(tuple(diff[8*k:8*k+8] for k in range(3))), 'selected contact exact cubic norm')
  selected.append(dict(point=[str(F(v,den)) for v in z],neighbours=sorted(item['neighbours']),witnesses=item['witnesses']))
  csp.attach(item['neighbours'])
 out=dict(status=status,source_sha256=source_hash,vertices=old['vertices'],grid_edges=len(E),candidate_count=len(candidates),
          candidates_obstructing_saved_colouring=len(fixed_bad),selected=selected,remaining_additive_base_colours=models,
          additive_constraint_count=csp.constraint_count,clauses=len(csp.clauses),checks=checks)
 csp.solver.delete();return out

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=128);a=p.parse_args()
 print(json.dumps(run(a.limit),sort_keys=True,separators=(',',':')))
