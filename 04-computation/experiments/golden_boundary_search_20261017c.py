#!/usr/bin/env python3
"""Bounded six-colour search on a residue graph with separate literal ports.
SAT witnesses are independently replayed; UNKNOWN/UNSAT are not forcing claims.
"""
import argparse,hashlib,json
from pathlib import Path
from itertools import combinations
from pysat.solvers import Solver
import golden_boundary_20261017c as G

def run(path,budget):
 raw=Path(path).read_bytes();old=json.loads(raw);points,base_edges=G.residue_graph();answers=[]
 for pattern in old['remaining_template_patterns']:
  attempts=[];answer=None
  for branch in old['branches']:
   edges=set(base_edges)
   for k,p in enumerate(branch['points']):
    for a,z in enumerate(points):
     if G.q(G.sub(z,p))==1:edges.add((a,361+k))
   edges.update((361+a,361+b) for a,b in G.PORT_EDGES)
   def v(a,c):return 6*a+c+1
   clauses=[]
   for a in range(369):
    clauses.append([v(a,c) for c in range(6)])
    clauses.extend([[-v(a,c),-v(a,d)] for c,d in combinations(range(6),2)])
   for a,b in edges:
    for c in range(6):clauses.append([-v(a,c),-v(b,c)])
   clauses.extend([[v(361+a,c)] for a,c in enumerate(pattern)])
   with Solver(name='g4',bootstrap_with=clauses) as solver:
    solver.set_phases([v(a,G.TABLE[(x+y)%19]) for a,(x,y) in enumerate(points)])
    solver.conf_budget(budget);res=solver.solve_limited(expect_interrupt=True)
    attempts.append(dict(branch=old['branches'].index(branch),answer=res,stats=solver.accum_stats()))
    if res:
     model=set(solver.get_model());col=[next(c for c in range(6) if v(a,c) in model) for a in range(369)]
     G.need(col[361:]==pattern,'all eight original port colours retained')
     G.need(all(col[a]!=col[b] for a,b in edges),'every quotient and literal edge replayed outside solver')
     answer=dict(branch=old['branches'].index(branch),colours=col,edges=len(edges));break
  answers.append(dict(pattern=pattern,attempts=attempts,witness=answer))
 return dict(status='FINITE-EXACT positive witnesses; bounded search failures do not imply geometric forcing',source_sha256=hashlib.sha256(raw).hexdigest(),conflict_budget=budget,answers=answers,checks=G.checks)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--budget',type=int,default=10000);a=p.parse_args()
 print(json.dumps(run(a.input,a.budget),sort_keys=True,separators=(',',':')))
