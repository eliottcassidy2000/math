#!/usr/bin/env python3
"""A 19-state quotient solves the five unsolved template patterns exactly.
Every SAT model is replayed on the full 369-vertex port-split graph.
"""
import argparse,hashlib,json
from itertools import combinations
from pathlib import Path
from pysat.solvers import Solver
import golden_boundary_20261017c as G

def main(path):
 raw=Path(path).read_bytes();d=json.loads(raw);points,edges=G.residue_graph()
 D={1,2,4,6,9,10,13,15,17,18};answers=[]
 for pat in d['remaining_template_patterns']:
  witness=None;hits=0;tested=0
  for branch_no,b in enumerate(d['branches']):
   for aa,bb,t in G.projection_params():
    if t:continue
    def v(i,c):return 6*i+c+1
    cl=[]
    for i in range(19):
     cl.append([v(i,c) for c in range(6)])
     for c,e in combinations(range(6),2):cl.append([-v(i,c),-v(i,e)])
    for i in range(19):
     for de in D:
      for c in range(6):cl.append([-v(i,c),-v((i+de)%19,c)])
    for c,(x,y) in zip(pat,b['points']):
     for de in D:cl.append([-v((aa*x+bb*y+de)%19,c)])
    with Solver(name='g3',bootstrap_with=cl) as so:
     tested+=1
     if not so.solve():continue
     hits+=1;model=set(so.get_model())
     circle=[next(c for c in range(6) if v(i,c) in model) for i in range(19)]
    G.need(all(circle[i]!=circle[(i+de)%19] for i in range(19) for de in D),'complete 19-state graph replay independent of solver')
    colour=[circle[(aa*x+bb*y)%19] for x,y in points]+pat
    G.need(all(colour[a]!=colour[b] for a,b in edges),'all 3610 residue edges replay')
    for k,p in enumerate(b['points']):
     G.need(all(colour[a]!=pat[k] for a,z in enumerate(points) if G.q(G.sub(z,p))==1),'all residue neighbours of each separate literal port replay')
    G.need(all(pat[a]!=pat[b] for a,b in G.PORT_EDGES),'literal eight-port edges replay')
    if witness is None:witness=dict(branch=branch_no,projection=[aa,bb],circle=circle,colours=colour)
  G.need(witness is not None,'every remaining pattern extends the complete field carrier')
  answers.append(dict(pattern=pat,projection_tests=tested,positive_templates=hits,witness=witness))
 G.need(len(d['cyclotomic_templates'])+len(answers)==570,'full boundary relation preserved')
 return dict(status='FINITE-EXACT witnesses for all five patterns; with companion field proof all 570 extend',source_sha256=hashlib.sha256(raw).hexdigest(),answers=answers,checks=G.checks)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--input',required=True);a=p.parse_args();print(json.dumps(main(a.input),sort_keys=True,separators=(',',':')))
