#!/usr/bin/env python3
import json
from pathlib import Path
from collections import Counter
from sympy import Matrix
ROOT=next(p for p in [*Path(__file__).resolve().parents, Path(__file__).resolve().parent.parent/'math'] if (p/'AGENTS.md').is_file() and (p/'01-canon').is_dir())
d=json.loads((ROOT/'05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020/census_n5.json').read_text())
B=Counter({(u,v):w for u,v,w in d['blue_edges']});Q=Counter({(u,v):w for u,v,w in d['paired_black_edges']});n=len(d['nodes'])
def L(es):
 a=[[0]*n for _ in range(n)]
 for (u,v),w in es.items():
  if u!=v:a[u][u]+=w;a[v][v]+=w;a[u][v]-=w;a[v][u]-=w
 return Matrix([row[1:] for row in a[1:]])
def degree(es):
 ds=[0]*n
 for (u,v),w in es.items():ds[u]+=w;ds[v]+=w
 return ds
q=Q.copy()
for e in [(1,4),(2,3)]:q[e]-=1
for e in [(1,3),(2,4)]:q[e]+=1
if degree(Q)!=degree(q):raise RuntimeError('Degree preservation failed')
true=int((L(B)+2*L(Q)).det());fake=int((L(B)+2*L(q)).det())
r={'n':5,'removed_paired_black_edges':[[1,4],[2,3]],'added_paired_black_edges':[[1,3],[2,4]],'paired_black_degrees':degree(Q),'literal_tree_count_actual':true,'literal_tree_count_rewired':fake,'same_blue_graph':True,'same_nodes_fibres_colours_loops_degrees_charges':True,'not_a_tiling_realizability_claim':True}
if true==fake:raise RuntimeError('Rewire failed to change tree count')
print(json.dumps(r,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n')
