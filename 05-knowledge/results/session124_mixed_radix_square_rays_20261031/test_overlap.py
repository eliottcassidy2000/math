#!/usr/bin/env python3
"""Exact finite controls for the piece-overlap parity universal property."""
import itertools,json
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise ValueError(msg)
def components(n,edges):
 parent=list(range(n))
 def find(x):
  while parent[x]!=x:x=parent[x]
  return x
 for a,b in edges:
  a,b=find(a),find(b)
  if a!=b:parent[b]=a
 return tuple(sorted(tuple(i for i in range(n) if find(i)==root) for root in set(find(i) for i in range(n))))
cases=0
for A in range(4):
 for B in range(4):
  for S in range(4):
   for left in itertools.product(range(A),repeat=S):
    for right in itertools.product(range(B),repeat=S):
     edges=list(zip(left,[A+x for x in right]));groups=components(A+B,edges);solutions=[]
     for colors in itertools.product([0,1],repeat=A+B):
      edge_agreement=all(colors[x]==colors[y] for x,y in edges)
      component_constant=all(len({colors[x] for x in g})==1 for g in groups)
      need(edge_agreement==component_constant,'edge equalizer versus pushout coloring')
      if edge_agreement:solutions.append(colors)
     need(len(solutions)==2**len(groups),'independent component bits');cases+=1
head_tail=components(4,[(0,2),(1,3)]);need(len(head_tail)==2,'head tail independent components')
need(len(components(4,[(0,2),(0,3),(1,3)]))==1,'overlapping nonprefix pieces force common shift')
need(len(components(4,[(0,2)]))==3,'unshared pieces isolated')
for mask in range(1<<6):
 edges=[e for j,e in enumerate([(0,1),(1,2),(2,3),(3,4),(0,4),(1,3)]) if mask>>j&1]
 need(components(5,edges)==components(5,list(reversed(edges))),'identification order independent')
 left=components(2,[(0,1)]);right=components(3,[(0,1)])
 combined=components(5,[(0,1),(2,3)])
 need(2**len(combined)==2**len(left)*2**len(right),'disjoint union product of assignments')
# Single-head corruption leaves every other abstract source-piece entry unchanged.
colors=[0,1,0,1];need(all(colors[x]==colors[y] for x,y in [(0,2),(1,3)]),'unmutated split profile')
colors[2]^=1;need(not all(colors[x]==colors[y] for x,y in [(0,2),(1,3)]),'one head color mutation rejected')
print(json.dumps({'status':'passed','checks':checks,'exhaustive_overlap_diagrams':cases,'max_pieces_per_side':3,'max_shared_incidences':3,'head_tail_components':len(head_tail)},indent=2))
