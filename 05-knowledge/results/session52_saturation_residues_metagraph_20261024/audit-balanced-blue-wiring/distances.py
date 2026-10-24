#!/usr/bin/env python3
"""Independent all-pairs unweighted blue-graph distance receipt."""
from pathlib import Path
from collections import Counter,deque
import json
P=Path(__file__).resolve().parent;d=json.loads((P/'producer/normal.json').read_text());nei={r['atlas_index']:set() for r in d['nodes']}
for a,b,m in d['adjacency_multiplicities']:
 if a!=b and m>0:nei[a].add(b);nei[b].add(a)
hist=Counter();exceptions=[]
for start in sorted(nei):
 distances={start:0};q=deque([start])
 while q:
  v=q.popleft()
  for w in nei[v]:
   if w not in distances:distances[w]=distances[v]+1;q.append(w)
 if len(distances)!=229:raise RuntimeError('disconnected')
 for end in sorted(nei):
  if start<end:
   value=distances[end];hist[value]+=1
   if value>2:exceptions.append([start,end,value])
if hist!={1:9219,2:16886,3:1} or exceptions!=[[63,86,3]]:raise RuntimeError('distance profile')
if len(nei[63])!=17 or len(nei[86])!=17 or nei[63]&nei[86]:raise RuntimeError('exception neighborhoods')
print(json.dumps({'status':'PASS','scope':'unweighted simple support of blue-only induced multigraph','diameter':3,'unordered_pair_distance_histogram':sorted(hist.items()),'unique_distance_three_pair':[63,86],'exception_node_records':[r for r in d['nodes'] if r['atlas_index'] in [63,86]],'exception_neighbor_lists':[[v,sorted(nei[v])] for v in [63,86]],'exception_simple_degrees':[17,17],'exception_common_neighbors':[],'pair_total':sum(hist.values())},indent=2))
