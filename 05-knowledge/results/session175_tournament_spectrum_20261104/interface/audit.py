#!/usr/bin/env python3
"""Independent literal-path audit; imports no author counting code."""
import itertools,json
from pathlib import Path
p=Path(__file__).parent
j=json.loads((p/'certificate.json').read_text());perms=list(itertools.permutations(range(7)));total=0
for f in j['families']:
 A=(6,4,0) if f['name'][0]=='T' else (2,4,1)
 B=(6,4,0) if f['name'][1]=='T' else (2,4,1)
 actual=[]
 for k in range(512):
  arcs={(i,t) for i in range(3) for t in range(3) if A[i]>>t&1}
  arcs|={(i+3,t+3) for i in range(3) for t in range(3) if B[i]>>t&1}
  arcs|={(6,i) for i in range(3)}|{(i,6) for i in range(3,6)}
  arcs|={(i,t+3) if k>>(3*i+t)&1 else (t+3,i) for i in range(3) for t in range(3)}
  h=sum(all((q[i],q[i+1]) in arcs for i in range(6)) for q in perms)
  if h!=f['H'][k]:raise RuntimeError(('literal count',f['name'],k,h))
  actual.append(h);total+=1
 for label in ['components','strong_components']:
  for c in f[label]:
   states=set(c['states']);start=next(iter(states));seen={start};todo=[start]
   while todo:
    k=todo.pop()
    for bit in range(9):
     t=k^(1<<bit)
     if t in states and t not in seen and abs(actual[t]-actual[k])<=2:seen.add(t);todo.append(t)
   if seen!=states:raise RuntimeError('disconnected certificate')
   if sorted({actual[k] for k in states})!=list(range(c['min'],c['max']+1,2)):raise RuntimeError('interval hole')
print(json.dumps({'independent_literal_order7_states':total,'permutations_per_state':5040,'all_component_connectivity_and_interval_checks':'PASS'},sort_keys=True))
