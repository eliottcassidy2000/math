import json, math
from collections import deque
checks=0
for p in range(3,30,2):
 states=[(a,q) for a in range(2*p) for q in range(p)]
 edges={s:[((2*s[0]+b)%(2*p),(s[1]+b*(s[0]+s[0]%2))%p) for b in (0,1)] for s in states}
 incoming={s:0 for s in states}
 for s in states:
  for t in edges[s]: incoming[t]+=1
 if set(incoming.values())!={2}:raise RuntimeError('stationarity')
 # A path from zero to every state and from every state to zero proves SCC.
 rev={s:[] for s in states}
 for s in states:
  for t in edges[s]:rev[t].append(s)
 for graph in (edges,rev):
  seen={(0,0)};todo=deque(seen)
  while todo:
   for t in graph[todo.popleft()]:
    if t not in seen:seen.add(t);todo.append(t)
  if len(seen)!=len(states):raise RuntimeError('irreducibility')
 if (0,0) not in edges[(0,0)]:raise RuntimeError('aperiodicity')
 checks+=len(states)
N=1000000;q=[0]*(N+1);counts={p:[0]*p for p in (3,5,7)}
for m in range(N+1):
 t=m//2;q[m]=q[t]+(t+t%2 if m%2 else 0)
 h=1+2*m-q[m]
 if h<=N:
  for p in counts:counts[p][h%p]+=1
print(json.dumps({'status':'PASS','odd_moduli_checked':14,'states_checked':checks,'finite_parameter_cutoff':N,'parameter_residue_counts':counts,'predicted_coefficient':{p:math.log(2)/p for p in counts}},indent=2))
