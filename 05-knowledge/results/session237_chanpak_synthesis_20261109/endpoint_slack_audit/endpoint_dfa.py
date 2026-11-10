"""Independent reachable automaton partition refinement for three observers."""
import json,sys
from collections import deque

def require(x,msg):
 if not x:raise RuntimeError(msg)
require('--fail' not in sys.argv,'deliberate checker failure')
out=[]
for M in range(1,41):
 states=[(0,0)];idx={(0,0):0};edges=[]
 for r,q in states:
  targets=[((2*r)%(2*M),q),((2*r+1)%(2*M),(q+(r+1)//2)%M)]
  row=[]
  for t in targets:
   if t not in idx:idx[t]=len(states);states.append(t)
   row.append(idx[t])
  edges.append(row)
 require(len(states)==2*M*M,'joint reachable')
 def minimize(outputs):
  part=outputs
  while True:
   ids={};new=[]
   for i,y in enumerate(outputs):
    key=(y,part[edges[i][0]],part[edges[i][1]])
    if key not in ids:ids[key]=len(ids)
    new.append(ids[key])
   if len(set(new))==len(set(part)):return len(ids)
   part=new
 E=minimize([(r-2*q)%M for r,q in states])
 F=minimize([(r-q)%M for r,q in states])
 H=minimize([(1+2*r-2*q)%M for r,q in states])
 expectedE=1 if M==1 else 2*M*M if M%2 else M*M//2
 expectedF=1 if M==1 else 2*M*M
 K=M if M%2 else M//2
 expectedH=1 if K==1 else 2*K*K
 require((E,F,H)==(expectedE,expectedF,expectedH),f'counts at M={M}')
 out.append({'M':M,'endpoint_E':E,'normalized_total':F,'actual_total':H})
print(json.dumps({'status':'PASS','method':'BFS and Moore partition refinement, M=1..40','counts':out},indent=2))
