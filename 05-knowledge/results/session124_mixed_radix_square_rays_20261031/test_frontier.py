"""Brute-force maximum antichain controls for the schedule-frontier counter."""
import json,itertools
from schedule_frontier import count_frontier,frontier
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def literal_prefix(word):
 p=[1,3]
 for a in reversed(word):
  q,c,h={'A':(25,12,[1,3]),'B':(25,8,[1,3,6,19]),'C':(49,4,[1,3,6,19])}[a]
  p=h+[q*x+(-1)**(i+1)*c for i,x in enumerate(p)]
 return p
small=[]
for N in [1062,2082,26562,52062,102042]:
 nodes=[('',1)];i=0
 while i<len(nodes):
  w,r=nodes[i];i+=1
  nodes.extend((w+a,r*q) for a,q in [('A',25),('B',25),('C',49)] if 85*r*q<=2*N+1)
 paths=[literal_prefix(w) for w,r in nodes];conflicts=[]
 for i,x in enumerate(paths):
  conflicts.append(sum(1<<j for j,y in enumerate(paths[:i]) if x[:min(len(x),len(y))]==y[:min(len(x),len(y))]))
 best=0
 for mask in range(1<<len(nodes)):
  ok=True
  for i in range(len(nodes)):
   if(mask>>i)&1 and mask&conflicts[i]:ok=False;break
  if ok:best=max(best,mask.bit_count())
  checks+=1
 result=count_frontier(N);need(best==result['distinct_cycles_lower_bound'],'literal-prefix maximum')
 small.append({'N':N,'tree_nodes':len(nodes),'maximum':best})
large=[]
for N in [16601562,10**9]:
 leaves=list(frontier(N));need(len(leaves)==count_frontier(N)['distinct_cycles_lower_bound'],'memoized versus listed')
 paths=[literal_prefix(w) for w,r in leaves]
 for w,r in leaves:need(bool(w) and 85*r<=2*N+1,'certified product budget')
 for i,x in enumerate(paths):
  for y in paths[:i]:need(x[:min(len(x),len(y))]!=y[:min(len(x),len(y))],'literal-prefix antichain')
 large.append({'N':N,'leaves':len(leaves)})
try:count_frontier(1061)
except ValueError:checks+=1
else:raise RuntimeError('invalid cycle-count domain accepted')
print(json.dumps({'status':'passed','checks':checks,'brute_controls':small,'large_controls':large},indent=2))
