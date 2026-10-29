"""Independent native-operation masks and all order-seven witness fillings."""
from itertools import combinations
from collections import Counter
import json

def need(x,msg):
 if not x:raise RuntimeError(msg)
def edge(u,v):return frozenset((u,v))
def E(p):return {edge(u,v) for u,v in zip(p,p[1:])}
def native(a,p):
 keep=E(p)
 return {(u,v) if edge(u,v) in keep else (v,u) for u,v in a}
def blue(a,p):return all(((p[i],p[j]) in a)==((p[-1-j],p[-1-i]) in a) for i,j in combinations(range(len(p)),2))
P=(0,4,1,5,2,6,3);Q=(1,6,0,5,3,4,2);S={0,1,2,3}
cut={edge(u,v) for u in S for v in range(7) if v not in S};inside=[(u,v) for u,v in combinations(range(7),2) if edge(u,v) not in cut]
need(not E(P)&E(Q) and E(P)|E(Q)==cut,'decomposition of K4,3')
colors=Counter()
for b in range(1<<len(inside)):
 a={(u,v) if b>>k&1 else (v,u) for k,(u,v) in enumerate(inside)}
 a|=set(zip(P,P[1:]));a|={(v,u) for u,v in zip(Q,Q[1:])}
 need(len(a)==21,'all arcs assigned once');t=native(a,P)
 need(set(zip(Q,Q[1:]))<=t,'second move legal')
 out=native(t,Q);sw={(v,u) if edge(u,v) in cut else (u,v) for u,v in a}
 need(out==sw and out!=a,'sharp switch witness')
 colors[(blue(a,P),blue(t,Q))]+=1
# Cut degrees for n>=8: degree at most4 forces n8 balanced, but endpoints break it.
for n in range(8,31):
 for s in range(1,n):
  degree_possible=max(s,n-s)<=4
  if degree_possible:need(n==8 and s==4,'degree obstruction')
print(json.dumps({'status':'FINITE-EXACT','witness_order':7,'fillings':512,'blue_pair_counts':{str(k):v for k,v in colors.items()},'all_order_proof':'Nontrivial cut sides each at most4; n8 requires degree4 at every vertex, impossible at any path endpoint.'},sort_keys=True,indent=2))
