"""Independent directed-arc-set verification of eight syndrome fibres and V4."""
from itertools import combinations,permutations
from collections import Counter
import json

def need(x,msg):
 if not x:raise RuntimeError(msg)
P=(0,4,1,5,2,6,3);Q=(1,6,0,5,3,4,2);L={0,1,2,3}
def edge(u,v):return frozenset((u,v))
def path(p):return frozenset(edge(u,v) for u,v in zip(p,p[1:]))
EP,EQ=path(P),path(Q);family={EP,EQ};free=[(i,j) for i,j in combinations(range(7),2) if (i in L)==(j in L)]
G=[g for g in permutations(range(7)) if {path(tuple(g[i] for i in P)),path(tuple(g[i] for i in Q))}==family]
need(len(G)==4,'V4 size')
I=tuple(range(7));need(all(tuple(g[g[i]] for i in range(7))==I for g in G),'V4 exponent two')
def action(A,g):
 B={(g[u],g[v]) for u,v in A}
 if (P[0],P[1]) not in B:B={(v,u) if (u in L)!=(v in L) else (u,v) for u,v in B}
 return B
def syndrome(A):return tuple(int((u,v) in A)^int((w,z) in A) for u,v,w,z in [(0,1,2,3),(0,2,1,3),(4,5,5,6)])
counts=Counter()
for b in range(512):
 A={(u,v) if b>>k&1 else (v,u) for k,(u,v) in enumerate(free)}|set(zip(P,P[1:]))|{(v,u) for u,v in zip(Q,Q[1:])}
 sig=syndrome(A);counts[str(sig)]+=1
 rev={P[i]:P[-1-i] for i in range(7)}
 anti={(rev[v],rev[u]) for u,v in A}==A
 need(anti==(sig==(0,0,0)),'blue exactly kernel')
 for g in G:
  B=action(A,g);need(syndrome(B)==sig,'trivial syndrome action')
  need(set(zip(P,P[1:]))<=B and {(v,u) for u,v in zip(Q,Q[1:])}<=B,'family preserved')
  for h in G:need(action(B,h)==action(A,tuple(h[g[i]] for i in range(7))),'corrected group law')
need(len(counts)==8 and set(counts.values())=={64},'eight equal fibres')
# Check claimed parity reduction in every n residue and all r,epsilon,s parities.
for n in range(3,24):
 for r in range(1,8):
  for eps in (0,1):
   for s in range(n+1):
    formula=(r*((n-1)*(n-2)//2)-eps*(n*(n-1)//2)-s*(n-s))%2==0
    reduced=[(r-s)%2==0,True,(eps-s)%2==0,(r-eps)%2==0][n%4]
    need(formula==reduced,'mod4 parity table')
print(json.dumps({'status':'FINITE-EXACT','fibres':dict(counts),'group':G,'syndrome_action':'trivial','parity_table_checked_n':'3..23'},indent=2,sort_keys=True))
