#!/usr/bin/env python3
import itertools as it,json
from collections import Counter
checks=0
def ck(x):
 global checks;checks+=1
 if not x:raise RuntimeError(checks)
P=(0,4,1,5,2,6,3);Q=(1,6,0,5,3,4,2);L={0,1,2,3}
def edges(p):return {tuple(sorted(e)) for e in zip(p,p[1:])}
EP,EQ=edges(P),edges(Q);free=sorted(set(it.combinations(range(7),2))-EP-EQ)
def syndrome(a):return (a[0][1]^a[2][3],a[0][2]^a[1][3],a[4][5]^a[5][6])
def blue(a,p):return all(a[p[i]][p[j]]==a[p[-1-j]][p[-1-i]] for i in range(7) for j in range(i+1,7))
def relabel(a,g):
 b=[[0]*7 for _ in range(7)]
 for u in range(7):
  for v in range(7):b[g[u]][g[v]]=a[u][v]
 return b
def normalize(b):
 # The skeleton group either preserves or reverses all forced cross directions.
 need=not b[P[0]][P[1]]
 if need:
  for u in L:
   for v in set(range(7))-L:b[u][v]^=1;b[v][u]^=1
 return b,need
G=[]
for g in it.permutations(range(7)):
 es={frozenset(edges(tuple(g[v] for v in p))) for p in (P,Q)}
 if es=={frozenset(EP),frozenset(EQ)}:G.append(g)
ck(len(G)==4)
counts=Counter();action=Counter();fix=Counter()
for bits in range(512):
 a=[[0]*7 for _ in range(7)]
 for k,(u,v) in enumerate(free):a[u][v]=(bits>>k)&1;a[v][u]=1-a[u][v]
 for u,v in zip(P,P[1:]):a[u][v]=1;a[v][u]=0
 for u,v in zip(Q,Q[1:]):a[u][v]=0;a[v][u]=1
 sig=syndrome(a);counts[''.join(map(str,sig))]+=1;ck(blue(a,P)==(sig==(0,0,0)))
 for g in G:
  b,sw=normalize(relabel(a,g));ck(all(b[u][v] for u,v in zip(P,P[1:])));ck(all(not b[u][v] for u,v in zip(Q,Q[1:])));ck(syndrome(b)==sig)
  action[(g,sw)]+=1
  if b==a:fix[str(g)]+=1
  for h in G:
   c,_=normalize(relabel(b,h));gh=tuple(h[g[v]] for v in range(7));d,_=normalize(relabel(a,gh));ck(c==d)
ck(set(counts.values())=={64})
print(json.dumps({'status':'PASS','checks':checks,'free_pairs':free,'syndrome_equations':['x01 XOR x23','x02 XOR x13','x45 XOR x56'],'fibres':dict(sorted(counts.items())),'skeleton_group':G,'syndrome_action':'identity for every element, after explicit cross-cut normalization','normalization_records':[{'permutation':g,'switch_L_needed':sw,'tested_fillings':c} for (g,sw),c in action.items()],'fixed_fillings':dict(fix)},indent=2))
