#!/usr/bin/env python3
"""Independent literal replay and prepend-style FAS DP for saved near miss."""
import json,math,hashlib
from pathlib import Path
from functools import lru_cache
root=Path(__file__).resolve().parent
r=json.loads((root/'receipt.json').read_text());b=r['best_literal_candidate'];n=11

def check(t,label):
 if not t:raise RuntimeError(label)
def fas(T):
 @lru_cache(None)
 def solve(S):
  if not S:return 0
  return min(solve(S-{v})+sum(T[u][v] for u in S if u!=v) for v in S)
 return solve(frozenset(range(n)))
def flip(T,P):
 check(set(P)==set(range(n)) and len(P)==n,'path not permutation')
 check(all(T[x][y] for x,y in zip(P,P[1:])),'path not legal')
 keep={frozenset(x) for x in zip(P,P[1:])}
 return [[False if i==j else T[i][j] if frozenset((i,j)) in keep else not T[i][j] for j in range(n)] for i in range(n)]
def rows(T):return [sum(1<<j for j in range(n) if T[i][j]) for i in range(n)]
squares={i*i%n for i in range(1,n)};T=[[i!=j and (j-i)%n in squares for j in range(n)] for i in range(n)]
A=flip(T,b['P']);F=flip(A,b['Q']);check(rows(A)==b['intermediate_rows'],'intermediate rows');check(rows(F)==b['final_rows'],'final rows')
check(flip(F,b['Q'])==A and flip(A,b['P'])==T,'reverse word')
values=[fas(X) for X in (T,A,F)];check(values==[20,b['intermediate_FAS'],13],'independent prepend FAS')
pos={v:i for i,v in enumerate(b['final_rank'])};E=[(i,j) for i in range(n) for j in range(i+1,n) if F[i][j]!=(pos[i]<pos[j])]
check(E==[tuple(e) for e in b['final_edges']] and len(E)==13,'literal defects')
adj=[set() for _ in range(n)]
for x,y in E:adj[x].add(y);adj[y].add(x)
unseen=set(range(n));components=0
while unseen:
 components+=1;Q=[unseen.pop()]
 for v in Q:
  for u in adj[v]:
   if u in unseen:unseen.remove(u);Q.append(u)
beta=len(E)-n+components;check(beta==b['final_rank_beta'],'literal beta')
check(r['final_rank_draws']==0 and r['minimum_observed_final_beta'] is None,'no sampled forest ranks')
check(r['status'].startswith('NO_WITNESS'),'scope')
check(r['first_rank_draws']==64 and r['first_path_draws']==512 and r['distinct_first_states']==337,'first search domain')
check(r['second_rank_draws']==2696 and r['second_path_draws']==10784 and r['distinct_second_states']==8517,'second search domain')
clocks=[]
for N in (7,21,63):
 f=N*(N-1)//6;lower=math.ceil((f-(N-1))/(N-3));clocks.append({'n':N,'FAS':f,'complement_lower_bound':lower})
check([x['complement_lower_bound'] for x in clocks]==[1,3,10],'family clocks')
check((root/'receipt.json').read_bytes()==(root/'receipt.optimized.json').read_bytes(),'normal optimized mismatch')
print(json.dumps({'status':'VERIFIED finite receipt and literal near miss; no two-complement impossibility claim','independent_FAS':values,'final_rank_cycle_debt':beta,'final_components':components,'P':b['P'],'Q':b['Q'],'whole_complement_clock':2,'family_clock_lower_bounds':clocks},indent=2,sort_keys=True))
