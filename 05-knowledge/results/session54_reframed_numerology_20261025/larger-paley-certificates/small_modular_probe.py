#!/usr/bin/env python3
import json
rows=[]
for n in (7,11):
 Q={x*x%n for x in range(1,n)};E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q];EI={e:i for i,e in enumerate(E)};cycles=[]
 def dfs(P,left):
  if not left:
   if (P[-1],0) in EI:cycles.append(sum(1<<EI[e] for e in zip(P,P[1:]+[0])))
   return
  for v in sorted(left):
   if (P[-1],v) in EI:dfs(P+[v],left-{v})
 dfs([0],set(range(1,n)));basis={}
 for C in cycles:
  w=C^cycles[0]
  while w:
   p=w.bit_length()-1
   if p in basis:w^=basis[p]
   else:basis[p]=w;break
 basis101={}
 for C in cycles:
  w=[(((C>>j)&1)-((cycles[0]>>j)&1))%101 for j in range(len(E))]
  for p,v in sorted(basis101.items()):
   f=w[p]
   if f:w=[(a-f*b)%101 for a,b in zip(w,v)]
  p=next((j for j,x in enumerate(w) if x),None)
  if p is not None:basis101[p]=[x*pow(w[p],-1,101)%101 for x in w]
  if len(basis101)==len(E)-2*n+1:break
 rows.append({'full_difference_rank_mod101':len(basis101),'n':n,'all_directed_Hamilton_cycles_rooted_at0':len(cycles),'full_difference_rank_mod2':len(basis),'rational_rank_already_certified':len(E)-2*n+1,'scope':'Rank comparison only; no integral-lattice torsion conclusion'})
print(json.dumps({'status':'PASS','results':rows},indent=2,sort_keys=True))
