#!/usr/bin/env python3
"""Independent exact family control for a=7,b=11,k=4 (n=22).
Enumerates all 32768 independent-large-block graphs on the six-node weighted
quotient, computes N-orbits and full weighted false-twin isomorphism/sign.
No producer imports; full graph automorphisms are handled by twin reduction.
"""
from itertools import combinations,permutations,product
from math import factorial
from pathlib import Path
import json

def need(x,msg):
 if not x:raise AssertionError(msg)

def main():
 a,b=7,11;k=b-a;s=k+2;weights=[a,b]+[1]*k;edges=list(combinations(range(s),2));idx={e:i for i,e in enumerate(edges)}
 roots=[(0,1)+p for p in permutations(range(2,s))]
 actions=[]
 for p in roots:
  actions.append([1<<idx[tuple(sorted((p[u],p[v]))) ] for u,v in edges])
 def apply(q,act):
  z=0
  while q:
   bit=q&-q;z|=act[bit.bit_length()-1];q-=bit
  return z
 seen=set();normalizer=[];globaltypes={};extraodd=[];counts={'presentations':1<<len(edges),'N_orbits':0,'N_even':0,'global_even':0};familyrecords=[]
 for q in range(1<<len(edges)):
  if q in seen:continue
  orb={apply(q,act) for act in actions};seen|=orb;counts['N_orbits']+=1
  Q=[[0]*s for _ in range(s)]
  for bit,(u,v) in enumerate(edges):Q[u][v]=Q[v][u]=(q>>bit)&1
  ne=True
  for p,act in zip(roots,actions):
   if apply(q,act)==q:
    parity=sum(weights[u]*weights[v] for u,v in edges if Q[u][v] and p[u]>p[v])%2
    if parity:ne=False
  if ne:counts['N_even']+=1
  # Intrinsic false-twin partition of the full graph, represented on quotient nodes.
  cells={}
  for i,row in enumerate(Q):cells.setdefault(tuple(row),[]).append(i)
  cells=list(cells.values());r=len(cells);W=[sum(weights[i] for i in c) for c in cells]
  R=[[Q[cells[i][0]][cells[j][0]] for j in range(r)] for i in range(r)]
  groups={}
  for i,w in enumerate(W):groups.setdefault(w,[]).append(i)
  orders=[]
  for parts in product(*(list(permutations(groups[w])) for w in sorted(groups))):orders.append(tuple(x for part in parts for x in part))
  canon=min(tuple(W[p[i]] for i in range(r))+tuple(R[p[i]][p[j]] for i,j in combinations(range(r),2)) for p in orders)
  # All equal-weight cell automorphisms; sorted equal-weight indexing restored.
  base=orders[0];even=True;oddcell=None
  for order in orders:
   p=list(range(r))
   for x,y in zip(base,order):p[x]=y
   if any(R[i][j]!=R[p[i]][p[j]] for i in range(r) for j in range(r)):continue
   parity=sum(W[i]*W[j] for i,j in combinations(range(r),2) if R[i][j] and p[i]>p[j])%2
   if parity:even=False;oddcell=p
  if canon in globaltypes:need(globaltypes[canon]==q,'no full graph fusion beyond N')
  globaltypes[canon]=q
  if even:counts['global_even']+=1
  if ne and not even:extraodd.append(q);familyrecords.append({'q':q,'false_twin_weights':W,'false_twin_adjacency':R,'odd_cell_permutation':oddcell})
 # Independent tournament quotient orbit enumeration, affine under root permutation.
 tacts=[]
 for p,act in zip(roots,actions):
  inv=sum(1<<idx[tuple(sorted((p[u],p[v])))] for u,v in edges if p[u]>p[v]);tacts.append((act,inv))
 tseen=set();tc=0
 for q in range(1<<len(edges)):
  if q in tseen:continue
  tseen|={apply(q,act)^inv for act,inv in tacts};tc+=1
 need(len(seen)==1<<len(edges) and len(tseen)==1<<len(edges),'complete universes')
 need(counts['N_even']==tc==1408,'N even equals tournament types')
 need(counts['global_even']==1407 and len(extraodd)==1,'exact gap one')
 need(familyrecords[0]['false_twin_weights']==[b,b],'exception Kbb')
 out={'a':a,'b':b,'k':k,'n':2*b,**counts,'tournament_types':tc,'full_graph_types':len(globaltypes),'extra_odd':familyrecords}
 Path('exact_gap_family_check.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2));print('ALL INDEPENDENT ORDER-22 FAMILY CHECKS PASSED')
if __name__=='__main__':main()
