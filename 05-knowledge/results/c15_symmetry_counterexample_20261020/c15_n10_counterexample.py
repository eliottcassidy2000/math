#!/usr/bin/env python3
"""Exact H=C15 on 3+5+1+1 census, no external packages.
Graphs: all 512 H-fixed graphs; explicit internal odd witnesses exclude448,
then intrinsic twin-class canonicalization and full automorphism sign.
Tournaments: all64 quotient presentations with fixed C3/R5 blocks;
independent exact directed-isomorphism backtracking, no normalizer assumption.
"""
from itertools import combinations,permutations,product
from collections import Counter
import json
checks=0
def ck(p,msg):
 global checks
 checks+=1
 if not p: raise AssertionError(msg)

def edges(n):return list(combinations(range(n),2))
def graph(n,E):
 A=[[0]*n for _ in range(n)]
 for u,v in E:A[u][v]=A[v][u]=1
 return A

def sign(A,p):return -1 if sum(A[u][v] for u,v in edges(len(A)) if p[u]>p[v])%2 else 1

def automorphism(A,p):return all(A[u][v]==A[p[u]][p[v]] for u in range(len(A)) for v in range(len(A)))

def twins(A):
 n=len(A); rem=set(range(n)); out=[]
 while rem:
  u=min(rem); C=tuple(v for v in sorted(rem) if all(A[u][w]==A[v][w] for w in range(n) if w not in (u,v)))
  out.append(C); rem.difference_update(C)
 return out

def full_graph_info(A):
 C=twins(A); m=len(C); n=len(A)
 types=tuple((len(B),A[B[0]][B[1]] if len(B)>1 else 0) for B in C)
 B=[[A[C[i][0]][C[j][0]] if i!=j else 0 for j in range(m)] for i in range(m)]
 # Canonical weighted intrinsic-twin quotient, all class orders.
 can=min((tuple(types[i] for i in p),tuple(B[p[i]][p[j]] for i,j in edges(m))) for p in permutations(range(m)))
 if any(size>1 and internal for size,internal in types):return can,False,{'reason':'true_twins','classes':C}
 signs=[]; odd_witness=None
 for p in permutations(range(m)):
  if any(types[i]!=types[p[i]] for i in range(m)):continue
  if any(B[i][j]!=B[p[i]][p[j]] for i,j in edges(m)):continue
  lift=list(range(n))
  for i in range(m):
   for u,v in zip(C[i],C[p[i]]):lift[u]=v
  ck(automorphism(A,lift),'quotient lift not automorphism')
  s=sign(A,lift)
  weighted=(-1)**sum(len(C[i])*len(C[j])*B[i][j] for i,j in edges(m) if p[i]>p[j])
  ck(s==weighted,'weighted parity mismatch')
  signs.append(s)
  if s<0:odd_witness=lift
 return can,all(s==1 for s in signs),{'classes':C,'quotient_automorphisms':len(signs),'odd_witness':odd_witness}

def tournament(qmask):
 n=10; A=[[0]*n for _ in range(n)];blocks=[list(range(3)),list(range(3,8)),[8],[9]]
 for C,S in ((blocks[0],{1}),(blocks[1],{1,2})):
  for i,u in enumerate(C):
   for j,v in enumerate(C):
    if (j-i)%len(C) in S:A[u][v]=1
 for z,(i,j) in enumerate(edges(4)):
  for u in blocks[i]:
   for v in blocks[j]:A[u][v]=1 if qmask>>z&1 else 0;A[v][u]=1-A[u][v]
 return A

def directed_iso(A,B):
 n=len(A)
 # Joint equitable coloring, respecting incoming and outgoing neighborhoods.
 ca=[sum(r) for r in A];cb=[sum(r) for r in B]
 while True:
  palette=sorted(set(ca+cb)); da=[(ca[i],tuple(sum(A[i][j] for j in range(n) if ca[j]==c) for c in palette)) for i in range(n)]
  db=[(cb[i],tuple(sum(B[i][j] for j in range(n) if cb[j]==c) for c in palette)) for i in range(n)]
  keys={x:k for k,x in enumerate(sorted(set(da+db)))};na=[keys[x] for x in da];nb=[keys[x] for x in db]
  if Counter(na)!=Counter(nb):return False
  if len(set(na))==len(set(ca)) and len(set(nb))==len(set(cb)):ca,cb=na,nb;break
  ca,cb=na,nb
 mapping={};used=set()
 def rec():
  if len(mapping)==n:return True
  best=None;choices=None
  for u in range(n):
   if u in mapping:continue
   opts=[v for v in range(n) if v not in used and ca[u]==cb[v] and all(A[u][x]==B[v][y] and A[x][u]==B[y][v] for x,y in mapping.items())]
   if not opts:return False
   if choices is None or len(opts)<len(choices):best,choices=u,opts
  for v in choices:
   mapping[best]=v;used.add(v)
   if rec():return True
   del mapping[best];used.remove(v)
  return False
 return rec()

def main():
 n=10;g=list(range(n))
 for C in ([0,1,2],[3,4,5,6,7]):
  for i,u in enumerate(C):g[u]=C[(i+1)%len(C)]
 H=[];p=list(range(n))
 for _ in range(15):H.append(p);p=[g[p[i]] for i in range(n)]
 ck(p==list(range(n)),'order');pairs=edges(n);remaining=set(pairs);orbits=[]
 while remaining:
  e=min(remaining);O={tuple(sorted((h[e[0]],h[e[1]]))) for h in H};orbits.append(O);remaining-=O
 ck(len(orbits)==9,'pair orbit count')
 excluded=0;relative_even_keys=set();even_keys=set();graph_keys=set();bad_relative=[];counter=Counter();records=[]
 root_swap=list(range(n));root_swap[8],root_swap[9]=9,8
 # A product group preserving the distinguished3/5 blocks, with root swap only.
 for bits in product((0,1),repeat=len(orbits)):
  E=set().union(*(O for b,O in zip(bits,orbits) if b))
  A=graph(n,E);ck(automorphism(A,g),'not H-fixed')
  internal3=sum(A[i][j] for i,j in combinations(range(3),2));internal5=sum(A[i][j] for i,j in combinations(range(3,8),2))
  if internal3 or internal5:
   witness=list(range(n))
   if internal3:witness[0],witness[1]=1,0
   elif internal5==10:witness[3],witness[4]=4,3
   else:
    for i in range(5):witness[3+i]=3+(-i)%5
   ck(automorphism(A,witness) and sign(A,witness)==-1,'internal odd witness')
   excluded+=1;continue
  # Internal empty blocks are fixed; relative quotient only exchanges two singleton roots.
  qedges=tuple(A[u][v] for u,v in ((0,3),(0,8),(0,9),(3,8),(3,9),(8,9)))
  swapq=tuple(A[root_swap[u]][root_swap[v]] for u,v in ((0,3),(0,8),(0,9),(3,8),(3,9),(8,9)))
  relative_key=min(qedges,swapq)
  relative_even=not automorphism(A,root_swap) or sign(A,root_swap)==1
  key,even,info=full_graph_info(A);graph_keys.add(key)
  if relative_even:relative_even_keys.add(relative_key)
  if even:even_keys.add(key);ck(relative_even,'full-even must be relative-even')
  if relative_even and not even:bad_relative.append({'qedges':qedges,'degrees':[sum(r) for r in A],'info':info})
  records.append((relative_key,key,relative_even,even))
 ck(excluded==448,'excluded universe')
 # Independently classify all 64 expanded tournaments by full directed isomorphism.
 treps=[];tclasses=[]
 for mask in range(64):
  A=tournament(mask)
  ck(automorphism(A,g),'tournament H invariance')
  i=next((i for i,B in enumerate(treps) if directed_iso(A,B)),None)
  if i is None:i=len(treps);treps.append(A)
  tclasses.append(i)
 # Directly enumerate all512 H-fixed tournament orientations independently.
 arc_orbits=[]
 for O in orbits:
  u,v=min(O);arcs={(h[u],h[v]) for h in H}
  ck(len(arcs)==len(O),'oriented orbit reverses pair')
  arc_orbits.append(arcs)
 full_tclasses=set()
 for choice in product((0,1),repeat=len(orbits)):
  A=[[0]*n for _ in range(n)]
  for bit,AO in zip(choice,arc_orbits):
   for u,v in AO:
    if bit:A[u][v]=1
    else:A[v][u]=1
  ck(automorphism(A,g),'direct tournament not H-fixed')
  idx=next((i for i,B in enumerate(treps) if directed_iso(A,B)),None)
  ck(idx is not None,'internal-orientation reduction failed')
  full_tclasses.add(idx)
 ck(len(full_tclasses)==len(treps),'direct universe missing type')
 # Quotient root swap must give exactly the same equivalence classes.
 def swap_mask(mask):
  ep=edges(4);ix={e:i for i,e in enumerate(ep)};p=[0,1,3,2];out=0
  for z,(u,v) in enumerate(ep):
   bit=(mask>>z)&1;a,b=p[u],p[v]
   if a>b:a,b=b,a;bit^=1
   out|=bit<<ix[a,b]
  return out
 for a in range(64):
  for b in range(64):ck((tclasses[a]==tclasses[b])==(b==a or b==swap_mask(a)),'unexpected tournament identification')
 result={'status':'PASS','checks':checks,'H':'C15, cycle type (5,3,1,1)','n':10,'pair_orbits':len(orbits),'H_fixed_graphs':512,'H_fixed_tournaments':512,'direct_full_tournament_types':len(full_tclasses),'explicit_internal_odd_exclusions':excluded,'remaining_labelled_graphs':64,'remaining_graph_types':len(graph_keys),'relative_even_quotient_types':len(relative_even_keys),'actual_sign_even_graph_types':len(even_keys),'actual_tournament_types':len(treps),'counterexample_difference':len(treps)-len(even_keys),'relative_even_but_globally_odd':bad_relative,'tournament_source_masks':64}
 print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
