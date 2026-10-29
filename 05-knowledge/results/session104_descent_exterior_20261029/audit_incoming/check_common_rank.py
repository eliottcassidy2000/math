#!/usr/bin/env python3
"""Literal native K4 path covers; independent W1/W2 audit, stdlib only."""
from itertools import permutations
from collections import Counter
from math import comb
import json
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def mul(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def add(A,B):return [[x+y for x,y in zip(a,b)] for a,b in zip(A,B)]
def scale(c,A):return [[c*x for x in r] for r in A]
def trace(A):return sum(A[i][i] for i in range(len(A)))
def topo(arcs,n):
 order=[];remaining=set(range(n))
 while remaining:
  ready=sorted(v for v in remaining if not any(u in remaining and w==v for u,w in arcs))
  if not ready:return None
  order+=ready;remaining.difference_update(ready)
 return order
def cycle(arcs,n):
 def walk(v,path):
  for u,w in arcs:
   if u!=v:continue
   if w in path:return path[path.index(w):]+[w]
   found=walk(w,path+[w])
   if found:return found
  return None
 for v in range(n):
  found=walk(v,[v])
  if found:return found
 return None
arcs={(0,1),(1,2),(2,3),(3,0),(0,2),(1,3)}
ports=list(permutations(range(4),2));index={p:i for i,p in enumerate(ports)};m=len(ports)
I=[[int(i==j) for j in range(m)] for i in range(m)]
M=[[0]*m for _ in range(m)];cover_count=0
# A permutation and one cut encode exactly one ordered pair of disjoint paths.
for word in permutations(range(4)):
 for cut in range(1,4):
  paths=(word[:cut],word[cut:])
  if all((u,v) in arcs for path in paths for u,v in zip(path,path[1:])):
   a=tuple(path[0] for path in paths);t=tuple(path[-1] for path in paths)
   M[index[a]][index[t]]+=1;cover_count+=1
phi={a:((a[1]-1)%4,(a[0]-1)%4) for a in ports}
need(len(set(phi.values()))==m,'canonical cut matching bijective')
B=[[M[i][index[phi[b]]] for b in ports] for i in range(m)]
N=add(B,scale(-1,I));need(all(x>=0 for row in N for x in row),'nonnegative N')
D={( (u+1)%4,v) for u,v in arcs if v!=(u+1)%4}
order=topo(D,4);need(order is not None,'seed auxiliary DAG')
f={v:i for i,v in enumerate(order)};w=[sum(f[v] for v in a) for a in ports]
for i in range(m):
 for j in range(m):need(not N[i][j] or w[j]>w[i],'weight ascent')
need(max(w)-min(w)==4,'k(n-k) range')
powers=[I]
for _ in range(5):powers.append(mul(powers[-1],N))
need(not any(x for row in powers[5] for x in row),'W1 N^5 zero')
actual_index=next(j for j in range(1,6) if not any(x for row in powers[j] for x in row))
inv=I
for j in range(1,5):inv=add(inv,scale((-1)**j,powers[j]))
need(mul(B,inv)==I and mul(inv,B)==I,'finite inverse')
Bd=I
for d in range(11):
 formula=I
 for j in range(1,5):formula=add(formula,scale(comb(d,j) if j<=d else 0,powers[j]))
 need(Bd==formula,'binomial power');Bd=mul(Bd,B)
hist=Counter();swap_trace=None;cyclic_trace12=0
support={(i,j) for i in range(m) for j in range(m) if N[i][j]}
for p in permutations(range(4)):
 perm=[index[tuple(p[v] for v in a)] for a in ports]
 C=[[B[perm[i]][perm[j]] for j in range(m)] for i in range(m)]
 NC=add(C,scale(-1,I));powC=I
 for _ in range(5):powC=mul(powC,NC)
 need(not any(x for row in powC for x in row),'remarked nilpotence')
 product=mul(B,C);t=trace(product);hist[t]+=1
 need(all(product[i][i]>=1 for i in range(m)),'product positive diagonal')
 union=support|{(i,j) for i in range(m) for j in range(m) if NC[i][j]}
 if topo(union,m) is not None:need(t==12,'DAG union trace')
 else:
  if t==12:cyclic_trace12+=1
  cyc=cycle(union,m);need(cyc is not None and len(cyc)-1<=m,'simple cycle witness')
  witness=I
  for u,v in zip(cyc,cyc[1:]):witness=mul(witness,B if N[u][v] else C)
  need(trace(witness)>m,'constructive nonunipotent word')
 if p==(0,3,2,1):swap_trace=t;need(t==22,'swap13 trace22');need(max(product[i][i] for i in range(m))>1,'exponential diagonal')
need(hist==Counter({12:17,16:6,22:1}),'all24 trace histogram')
# Opposite shears independently expose failure of individual unipotence closure.
need(mul([[1,1],[0,1]],[[1,0],[1,1]])==[[2,1],[1,1]],'opposite shears')
print(json.dumps({'status':'PASS','checks':checks,'native_ordered_covers':cover_count,'port_dimension':m,'N_nilpotence_index':actual_index,'all_vertex_remarkings':24,'trace_histogram':dict(hist),'swap_1_3_trace':swap_trace,'cyclic_union_with_trace12':cyclic_trace12,'producer_imports':False,'uses_assert':False,'scope':'K4 W1/W2 only; no n6 k3 sharpness claim'},sort_keys=True,indent=2))
