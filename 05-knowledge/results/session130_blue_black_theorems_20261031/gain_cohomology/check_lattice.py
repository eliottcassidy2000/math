#!/usr/bin/env python3
"""Exact cycle-lattice certificate for the pinned order-five colored atlas."""
from fractions import Fraction
from collections import deque
import json
# Owner display-node convention of Session129, not least-mask node ranks.
BLUE=[(0,5),(3,8),(4,6),(5,6),(5,8),(6,9),(7,9),(8,9)]
BLACK=[(1,4),(2,4),(2,5),(2,6),(2,7),(2,8),(4,5),(4,8),(5,7),(6,6),(6,8),(8,8)]
def req(c,m):
 if not c:raise RuntimeError(m)
def comps(n,es):
 p=list(range(n))
 def root(a):
  while p[a]!=a:a=p[a]
  return a
 for u,v in es:
  a,b=root(u),root(v)
  if a!=b:p[a]=b
 out={}
 for u in range(n):out.setdefault(root(u),[]).append(u)
 return sorted(out.values())
def cycles(n,es):
 p=list(range(n));adj=[[] for _ in range(n)];tree=[];ch=[]
 def root(a):
  while p[a]!=a:a=p[a]
  return a
 for j,(u,v) in enumerate(es):
  a,b=root(u),root(v)
  if a==b:ch.append(j)
  else:p[a]=b;tree.append(j);adj[u].append((v,j,1));adj[v].append((u,j,-1))
 basis=[]
 for e in ch:
  u,v=es[e];z=[0]*len(es);z[e]=1
  seen={u:None};q=deque([u])
  while v not in seen:
   a=q.popleft()
   for b,j,sign in adj[a]:
    if b not in seen:seen[b]=(a,j,sign);q.append(b)
  b=v
  while b!=u:
   a,j,sign=seen[b];z[j]-=sign;b=a
  req(boundary(n,es,z)==[0]*n,'basis not cycle');basis.append(z)
 return basis,ch,tree
def boundary(n,es,z):
 b=[0]*n
 for (u,v),x in zip(es,z):b[u]-=x;b[v]+=x
 return b
def det(a):
 a=[list(map(int,row)) for row in a];n=len(a);prev=1;sgn=1
 for k in range(n-1):
  if a[k][k]==0:
   j=next((j for j in range(k+1,n) if a[j][k]),None)
   if j is None:return 0
   a[k],a[j]=a[j],a[k];sgn=-sgn
  piv=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    x=a[i][j]*piv-a[i][k]*a[k][j];req(x%prev==0,'Bareiss division');a[i][j]=x//prev
   a[i][k]=0
  prev=piv
 return sgn*a[-1][-1] if n else 1
def rref2(a):
 a=[[x%2 for x in row] for row in a];r=0;ps=[]
 for c in range(len(a[0])):
  j=next((j for j in range(r,len(a)) if a[j][c]),None)
  if j is None:continue
  a[r],a[j]=a[j],a[r]
  for i in range(len(a)):
   if i!=r and a[i][c]:a[i]=[x^y for x,y in zip(a[i],a[r])]
  ps.append(c);r+=1
 return a[:r],ps
Q=BLUE+BLACK;G=BLUE+[e for e in BLACK for _ in (0,1)]
qbas,qch,qt=cycles(10,Q);gbas,gch,gt=cycles(10,G)
req(len(qbas)==11 and len(gbas)==23,'ranks')
projection=[[z[i] for z in qbas] for i in range(8,20)]
a,pivots=rref2(projection);r=len(pivots);req(r==10,'black restriction rank')
K=[]
for p in pivots:
 z=[0]*11;z[p]=2;K.append(z)
for f in range(11):
 if f in pivots:continue
 z=[0]*11;z[f]=1
 for row,p in zip(a,pivots):z[p]=row[f]
 K.append(z)
req(abs(det(list(map(list,zip(*K)))))==2**10,'even quotient lattice index')
plus=[]
for k in K:
 z=[sum(x*b[i] for x,b in zip(k,qbas)) for i in range(20)]
 req(all(x%2==0 for x in z[8:]),'black coefficient parity')
 w=z[:8]+[y for x in z[8:] for y in (x//2,x//2)]
 req(boundary(10,G,w)==[0]*10,'even lift cycle');plus.append(w)
minus=[]
for j in range(12):
 z=[0]*32;z[8+2*j]=1;z[9+2*j]=-1;minus.append(z)
R=list(range(8))+[i^1 for i in range(8,32)]
for z in plus:req([z[i] for i in R]==z,'not even')
for z in minus:req([z[i] for i in R]==[-x for x in z],'not odd')
M=[[z[c] for z in plus+minus] for c in gch]
index=abs(det(M));req(index==1024,'integral eigensplitting index')
# Cycle lift plus twin-difference is an independent integral splitting.
lifts=[z[:8]+[y for x in z[8:] for y in (x,0)] for z in qbas]
req(abs(det([[z[c] for z in lifts+minus] for c in gch]))==1,'split sequence not primitive')
cb=comps(10,BLUE);ck=comps(10,BLACK);cq=comps(10,Q)
bb=len(BLUE)-10+len(cb);bk=len(BLACK)-10+len(ck)
req((bb,bk,len(qbas)-bb-bk)==(1,6,4),'color cycle decomposition')
# Component incidence graph; the five mixed nodes are five parallel edges.
bi={v:i for i,c in enumerate(cb) for v in c};ki={v:len(cb)+i for i,c in enumerate(ck) for v in c}
component_edges=[(bi[v],ki[v]) for v in range(10)]
req(len(component_edges)-(len(cb)+len(ck))+len(comps(len(cb)+len(ck),component_edges))==4,'mixed rank')
P,Qg,Rg=Fraction(55,637),Fraction(539,153),Fraction(585,34)
req(P*Rg/Qg==Fraction(2025,4802),'holonomy')
# Concrete integral splitting hostile: one of the black self-loop twins.
loop=G.index((6,6));e=[0]*32;e[loop]=1
req(boundary(10,G,e)==[0]*10 and e[R[loop]]==0,'loop hostile')
if '--test-failing-gate' in __import__('sys').argv:req(False,'intentional failing gate')
print(json.dumps({'status':'PASS','vertices':10,'blue_edges':BLUE,'black_orbit_edges':BLACK,'literal_edges':32,'quotient_edges':20,'literal_cycle_rank':23,'quotient_cycle_rank':11,'anti_cycle_rank':12,'blue_components':cb,'black_components':ck,'full_components':cq,'blue_cycle_rank':bb,'half_black_cycle_rank':bk,'mixed_cycle_rank':4,'component_incidence_edges':component_edges,'black_mod2_rank':r,'eigensublattice_index':index,'quotient_torsion':'(Z/2)^10','quotient_cycle_basis':qbas,'even_quotient_coordinate_basis':K,'eigenbasis_determinant':det(M),'control_triangle_holonomy':str(P*Rg/Qg)},indent=2,sort_keys=True))
