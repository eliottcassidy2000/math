#!/usr/bin/env python3
"""Exact golden-cell quotient and flat-torus Paley action controls. Stdlib only."""
from fractions import Fraction as F
from itertools import product
import json
CHECKS=0
def ck(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise AssertionError(msg)
def mul(A,B):return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def add(A,B):return tuple(tuple(A[i][j]+B[i][j] for j in range(2)) for i in range(2))
def scale(t,A):return tuple(tuple(t*x for x in r) for r in A)
def trans(A):return tuple(zip(*A))
def det(A):return A[0][0]*A[1][1]-A[0][1]*A[1][0]
def power(A,k):
 B=I
 for _ in range(k):B=mul(B,A)
 return B
def mv(A,v):return tuple(sum(A[i][j]*v[j] for j in range(2)) for i in range(2))
I=((1,0),(0,1));Q=((0,1),(1,1));C=((4,-1),(-1,3))
def main():
 ck(mul(Q,Q)==add(Q,I),'golden relation');ck(det(C)==11,'eleven cells')
 ck(add(power(Q,5),scale(-1,I))==mul(power(Q,3),C),'five-step cell kernel')
 ck(power(Q,5)!=I,'lift not periodic')
 for a,b in product(range(-12,13),repeat=2):
  label=lambda x:(x[0]+4*x[1])%11
  ck(label(mv(Q,(a,b)))==4*label((a,b))%11,'clock label')
  for z in ((1,0),(0,1)):
   c=mv(C,z);ck(label((a+c[0],b+c[1]))==label((a,b)),'quotient relation')
 ck([pow(4,k,11) for k in range(5)]==[1,4,5,9,3],'order five')
 for den in range(1,16):
  for num in range(-20,21):
   w=F(num,den);M=add(I,scale(w,Q));B=mul(M,C)
   ck(B==add(scale(4-w,I),scale(3*w-1,Q)),'lattice polynomial')
   ck(det(B)==11*(1+w-w*w),'covolume polynomial')
   ck(mul(B,Q)==mul(Q,B),'clock descends')
   if det(B):
    G=mul(trans(B),B)
    ck(det(G)>0 and G[0][0]>0,'positive Gram')
    ck(mul(mul(trans(Q),G),Q)!=G,'not an isometry')
 ck(mul(add(I,scale(F(1,3),Q)),C)==scale(F(11,3),I),'square lattice special value')
 # The linear invariance equations Q^t [[a,b],[b,c]] Q=G force zero:
 for a,b,c in product(range(-5,6),repeat=3):
  G=((a,b),(b,c));ck((mul(mul(trans(Q),G),Q)==G)==(a==b==c==0),'no nonzero invariant bilinear form')
 # Explicit hexagonal torus, translations represented modulo 7.
 R=((0,-1),(1,-1));G=((2,-1),(-1,2));v=(3,1)
 ck(power(R,3)==I,'third turn');ck(mul(mul(trans(R),G),R)==G,'hexagonal metric')
 ck(tuple(x%7 for x in mv(R,v))==tuple(2*x%7 for x in v),'multiplier two')
 elems=[(j,b) for j in range(3) for b in range(7)]
 def comp(x,y):return ((x[0]+y[0])%3,(x[1]+pow(2,x[0],7)*y[1])%7)
 def act(g,x):return (pow(2,g[0],7)*x+g[1])%7
 perms={tuple(act(g,x) for x in range(7)) for g in elems};ck(len(perms)==21,'faithful action')
 for g,h in product(elems,repeat=2):
  for x in range(7):ck(act(comp(g,h),x)==act(g,act(h,x)),'affine law')
 for g in elems:
  for x,y in product(range(7),repeat=2):
   if x!=y:ck(((y-x)%7 in {1,2,4})==((act(g,y)-act(g,x))%7 in {1,2,4}),'Paley arcs')
 # Short displacement check on the seven torus points: all six differences
 # have equal minimal squared distance 2/7, unlike the planar seven-point bound.
 distances=[]
 for d in range(1,7):
  q=tuple(F(d*x,7) for x in v)
  vals=[]
  for z in product(range(-6,7),repeat=2):
   u=tuple(q[i]-z[i] for i in range(2));vals.append(sum(u[i]*G[i][j]*u[j] for i in range(2) for j in range(2)))
  distances.append(min(vals))
 ck(distances==[F(2,7)]*6,'seven torus points equidistant')
 # The torus quotient distance is not the same as distance between fixed planar lifts.
 lifts=[tuple(F(x*a,7) for a in v) for x in range(7)]
 direct=sum(lifts[6][i]*G[i][j]*lifts[6][j] for i in range(2) for j in range(2))
 ck(direct>F(2,7),'planar lifts lose wrapped distance')
 # Signed norm, triangular torus and two Fano systems.
 N=((2,1),(1,-2))
 ck(mul(mul(trans(Q),N),Q)==scale(-1,N),'golden norm anti-invariance')
 faces=[(r,(r+1)%7,(r+3)%7) for r in range(7)]+[(r,(r+3)%7,(r+2)%7) for r in range(7)]
 from collections import Counter
 edgecount=Counter(tuple(sorted((f[i],f[(i+1)%3]))) for f in faces for i in range(3))
 directed=Counter((f[i],f[(i+1)%3]) for f in faces for i in range(3))
 ck(len({frozenset(f) for f in faces})==14,'distinct faces')
 ck(len(edgecount)==21 and set(edgecount.values())=={2},'two faces per edge')
 ck(all(directed[(x,y)]==directed[(y,x)]==1 for x,y in edgecount),'oriented gluing')
 ck(7-len(edgecount)+len(faces)==0,'Euler characteristic')
 for vertex in range(7):
  link=[tuple(x for x in f if x!=vertex) for f in faces if vertex in f]
  adj={x:set() for x in range(7) if x!=vertex}
  for x,y in link:adj[x].add(y);adj[y].add(x)
  seen=set();todo=[next(iter(adj))]
  while todo:
   x=todo.pop()
   if x not in seen:seen.add(x);todo.extend(adj[x]-seen)
  ck(len(seen)==6 and all(len(a)==2 for a in adj.values()),'vertex link is six-cycle')
 def span(rows):
  result={0}
  for row in rows:result|={x^row for x in result}
  return result
 codes=[]
 for family in (faces[:7],faces[7:]):
  pairs=Counter(tuple(sorted((f[i],f[(i+1)%3]))) for f in family for i in range(3))
  ck(len(pairs)==21 and set(pairs.values())=={1},'Fano pair incidence')
  code=span([sum(1<<x for x in f) for f in family]);codes.append(code)
  ck(len(code)==16 and min(x.bit_count() for x in code if x)==3,'Hamming parameters')
  ck(Counter(x.bit_count() for x in code)=={0:1,3:7,4:7,7:1},'Hamming weight enumerator')
 for f in faces:
  signs=[(f[(i+1)%3]-f[i])%7 in {1,2,4} for i in range(3)]
  ck(len(set(signs))==1,'face is Paley cyclic triple')
 ck(codes[0]&codes[1]=={0,127},'two Hamming codes intersection')
 ck(tuple(sum(t[i] for t in ((3,1),(3,1),(1,-2))) for i in range(2))==(7,0),'nonzero cycle lift defect')
 print(json.dumps({'status':'PASS','checks':CHECKS,'Q5':power(Q,5),'paley7_group':len(perms),'torus_squared_distances':[str(x) for x in distances],'scope':'exact finite controls; all-order statements proved in note'},indent=2))
if __name__=='__main__':main()
