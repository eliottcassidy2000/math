#!/usr/bin/env python3
"""Independent exact metric, quotient and shortest-displacement controls."""
from itertools import product,combinations
from fractions import Fraction as F
import json
N=0
def ck(x,m):
 global N
 N+=1
 if not x:raise RuntimeError(m)
def matmul(A,B):return tuple(tuple(sum(A[i][k]*B[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def transpose(A):return tuple(zip(*A))
def vect(A,v):return tuple(sum(A[i][j]*v[j] for j in range(2)) for i in range(2))
def norm(v):a,b=v;return a*a-a*b+b*b

def main():
 Q=((0,1),(1,1));Nform=((2,1),(1,-2));G=((2,-1),(-1,2));R=((0,-1),(1,-1));I=((1,0),(0,1))
 ck(matmul(transpose(Q),matmul(Nform,Q))==tuple(tuple(-x for x in r) for r in Nform),'split norm sign')
 ck(matmul(transpose(R),matmul(G,R))==G,'metric rotation')
 ck(matmul(R,matmul(R,R))==I,'order3')
 reps=[(3,1),(-1,2),(2,3),(-2,-3),(1,-2),(-3,-1)]
 for d,z in enumerate(reps,1):
  ck(tuple(t%7 for t in z)==((3*d)%7,d),'short lift residue')
  ck(norm(z)==7,'attained norm lower bound')
  for a,b in product(range(-15,16),repeat=2):
   u=(3*d+7*a,d+7*b)
   ck(norm(u)%7==0 and norm(u)>=7,'all tested lifts obey exact modular lower bound')
  rv=vect(R,z)
  ck(tuple(t%7 for t in rv)==((6*d)%7,(2*d)%7),'rotation multiplication2')
 ck(tuple(reps[0][i]+reps[0][i]+reps[4][i] for i in range(2))==(7,0),'nonzero cycle lift debt')
 # Verify all21 affine label permutations preserve the complete torus metric
 # and the proper Paley orientation, not identify their planar lifts.
 perms=[]
 for a in (1,2,4):
  for b in range(7):
   p=tuple((a*x+b)%7 for x in range(7));perms.append(p)
   for x,y in combinations(range(7),2):
    d=(y-x)%7;e=(p[y]-p[x])%7
    ck(norm(reps[d-1])==norm(reps[e-1]),'metric permutation')
    ck((d in(1,2,4))==(e in(1,2,4)),'Paley permutation')
 ck(len(set(perms))==21,'faithfulness')
 # Congruence labels mod11 plus explicit hyperbolic fifth iterate.
 Q5=I
 for _ in range(5):Q5=matmul(Q,Q5)
 ck(Q5==((3,5),(5,8)),'Q5')
 for a,b in product(range(11),repeat=2):
  z=vect(Q,(a,b));ck((z[0]+4*z[1])%11==4*(a+4*b)%11,'label clock')
 ck(pow(4,5,11)==1 and all(pow(4,j,11)!=1 for j in range(1,5)),'label order5')
 # Integertrace restriction: odd finite planar rotation orders divide3.
 permitted=[]
 for t in range(-2,3):
  A=((0,-1),(1,t));P=I;order=None
  for k in range(1,13):
   P=matmul(A,P)
   if P==I:order=k;break
  if order:permitted.append(order)
 ck(permitted==[3,4,6],'nontrivial finite companion rotations')
 # The two oriented face systems of the seven-vertex torus triangulation.
 faces1=[(r,(r+1)%7,(r+3)%7) for r in range(7)]
 faces2=[(r,(r+3)%7,(r+2)%7) for r in range(7)]
 from collections import Counter
 directed=Counter((f[i],f[(i+1)%3]) for f in faces1+faces2 for i in range(3))
 ck(len(directed)==42 and all(x==1 for x in directed.values()),'each directed edge once')
 ck(all(directed[u,v]==directed[v,u] for u,v in directed),'oriented face gluing')
 ck(7-21+14==0,'Euler torus')
 for v in range(7):
  link={u:set() for u in range(7) if u!=v}
  for f in faces1+faces2:
   if v in f:
    u,w=[x for x in f if x!=v];link[u].add(w);link[w].add(u)
  ck(all(len(x)==2 for x in link.values()),'link degree2')
  seen={next(iter(link))};todo=list(seen)
  while todo:
   for w in link[todo.pop()]:
    if w not in seen:seen.add(w);todo.append(w)
  ck(len(seen)==6,'link connected cycle')
 codes=[]
 for faces in(faces1,faces2):
  pairs=Counter(tuple(sorted(e)) for f in faces for e in combinations(f,2))
  ck(len(pairs)==21 and all(x==1 for x in pairs.values()),'Fano pair incidence')
  code={0}
  for f in faces:
   z=sum(1<<i for i in f);code|={x^z for x in list(code)}
  ck(Counter(x.bit_count() for x in code)==Counter({0:1,3:7,4:7,7:1}),'Hamming weights')
  codes.append(code)
 ck(codes[0]&codes[1]=={0,127},'two code intersection')
 print(json.dumps(dict(status='PASS',checks=N,group_order=21,torus_squared_distance='2/7',lift_cycle=[1,0],finite_companion_orders=permitted,torus_faces=14,hamming_intersection=[0,127]),indent=2))
if __name__=='__main__':main()
