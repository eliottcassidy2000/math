#!/usr/bin/env python3
"""Standalone exact S93 check: no imports from repository; checks survive -O."""
from itertools import product, permutations
import json
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x: raise RuntimeError(msg)
def mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inv(a):return tuple(a.index(i) for i in range(len(a)))
def conj(g,s):return mul(mul(g,s),inv(g))
def root(a):
 out=[-1]*len(a)
 for i in range(len(a)):
  if out[i]>=0:continue
  c=[i];j=a[i]
  while j!=i:c.append(j);j=a[j]
  need(len(c)%2==1,'odd cycles')
  for k,v in enumerate(c):out[v]=c[(k+(len(c)+1)//2)%len(c)]
 return tuple(out)
def transport(t,s):return root(mul(t,s))
def action(g,P):return tuple(g[v] for v in P)
def anti(T,s):return all(T[i][j]==T[s[j]][s[i]] for i in range(len(T)) for j in range(len(T)))
def auto(T,g):return all(T[i][j]==T[g[i]][g[j]] for i in range(len(T)) for j in range(len(T)))
def paths(T):
 n=len(T);out=[]
 def visit(P,rest):
  if not rest:out.append(P);return
  for v in rest:
   if not P or T[P[-1]][v]:visit(P+(v,),rest-{v})
 visit((),set(range(n)));return out
def reversal(P):
 r=[0]*len(P)
 for a,b in zip(P,reversed(P)):r[a]=b
 return tuple(r)
def partner(T,P):
 edges=set(zip(P,P[1:]));n=len(T)
 return tuple(tuple(0 if i==j else T[i][j] if (i,j) in edges or (j,i) in edges else 1-T[i][j] for j in range(n)) for i in range(n))
def relabel(T,g):
 gi=inv(g);return tuple(tuple(T[gi[i]][gi[j]] for j in range(len(g))) for i in range(len(g)))
def run(T,A,s,name,expected):
 n=len(T);e=tuple(range(n));need(all(auto(T,g) for g in A),'automorphisms')
 C=tuple(g for g in A if mul(g,s)==mul(s,g));I=tuple(sorted({conj(g,s) for g in A}))
 need(all(anti(T,t) and mul(t,t)==e for t in I),'marks')
 P=paths(T);X=tuple(p for p in P if reversal(p)==s)
 key=lambda p:min(action(c,p) for c in C)
 B={key(p) for p in X};H=set();nontrivial=0
 for t,u in product(I,repeat=2):
  h=mul(inv(transport(u,s)),mul(transport(u,t),transport(t,s)));H.add(h)
  need(h in C,'centralizer triangle')
  if h!=e:nontrivial+=1
  for p in X:
   hp=action(h,p)
   need(key(hp)==key(p),'same blue presentation')
   need((hp==p)==(h==e),'free path action')
 # Here triangle values themselves form the subgroup, verified rather than assumed.
 need(all(mul(a,b) in H for a,b in product(H,repeat=2)),'triangle value closure')
 need(len({min(action(h,p) for h in H) for p in X})==len(B)*len(C)//len(H),'holonomy fibre index')
 for t in I:
  g=transport(t,s);Ct=tuple(a for a in A if mul(a,t)==mul(t,a))
  Xt={p for p in P if reversal(p)==t}
  need({action(g,p) for p in X}==Xt,'fixed path transport')
  need({conj(g,c) for c in C}==set(Ct),'centralizer transport')
 # Check actual partner equivariance and colour for all paths, one nonidentity automorphism.
 a=next(g for g in A if g!=e)
 blue=0
 for p in P:
  ap=action(a,p);isblue=anti(T,reversal(p));blue+=isblue
  need(anti(T,reversal(ap))==isblue,'black/blue equivariance')
  need(partner(T,ap)==relabel(partner(T,p),a),'native partner equivariance')
  if isblue:need(anti(partner(T,p),reversal(p)),'blue preserved across native line')
 # A four-step actual native word, with simultaneous transport of all guards.
 pole=next(v for v in range(n) if s[v]==v)
 for c in C:
  need(c[pole]==pole,'centralizer fixes named pole')
  current=T;transported=T
  for step in range(4):
   path=paths(current)[-1]
   image=action(c,path)
   need(anti(current,reversal(path))==anti(transported,reversal(image)),'word colours')
   following=partner(current,path)
   need(partner(following,path)==current,'native involution')
   transported=partner(transported,image)
   need(transported==relabel(following,c),'whole word equivariance')
   current=following
 need((len(A),len(C),len(I),len(P),len(X),len(B),len(H))==expected,'known exact totals')
 need(blue==len(I)*len(X),'all blue labelled paths')
 return dict(name=name,aut=len(A),centralizer=len(C),marks=len(I),paths=len(P),fixed_paths=len(X),blue_presentations=len(B),triangle_subgroup=len(H),nonidentity_base_triangles=nontrivial,black_presentations=(len(P)-blue)//len(A))
def main():
 W=tuple(tuple(int((j//3-i//3)%3==1 if i//3!=j//3 else (j-i)%3==1) if i!=j else 0 for j in range(9)) for i in range(9))
 A=tuple(tuple(3*((i//3+a)%3)+(i%3+b[i//3])%3 for i in range(9)) for a in range(3) for b in product(range(3),repeat=3))
 s=tuple(3*(-(i//3)%3)+(-i%3) for i in range(9))
 out=[run(W,A,s,'W9',(81,3,27,3159,33,11,3))]
 T=tuple(tuple(int((j-i)%7 in {1,2,4}) for j in range(7)) for i in range(7))
 # Independently exhaustive full symmetric-group automorphism enumeration.
 A=tuple(g for g in permutations(range(7)) if auto(T,g));s=tuple(-i%7 for i in range(7))
 out.append(run(T,A,s,'Paley7',(21,3,7,189,9,3,1)))
 print(json.dumps(dict(checks=checks,objects=out),sort_keys=True,indent=2))
if __name__=='__main__':main()
