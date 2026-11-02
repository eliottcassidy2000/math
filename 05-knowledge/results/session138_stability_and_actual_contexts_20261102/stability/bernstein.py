"""Exact tensor degree-two Bernstein certificates in orientation orthants.
Build dense coefficients via determinant and Hamilton permutation expansions,
no imported research kernel. Orthant z∈[0,1]^m maps p=(1+s*z)/2.
"""
import numpy as np
from itertools import permutations,combinations
from fractions import Fraction as F
import math,json,sys

def polys(n):
 edges=list(combinations(range(n),2));m=len(edges);idx={p:k for k,p in enumerate(edges)};pw=3**np.arange(m)
 h=np.zeros(3**m,dtype=np.int64);u=h.copy();v=h.copy();var=h.copy()
 def add_terms(out,arcs,coef):
  terms={0:coef}
  for i,j in arcs:
   k=idx[min(i,j),max(i,j)];s=1 if i<j else -1
   new={}
   for exp,c in terms.items():
    new[exp]=new.get(exp,0)+c
    new[exp+pw[k]]=new.get(exp+pw[k],0)+s*c
   terms=new
  for exp,c in terms.items():out[exp]+=c
 for p in permutations(range(n)):
  add_terms(h,list(zip(p,p[1:])),2)
  inversions=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n));sgn=(-1)**inversions
  arcs=[(i,p[i]) for i in range(n) if i!=p[i]]
  add_terms(u,arcs,sgn*2**(n-len(arcs)))
  if len(arcs)==n:add_terms(v,arcs,sgn*(-1)**n)
 for k in range(m):var[0]+=2**(n-2);var[2*pw[k]]-=2**(n-2)
 return h,u,v,var

def bernstein(a,m,orientation):
 # every axis transform doubled (integral), final Bernstein coefficients /2**m
 b=a.reshape((3,)*m).copy()
 for k in range(m):
  axis=m-1-k;s=1 if orientation>>k&1 else -1
  b=np.moveaxis(b,axis,0);c=b.copy();b[0]=2*c[0];b[1]=2*c[0]+s*c[1];b[2]=2*c[0]+2*s*c[1]+2*c[2];b=np.moveaxis(b,0,axis)
 return b.ravel()

def representatives(n):
 edges=list(combinations(range(n),2));m=len(edges);unseen=set(range(1<<m));out=[]
 while unseen:
  b=min(unseen);out.append(b)
  for p in permutations(range(n)):
   c=0
   for k,(i,j) in enumerate(edges):
    x,y=p[i],p[j];kk=edges.index((min(x,y),max(x,y)));bit=(b>>kk&1)^(x>y);c|=bit<<k
   unseen.discard(c)
 return out
if __name__=='__main__':
 n=int(sys.argv[1]) if len(sys.argv)>1 else 5;k=F(sys.argv[2]) if len(sys.argv)>2 else F(7,8)
 h,u,v,V=polys(n);m=n*(n-1)//2
 for o in representatives(n):
  for s in (-1,1):
   a=(h-u+s*v)*k.denominator-k.numerator*V;b=bernstein(a,m,o)
   print(json.dumps({'n':n,'k':str(k),'orientation':o,'sign':s,'min':int(b.min()),'negatives':int((b<0).sum()),'zero':int((b==0).sum()),'scale':2**(n+m)*k.denominator,'argmin':int(b.argmin())}),flush=True)
