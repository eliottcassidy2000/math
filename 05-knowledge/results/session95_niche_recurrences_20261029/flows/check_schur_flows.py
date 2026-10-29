#!/usr/bin/env python3
"""Exact small checks for Schur boundary circulation recurrences (stdlib only)."""
from collections import defaultdict
from functools import lru_cache
from itertools import product
import json

CHECKS=0
def check(a,b,tag):
 global CHECKS
 CHECKS+=1
 if a!=b:raise RuntimeError((tag,a,b))
@lru_cache(None)
def parts(t,n=100):
 if t==0:return ((),)
 out=[]
 for first in range(min(t,n),0,-1):
  out.extend((first,)+p for p in parts(t-first,first))
 return tuple(out)
@lru_cache(None)
def schur(shape,n):
 cells=[(i,j) for i,row in enumerate(shape) for j in range(row)]
 out=defaultdict(int);values={};counts=[0]*n
 def visit(pos):
  if pos==len(cells):out[tuple(counts)]+=1;return
  i,j=cells[pos];lo=max(values.get((i,j-1),0),values.get((i-1,j),-1)+1)
  for v in range(lo,n):
   values[i,j]=v;counts[v]+=1;visit(pos+1);counts[v]-=1
  values.pop((i,j),None)
 visit(0)
 return dict(out)
def mul(a,b):
 out=defaultdict(int)
 for u,c in a.items():
  for v,d in b.items():out[tuple(x+y for x,y in zip(u,v))]+=c*d
 return {u:c for u,c in out.items() if c}
@lru_cache(None)
def lr(a,b):
 n=max(sum(a)+sum(b),1);rest=mul(schur(a,n),schur(b,n));out={}
 while rest:
  lead=max(rest);shape=tuple(x for x in lead if x);c=rest[lead]
  if tuple(sorted(lead,reverse=True))!=lead or c<0:raise RuntimeError('Schur triangularity')
  out[shape]=c
  for x,z in schur(shape,n).items():
   rest[x]=rest.get(x,0)-c*z
   if not rest[x]:del rest[x]
 return out
@lru_cache(None)
def splits3(shape):
 out=[];r=sum(shape)
 for x in range(r+1):
  for y in range(r-x+1):
   for a,b,c in product(parts(x),parts(y),parts(r-x-y)):
    coef=sum(v*lr(d,c).get(shape,0) for d,v in lr(a,b).items())
    if coef:out.append(((a,b,c),coef))
 return tuple(out)
def conv(a,b,D):
 out=[0]*(D+1)
 for i,x in enumerate(a[:D+1]):
  for j,y in enumerate(b[:D+1-i]):out[i+j]+=x*y
 return out
def flow_kernel(n,arcs,D):
 a=[{(0,)*n:1}]+[{} for _ in range(D)]
 for u,v in arcs:
  e=[0]*n;e[u]=1;e[v]-=1
  for d in range(1,D+1):
   for b,c in list(a[d-1].items()):
    k=tuple(x+y for x,y in zip(b,e));a[d][k]=a[d].get(k,0)+c
 return a
def boundary(kernel,n,lam,mu):
 if sum(lam)!=sum(mu):return [0]*len(kernel)
 a=schur(lam,n);b=schur(mu,n)
 return [sum(c*e*layer.get(tuple(v-u for u,v in zip(x,y)),0) for x,c in a.items() for y,e in b.items()) for layer in kernel]
def inflate(children):
 n=sum(x[0] for x in children);offset=[0,children[0][0],children[0][0]+children[1][0]]
 arcs=[]
 for i,(size,inside) in enumerate(children):
  arcs.extend((u+offset[i],v+offset[i]) for u,v in inside)
  j=(i+1)%3
  arcs.extend((u+offset[i],v+offset[j]) for u in range(size) for v in range(children[j][0]))
 return n,arcs
def trace_formula(child,D,fixed=False):
 n,arcs=child;kernel=flow_kernel(n,arcs,D);out=[0]*(D+1)
 for t in range(D//3+1):
  bank=[p for p in parts(t) if len(p)<=n]
  M={(a,b):boundary(kernel,n,a,b) for a,b in product(bank,repeat=2)}
  if fixed:
   for a in bank:
    for j,c in enumerate(M[a,a]):
     if 3*(t+j)<=D:out[3*(t+j)]+=c
  else:
   for a,b,c in product(bank,repeat=3):
    v=conv(conv(M[a,b],M[b,c],D-3*t),M[c,a],D-3*t)
    for j,z in enumerate(v):out[3*t+j]+=z
 return out
def closed_boundary(child,lam,mu,D):
 n,arcs=child;kernel=flow_kernel(n,arcs,D);out=[0]*(D+1)
 @lru_cache(None)
 def M(a,b):return boundary(kernel,n,a,b)
 for aa,ca in splits3(lam):
  for bb,cb in splits3(mu):
   for t0 in range(D+1):
    for t1 in range(D-t0+1):
     for t2 in range(D-t0-t1+1):
      ts=(t0,t1,t2);shift=sum(ts)
      if any(sum(aa[i])+ts[i]!=sum(bb[i])+ts[(i-1)%3] for i in range(3)):continue
      for al in product(*(tuple(p for p in parts(t) if len(p)<=n) for t in ts)):
       value=[ca*cb]+[0]*(D-shift)
       for i in range(3):
        local=[0]*(D-shift+1)
        for rho,cr in lr(aa[i],al[i]).items():
         for sig,cs in lr(bb[i],al[(i-1)%3]).items():
          for j,z in enumerate(M(rho,sig)[:D-shift+1]):local[j]+=cr*cs*z
        value=conv(value,local,D-shift)
       for j,z in enumerate(value):out[shift+j]+=z
 return out

def main():
 # Independently materialize Cauchy kernels by geometric factors.
 for n in (1,2):
  D=4;kernel=flow_kernel(2*n,[(i,n+j) for i in range(n) for j in range(n)],D)
  for t in range(D+1):
   rhs=defaultdict(int)
   for p in parts(t):
    for a,c in schur(p,n).items():
     for b,d in schur(p,n).items():rhs[a+tuple(-x for x in b)]+=c*d
   check(kernel[t],dict(rhs),('Cauchy',n,t))
 children=[(1,[]),(2,[(0,1)]),(3,[(0,1),(1,2),(2,0)])]
 raw={}
 for child,D in zip(children,(9,6,5)):
  n,arcs=inflate([child]*3);k=flow_kernel(n,arcs,D)
  direct=[a.get((0,)*n,0) for a in k];formula=trace_formula(child,D)
  check(formula,direct,('trace',child[0]));raw[child[0]]=direct
  for degree in range(D+1):check(formula[degree],direct[degree],('trace coefficient',child[0],degree))
 # Full Littlewood-Richardson closed carrier versus direct constant terms.
 for child,rmax,D in [(children[0],2,4),(children[1],1,3)]:
  n,arcs=inflate([child]*3);k=flow_kernel(n,arcs,D)
  for r in range(rmax+1):
   for a,b in product(parts(r),repeat=2):
    check(closed_boundary(child,a,b,D),boundary(k,n,a,b),('LR carrier',child[0],a,b))
 # Literal root-rotation fixed arc-orbit kernel: each orbit has three arcs.
 for child in children:
  n,arcs=child;D=12
  # Collapsing vertex balance leaves internal child arcs and all ordered pairs,
  # including diagonal cross-block edges, each with total cost 3.
  collapsed=flow_kernel(n,arcs+[(i,j) for i in range(n) for j in range(n)],D//3)
  direct=[0]*(D+1)
  for j,layer in enumerate(collapsed):direct[3*j]=layer.get((0,)*n,0)
  check(trace_formula(child,D,True),direct,('fixed root rotation',n))
 # First literal quotient universe: C3[T2,T2,T2], all degree<=6 flows.
 # Burnside integrality at all coefficients and the triangle port-loss witness.
 fixed=trace_formula(children[1],6,True)
 check(raw[2][3],8,'eight cross triangles')
 check(fixed[3],2,'two root-rotation-fixed cross triangles')
 check((raw[2][3]+2*fixed[3])//3,4,'four root-C3 triangle orbits')
 for a,b in zip(raw[2],fixed):check((a+2*b)%3,0,'Burnside divisibility')
 # Same-size, same scalar R, distinct charged responses and cyclic contexts.
 T=[(1,0),(2,0),(3,0),(4,0),(2,1),(3,1),(4,1),(3,2),(2,4),(4,3)]
 K=flow_kernel(5,T,12)
 check(boundary(K,5,(3,),(2,1))[3],878,'same-scalar forward boundary')
 check(boundary(K,5,(2,1),(3,))[3],874,'same-scalar reversed boundary')
 context=[]
 B=flow_kernel(2,[(0,1)],12)
 for arcs in [T,[(v,u) for u,v in T]]:
  K=flow_kernel(5,arcs,12);z=[0]*13
  for t in range(5):
   row=(t,) if t else ()
   for a in parts(t):
    if len(a)>2:continue
    v=conv(boundary(K,5,row,a),boundary(B,2,a,row),12-3*t)
    for j,c in enumerate(v):z[3*t+j]+=c
  context.append(z)
 check(context[0], [1,0,0,11,25,34,96,270,573,1127,2344,4843,9470],'forward cyclic context')
 check(context[1], [1,0,0,11,25,34,96,270,573,1127,2344,4843,9462],'reversed cyclic context')
 # Exact first-separation attribution: only t=3, shape(2,1) contributes.
 K=flow_kernel(5,T,12);separators=[]
 for t in range(5):
  row=(t,) if t else ()
  for a in parts(t):
   if len(a)>2:continue
   x=boundary(K,5,row,a);xr=boundary(K,5,a,row)
   difference=conv([v-w for v,w in zip(x,xr)],boundary(B,2,a,row),12-3*t)
   if any(difference):separators.append((t,a,difference))
 check(separators,[(3,(2,1),[0,0,0,8])],'first separator attribution')
 # Closed low-degree recurrence on recursive ternary tournaments.
 c3=c4=c5=0;n=1;low=[]
 for k in range(1,7):
  m=n*(n-1)//2;two=n*((n-1)//2)**2
  c3,c4,c5=3*c3+n**3,3*c4+3*n*n*m,3*c5+3*n*n*two+3*n*m*m
  n*=3
  check(c3,n*(n*n-1)//24,('triangle closed form',k))
  check(c4,(4*n**4-13*n**3+9*n)//208,('four-cycle closed form',k))
  from fractions import Fraction
  check(c5,Fraction(n**5-n,160)-Fraction(n**4-n,26)+Fraction(n**3-n,16),('five-cycle closed form',k))
  low.append([k,n,c3,c4,c5])
 check(raw[3][3:6],[30,81,162],'nine-vertex low coefficients')
 check(boundary(flow_kernel(3,children[2][1],12),3,(1,),(1,)),[3]*13,'charged reciprocity hostile')
 print(json.dumps({'status':'PASS','checks':CHECKS,'raw_prefixes':raw,'root_rotation_fixed_T2':fixed,'hostile_context_prefixes':context,'low_cycle_recurrences':low},sort_keys=True))
if __name__=='__main__':main()
