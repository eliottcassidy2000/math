#!/usr/bin/env python3
"""Exact all-pair edges and explicit global-colouring certificates. Standard library only."""
from fractions import Fraction as F
from itertools import product,combinations
import json
from math import lcm
B=[1,3,11,33,47,141,517,1551]
Z=(F(0),)*8
ONE=(F(1),)+(F(0),)*7

def scalar(q):return (F(q),)+(F(0),)*7

def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def mul(a,b):
 out=[F(0)]*8
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:out[i^j]+=x*y*B[i&j]
 return tuple(out)
def cmul(u,v):
 a,b=u;c,d=v
 return add(mul(a,c),neg(mul(b,d))),add(mul(a,d),mul(b,c))
def cpow(u,e):
 if e<0:u=(u[0],neg(u[1]));e=-e
 out=(ONE,Z)
 for _ in range(e):out=cmul(out,u)
 return out
r1=(scalar(F(5,6)),tuple(F(1,6) if k==2 else F(0) for k in range(8)))
r2=(scalar(F(23,24)),tuple(F(1,24) if k==4 else F(0) for k in range(8)))

def vertices(m):
 points={}
 for a,b in product(range(-m,m+1),repeat=2):
  rot=cmul(cpow(r1,a),cpow(r2,b))
  for i,j in product(range(-m,m+1),repeat=2):
   v=(scalar(F(2*i+j,2)),tuple(F(j,2) if k==1 else F(0) for k in range(8)))
   p=cmul(rot,v)
   points.setdefault(p,[]).append((a,b,i,j))
 den=lcm(*(q.denominator for p in points for co in p for q in co))
 return [[tuple(int(q*den) for q in co) for co in p] for p in points],list(points.values()),den

def unit(u,v,D2):
 dx=tuple(a-b for a,b in zip(u[0],v[0]));dy=tuple(a-b for a,b in zip(u[1],v[1]))
 if sum((a*a+b*b)*k for a,b,k in zip(dx,dy,B))!=D2:return False
 other=[0]*8
 for d in [dx,dy]:
  for i in range(8):
   for j in range(i+1,8):
    other[i^j]+=2*d[i]*d[j]*B[i&j]
 return not any(other[1:])


def power_exponent(q,d):
    if q<d:return None
    k=0
    while q%d==0:q//=d;k+=1
    return k if q==1 else None

def colour(a,b,i,j):
    c=(i-j)%3
    if c:return c
    q=i*i+i*j+j*j
    parity=0
    for exponent,d in [(a,3),(b,12)]:
        k=power_exponent(q,d)
        if k is not None:parity+=exponent//k
    return 3*(parity%2)

def verify(m):
    p,labels,D=vertices(m)
    edges=[(i,j) for i,j in combinations(range(len(p)),2) if unit(p[i],p[j],D*D)]
    cols=[]
    for lab in labels:
        values={colour(*a) for a in lab}
        if len(values)!=1:raise RuntimeError('ill-defined formula')
        cols.append(next(iter(values)))
        if len(lab)>1 and any(a[2:]!=(0,0) for a in lab):raise RuntimeError('nonzero representation collision')
    if any(cols[a]==cols[b] for a,b in edges):raise RuntimeError('improper colouring')
    nonzero_seams=[]
    for u,v in edges:
        a,b,i,j=labels[u][0]; c,d,k,l=labels[v][0]
        if (i,j)==(0,0) or (k,l)==(0,0):continue
        if (a,b)==(c,d):continue
        if (a!=c)+(b!=d)!=1:raise RuntimeError('mixed-index seam')
        if (i,j) not in [(k,l),(-k,-l)]:raise RuntimeError('different radial vector seam')
        q=i*i+i*j+j*j
        base,gap=(3,abs(a-c)) if a!=c else (12,abs(b-d))
        if q!=base**gap:raise RuntimeError('wrong seam norm')
        if (i-j)%3:raise RuntimeError('seam outside zero colour class')
        nonzero_seams.append([u,v,base,gap,q])
    return {'m':m,'vertices':len(p),'edges':len(edges),'common_denominator':D,
            'all_pair_exact_tests':len(p)*(len(p)-1)//2,
            'formula_colouring':cols,'labels':[lab[0] for lab in labels],
            'edge_list':edges,'nonzero_seams':nonzero_seams,
            'colours_used':sorted(set(cols)),'verified':True}

if __name__=='__main__':
    print(json.dumps({'status':'FINITE-EXACT certificates for the global theorem; H-G is REFUTED by the proof',
      'arithmetic_basis':[1,3,11,33,47,141,517,1551],
      'tests':[verify(1),verify(2)],
      'no_float_filters':True,'no_sat_dependency':True},indent=2,sort_keys=True))
