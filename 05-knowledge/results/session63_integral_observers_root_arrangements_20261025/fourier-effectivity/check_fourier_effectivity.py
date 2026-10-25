#!/usr/bin/env python3
from itertools import product
from math import comb
import json
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
checks=0
def need(c,m):
 global checks
 checks+=1
 if not c:raise RuntimeError(m)
def H(r):return [[-1 if (a&x).bit_count()%2 else 1 for x in range(2**r)] for a in range(2**r)]
def transform(h,x):return tuple(sum(a*b for a,b in zip(row,x)) for row in h)
smith=[]
for r in range(1,7):
 h=sp.Matrix(H(r));N=2**r
 need(h*h==N*sp.eye(N),'orthogonality')
 d=smith_normal_form(h,domain=sp.ZZ)
 expected=sorted(2**j for j in range(r+1) for _ in range(comb(r,j)))
 need(sorted(abs(int(d[i,i])) for i in range(N))==expected,'Smith invariants')
 smith.append({'r':r,'multiplicities':{str(2**j):comb(r,j) for j in range(r+1)}})
h=H(2);N=4;tables=list(product(range(3),repeat=4));T={f:transform(h,f) for f in tables}
for f in tables:
 a=T[f];c=[[a[i^j] for j in range(N)] for i in range(N)]
 for x in range(N):
  need(all(sum(c[i][j]*h[x][j] for j in range(N))==N*f[x]*h[x][i] for i in range(N)),'circulant spectrum')
 for s in range(N):
  b=T[tuple(f[x^s] for x in range(N))]
  need(all(b[i]==h[i][s]*a[i] for i in range(N)),'phase shift')
  need(all(b[i]*b[i]==a[i]*a[i] for i in range(N)),'energy blindness')
 for g in tables:
  b=T[g];expected=transform(h,tuple(f[i]*g[i] for i in range(N)))
  need(all(sum(a[i]*b[k^i] for i in range(N))==N*expected[k] for k in range(N)),'integer convolution')
bad=(1,1,1,-1);a=transform(h,bad)
need(a==(2,2,2,-2) and all(abs(x)<=a[0] for x in a),'scalar-bound hostile')
need(min(transform(h,a))<0,'negative inverse')
controls=[]
for q,l in [(7,2),(7,3),(5,2),(5,3)]:
 U=[x for x in range(q) if x not in (0,1,l)]
 nonsquare=next(a for a in range(1,q) if pow(a,(q-1)//2,q)==q-1)
 counts=[sum(y*y%q==e*x*(x-1)*(x-l)%q for x in U for y in range(1,q)) for e in (1,nonsquare)]
 S=sum(1 if pow(x*(x-1)*(x-l)%q,(q-1)//2,q)==1 else -1 for x in U)
 need(counts==[q-3+S,q-3-S],'direct Legendre roots')
 controls.append({'q':q,'lambda':l,'S':S,'twists':counts})
need([c['twists'] for c in controls]==[[4,4],[0,8],[4,0],[0,4]],'specified controls')
print(json.dumps({'status':'PASS','checks':checks,'smith':smith,'universe':'r1..6 Smith matrices; all81 ternary tables onF2², all6561 pair products and all4 translations; four exact Legendre controls','legendre':controls},indent=2))
