#!/usr/bin/env python3
from fractions import Fraction as Q
from itertools import product, combinations
from math import factorial, comb
from functools import lru_cache
from pathlib import Path
import json, hashlib, subprocess, sys
ROOT=Path(__file__).resolve().parent

def need(ok, msg):
    if not ok: raise RuntimeError(msg)

def poly(n, edges):
    out=[0]*n
    for bits in product((0,1),repeat=len(edges)):
        blocks=[{v} for v in range(n)]
        good=True
        for bit,(u,v) in zip(bits,edges):
            if not bit: continue
            x=next(b for b in blocks if u in b); y=next(b for b in blocks if v in b)
            if x is y: good=False; break
            blocks.remove(x); blocks.remove(y); blocks.append(x|y)
        if good: out[sum(bits)]+=1
    return out

def points(n,edges,t):
    # Integer endpoint allocation in each dilation segment; cycle rounding
    # establishes these equal lattice points, independently audited in REPORT.
    out={(0,)*n}
    for u,v in edges:
        nxt=set()
        for d in out:
            for k in range(t+1):
                e=list(d); e[u]+=k; e[v]+=t-k; nxt.add(tuple(e))
        out=nxt
    return out

def contracted(n,es,u,v):
    labels=[i for i in range(n) if i!=v]
    def f(x): return labels.index(u if x==v else x)
    return [(f(a),f(b)) for a,b in es]

checks=0
for mult in product(range(3),repeat=3):
 for loop in range(2):
    es=[e for e,m in zip([(0,1),(0,2),(1,2)],mult) for _ in range(m)]+[(0,0)]*loop
    p=poly(3,es)
    for t in range(4):
        count=len(points(3,es,t))
        need(count==sum(c*t**k for k,c in enumerate(p)),('multigraph Ehrhart',es,t))
        for i,(u,v) in enumerate(es):
            if u==v: continue
            deletion=es[:i]+es[i+1:]
            rhs=len(points(3,deletion,t))+t*len(points(2,contracted(3,deletion,u,v),t))
            need(count==rhs,('deletion contraction',es,t)); checks+=1
# Explicit no-integer fibre countercheck: formula still gains exactly t.
for lo,hi in [(Q(1,3),Q(2,3)),(Q(-5,3),Q(-4,3)),(Q(2),Q(2)),(Q(1,7),Q(15,7))]:
 for t in range(5):
    ceil=lambda x:-((-x.numerator)//x.denominator)
    floor=lambda x:x.numerator//x.denominator
    need(floor(hi+t)-ceil(lo)+1-(floor(hi)-ceil(lo)+1)==t,'interval gain')
rows=json.loads((ROOT/'inputs/forest/certificate.json').read_text())['rows']
for row in rows:
    n=row['n']; p=row['forest_coefficients']; m=(n-1)*(n-2)//2
    need(p[0]==1,'constant')
    if n>=3: need(p[1]==m,'linear')
    if n>=4: need(p[2]==comb(m,2),'quadratic')
    if n>=5:
        es=set((u,v) for u in range(n) for v in range(u+2,n))
        triangles=sum(all(e in es for e in combinations(vs,2)) for vs in combinations(range(n),3))
        need(p[3]==comb(m,3)-triangles,'cubic')
    if n<=7: need(poly(n,[(u,v) for u in range(n) for v in range(u+2,n)])==p+[0]*(n-len(p)),'independent edge subsets')
# Independent formal fixed-point iteration, not exponential differentiation.
N=32
A=[Q(1)]+[Q((m+1)**(m-1),factorial(m)) for m in range(1,N+1)]
def mul(a,b): return [sum(a[j]*b[k-j] for j in range(k+1)) for k in range(N+1)]
def exponential(a):
    out=[Q(1)]+[Q(0)]*N; power=out[:]
    for k in range(1,N+1):
        power=mul(power,a)
        out=[x+y/factorial(k) for x,y in zip(out,power)]
    return out
B=[Q(1)]+[Q(0)]*N
for iteration in range(N): B=exponential([Q(0)]+B[:-1])
need(A==B,'formal fixed point')
# Exception-safe negative control is active under -O.
try: need(False,'hostile sentinel')
except RuntimeError as err: need(str(err)=='hostile sentinel','wrong exception')
else: raise RuntimeError('negative check disabled')
result={'status':'PASS','multigraph_deletion_contraction_checks':checks,'multigraph_universe':'3 labelled vertices; each nonloop multiplicity 0..2; loop multiplicity 0..1; dilation t=0..3','independent_full_forest_rows_through':7,'coefficient_structure_checks_through':14,'formal_tree_fixed_point_degree':N,'exception_safe_negative_control':'PASS'}
print(json.dumps(result,indent=2))
