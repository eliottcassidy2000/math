#!/usr/bin/env python3
"""Independent direct encoding of Q_inv(n,t) (fully gap-determined) for small n,t: explicit triangles and
explicit binary subgrids of [t]^n.  Also the row-invariant (THM-453 F1) variant for n=2."""
import sys, itertools
from pysat.solvers import Solver
n=int(sys.argv[1]); t=int(sys.argv[2])
pts=sorted(itertools.product(range(t),repeat=n))
def gap(x,y): return tuple(y[i]-x[i] for i in range(n))
gid={}
def g(x,y):
    d=gap(x,y)
    if d not in gid: gid[d]=len(gid)+1
    return gid[d]
s=Solver(name='cadical195')
for x,y,z in itertools.combinations(pts,3):   # pts sorted lex => x<y<z
    s.add_clause([-g(x,y),-g(y,z),-g(x,z)])
def subgrids(h):
    if h==0:
        yield ((),); return
    for b0,b1 in itertools.combinations(range(t),2):
        for X0 in subgrids(h-1):
            for X1 in subgrids(h-1):
                yield tuple((b0,)+x for x in X0)+tuple((b1,)+x for x in X1)
cnt=0
for L in subgrids(n):
    s.add_clause(sorted({g(L[i],L[j]) for i in range(len(L)) for j in range(i+1,len(L))})); cnt+=1
ok=s.solve()
print(f"fully gap-determined Q_inv({n},{t}): {'SAT' if ok else 'UNSAT'}  ({len(gid)} gap classes, {cnt} subgrids)")
if n==2:
    # row-invariant: within-row graph R on columns (pairs b<b'), cross relation B_h(b,c) for root gap h, arbitrary
    vid={}
    def var(key):
        if key not in vid: vid[key]=len(vid)+1
        return vid[key]
    def e(x,y):  # x<lex y
        if x[0]==y[0]: return var(('R',x[1],y[1]))
        return var(('B',y[0]-x[0],x[1],y[1]))
    s2=Solver(name='cadical195')
    for x,y,z in itertools.combinations(pts,3):
        s2.add_clause([-e(x,y),-e(y,z),-e(x,z)])
    for L in subgrids(2):
        s2.add_clause(sorted({e(L[i],L[j]) for i in range(4) for j in range(i+1,4)}))
    ok2=s2.solve()
    print(f"row-invariant (THM-453 F1) invQ(2,{t}): {'SAT' if ok2 else 'UNSAT'}")
