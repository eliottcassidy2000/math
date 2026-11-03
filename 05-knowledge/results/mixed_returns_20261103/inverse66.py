"""Inherited exact inverse, extended only by the verified s78 matrix."""
from functools import lru_cache
import sympy as sp
from support import *
from inverse import physical_cone,norm,INVERSES as OLD_INVERSES
INVERSES=dict(OLD_INVERSES)
p=sp.Matrix(MATS['s78']);d=int(p.det())
INVERSES['s78']=(d,tuple(tuple(int(x) for x in row) for row in p.adjugate().tolist()))
need(d!=0,'invertible extra atom')
weight=(6,9,5)
need(all(5*sum(MATS['s78'][i][j]*weight[j] for j in range(3))>=9*weight[i] for i in range(3)),'norm contraction remains valid')
need(all(MATS['s78'][i][j]<=16*MATS['s78'][i][0] for i in range(3) for j in (1,2)),'bank cone remains valid')

@lru_cache(None)
def decode(v):
    if v==(1,0,0):return ()
    if not physical_cone(v):return None
    opts=[]
    for l,(det,adj) in INVERSES.items():
        z=step(v,adj)
        if any(x%det for x in z):continue
        u=tuple(x//det for x in z)
        if not physical_cone(u) or u[0]>v[0] or 9*norm(u)>5*norm(v):continue
        opts.append((norm(u),COSTS[l],l,u))
    for _,__,l,u in sorted(opts):
        w=decode(u)
        if w is not None:return w+(l,)
    return None
