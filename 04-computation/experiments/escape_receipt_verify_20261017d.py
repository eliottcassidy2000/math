#!/usr/bin/env python3
"""Independent combinatorial replay and small involution controls."""
import argparse
import json
from itertools import combinations,product
from pathlib import Path
from functools import lru_cache
from pysat.solvers import Solver


def need(ok,msg):
    if not ok:raise RuntimeError(msg)


def verify(path):
    d=json.loads(Path(path).read_text());n=d['vertices'];edges=d['edges']
    need(len(set(map(tuple,edges)))==len(edges),'distinct edges')
    need(all(0<=a<b<n for a,b in edges),'edge universe')
    port_edges=[(a,b) for a,b in edges if b<8]
    actual=set()
    def grow(a):
        if len(a)==8:actual.add(tuple(a));return
        for c in range(min(6,2+max(a,default=-1))):
            b=a+[c]
            if all(b[i]!=b[j] for i,j in port_edges if j<len(b)):grow(b)
    grow([])
    need(len(actual)==570 and actual==set(map(tuple,d['patterns'])),'independent exhaustive boundary universe')
    checks=0
    for pat,w in zip(d['patterns'],d['witnesses']):
        need(len(w)==n and set(w)<=set('012345'),'complete six-colouring')
        need(tuple(map(int,w[:8]))==tuple(pat),'literal original ports')
        for a,b in edges:
            need(w[a]!=w[b],'all edges in each boundary witness');checks+=1
    adj=[set() for _ in range(n)]
    for a,b in edges:adj[a].add(b);adj[b].add(a)
    alive=set(range(n))
    for a in d['peeling_order']:
        need(a>=8 and a in alive and len(adj[a]&alive)<=5,'certified relative peeling')
        alive.remove(a)
    need(len(alive)-8==d['relative_six_core'],'relative six-core size')
    with Solver(name='g3') as solver:
        for a in range(n):
            solver.add_clause([4*a+c+1 for c in range(4)])
            for c in range(4):
                for b in range(c):solver.add_clause([-4*a-b-1,-4*a-c-1])
        for a,b in edges:
            for c in range(4):solver.add_clause([-4*a-c-1,-4*b-c-1])
        ok=solver.solve()
        witness=None
        if ok:
            positive={a for a in solver.get_model() if a>0}
            witness=[next(c for c in range(4) if 4*a+c+1 in positive) for a in range(n)]
            need(all(witness[a]!=witness[b] for a,b in edges),'four-colour upper witness')
    return dict(file=path,vertices=n,edges=len(edges),boundary_checks=checks,
                relative_six_core=d['relative_six_core'],four_colourable=ok,four_colour_witness=witness)


def cd():
    # Baez's displayed convention: (a,b)(c,d)=(ac-d b*, a* d+c b).
    @lru_cache(None)
    def mult(n,i,j):
        if n==1:return (1,0)
        h=n//2
        if i<h and j<h:
            return mult(h,i,j)
        if i<h:
            sign,k=mult(h,i,j-h);return (sign*(1 if i==0 else -1),h+k)
        if j<h:
            sign,k=mult(h,j,i-h);return (sign,h+k)
        sign,k=mult(h,j-h,i-h)
        return (-sign*(1 if i-h==0 else -1),k)
    found=None
    for a,b in combinations(range(1,16),2):
        for c,d in combinations(range(1,16),2):
            for sign in (-1,1):
                v=[0]*16
                for i,ci in ((a,1),(b,1)):
                    for j,cj in ((c,1),(d,sign)):
                        eps,k=mult(16,i,j);v[k]+=ci*cj*eps
                if not any(v):found=[a,b,c,d,sign];break
            if found:break
        if found:break
    need(found is not None,'sedenion zero-divisor hostile')
    need(all(mult(16,i,j)[1]==i^j for i in range(16) for j in range(16)),'same binary address, signed multiplication')
    # All rotational tournaments are self-converse by v -> -v.
    sets=[]
    for bits in product((0,1),repeat=8):sets.append(frozenset((i+1 if b else 16-i) for i,b in enumerate(bits)))
    orbits=set()
    for A in sets:
        orbit={tuple(sorted((m*a)%17 for a in A)) for m in range(1,17)}
        need(len(orbit)==16,'free multiplier action on presentations')
        need(tuple(sorted((-a)%17 for a in A)) in orbit,'converse lies in the same rotational class')
        orbits.add(min(orbit))
    need(len(orbits)==16,'sixteen multiplier orbits, each a torsor; quotient action trivial')
    return dict(cayley_dickson_convention='Baez (ac-db*,a*d+cb)',zero_divisors=found,
                zero_divisor_input_squared_norms=[2,2],product_squared_norm=0,
                rotational_presentations=256,rotational_classes=16,self_converse_rotational_classes=16)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('files',nargs='+');a=p.parse_args()
    print(json.dumps(dict(status='FINITE-EXACT independent positive certificate replay',
                         geometry=[verify(p) for p in a.files],involutions=cd()),sort_keys=True,separators=(',',':')))
