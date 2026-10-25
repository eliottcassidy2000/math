#!/usr/bin/env python3
"""Direct integer matrices versus independently computed path invariants.

Uses SymPy for exact Smith forms, not the producer determinant routine.
"""
import hashlib
import itertools as it
import json
import random
from collections import Counter
from math import gcd
from pathlib import Path
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def invariant(n, subset):
    h = {j for j in subset if (j+9) % n in subset}
    adj = {j:set() for j in range(n)}
    for j in set(range(n))-subset:
        a,b=(j-6)%n,(j-3)%n
        adj[a].add(b);adj[b].add(a)
    seen=set();k=0;o=False;b=True;beta=0
    for j in range(n):
        if j in seen:
            continue
        color={j:0};stack=[j];seen.add(j)
        while stack:
            v=stack.pop()
            for w in adj[v]:
                if w not in color:
                    color[w]=1-color[v];seen.add(w);stack.append(w)
                elif color[w] == color[v]:
                    raise RuntimeError('expected path forest')
        anchors={color[v] for v in color if v not in h}
        if len(anchors)>1:
            b=False
        if not anchors:
            k+=1;o |= bool(len(color)%2)
        zero=next(iter(anchors)) if anchors else 0
        beta+=sum(x != zero for x in color.values())
    d=len(subset)-1+k-int(o)
    torsion=1 if not b else n if o else gcd(n,beta)
    return d,torsion,{'free_paths':k,'odd_free':o,'anchor_compatible':b,'beta':beta}


def physical_matrix(n, subset):
    points=[((-a-3*b)%n,(j-9*b)%n)
            for a,b,j in it.product((0,1),(0,1),sorted(subset))]
    if len(set(points)) != 4*len(subset):
        raise RuntimeError('physical point collision')
    # Direct cross products; the theory instead used divided sloped labels.
    labels=[(y,x,(3*y-9*x)%n) for x,y in points]
    return sp.Matrix([[int(labels[j][axis]==r) for j in range(len(points))]
                      for axis in range(3) for r in range(n)])


def main():
    universe=set()
    for n in (5,7,8,10):
        for mask in range(1,1<<n):
            universe.add((n,tuple(j for j in range(n) if mask>>j&1)))
    for n in (11,13,14,16,17,19,20,23,25,29,31,35,37,41,43,49):
        for length in range(1,n+1):
            universe.add((n,tuple(range(length))))
        rng=random.Random(701+n)
        for _ in range(12):
            mask=rng.randrange(1,1<<n)
            universe.add((n,tuple(j for j in range(n) if mask>>j&1)))
    hostile_keys=[(8,(0,1,2,5)),(16,(0,1,2,4,5,7,8,9)),(49,tuple(range(49)))]
    universe.update(hostile_keys)
    histogram=Counter();hostiles=[];checks=0
    for n,sub in sorted(universe):
        d,D,info=invariant(n,set(sub))
        matrix=physical_matrix(n,set(sub))
        smith=smith_normal_form(matrix,domain=ZZ)
        diagonal=[abs(int(v)) for v in smith.diagonal() if v]
        nonunits=[v for v in diagonal if v!=1]
        if len(diagonal)!=4*len(sub)-d or nonunits!=([D] if D>1 else []):
            raise RuntimeError(('Smith mismatch',n,sub,d,D,diagonal))
        checks+=1;histogram[D]+=1
        # Prime-power kernels detect torsion exponent, rather than only rank mod p.
        for ell in (2,3,5,7):
            for e in (1,2,3):
                q=ell**e
                actual=q**d
                for v in diagonal:
                    actual*=gcd(q,v)
                predicted=q**d*gcd(q,D)
                if actual!=predicted:
                    raise RuntimeError('prime-power kernel order mismatch')
                checks+=1
        if (n,sub) in hostile_keys:
            hostiles.append({'N':n,'I':sub,'d':d,'D':D,**info})
    packet=Path(__file__).parent
    receipt=json.loads((packet/'normal.json').read_text())
    matrix_controls={}
    for name in ('C0','C1','C4'):
        m=sp.Matrix(receipt[name]);snf=smith_normal_form(m,domain=ZZ)
        diag=[abs(int(v)) for v in snf.diagonal()]
        if abs(int(m.det()))!= (1 if name=='C1' else 3):
            raise RuntimeError('independent determinant')
        matrix_controls[name]={'det':int(m.det()),'smith':diag}
        checks+=1
    if sp.Matrix(receipt['C1'])*sp.Matrix(receipt['C1_inverse'])!=sp.eye(6):
        raise RuntimeError('independent printed inverse')
    checks+=1
    for detv,c in receipt['determinant_representatives'].items():
        if abs(sp.Matrix(c).det())!=int(detv):
            raise RuntimeError('independent Bareiss audit')
        checks+=1
    result={'status':'FINITE-EXACT','checks':checks,'physical_matrices':len(universe),
            'torsion_histogram':dict(sorted(histogram.items())), 'hostiles':hostiles,
            'matrix_controls':matrix_controls,'sympy_version':sp.__version__,
            'universe':'All nonempty subsets at N=5,7,8,10; all intervals and 12 seeded nonempty masks for N=11,13,14,16,17,19,20,23,25,29,31,35,37,41,43,49; printed hostiles. Duplicate supports removed. No symmetry filtering.',
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
