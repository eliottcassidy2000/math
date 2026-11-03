#!/usr/bin/env python3
"""Independent exact non-odd-cubic obstruction; Python standard library only.
Universe f(k)=A*k**3 + k**2 + C*k + D, k=-4,...,4.
B!=0 is normalized to 1 by scaling entries; D cancels.
Every possible center index and opposite-pair matching is covered.
"""
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json

def require(v, message):
    if not v:
        raise RuntimeError(message)

def perfect_matchings(xs, allowed):
    if not xs:
        yield []
        return
    a=xs[0]
    for b in xs[1:]:
        if (a,b) in allowed:
            for rest in perfect_matchings([x for x in xs[1:] if x!=b],allowed):
                yield [(a,b)]+rest

def quartet_is_magic(p):
    for pp in permutations(p):
        for signs in product((-1,1),repeat=4):
            u,v,z,w=[x*s for x,s in zip(pp,signs)]
            if u+v==z and u-v==w:
                return True
    return False

out=[]
for j in range(-4,5):
    if j==0:
        # Center=mean requires 0=20/3 after B=1 normalization.
        continue
    xs=[k for k in range(-4,5) if k!=j]
    edges=[]
    # Center condition forces C = (20/3-j^2)/j - A*j^2.
    cc=F(20-3*j*j,3*j)
    for u,v in combinations(xs,2):
        aa=u**3+v**3-j*j*(u+v)
        bb=F(u*u+v*v)+cc*(u+v)-F(40,3)
        edges.append((u,v,aa,bb))
    identities=[(u,v) for u,v,aa,bb in edges if aa==0 and bb==0]
    require(not identities,'Unexpected parameter-free pair edge')
    candidates=sorted({-bb/aa for u,v,aa,bb in edges if aa})
    survivors=[]
    for a in candidates:
        allowed={(u,v) for u,v,aa,bb in edges if aa*a+bb==0}
        matching=list(perfect_matchings(xs,allowed))
        if not matching:
            continue
        c=cc-a*j*j
        vals=[a*k**3+k*k+c*k-F(20,3) for k in range(-4,5)]
        pp=sorted(v for v in vals if v>0)
        # Here surviving values have only one zero; verify this explicitly.
        require(vals.count(0)==1 and len(pp)==4,'Unexpected zero degeneracy')
        require(sorted(vals)==sorted([F(0)]+pp+[-p for p in pp]),'Not symmetric')
        magic=quartet_is_magic(pp)
        require(not magic,'Nonodd cubic magic counterexample')
        survivors.append({'A':str(a),'B':'1','C':str(c),
            'positive_centered_quartet':[str(p) for p in pp],
            'perfect_matchings':matching,'magic':magic})
    out.append({'center_index':j,'candidate_A_count':len(candidates),
        'parameter_free_edges':identities,'survivors':survivors})
expected={-4:(F(-2,3),F(13)),-1:(F(1,3),F(-6)),
          1:(F(-1,3),F(6)),4:(F(2,3),F(-13))}
actual={r['center_index']:(F(s['A']),F(s['C'])) for r in out for s in r['survivors']}
require(actual==expected,'Unexpected nonodd symmetric cubic classification')
require(sum(len(r['survivors']) for r in out)==4,'Duplicate or missing survivor')
print(json.dumps({'status':'PASS','scope':'All real cubic coefficient tuples with quadratic coefficient B!=0; repeated entries allowed',
    'normalization':'B=1, D=0 by affine entry transformation',
    'classification':'No such nonodd cubic sample can be a 3x3 magic square',
    'centers':out},indent=2))
