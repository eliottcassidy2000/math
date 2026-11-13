#!/usr/bin/env python3
"""Finite projective/affine observer controls; no imports from prior checkers."""
import hashlib
import itertools
import json
from pathlib import Path

checks = 0

def require(test, label):
    global checks
    checks += 1
    if not test:
        raise RuntimeError(label)


class Field:
    def __init__(self, q):
        if q not in (2, 3, 4, 5, 7):
            raise ValueError('declared fields only')
        self.q = q
        self.p = 2 if q == 4 else q

    def add(self, a, b):
        return a ^ b if self.q == 4 else (a+b) % self.q

    def mul(self, a, b):
        if self.q != 4:
            return a*b % self.q
        z = 0
        while b:
            if b & 1:
                z ^= a
            a <<= 1
            if a & 4:
                a ^= 7  # X^2+X+1
            b >>= 1
        return z

    def inv(self, a):
        if not a:
            raise ValueError('inverse of zero')
        return next(b for b in range(1,self.q) if self.mul(a,b) == 1)

    def sum(self, values):
        z=0
        for a in values:
            z=self.add(z,a)
        return z


def labels(keys):
    seen = {}
    return [seen.setdefault(k, len(seen)) for k in keys]


def projective_key(x, F):
    lead = next((a for a in x if a), 0)
    if not lead:
        return x
    inverse = F.inv(lead)
    return tuple(F.mul(a,inverse) for a in x)


def refine(partition, operations):
    return labels([tuple([partition[x]]+[partition[f[x]] for f in operations])
                   for x in range(len(partition))])


def chain(partition, operations):
    result = [partition]
    while True:
        nxt = refine(result[-1],operations)
        if nxt == result[-1]:
            return result
        result.append(nxt)


def apply_matrix(matrix, x, F):
    return tuple(F.sum(F.mul(a,b) for a,b in zip(row,x)) for row in matrix)


def states(q, d):
    F=Field(q)
    V=list(itertools.product(range(q),repeat=d))
    index={x:i for i,x in enumerate(V)}
    P=labels([projective_key(x,F) for x in V])
    require(len(set(P))==1+(q**d-1)//(q-1),'projective cells plus zero')
    for x in V:
        orbit={tuple(F.mul(c,a) for a in x) for c in range(1,q)}
        observed={y for y in V if projective_key(y,F)==projective_key(x,F)}
        require(orbit==observed,'canonical labels exactly scalar orbits')
    return F,V,index,P


def operations(F,V,index,matrix,v):
    B=[index[apply_matrix(matrix,x,F)] for x in V]
    A=[index[tuple(F.add(a,b) for a,b in zip(x,v))] for x in V]
    return B,A


summary=[]
for q in (2,3,4,5):
    F,V,index,P=states(q,2)
    invertible=positive=preserving=0
    for entries in itertools.product(range(q),repeat=4):
        matrix=(entries[:2],entries[2:])
        image=[apply_matrix(matrix,x,F) for x in V]
        if len(set(image))!=q*q:
            continue
        invertible+=1
        for v in V[1:]:
            B,A=operations(F,V,index,matrix,v)
            require(refine(P,[B])==P,'every invertible linear map descends to projective quotient')
            if projective_key(apply_matrix(matrix,v,F),F)==projective_key(v,F):
                preserving+=1
                continue
            positive+=1
            trajectory=chain(P,[B,A])
            got=[len(set(r)) for r in trajectory]
            expected=([q*q] if q==2 else [q+2,q*q] if q==3 else [q+2,q*q-q+3,q*q])
            require(got==expected,'line-moving theorem counts and strict depths')
            if q>=4:
                first=trajectory[1]
                large=[tuple(V[i] for i,a in enumerate(first) if a==label)
                       for label in set(first) if first.count(label)>1]
                minus_v=tuple(F.mul(q-1,a) for a in v) if q!=4 else v
                expected_residual={tuple(F.mul(c,a) for a in v) for c in range(q)}-{(0,0),minus_v}
                require(len(large)==1 and set(large[0])==expected_residual,'unique first-round residual cell')
    require(invertible==(q*q-1)*(q*q-q),'GL2 complete enumeration')
    require(positive+preserving==invertible*(q*q-1),'every nonzero marked vector covered')
    summary.append({'q':q,'dimension':2,'invertible_matrices':invertible,
                    'line_moving_cases':positive,'line_preserving_cases':preserving,
                    'line_moving_counts':([q*q] if q==2 else [q+2,q*q] if q==3 else [q+2,q*q-q+3,q*q])})

hostiles=[]
for q in (2,3,4,5,7):
    F,V,index,P=states(q,2)
    v=(1,0)
    B,A=operations(F,V,index,((1,0),(0,1)),v)
    result=chain(P,[B,A])
    counts=[len(set(r)) for r in result]
    expected={2:[4],3:[5,9],4:[6,15],5:[7,23,24,25],7:[9,45,46,47,48,49]}[q]
    require(counts==expected,'line-preserving hostile')
    # Direct word observations for identity B are successive translations.
    funcs=[list(range(len(V)))]
    for depth,part in enumerate(result):
        direct=labels([tuple(P[f[x]] for f in funcs) for x in range(len(V))])
        require(direct==part,'translation word description of hostile depth')
        funcs.append([A[funcs[-1][x]] for x in range(len(V))])
    hostiles.append({'q':q,'dimension':2,'B':'identity','counts':counts,
                     'remaining_nonsingleton_sizes':sorted(result[-1].count(k) for k in set(result[-1]) if result[-1].count(k)>1)})

# The earlier 64-state norm model has q=4,d=3. This control uses an actual
# coordinate-cycle matrix, not multiplication in the prior sextic encoding.
F,V,index,P=states(4,3)
B,A=operations(F,V,index,((0,0,1),(1,0,0),(0,1,0)),(1,0,0))
trajectory=chain(P,[B,A])
require([len(set(r)) for r in trajectory]==[22,63,64],'dimension3 independent model')
words=[list(range(64))]
for part in trajectory:
    direct=labels([tuple(P[f[x]] for f in words) for x in range(64)])
    require(direct==part,'all two-letter bounded words independently recover refinement')
    words += [[op[f[x]] for x in range(64)] for f in list(words) for op in (B,A)]

print(json.dumps({'status':'FINITE-EXACT','checks':checks,
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'universe':'all GL2(Fq) and all nonzero v for q=2,3,4,5; explicit identity controls q=2,3,4,5,7; coordinate-cycle q=4,d=3',
 'field4':'F2[X]/(X^2+X+1), polynomial coefficient bits',
 'complete_GL2_controls':summary,'line_preserving_controls':hostiles,
 'independent_dimension3_counts':[22,63,64],
 'scope':'finite controls accompany a general proof candidate; no native tournament-action claim'},indent=2,sort_keys=True))
