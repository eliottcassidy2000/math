#!/usr/bin/env python3
"""Exact colored-lattice descent controls and complete three-ray proof head.

Reuses the earlier congruence-row and literal-sheet engines transparently;
the subgroup-body tests and lattice-basis/classification logic are new.
"""

import importlib.util
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path


CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def det(u,v):
    return u[0]*v[1]-u[1]*v[0]


def primitive2(u):
    g=gcd(abs(u[0]),abs(u[1]))
    v=tuple(x//g for x in u)
    return v if next(x for x in v if x)>0 else tuple(-x for x in v)


def third_coordinates(vectors, color, unit_det):
    basis=next((u,v) for u,v in combinations(vectors,2) if abs(det(u,v))==unit_det)
    u,v=(tuple(-x for x in t) if color(t)==2 else t for t in basis)
    z=next(t for t in vectors if t!=u and t!=tuple(-x for x in u)
           and t!=v and t!=tuple(-x for x in v))
    denominator=det(u,v)
    m,n=Q(det(z,v),denominator),Q(det(u,z),denominator)
    need(m.denominator==n.denominator==1,('basis integral coordinates',vectors,m,n))
    m,n=int(m),int(n)
    need(gcd(abs(m),abs(n))==1 and (m+n)%3!=0,('primitive live third vector',m,n))
    need(abs(m)==1 or abs(n)==1,('three-direction geometry',vectors,m,n))
    if abs(n)!=1:
        m,n=n,m
    if n<0:
        m,n=-m,-n
    need(n==1 and m in (1,-2),('complete three-ray normal forms',vectors,m,n))
    return m


def lattice_body_controls():
    kinds=Counter()
    three_types=Counter()
    cap_count=0
    body_count=0
    characters=[(mod,a,b) for mod in (2,3,4,5) for a in range(mod) for b in range(mod)
                if (a,b)!=(0,0)]
    characters.append(('noncyclic',0,0))
    for A,B,(p,q),C in product((3,5,7),(3,5,7),((1,1),(1,2),(2,-1),(3,1)),(3,5,9)):
        def inside(x):
            return 2*abs(x[0])<A and 2*abs(x[1])<B and 2*abs(p*x[0]+q*x[1])<C
        points=[x for x in product(range(-3,4),repeat=2) if inside(x)]
        for mod,a,b in characters:
            def dead(x):
                return (x[0]%2==x[1]%2==0) if mod=='noncyclic' else (a*x[0]+b*x[1])%mod==0
            live=[x for x in points if not dead(x)]
            vectors=sorted({primitive2(x) for x in live})
            need(all(inside(x) and not dead(x) for x in vectors),'primitive reduction preserves live body')
            if len(vectors)<2:
                continue
            body_count+=1
            pairs=[(abs(det(u,v)),u,v) for u,v in combinations(vectors,2)]
            need(min(t[0] for t in pairs)==1,('live basis',A,B,p,q,C,mod,a,b,vectors))
            if mod==3 and len(vectors)==3:
                shape=third_coordinates(vectors,lambda x:(a*x[0]+b*x[1])%3,1)
                three_types[shape]+=1
            if mod==3:
                cap=all(det((v[0]-u[0],v[1]-u[1]),(z[0]-u[0],z[1]-u[1]))!=0
                        for u,v,z in combinations(live,3))
                if cap:
                    cap_count+=1
                    need(len(live)<=6 and not any(x[0]%2==x[1]%2==0 for x in live),
                         'full cap excludes zero parity bucket')
                    if len(vectors)==3:
                        need(shape==1,'three-ray cap has additive circuit only')
            hostile_pair=next(((u,v) for index,u,v in pairs if index>1),None)
            if hostile_pair is None:
                continue
            u,v=hostile_pair
            denominator=det(u,v)
            # Every centered-parallelogram point lies in this fixed box.
            found=None
            for x in product(range(-3,4),repeat=2):
                t,s=Q(det(x,v),denominator),Q(det(u,x),denominator)
                if abs(t)<=Q(1,2) and abs(s)<=Q(1,2) and (t.denominator>1 or s.denominator>1):
                    found=x,t,s
                    break
            need(found is not None,('centered quotient witness',u,v))
            x,t,s=found
            need(t and s and inside(x),('centered witness primitive and convex',u,v,x,t,s))
            if not dead(x):
                newpair=(u,x)
                kinds['live representative']+=1
            elif abs(s)<=abs(t):
                y=tuple(u[i]-(1 if t>0 else -1)*x[i] for i in range(2))
                newpair=(u,y)
                kinds['dead representative repair']+=1
            else:
                y=tuple(v[i]-(1 if s>0 else -1)*x[i] for i in range(2))
                newpair=(v,y)
                kinds['dead representative repair']+=1
            need(all(inside(y) and not dead(y) for y in newpair),'repair stays in native live body')
            need(0<abs(det(*newpair))<abs(denominator),'strict determinant descent')
    need(all(kinds.values()) and len(kinds)==2,'both proof branches exercised')
    # An arbitrary deleted set is not a subgroup: it can kill every basis.
    nongroup_live=((1,1),(1,-1),(-1,1),(-1,-1))
    need(min(abs(det(u,v)) for u,v in combinations(nongroup_live,2) if det(u,v))==2,
         'non-subgroup deletion hostile')
    # Convexity without symmetry around zero does not suffice.
    nonsymmetric_live=((1,0),(1,2))
    need(abs(det(*nonsymmetric_live))==2,'noncentral rectangle / parity subgroup hostile')
    print('GENERIC LATTICE BODIES',body_count,'descent controls',dict(kinds),'three-ray forms',dict(three_types))
    print('GENERIC NONCOLLINEAR INDEX-THREE CAPS',cap_count,'all N<=6; zero parity absent')


def load(filename,name):
    path=Path(__file__).with_name(filename)
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def lrc_controls():
    rows=load('lrc14_two_ray_overnight_hexagon_sep05.py','rows')
    literal=load('lrc14_one_ray_overnight_hexagon_sep05.py','literal')
    types=Counter()
    total=selected=basis_rows=0
    leader=(Q(0),None,None)
    digest=sha256()
    def analyze(w,live,check_projection):
        vectors=sorted({rows.primitive(C) for C in live})
        need(len(vectors)>=2,'multi-ray control')
        need(any(abs(det(u,v))==3*w[2] for u,v in combinations(vectors,2)),
             ('native owner-lattice basis',w,vectors))
        if len(vectors)!=3:
            return None
        shape=third_coordinates(vectors,lambda v:(w[0]*v[0])%3,3*w[2])
        ms=[max(map(abs,v)) for v in vectors]
        need(all(sum(map(abs,v))>=16 and m>=7 for v,m in zip(vectors,ms)),
             ('multi-ray coefficient floor',w,vectors))
        P=Q(3*w[2],2)
        need(all(x*y>=P for x,y in combinations(ms,2)),('determinant product',w,ms))
        S=sum((Q(1,m) for m in ms),Q(0))
        need(S<=Q(1,7)+14/P or S*S<=9/P,('three reciprocal envelope',w,ms))
        projections=rows.projections(w,live)
        need(max(projections)<Q(12,49)*S+Q(12,7*w[2]),('three-ray counting envelope',w))
        if check_projection:
            need(projections==literal.literal_projection_data(w)[0],('physical projection replay',w))
        need(max(projections)<Q(6,77),('three-ray network target',w,projections))
        return shape,projections
    eligible=[v for v in range(1,99,2) if v%3]
    for w in combinations(eligible,3):
        if gcd(gcd(w[0],w[1]),w[2])!=1:
            continue
        total+=1
        live=rows.row_carriers(w)
        vectors={rows.primitive(C) for C in live}
        if len(vectors)<2:
            continue
        basis_rows+=1
        result=analyze(w,live,True)
        if result is None:
            continue
        selected+=1
        shape,projections=result
        types[shape]+=1
        if max(projections)>leader[0]:
            leader=max(projections),w,projections
        digest.update(repr((w,tuple(sorted(live)),shape,projections)).encode())
    need(selected==1791,('complete three-ray head',selected))
    need(leader[0]==Q(18,301) and leader[1]==(5,37,43),('finite head leader',leader))
    need(Q(96,26411)<Q(4,1089),'radical tail at c99 after positive squaring')
    need(Q(12,343)+Q(4,99)==Q(2560,33957)<Q(6,77),'linear tail at c99')
    print('LRC PROOF HEAD c<99',total,'primitive eligible;',basis_rows,'multi-ray basis controls;',
          selected,'three-ray;',dict(types))
    print('THREE-RAY HEAD LEADER',leader)
    for w in ((19,23,29),(41,47,49),(5,191,199),(7,611,613)):
        live=rows.row_carriers(w)
        need(live==literal.carriers(w),('wide row/box check',w))
        result=analyze(w,live,True)
        print('WIDE CONTROL',w,'N',len(live),'directions',len({rows.primitive(C) for C in live}),
              'three-ray result',result)
        if w==(41,47,49):
            need(result[0]==-2,'incoming non-A2 circuit hostile retained')
    need(len(rows.row_carriers((5,191,199)))==16,'three-ray multiplicity hostile to N<=8')
    print('SEMANTIC SHA256',digest.hexdigest())
    print('HELPER CHECKS',rows.CHECKS,literal.CHECKS)


def main():
    lattice_body_controls()
    lrc_controls()
    print('PROVED: native live lattice basis; exactly-three-ray projective classification; every E_i<6/77')
    print('CHECKS',CHECKS)
    print('OPEN: four or more primitive directions; weighted arbitrary-support geometry; entry; LRC14')


if __name__=='__main__':
    main()
