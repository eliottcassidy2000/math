#!/usr/bin/env python3
"""Six-coordinate unit circuits, with exact support-set contact compilation.

For a Cartesian contact family, a factor's attainable colour differences form
a subset of Z/6. The edge condition is exactly that zero is absent from their
Minkowski sum. This avoids enumerating the Cartesian product of all edge pairs;
it does not weaken the selected graph's boundary extension predicate.
"""
import argparse
import json
import math
from fractions import Fraction as Q
from pathlib import Path
from pysat.solvers import Solver
from golden_forcing_search_20261017 import Geo, recipe
from golden_boundary_20261017c import PORT_EDGES, patterns
from mixed_quad_labels_20261018k import mixed_quad


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--six',action='store_true')
    ap.add_argument('--dense-six',action='store_true')
    ap.add_argument('--mixed-quads',type=int,default=2,choices=[0,1,2])
    ap.add_argument('--patterns',type=int,default=570)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    assert not args.dense_six or args.six
    geo=Geo();f=geo.f;bar0=geo.bar
    A,C=recipe(geo,'base');d=geo.d;a=bar0(d)
    h=f.add(f.gen(2),f.scalar(3));u=f.mul(f.sub(h,geo.j),f.inv(f.add(h,geo.j)))
    Z=f.add(d,u);Y=f.scale(f.mul(Z,f.inv(bar0(Z))),-1)
    m2=f.scale(f.add(f.scalar(20),f.scale(geo.j,11)),Q(-1,19))
    m1=f.scale(f.sub(f.scalar(20),f.scale(geo.j,11)),Q(1,19))
    a2=f.power(a,2)
    c=f.mul(f.add(f.sub(f.add(f.mul(a,m2),f.one),f.mul(a2,m1)),f.mul(a2,a)),f.inv(f.sub(f.mul(d,a),f.one)))
    b=f.add(f.mul(c,d),f.sub(f.mul(a,m1),a2))
    assert f.power(u,2)==f.add(f.add(a,f.mul(b,u)),f.mul(c,Y))
    M=list(geo.M);quads=[]
    for name,k in list(zip(('b','c'),(b,c)))[:args.mixed_quads]:
        pts=mixed_quad(geo,k);quads.append(dict(coefficient=name,points=[[str(t) for t in p] for p in pts]))
        M=list(dict.fromkeys(M+pts))
    sets=[A,C,C,C,C]+([C] if args.six else [])
    def enlarge(axis,co):
        base=A if axis==0 else C
        sets[axis]=list(dict.fromkeys(sets[axis]+[f.add(x,f.mul(co,m)) for x in base for m in M]))
    for axis,co in ((0,a),(1,b),(4,c)):enlarge(axis,co)
    zero=(f.zero,f.zero);one=(f.one,f.zero)
    def plus(x,y):return f.add(x[0],y[0]),f.add(x[1],y[1])
    def scale(x,k):return f.mul(x[0],k),f.mul(x[1],k)
    def times(x,y):
        bb=f.mul(x[1],y[1])
        return f.add(f.mul(x[0],y[0]),f.mul(Y,bb)),f.sub(f.add(f.mul(x[0],y[1]),f.mul(x[1],y[0])),f.mul(Z,bb))
    def bar(x):
        t=f.mul(bar0(x[1]),f.inv(Y))
        return f.add(bar0(x[0]),f.mul(t,Z)),t
    v=(f.zero,f.one)
    dirs=[one,(u,f.zero),v,(f.zero,u),(Y,f.zero)]+([(f.zero,Y)] if args.six else [])
    families=[('golden',[(0,d),(1,f.one),(2,f.one)]),
              ('reflected',[(1,f.one),(2,f.one),(3,a)]),
              ('Y_contact',[(0,d),(1,f.one),(4,a)]),
              ('v_square',[(2,f.scale(d,-1)),(3,f.scalar(-1)),(4,f.one)]),
              ('u_square',[(0,a),(1,b),(4,c)])]
    if args.six:
        for axis,co in ((2,a),(3,b),(5,c)):enlarge(axis,co)
        families.append(('u_square_times_v',[(2,a),(3,b),(5,c)]))
        # uY = r0 + r1*u + r2*Y, followed by uv^2 = uY-d*uv-u^2v.
        r2=f.sub(m2,b)
        r0=f.mul(f.add(f.mul(r2,a),f.one),f.inv(c))
        r1=f.mul(f.add(f.mul(r2,b),f.sub(m1,a)),f.inv(c))
        assert f.mul(u,Y)==f.add(f.add(r0,f.mul(r1,u)),f.mul(r2,Y))
        six_terms=[(0,r0),(1,r1),(2,f.scale(a,-1)),(3,f.scale(f.add(d,b),-1)),(4,r2),(5,f.scale(c,-1))]
        if args.dense_six:
            for axis,co in six_terms:enlarge(axis,co)
        families.append(('u_times_v_square',six_terms))
        # The v-multiples of three inherited unit directions use the new axis.
        families.append(('Y_contact_times_v',[(2,d),(3,f.one),(5,a)]))
    for name,terms in families:
        z=zero
        assert all(co!=f.zero for _,co in terms)
        for axis,co in terms:z=plus(z,scale(dirs[axis],co))
        assert times(z,bar(z))==one, name
        if name=='u_times_v_square':assert z==times((u,f.zero),times(v,v))
    axes=[geo.unit_edges(s) for s in sets]
    print('sizes',list(map(len,sets)),'axis edges',list(map(len,axes)),flush=True)
    maps=[{x:i for i,x in enumerate(s)} for s in sets]
    cache={}
    def pairs(axis,delta):
        key=axis,delta
        if key not in cache:cache[key]=tuple((i,maps[axis][f.add(x,delta)]) for i,x in enumerate(sets[axis]) if f.add(x,delta) in maps[axis])
        return cache[key]
    groups=[];group_names=[];counts={}
    for name,terms in families:
        axis,co=min(terms,key=lambda t:len(sets[t[0]])); inv=f.inv(co)
        scaled=[f.mul(x,inv) for x in sets[axis]];edges=geo.unit_edges(scaled)
        # One representative per q/-q suffices: every edge is unoriented.
        qs=sorted({min(f.sub(scaled[j],scaled[i]),f.sub(scaled[i],scaled[j])) for i,j in edges})
        total=0;ng=0
        for q in qs:
            group=tuple((axis,pairs(axis,f.mul(co,q))) for axis,co in terms)
            if any(not ps for _,ps in group):continue
            total+=math.prod(len(ps) for _,ps in group);ng+=1;groups.append(group);group_names.append(name)
        counts[name]=dict(cartesian_constraints=total,support_groups=ng,arity=len(terms))
        print(name,counts[name],flush=True)
    top=0
    def new():
        nonlocal top
        top+=1;return top
    colours=[[[new() for _ in range(6)] for x in s] for s in sets]
    clauses=[]
    for axis,rows in enumerate(colours):
        for row in rows:
            clauses.append(row)
            for x in range(6):
                for y in range(x):clauses.append([-row[x],-row[y]])
        for i,j in axes[axis]:
            for q in range(6):clauses.append([-rows[i][q],-rows[j][q]])
        if axis:
            assert sets[axis][0]==f.zero;clauses.append([rows[0][0]])
    supports={};sums={}
    def support(axis,ps):
        key=axis,ps
        if key not in supports:
            S=tuple(new() for _ in range(6));supports[key]=S
            for i,j in ps:
                for x in range(6):
                    for y in range(6):clauses.append([-colours[axis][j][x],-colours[axis][i][y],S[(x-y)%6]])
        return supports[key]
    def sumset(A,B):
        key=A,B
        if key not in sums:
            S=tuple(new() for _ in range(6));sums[key]=S
            for x in range(6):
                for y in range(6):clauses.append([-A[x],-B[y],S[(x+y)%6]])
        return sums[key]
    for group in groups:
        factors=[support(axis,ps) for axis,ps in group]
        acc=factors[0]
        for B in factors[1:-1]:acc=sumset(acc,B)
        for x in range(6):clauses.append([-acc[x],-factors[-1][-x%6]])
    print('SAT variables',top,'clauses',len(clauses),'support groups',len(groups),flush=True)
    blocked=[];witnesses=[];checks=0;implied=0;rich_profiles={}
    pats=patterns(8,PORT_EDGES)[:args.patterns]
    with Solver(name='g3',bootstrap_with=clauses) as solver:
        for pat in pats:
            if not solver.solve(assumptions=[colours[0][i][q] for i,q in enumerate(pat)]):
                blocked.append(pat);continue
            model={x for x in solver.get_model() if x>0}
            col=[[next(q for q in range(6) if row[q] in model) for row in rows] for rows in colours]
            assert col[0][:8]==list(pat)
            for axis,E in enumerate(axes):
                for i,j in E:assert col[axis][i]!=col[axis][j];checks+=1
            for name,group in zip(group_names,groups):
                acc={0};rich=0;sizes=[]
                for axis,ps in group:
                    D={(col[axis][j]-col[axis][i])%6 for i,j in ps}
                    rich+=any((d+1)%6 in D for d in D);sizes.append(len(D))
                    acc={(a+b)%6 for a in acc for b in D}
                assert 0 not in acc;checks+=1;implied+=math.prod(len(ps) for _,ps in group)
                assert rich<=4
                profile=','.join(map(str, [rich]+sizes))
                rich_profiles.setdefault(name,{})[profile]=rich_profiles.setdefault(name,{}).get(profile,0)+1
            witnesses.append(dict(pattern=pat,colours=col))
            if len(witnesses)%100==0:print('verified',len(witnesses),'blocked',len(blocked),flush=True)
    out=dict(status='FINITE-EXACT',six=args.six,dense_six=args.dense_six,mixed_quads=args.mixed_quads,
             coordinate_sizes=list(map(len,sets)),implicit_vertices=math.prod(map(len,sets)),
             axis_edges=list(map(len,axes)),families=counts,patterns=len(pats),blocked=blocked,
             variables=top,clauses=len(clauses),witness_checks=checks,implied_contact_checks=implied,
             support_groups=[dict(family=name,factors=group) for name,group in zip(group_names,groups)],
             rich_profiles=rich_profiles,
             family_coefficients={name:[(axis,[str(t) for t in co]) for axis,co in terms] for name,terms in families},
             quadrilaterals=quads,witnesses=witnesses,
             scope='Exact extensions of axis and named-family edges only. No complete six-coordinate unit census.')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print('blocked',len(blocked),'checks',checks,'implied contacts',implied,flush=True)


if __name__=='__main__':main()
