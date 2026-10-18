#!/usr/bin/env python3
"""Reflection-linked coefficient copies, retaining the five-coordinate circuit.

Graph: A+uC+vC+(uv)C+Y C, Y=e^2, with axis edges and four
authenticated contact families. This is a selected unit-distance graph,
not a census of every unit difference in the five-coordinate space.
Every SAT witness is checked directly against all factor constraints.
"""
import argparse
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path

from pysat.solvers import Solver
from golden_forcing_search_20261017 import Geo, recipe
from golden_boundary_20261017c import PORT_EDGES, patterns


def require(p, text):
    if not p:
        raise RuntimeError(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--recipe', default='base', choices=['base','dM','MM'])
    ap.add_argument('--cubic', action='store_true',help='realize the u^2 circuit with 19-denominator coefficient translations')
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--links', type=int, default=2, help='new coefficient points obtained by b/c reflection links')
    ap.add_argument('--patterns', type=int, default=570)
    ap.add_argument('--copy-blind',action='store_true',help='force b/c colour factors to ignore the copy label')
    ap.add_argument('--overlap',type=int,default=0,choices=[0,1,2],help='add coefficient labels 1/b,1/b+1 and optionally 1/c,1/c+1')
    ap.add_argument('--inverse',action='store_true',help='add u^-1 circuit sharing its Y coefficient with u^2')
    args = ap.parse_args()
    geo = Geo(); f = geo.f
    A,C = recipe(geo,args.recipe)
    sets = [A,C,C,C,C]
    d,db = geo.d,geo.bar(geo.d)
    w = f.gen(2)
    h = f.add(w,f.scalar(3))
    u = f.mul(f.sub(h,geo.j),f.inv(f.add(h,geo.j)))
    Z = f.add(d,u)
    Y = f.scale(f.mul(Z,f.inv(geo.bar(Z))),-1)
    cubic_coefficients=None
    reflection_receipts=[]
    if args.cubic:
        m2=f.scale(f.add(f.scalar(20),f.scale(geo.j,11)),Q(-1,19))
        m1=f.scale(f.sub(f.scalar(20),f.scale(geo.j,11)),Q(1,19))
        require(f.power(u,3)==f.add(f.add(f.mul(m2,f.power(u,2)),f.mul(m1,u)),f.one),'certified cubic')
        db2=f.mul(db,db)
        cc=f.mul(f.add(f.sub(f.add(f.mul(db,m2),f.one),f.mul(db2,m1)),f.mul(db2,db)),
                 f.inv(f.sub(f.mul(d,db),f.one)))
        cb=f.add(f.mul(cc,d),f.sub(f.mul(db,m1),db2))
        cubic_coefficients=(db,cb,cc)
        require(f.power(u,2)==f.add(f.add(db,f.mul(cb,u)),f.mul(cc,Y)),'u^2 cross-coordinate circuit')
        oldM=list(geo.M)
        for co in (cb,cc)[:args.overlap]:
            inv=f.inv(co)
            geo.M.extend([inv,f.add(inv,f.one)])
        for i in range(args.links):
            co=(cb,cc)[i%2]
            require(not geo.is_unit(co),'nondegenerate Mobius generator')
            q=geo.M[1+i//2] if 1+i//2<len(oldM) else f.one
            require(geo.is_unit(q),'reflection starts on unit circle')
            s=f.one
            # R_s(q)=(q-s/co)/(1-conj(s/co)*q), inverse R_-s.
            v=f.mul(s,f.inv(co))
            qp=f.mul(f.sub(q,v),f.inv(f.sub(f.one,f.mul(geo.bar(v),q))))
            require(geo.is_unit(qp),'reflected copy label is unit')
            require(geo.is_unit(f.add(s,f.mul(co,f.sub(qp,q)))),'authenticated cross-copy edge')
            require(qp!=q and qp not in geo.M,'new copy')
            geo.M.append(qp)
            reflection_receipts.append(dict(axis=(1,4)[i%2],source_copy=oldM.index(q),new_copy=len(geo.M)-1,
                                            qprime=[str(c) for c in qp]))
        for axis,co in zip((0,1,4),cubic_coefficients):
            # dict order keeps the old ports and zero coefficient at their old indices.
            sets[axis]=list(dict.fromkeys(sets[axis]+[f.add(x,f.mul(co,m)) for x in sets[axis] for m in geo.M]))
        if args.inverse:
            inverse_coefficients=(f.sub(db,m1),f.sub(cb,m2),cc)
            require(f.inv(u)==f.add(f.add(inverse_coefficients[0],f.mul(inverse_coefficients[1],u)),f.mul(cc,Y)),
                    'inverse cubic circuit')
            for axis,co,base in ((0,inverse_coefficients[0],A),(1,inverse_coefficients[1],C)):
                sets[axis]=list(dict.fromkeys(sets[axis]+[f.add(x,f.mul(co,m)) for x in base for m in geo.M]))
    edge_sets=[geo.unit_edges(s) for s in sets]
    print('coordinate sizes',list(map(len,sets)),'axis edges',list(map(len,edge_sets)),flush=True)
    yi = f.inv(Y)
    zero = (f.zero,f.zero)
    one = (f.one,f.zero)
    def add(x,y): return (f.add(x[0],y[0]),f.add(x[1],y[1]))
    def scale(x,a): return (f.mul(x[0],a),f.mul(x[1],a))
    def mul(x,y):
        bb = f.mul(x[1],y[1])
        return (f.add(f.mul(x[0],y[0]),f.mul(Y,bb)),
                f.sub(f.add(f.mul(x[0],y[1]),f.mul(x[1],y[0])),f.mul(Z,bb)))
    def bar(x):
        z = f.mul(geo.bar(x[1]),yi)
        return (f.add(geo.bar(x[0]),f.mul(z,Z)),z)
    dirs = [one,(u,f.zero),(f.zero,f.one),(f.zero,u),(Y,f.zero)]
    families = [
        ('golden',[(0,d),(1,f.one),(2,f.one)]),
        ('reflected',[(1,f.one),(2,f.one),(3,db)]),
        ('Y_contact',[(0,d),(1,f.one),(4,db)]),
        ('v_square_circuit',[(2,f.scale(d,-1)),(3,f.scalar(-1)),(4,f.one)]),
    ]
    if args.cubic:
        families.append(('u_square_circuit',list(zip((0,1,4),cubic_coefficients))))
        if args.inverse:families.append(('u_inverse_circuit',list(zip((0,1,4),inverse_coefficients))))
    for name,terms in families:
        value = zero
        for axis,co in terms:
            value = add(value,scale(dirs[axis],co))
        require(mul(value,bar(value)) == one,'unit family '+name)
    base_edges=geo.unit_edges(C) if args.cubic else edge_sets[1]
    units = sorted({f.sub(C[b],C[a]) for a,b in base_edges} |
                   {f.sub(C[a],C[b]) for a,b in base_edges})
    maps = [{x:i for i,x in enumerate(s)} for s in sets]
    cache = {}
    def pairs(axis,delta):
        key = (axis,delta)
        if key not in cache:
            cache[key] = [(i,maps[axis][f.add(x,delta)]) for i,x in enumerate(sets[axis])
                          if f.add(x,delta) in maps[axis]]
        return cache[key]
    constraints = []
    family_counts = {}
    for name,terms in families:
        seen = set()
        family_units=units
        if args.cubic:
            axis,co=min(terms,key=lambda term: len(sets[term[0]]))
            inverse=f.inv(co)
            scaled=[f.mul(x,inverse) for x in sets[axis]]
            scaled_edges=geo.unit_edges(scaled)
            family_units=sorted({f.sub(scaled[j],scaled[i]) for i,j in scaled_edges} |
                                {f.sub(scaled[i],scaled[j]) for i,j in scaled_edges})
        for q in family_units:
            for choices in itertools.product(*(pairs(axis,f.mul(co,q)) for axis,co in terms)):
                row = tuple((axis,i,j) for (axis,_),(i,j) in zip(terms,choices))
                reverse = tuple((axis,j,i) for axis,i,j in row)
                seen.add(min(row,reverse))
        family_counts[name] = len(seen)
        constraints.extend(sorted(seen))
        print('family',name,len(seen),flush=True)
    top = 0
    def new():
        nonlocal top
        top += 1
        return top
    V = [[[new() for _ in range(6)] for _ in s] for s in sets]
    clauses = []
    for axis,s in enumerate(sets):
        for row in V[axis]:
            clauses.append(row)
            clauses.extend([[-row[a],-row[b]] for a in range(6) for b in range(a)])
        for i,j in edge_sets[axis]:
            clauses.extend([[-V[axis][i][c],-V[axis][j][c]] for c in range(6)])
        if axis:
            require(s[0] == f.zero,'coefficient zero port')
            clauses.append([V[axis][0][0]])
    if args.copy_blind:
        require(args.cubic,'copy-blind test needs cubic recipe')
        for axis,co in ((1,cb),(4,cc)):
            for m in geo.M:
                for x in C:
                    i=maps[axis][x]; j=maps[axis][f.add(x,f.mul(co,m))]
                    for c in range(6):clauses.append([-V[axis][i][c],V[axis][j][c]])
    diffs = {}
    def diff(axis,i,j):
        key = (axis,i,j)
        if key not in diffs:
            D = [new() for _ in range(6)]
            for x in range(6):
                for y in range(6):
                    clauses.append([-V[axis][j][x],-V[axis][i][y],D[(x-y)%6]])
            diffs[key] = D
        return diffs[key]
    for row in constraints:
        a,b,c = [diff(*entry) for entry in row]
        for x in range(6):
            for y in range(6):
                clauses.append([-a[x],-b[y],-c[(-x-y)%6]])
    pats = patterns(8,PORT_EDGES)[:args.patterns]
    blocked,witnesses = [],[]
    checks = 0
    with Solver(name='g3',bootstrap_with=clauses) as sol:
        for pat in pats:
            if not sol.solve(assumptions=[V[0][i][c] for i,c in enumerate(pat)]):
                blocked.append(pat)
                continue
            model = set(x for x in sol.get_model() if x>0)
            col = [[next(c for c in range(6) if row[c] in model) for row in vs] for vs in V]
            require(col[0][:8] == list(pat),'literal ports')
            for axis,edges in enumerate(edge_sets):
                for i,j in edges:
                    require(col[axis][i] != col[axis][j],'axis edge')
                    checks += 1
            for row in constraints:
                require(sum(col[axis][j]-col[axis][i] for axis,i,j in row)%6 != 0,'contact edge')
                checks += 1
            witnesses.append(dict(pattern=pat,colours=col))
            if len(witnesses)%100==0:print('verified patterns',len(witnesses),'blocked',len(blocked),flush=True)
    out = dict(status='FINITE-EXACT',recipe=args.recipe,coordinate_sizes=list(map(len,sets)),
               cubic=args.cubic,cubic_coefficients=None if cubic_coefficients is None else
                    [[str(x) for x in co] for co in cubic_coefficients],
               reflection_receipts=reflection_receipts,
               copy_blind=args.copy_blind,
               overlap=args.overlap,
               inverse=args.inverse,
               vertex_count=__import__('math').prod(map(len,sets)),
               selected_families=family_counts,axis_edge_counts=list(map(len,edge_sets)),
               blocked=blocked,patterns=len(pats),witness_checks=checks,witnesses=witnesses,
               scope='axis and named contact families only; not a complete unit census')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,separators=(',',':'))+'\n')
    summary={k:v for k,v in out.items() if k not in ('witnesses','reflection_receipts','cubic_coefficients','blocked')}
    summary['blocked_count']=len(blocked)
    print(json.dumps(summary,indent=2),flush=True)


if __name__ == '__main__':
    main()
