"""Exact unit-rhombus states and marked two-lattice six-colour certificates.

Local orientation choices, unit-edge presence, and palette equalities are
different types. This supplies a restricted positive colouring family.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations, product

import finite_algebraic_witness_20261007 as a
import boundary_palette_compression_20261011 as p

CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def plus(x,y):
    return tuple(a.add(u,v) for u,v in zip(x,y))


def times(k,x):
    return tuple(a.scale(u,k) for u in x)


def dot(x,y,field):
    return a.add(field.times(x[0],y[0]),field.times(x[1],y[1]))


def point_checked(x):
    need(type(x) is tuple and len(x)==2, 'plane point')
    need(all(type(t) is tuple and t and all(type(v) in (int,F) for v in t)
             for t in x), 'exact rational polynomial coordinates')


@dataclass(frozen=True)
class Carrier:
    origin: tuple
    u: tuple
    v: tuple


@dataclass(frozen=True)
class Cover:
    carriers: tuple
    addresses: tuple


def lift_cover(points, field, cover):
    need(type(points) is tuple and type(field) is a.AlgebraicReal, 'exact source points and field')
    for x in points:
        point_checked(x)
    need(type(cover) is Cover and type(cover.carriers) is tuple
         and 1 <= len(cover.carriers) <= 2, 'one or two marked carriers')
    for carrier in cover.carriers:
        need(type(carrier) is Carrier, 'typed carrier')
        for x in (carrier.origin,carrier.u,carrier.v):
            point_checked(x)
        need(field.equal(dot(carrier.u,carrier.u,field),a.ONE)
             and field.equal(dot(carrier.v,carrier.v,field),a.ONE)
             and field.equal(dot(carrier.u,carrier.v,field),a.poly((F(1,2),))),
             'unit triangular Gram matrix')
    need(type(cover.addresses) is tuple and len(cover.addresses)==len(points), 'every source point addressed')
    colours=[]
    for point,address in zip(points,cover.addresses):
        need(type(address) is tuple and len(address)==3
             and all(type(v) is int for v in address), 'exact layer and integer coordinates')
        layer,x,y=address
        need(0 <= layer < len(cover.carriers), 'declared layer')
        c=cover.carriers[layer]
        decoded=plus(c.origin,plus(times(x,c.u),times(y,c.v)))
        need(a.same_point(point,decoded,field), 'address decodes actual source point')
        colours.append(3*layer+(x-y)%3)
    # Reauthenticate all actual unit edges, including cross-carrier edges.
    edges=a.unit_graph(points,field)
    need(all(colours[u]!=colours[v] for u,v in edges), 'full actual colouring lift')
    return edges,tuple(colours)


def rhombi():
    field=a.AlgebraicReal((-3,0,1),1,2)
    zero=(a.ZERO,a.ZERO); one=(a.ONE,a.ZERO)
    high=a.poly((0,F(1,2)))
    return field,(
        (zero,one,(a.ONE,a.ONE),(a.ZERO,a.ONE)),
        (zero,one,(a.poly((F(3,2),)),high),(a.poly((F(1,2),)),high)),
        (zero,one,(a.poly((F(1,2),)),high),(a.poly((F(-1,2),)),high)),
    )


def marked_spindle():
    field=a.AlgebraicReal((64,0,-28,0,1),5,6)
    s=a.scale(a.poly((0,-20,0,1)),F(1,16))
    t=a.scale(a.poly((0,36,0,-1)),F(1,16))
    st=field.times(s,t)
    zero=(a.ZERO,a.ZERO)
    u=(a.scale(s,F(1,2)),a.poly((F(1,2),)))
    v=(a.scale(s,F(1,2)),a.poly((F(-1,2),)))
    ru=(a.scale(a.sub(a.scale(s,5),t),F(1,12)),a.scale(a.add(st,a.poly((5,))),F(1,12)))
    rv=(a.scale(a.add(a.scale(s,5),t),F(1,12)),a.scale(a.sub(st,a.poly((5,))),F(1,12)))
    points=(zero,plus(u,v),u,v,plus(ru,rv),ru,rv)
    cover=Cover((Carrier(zero,u,v),Carrier(zero,ru,rv)),
                ((0,0,0),(0,1,1),(0,1,0),(0,0,1),(1,1,1),(1,1,0),(1,0,1)))
    return points,field,cover


def labelled_colour_count(n,edges,k):
    return sum(all(c[u]!=c[v] for u,v in edges) for c in product(range(k),repeat=n))


def rejects(fn):
    try:
        fn()
    except ValueError:
        check(True,'hostile rejected')
    else:
        raise RuntimeError('invalid certificate accepted')


def main():
    field,examples=rhombi()
    ring=frozenset(((0,1),(1,2),(2,3),(0,3)))
    states=set()
    for points in examples:
        edges=a.unit_graph(points,field)
        check(ring <= edges,'four unit sides')
        d0=a.distance_squared(points[0],points[2],field)
        d1=a.distance_squared(points[1],points[3],field)
        check(field.equal(a.add(d0,d1),a.poly((4,))),'coupled diagonal squares')
        states.add(((0,2) in edges,(1,3) in edges))
    check(states=={(False,False),(False,True),(True,False)}, 'three geometric edge-presence states')
    rejects(lambda:a.unit_graph((examples[0][0],examples[0][1],examples[0][0],examples[0][3]),field))
    close=((F(0),F(0)),(F(1,3),F(0)),(F(0),F(1,3)))
    check(all(0 < (u[0]-v[0])**2+(u[1]-v[1])**2 < 1
              for u,v in combinations(close,2)),
          'three close points cannot lie in two unit triangular carriers')

    for k in range(2,8):
        signatures=set()
        for c in product(range(k),repeat=4):
            if all(c[u]!=c[v] for u,v in ring):
                signatures.add((c[0]!=c[2],c[1]!=c[3]))
        expected={(False,False)}
        if k>=3: expected |= {(False,True),(True,False)}
        if k>=4: expected.add((True,True))
        check(signatures==expected,'palette-dependent diagonal equality interface')
    # Every residual same-label subgraph of C4 is bipartite, unlike an odd cycle.
    for labels in product(range(3),repeat=4):
        residual=frozenset((u,v) for u,v in ring if labels[u]==labels[v])
        check(labelled_colour_count(4,residual,2)>0,'C4 imposes no ternary-channel obstruction')

    # Finite strip of edge-sharing triangles in an integer triangular lattice.
    for n in range(3,8):
        lattice=tuple((i//2,i%2) for i in range(n))
        edges=frozenset((i,j) for i,j in combinations(range(n),2)
                        if (lattice[i][0]-lattice[j][0])**2
                        +(lattice[i][0]-lattice[j][0])*(lattice[i][1]-lattice[j][1])
                        +(lattice[i][1]-lattice[j][1])**2 == 1)
        check(labelled_colour_count(n,edges,3)==6,'triangle-connected patch has one S3 colour orbit')
        check(all((x-y)%3!=(u-v)%3 for i,j in edges for x,y in (lattice[i],) for u,v in (lattice[j],)),
              'integer lattice phase colours every actual unit edge')
    diamond=ring|{(1,3)}
    for k in (3,4,5,6):
        count=labelled_colour_count(4,diamond,k)
        check(count==k*(k-1)*(k-2)**2,'diamond palette freedom')
        pairs={(c[0]==c[2]) for c in product(range(k),repeat=4)
               if all(c[u]!=c[v] for u,v in diamond)}
        check(pairs==({True} if k==3 else {True,False}), 'three-colour equality is not a six-colour rule')

    points,field,cover=marked_spindle()
    edges,colours=lift_cover(points,field,cover)
    check(edges==p.spindle().edges,'two-carrier coordinates reproduce exact unit spindle')
    check(labelled_colour_count(7,edges,3)==0,'single-lattice three-colour propagation cannot cover the spindle')
    check(all(colours[u]!=colours[v] for u,v in edges),'two-lattice certificate gives actual six-colouring')
    c=cover.carriers[0]
    for x,y in product(range(-5,6),repeat=2):
        pt=plus(c.origin,plus(times(x,c.u),times(y,c.v)))
        es,cs=lift_cover((pt,),field,Cover((c,),((0,x,y),)))
        check(cs==((x-y)%3,), 'signed integer coordinates decode')
    rejects(lambda:lift_cover(points,field,replace(cover,addresses=cover.addresses[:-1])))
    rejects(lambda:lift_cover(points,field,replace(cover,addresses=((False,0,0),)+cover.addresses[1:])))
    rejects(lambda:lift_cover(points,field,replace(cover,addresses=((0,1,0),)+cover.addresses[1:])))
    rejects(lambda:lift_cover(points,field,replace(cover,carriers=(replace(c,u=times(2,c.u)),cover.carriers[1]))))
    rejects(lambda:lift_cover(points,field,replace(cover,carriers=(c,)*3)))

    print('Unit C4: diagonal squared lengths sum4; edge-presence states00,01,10, never11.')
    print('Proper-colour C4 diagonal-difference states: k2 one, k3 three, k>=4 four.')
    print('Every ternary label assignment on C4 lifts to six colours; no odd-cycle obstruction.')
    print('Edge-sharing triangular strips orders3..7 have exactly6 labelled3-colourings.')
    print('Unit diamond: k(k-1)(k-2)^2 colourings; opposite equality forced only at k3.')
    print('Marked two-triangular-lattice cover exports six-colouring; exact spindle control:',colours)
    print('No universal two-carrier cover or plane six/seven decision claimed.')
    print('PASS:',CHECKS,'exact checks')


if __name__=='__main__':
    main()
