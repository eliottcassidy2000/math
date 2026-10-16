#!/usr/bin/env python3
"""Exact marked elliptic/contact controls; stdout is a deterministic receipt.

Dependencies: SymPy, PySAT, and the two explicitly imported repository files.
The graph universe, all coordinates, all edges, and colour witnesses are saved.
No geometric test uses floating point or assumes Collatz convergence.
"""
import os
os.environ['SYMPY_GROUND_TYPES'] = 'python'
import json
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import sympy as s
from sympy.polys.matrices import DomainMatrix
from pysat.solvers import Solver
import charge_ports_eleven_20261011 as G
import creation_elliptic_20260925 as E

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def sub(a, b):
    return G.add(a, G.sc(b, -1))


def power(a, n):
    if n < 0:
        return power(G.conj(a), -n)  # only used on norm-one elements
    out = G.I
    while n:
        if n & 1:
            out = G.mul(out, a)
        a = G.mul(a, a)
        n //= 2
    return out


def unit(a):
    return G.mul(a, G.conj(a)) == G.I


def colouring(n, edges, q):
    clauses = [[q*i+c+1 for c in range(q)] for i in range(n)]
    clauses += [[-q*i-c-1, -q*j-c-1] for i, j in edges for c in range(q)]
    with Solver(name='g3', bootstrap_with=clauses) as solver:
        if not solver.solve():
            return None
        model = set(solver.get_model())
    out = [next(c for c in range(q) if q*i+c+1 in model) for i in range(n)]
    check(all(out[i] != out[j] for i, j in edges), 'every actual edge coloured')
    return out


def graph_receipt(points):
    edges, survivors, den = G.graph(points)
    colours = colouring(len(points), edges, 3)
    return {'vertices': len(points), 'edges': edges,
            'modular_survivors': survivors, 'denominator': str(den),
            'three_colouring': colours,
            'coordinates': [[str(v) for v in p] for p in points]}


def contact_tree():
    a = G.sc(G.add(G.I, G.b), F(1, 2))
    gamma = G.sc(G.add(G.sc(G.I, 5), G.b), F(1, 6))
    parameters = [(F(1, 5), F(8, 5)), (F(91, 135), F(10552, 3645))]
    families, traces = [], []
    for r, y in parameters:
        radius = (1+55*r*r)/4
        delta = G.sc(G.add(G.sc(G.I, 55*r*r-1), G.sc(G.w, 2*r)),
                     1/(1+55*r*r))
        check(unit(delta), 'secondary rotation unit')
        for sign in (1, -1):
            b = G.sc(G.mul(G.add(G.sc(G.I, radius+2),
                                G.sc(G.b, sign*y/4)), a), F(1, 6))
            check(unit(sub(a, b)), 'elliptic contact with marked spine')
            check(G.mul(b, G.conj(b)) == G.sc(G.I, radius), 'contact radius')
            check(unit(sub(b, G.mul(delta, b))), 'secondary circle chord')
            families.append((r, sign, b, delta))
        u, v = (radius+2)/2, y/8
        trace = 2*(u*u-11*v*v)/(u*u+11*v*v)
        check(trace.denominator % 5 == 0, 'other contact has a 5-adic pole')
        traces.append(str(trace))
    vertices, labels, expected = [], [], []
    index = {}
    def put(label, point):
        check(point not in vertices, 'finite tree has no coordinate collisions')
        index[label] = len(vertices)
        labels.append(label)
        vertices.append(point)
    for i in range(-2, 3):
        put(('spine', i), G.mul(power(gamma, i), a))
        for f, (_, _, b, delta) in enumerate(families):
            for j in range(-2, 3):
                put(('branch', f, i, j), G.mul(power(gamma, i),
                                            G.mul(power(delta, j), b)))
    for label, idx in index.items():
        if label[0] == 'spine':
            i = label[1]
            if i < 2:
                expected.append(tuple(sorted((idx, index['spine', i+1]))))
            for f in range(4):
                expected.append(tuple(sorted((idx, index['branch', f, i, 0]))))
        else:
            _, f, i, j = label
            if j < 2:
                expected.append(tuple(sorted((idx, index['branch', f, i, j+1]))))
    receipt = graph_receipt(vertices)
    check(sorted(receipt['edges']) == sorted(expected), 'all tree edges and no others')
    two = colouring(len(vertices), receipt['edges'], 2)
    check(two is not None and len(expected) == len(vertices)-1, 'finite tree control')
    receipt.update({'labels': labels, 'two_colouring': two,
                    'universe': 'i,j=-2..2; both signs at r=1/5 and 91/135',
                    'other_contact_traces': traces})
    # A smaller source set for the five-conjugate linkage control.
    seed = list(G.ports)
    for i in range(4):
        seed.append(G.sc(gamma, i))
    seed += [sub(G.sc(gamma, 3), G.I), a]
    for i in range(-1, 2):
        seed.append(G.mul(power(gamma, i), a))
        for _, _, b, delta in families:
            for j in range(-1, 2):
                seed.append(G.mul(power(gamma, i), G.mul(power(delta, j), b)))
    seed = list(dict.fromkeys(seed))
    check(seed[:8] == G.ports, 'original eight ports fixed in order')
    return receipt, a, seed


def inverse_a_minus_root(a, root):
    # (a-root)^-1 = sum a^(10-k) root^k / (a^11-1), root^11=1.
    ap = [G.I]
    rp = [G.I]
    for _ in range(11):
        ap.append(G.mul(ap[-1], a))
        rp.append(G.mul(rp[-1], root))
    den = sub(ap[11], G.I)
    norm = G.mul(den, G.conj(den))
    check(all(t == 0 for t in norm[1:]) and norm[0] != 0,
          'quadratic denominator has nonzero rational norm')
    invden = G.sc(G.conj(den), F(1)/norm[0])
    total = G.O
    for k in range(11):
        total = G.add(total, G.mul(ap[10-k], rp[k]))
    inverse = G.mul(total, invden)
    check(G.mul(sub(a, root), inverse) == G.I, 'exact geometric-series inverse')
    return inverse


def linkages(anchors, seed, direct_all=False):
    pieces = [seed]
    rotations, specifications = [], []
    for anchor_id, a in enumerate(anchors):
        for h in (1, 3, 4, 5, 9):
            root = G.powers[5*h % 55]
            u = G.mul(sub(a, root), inverse_a_minus_root(G.conj(a), G.conj(root)))
            check(unit(u), 'linkage rotation norm')
            piece = [G.add(root, G.mul(u, G.conj(q))) for q in seed]
            check(not set(seed) & set(piece), 'linkage source and target disjoint')
            if direct_all or h == 1:
                # Direct norms independently of the modular edge filter.
                cross = [(i,j) for i,p in enumerate(seed) for j,q in enumerate(piece)
                         if unit(sub(p,q))]
                want = [(seed.index(G.O), seed.index(G.O)), (seed.index(a), seed.index(a))]
                check(sorted(cross) == sorted(want), 'only the two marked seam edges')
            pieces.append(piece)
            rotations.append([str(v) for v in u])
            specifications.append([anchor_id,h])
    points = list(dict.fromkeys(p for piece in pieces for p in piece))
    receipt = graph_receipt(points)
    lookup = {p:i for i,p in enumerate(points)}
    memberships = [{k for k,piece in enumerate(pieces) if p in piece} for p in points]
    new_cross = [(i,j) for i,j in receipt['edges']
                 if memberships[i].isdisjoint(memberships[j])
                 and 0 not in memberships[i] and 0 not in memberships[j]]
    check(points[:8] == G.ports, 'all marked ports remain literal coordinates')
    receipt.update({'seed_vertices': len(seed), 'ports': list(range(8)),
                    'conjugate_exponents': [1,3,4,5,9], 'rotations': rotations,
                    'piece_specifications':specifications,
                    'anchor_indices_in_source':[seed.index(a) for a in anchors],
                    'piece_indices': [[lookup[p] for p in piece] for piece in pieces],
                    'cross_edges_between_new_pieces': new_cross,
                    'source_copy_seams_counted_with_piece_multiplicity': 10*len(anchors)})
    return receipt


def specialization(anchors, seed):
    def build(root, foreign):
        pieces=[seed]
        for a in anchors:
            if foreign:
                inverse=inverse_a_minus_root(G.conj(a),G.conj(root))
            else:
                den=sub(G.conj(a),G.conj(root))
                norm=G.mul(den,G.conj(den))
                check(all(v==0 for v in norm[1:]) and norm[0]!=0,
                      'specialized denominator has rational nonzero norm')
                inverse=G.sc(G.conj(den),F(1)/norm[0])
            u=G.mul(sub(a,root),inverse)
            check(unit(u),'specialized linkage unit')
            pieces.append([G.add(root,G.mul(u,G.conj(q))) for q in seed])
        points=list(dict.fromkeys(p for piece in pieces for p in piece))
        lookup={p:i for i,p in enumerate(points)}
        address=[lookup[p] for piece in pieces for p in piece]
        edges,_,_=G.graph(points)
        return points,address,edges
    points,address,edges=build(G.g,True)
    trials=[]
    for t in (1,2,3,4,5):
        # (1+t*j)/(1-t*j) = (1-11*t^2+2*t*j)/(1+11*t^2).
        root=G.sc(G.add(G.sc(G.I,1-11*t*t),G.sc(G.b,2*t)),F(1,1+11*t*t))
        pts,addr,edg=build(root,False)
        matched=address==addr and edges==edg
        trials.append({'t':t,'vertices':len(pts),'edges':len(edg),'same_labelled_graph':matched})
        if matched:
            check(points[:8]==pts[:8]==G.ports,'specialization fixes every port')
            return {'foreign_vertices':len(points),'foreign_edges':len(edges),
                    'trials':trials,'selected_t':t,'address_map':address,
                    'edges':edges,'specialized_coordinates':[[str(v) for v in p] for p in pts]}
    raise RuntimeError('finite specialization bank exhausted; theorem is not refuted')


def symbolic_controls():
    z, j, aa, bb, cc, dd = s.symbols('z j aa bb cc dd')
    f = z**5+(1-j)*z**4/2-z**3+z**2-(1+j)*z/2-1
    barrev = s.expand(z**5*f.subs({z:1/z,j:-j}, simultaneous=True))
    check(s.expand(barrev+f) == 0, 'quintic conjugate reciprocal relation')
    # Verify the quintic at the actual root and selected Gauss square root.
    value = G.O
    for k in range(6):
        coeff = s.expand(f).coeff(z,k)
        c = G.add(G.sc(G.I,F(coeff.subs(j,0))),G.sc(G.b,F(coeff.coeff(j))))
        value = G.add(value,G.mul(c,G.powers[5*k % 55]))
    check(value == G.O, 'Paley quintic uses the actual embedding')
    # Real part of a product in the biquadratic field keeps a mixed coordinate.
    k = s.Symbol('k')
    h = (aa+bb*j)*(cc+dd*k)
    trace = s.expand(h+h.subs({j:-j,k:-k}, simultaneous=True))
    check(trace == 2*aa*cc+2*bb*dd*j*k, 'two-phase mixed term')
    # Positive control beyond width two: a genuine nonmonomial norm-one unit.
    u = G.unit(2)
    check(unit(u), 'width-four Hilbert90 positive control')
    # Formula for coefficient four in the quintic basis, t=2.
    t = s.Symbol('t')
    cs = s.symbols('c:5')
    v = sum(cs[i]*z**i for i in range(5))
    # Solving at one exact parameter avoids symbolic rational-expression swell.
    equations = s.Poly(s.rem((2*z-1)*v-z*(2-z),f,z),z).all_coeffs()
    sol = s.solve(equations,cs)
    coeff4 = s.cancel(sol[cs[4]])
    check(coeff4 != 0, 'first tested genuinely mixed unit has width four')
    # All real quadratic Hilbert90 seeds: exact Cramer numerator, no scan.
    a,b=s.symbols('a b')
    rem=lambda expr:s.Poly(s.rem(s.rem(expr,f,z),j*j+11,j),z)
    den=a*z*z+b*z+1
    matrix=s.Matrix(5,5,lambda i,k:rem(den*z**k).coeff_monomial(z**i))
    matrix[:,4]=s.Matrix([rem(z*z*(a+b*z+z*z)).coeff_monomial(z**i) for i in range(5)])
    dm=DomainMatrix.from_Matrix(matrix)
    numerator=s.rem(dm.domain.to_sympy(dm.det()),j*j+11,j).expand()
    real=-(a-1)*(a+b+1)*(2*a*a*b-a*a-2*a*b*b+3*a*b+2*a-2*b+2)/2
    imag=a*(a-1)*(a+b)*(a-b+1)/2
    check(s.expand(numerator-real-j*imag)==0,'real quadratic seed Cramer numerator')
    cubic=2*a**3+2*a*a-2*a-1
    check(all(cubic.subs(a,r)!=0 for r in [1,-1,s.Rational(1,2),s.Rational(-1,2)]),
          'rational root test excludes cubic roots')
    # The norm-one anchor is a genuine excluded case of the seam theorem.
    root=G.g
    exceptional_u=G.sc(root,-1)
    check(G.mul(exceptional_u,sub(G.I,G.conj(root)))==sub(G.I,root),
          'unit-anchor collapse by cross multiplication')
    check(G.add(root,G.mul(exceptional_u,G.I))==G.O,'unit-anchor collision hostile')
    return {'paley_quintic':str(f),'mixed_trace':str(trace),
            'hilbert90_parameter':2, 'quintic_coefficient_four':str(coeff4),
            'real_quadratic_seed_cramer_real':str(real),
            'real_quadratic_seed_cramer_imaginary':str(imag),
            'irreducible_cubic':str(cubic), 'unit_anchor_collision':True,
            'hilbert90_unit_coordinates':[str(v) for v in u]}


def fruit_and_clocks():
    p = E.mul(9,E.G)
    triple = E.fruit(p)
    check([len(str(v)) for v in triple] == [81,80,79], 'correct positive triple lengths')
    check(sum(F(a,b+c) for a,b,c in [triple, (triple[1],triple[2],triple[0]),
                                    (triple[2],triple[0],triple[1])]) == 4,
          'fruit equality exactly four')
    # The same elliptic point has different creation and native Collatz successors.
    decoded, r, t = E.decode_step(p)
    check(decoded == E.mul(4,E.G), 'creation successor')
    lifts=[]
    for n in range(1,33):
        m = 1 if n==1 else (n//2 if n%2==0 else (3*n+1)//2)
        point=E.mul(n,E.G)
        if n==1:
            dest=point
        else:
            doubled=point if n%2==0 else E.add(E.mul(3,point),E.G)
            candidates=E.halves(doubled)
            dest=next(q for q in candidates if E.alpha(q) in (1,-1))
        check(dest==E.mul(m,E.G), 'native clock-one elliptic lift')
        lifts.append([n,m,3*(n%3==2)-1,1 if n!=1 else 0])
    # Clock variable distinguishes histories that parity evaluation loses.
    # Difference of two legitimate singleton source rows: w(4)z^2-w(1).
    clock_polynomial=[1,0,-1]
    check(sum(clock_polynomial)==0 and sum(c*(-1)**i for i,c in enumerate(clock_polynomial))==0,
          'both parity specializations lose the different clocks')
    return {'primitive_positive_triple':[str(v) for v in triple],
            'digits':[81,80,79], 'creation_step':[9,4], 'native_step':[9,14],
            'native_lifts':lifts, 'clock_loss_hostile_coefficients':clock_polynomial}


def main():
    tree,a,seed=contact_tree()
    b=G.sc(G.mul(G.add(G.sc(G.I,F(14,5)),G.sc(G.b,F(2,5))),a),F(1,6))
    out={'status':'PROVED lemmas in companion note; finite controls only here',
         'contact_tree':tree,'five_conjugate_linkages':linkages([a],seed,direct_all=True),
         'two_anchor_linkages':linkages([a,b],seed),
         'common_parameter_specialization':specialization([a,b],seed),
         'symbolic':symbolic_controls(),'fruit_and_clocks':fruit_and_clocks()}
    out['always_active_checks']=CHECKS
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
