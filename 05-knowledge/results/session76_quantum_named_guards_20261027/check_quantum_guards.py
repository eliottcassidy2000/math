#!/usr/bin/env python3
"""Exact core checks for quantum_named_guards.md. Requires SymPy.

Constructs bases and response matrices directly; reads/writes no data files.
No external maximum-MUB result or QCD assertion is used.
"""
from collections import Counter
from itertools import combinations, product
import sympy as s


def require(ok, message):
    if not ok:
        raise ValueError(message)


def simp(x):
    return s.simplify(s.expand_complex(x))


def zero(a):
    return all(simp(x) == 0 for x in a)


def main():
    w = (-1+s.I*s.sqrt(3))/2
    f = s.Matrix(3,3,lambda y,c: simp(w**(y*c)))/s.sqrt(3)
    x = s.Matrix([[1,1],[1,-1]])/s.sqrt(2)
    y = s.Matrix([[1,1],[s.I,-s.I]])/s.sqrt(2)
    bases = [s.eye(6),s.kronecker_product(x,f),s.kronecker_product(y,s.diag(1,w,w)*f)]
    for b in bases:
        require(zero(b.H*b-s.eye(6)), 'basis not unitary')
    for b,c in combinations(bases,2):
        require(all(simp(z*s.conjugate(z)) == s.Rational(1,6) for z in b.H*c), 'not mutually unbiased')
    edges = list(combinations(range(6),2))
    response = s.Matrix([[simp(s.I*(s.conjugate(b[i,c])*b[j,c]-s.conjugate(b[j,c])*b[i,c]))
                          for i,j in edges] for b in bases for c in range(6)])
    def row(terms):
        return [terms.get(e,0) for e in edges]
    stats = s.Matrix([
        row({(0,1):1,(3,4):1}), row({(0,2):1,(3,5):1}), row({(1,2):1,(4,5):1}),
        row({(0,3):1,(1,4):1,(2,5):1}),
        row({(0,4):1,(1,3):1,(0,5):-1,(2,3):-1}),
        row({(0,5):1,(2,3):1,(1,5):-1,(2,4):-1}),
        row({(0,4):1,(1,3):-1,(0,5):-1,(2,3):1,(1,5):1,(2,4):-1})])
    require(response.rank() == stats.rank() == response.col_join(stats).rank() == 7, 'row-space equivalence')
    expanded = []
    for r in response.tolist():
        radical = [s.expand(v).coeff(s.sqrt(3)) for v in r]
        rational = [simp(v-s.sqrt(3)*b) for v,b in zip(r,radical)]
        require(all(v.is_Rational for v in rational+radical), 'coefficient field')
        expanded.extend((rational,radical))
    require(s.Matrix(expanded).rank() == 7, 'rational-expanded rank')

    k0 = s.Matrix(6,6,lambda i,j: 1 if i < j else -1 if i > j else 0)
    kt,ku = k0.copy(),k0.copy()
    kt[1,4],kt[4,1] = -1,1
    ku[0,3],ku[3,0] = -1,1
    eps = s.Rational(1,60)
    rt,ru = s.eye(6)/6+eps*s.I*kt, s.eye(6)/6+eps*s.I*ku
    require(s.Rational(1,6)-5*eps == s.Rational(1,12), 'positivity tariff')
    require(all(zero((b.H*(rt-ru)*b).diagonal()) for b in bases), 'collision failed')
    degrees = lambda k: [sum(k[i,j] == 1 for j in range(6)) for i in range(6)]
    require(5 in degrees(kt) and 5 not in degrees(ku), 'nonisomorphism witness failed')
    probe = s.Matrix([1,0,0,s.I,0,0])/s.sqrt(2)
    require(simp((probe.H*rt*probe)[0]) == s.Rational(3,20), 'T probe')
    require(simp((probe.H*ru*probe)[0]) == s.Rational(11,60), 'U probe')
    require(s.trace((rt-ru)**2) == s.Rational(1,225), 'HS distance')

    def cross(v):
        a,c,b,e,d,f0 = v
        return a+c-b-e,b+e-d-f0,a-c-b+e+d-f0
    counts = Counter(cross(v) for v in product((-1,1),repeat=6))
    require(Counter(counts.values()) == {1:44,2:7,3:2}, 'cross64 histogram')
    supports = {tuple(i for i,a in enumerate(v) if a) for v in product((-1,0,1),repeat=6)
                if any(v) and cross(v) == (0,0,0)}
    require(supports == {(0,1,2,3),(0,1,4,5),(2,3,4,5),tuple(range(6))}, 'cross collision proof')
    guards = {(0,1),(0,2),(0,3),(0,4),(0,5),(1,2),(1,4)}
    unknown = [e for e in edges if e not in guards]
    restricted = stats[:,[edges.index(e) for e in unknown]]
    witness = s.Matrix([2,1,2,1,0,0,0,0])
    require(unknown == [(1,3),(1,5),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)], 'unknown order')
    require(restricted.rank() == 7 and restricted*witness == s.zeros(7,1), 'primitive kernel')
    # Direct full sign-cube control, using the already-proved seven statistics.
    integer_rows = [list(map(int,r)) for r in stats.tolist()]
    hist = Counter(); augmented = set()
    for signs in product((-1,1),repeat=15):
        sig = tuple(sum(a*b for a,b in zip(r,signs)) for r in integer_rows)
        hist[sig] += 1
        augmented.add((sig,tuple(signs[edges.index(e)] for e in sorted(guards))))
    require(len(hist) == 5724 and max(hist.values()) == 72 and Counter(hist.values())[72] == 4, 'full signature count')
    require(len(augmented) == 32768, 'seven guards do not separate')

    # The exact colour-algebra coordinate: realification of -i*diag(2,-2,0).
    d = s.diag(2,-2,0); o = s.zeros(3); eye = s.eye(3)
    j = o.row_join(-eye).col_join(eye.row_join(o))
    kd = o.row_join(d).col_join((-d).row_join(o))
    require(kd == kt-ku and kd*j == j*kd and s.trace(d) == 0, 'Cartan realification')
    c3 = s.Matrix([[1,-2],[2,-1]])
    require(c3*c3 == -3*s.eye(2), 'Eisenstein operator square')
    # A deliberate false equality must fail even under optimized Python.
    try:
        require(s.Rational(1,6) == s.Rational(1,5), 'intentional bad overlap')
    except ValueError:
        pass
    else:
        raise ValueError('failure gate vanished')
    print('PASS: 3 MUBs, exact collision, rank7 statistics, 32768 -> 5724 signatures, seven guards, SU(3)/A2 identities; hostile control rejected')


if __name__ == '__main__':
    main()
