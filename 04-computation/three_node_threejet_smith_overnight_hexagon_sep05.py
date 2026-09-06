#!/usr/bin/env python3
"""Finite symbolic certificate for every dyadic three-node three-jet Smith form.

Exact universe: ALL886 minors of ranks1..4 of the6x6 cleared residual,
not sampled minors. Every term is retained as a polynomial in the unit-
scaled close-node parameter a, with separate global dilation weight.
Coordinatewise dominance after d=1+d' certifies the unbounded valuation
lower envelope. Eight explicit minors attain every surviving affine cost.
The largest factor is inherited from audited THM-4443, not guessed here.
"""

from collections import defaultdict
import hashlib
import itertools
import math

from sympy import Matrix, ZZ, binomial
from sympy.matrices.normalforms import smith_normal_form


ROWS = [(0,r) for r in range(3)] + [(1,r) for r in range(3)]
GATES = 0
DIGEST = hashlib.sha256()


def check(condition, data):
    global GATES
    if not condition:
        raise RuntimeError(repr(data))
    GATES += 1
    DIGEST.update(repr(data).encode())
    DIGEST.update(b'\n')


def valuation(a):
    if not a:
        return math.inf
    a = abs(int(a))
    out = 0
    while a % 2 == 0:
        a //= 2
        out += 1
    return out


def symbolic_minor(rr, cc):
    coefficients = defaultdict(int)
    size = len(rr)
    for perm in itertools.permutations(range(size)):
        inversions = sum(perm[i] > perm[j] for i in range(size)
                         for j in range(i+1,size))
        term = (-1)**inversions
        power = 0
        for i, j in enumerate(perm):
            node, r = ROWS[rr[i]]
            degree = cc[j]
            term *= math.comb(degree, r)
            power += node*(degree-r)
        coefficients[power] += term
    weight = sum(cc)-sum(ROWS[i][1] for i in rr)
    return weight, tuple(sorted((j,a) for j,a in coefficients.items() if a))


def polynomial_product(p,q):
    out = defaultdict(int)
    for i,a in p.items():
        for j,b in q.items():
            out[i+j] += a*b
    return {i:a for i,a in out.items() if a}


def primary_certificate():
    expected = {
        1: {(1,0,0)},
        2: {(5,0,0),(4,0,1),(3,1,2)},
        3: {(9,0,0),(7,1,2)},
        4: {(13,1,1),(12,4,5)},
    }
    counts = []
    witnesses = {}
    for size in range(1,5):
        costs = {}
        count = 0
        for rr in itertools.combinations(range(6),size):
            for cc in itertools.combinations(range(3,9),size):
                weight, polynomial = symbolic_minor(rr,cc)
                count += 1
                check(bool(polynomial), ('nonzero-minor',size,rr,cc))
                for exponent, coefficient in polynomial:
                    cost = weight,exponent,valuation(coefficient)+exponent
                    costs.setdefault(cost,(rr,cc,polynomial))
                for a in (2,4,6):
                    literal = Matrix([[binomial(j,ROWS[i][1])
                                       *a**(ROWS[i][0]*(j-ROWS[i][1]))
                                       for j in cc] for i in rr]).det()
                    value = sum(coefficient*a**exponent
                                for exponent,coefficient in polynomial)
                    check(literal == value, ('minor-evaluation',rr,cc,a,value))
                DIGEST.update(repr((rr,cc,weight,polynomial)).encode())
        pareto = {c for c in costs if not any(
            d != c and all(x <= y for x,y in zip(d,c)) for d in costs)}
        check(pareto == expected[size], ('complete-cost-frontier',size,sorted(pareto)))
        for cost in costs:
            check(any(all(x <= y for x,y in zip(p,cost)) for p in pareto),
                  ('unbounded-dominance',size,cost))
        for cost in sorted(pareto):
            witnesses[size,cost] = costs[cost]
        counts.append(count)
        print('RANK',size,'ALL MINORS',count,'PARETO (e,d-1,constant)',sorted(pareto))

    # Independently expand the displayed attained factorizations.
    a = {1:1}
    am1 = {1:1,0:-1}
    explicit = [
        ((2,),(3,),{0:3}),
        ((0,2),(3,4),{0:3}),
        ((1,2),(3,4),{0:6}),
        ((2,5),(3,4),polynomial_product({1:18},am1)),
        ((0,1,2),(3,4,5),{0:1}),
        ((1,2,5),(3,4,5),polynomial_product(
            polynomial_product({1:30},am1),{1:2,0:-1})),
        ((0,1,2,5),(3,4,5,6),polynomial_product(
            polynomial_product({1:3},am1),{2:5,1:-5,0:1})),
        ((1,2,4,5),(3,4,5,6),{4:90,5:-360,6:540,7:-360,8:90}),
    ]
    for rr,cc,p in explicit:
        weight,poly = symbolic_minor(rr,cc)
        check(dict(poly) == p, ('attaining-factorization',rr,cc,p))
        print('ATTAIN',rr,cc,'weight',weight,'polynomial',poly)
    check(sum(counts)==886,('complete-universe',counts))
    return sum(counts)


def predicted(e,d,gamma):
    largest = 8*e+5*d-1 if d>=2 else max(7*e+4,8*e+gamma)
    determinantal = [0,0,0,0,e,min(5*e,4*e+1,3*e+d+1),
                    min(9*e,7*e+d+1),min(13*e+d,12*e+4*d+1),
                    27*e+9*d-largest,27*e+9*d]
    return tuple(determinantal[j+1]-determinantal[j] for j in range(9))


def literal_controls():
    count = 0
    for e,d,t in itertools.product(range(9),range(1,6),(1,3,5,9)):
        nodes = (0,2**e*t,2**(e+d))
        matrix = Matrix([[binomial(j,r)*x**(j-r) if j>=r else 0
                          for j in range(9)] for x in nodes for r in range(3)])
        snf = smith_normal_form(matrix,domain=ZZ)
        actual = tuple(valuation(snf[i,i]) for i in range(9))
        gamma = {1:1,3:3,5:2,9:0}[t]
        expected = predicted(e,d,gamma)
        check(actual == expected, ('full-smith',e,d,t,actual,expected))
        count += 1
    print('LITERAL FULL9x9 CONTROLS',count,'e0..8,d1..5,t=1,3,5,9')
    print('SHALLOW-EXTRAPOLATION HOSTILE e5,d1: D7=65, not13e+d=66')
    return count


def main():
    minors = primary_certificate()
    matrices = literal_controls()
    print('ALL SYMBOLIC MINORS',minors,'LITERAL MATRICES',matrices)
    print('GATES',GATES,'SEMANTIC',DIGEST.hexdigest())
    print('PASS finite symbolic lower-envelope certificate and unbounded attained factors')


if __name__ == '__main__':
    main()
