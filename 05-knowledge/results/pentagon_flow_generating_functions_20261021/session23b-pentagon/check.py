#!/usr/bin/env python3
"""Exact, standalone checks for the R5 anti-fixed circulation series.

Uses only Python's standard library. No coefficient fitting or external data.
Run with python check.py and python -O check.py; outputs must agree.
"""
from fractions import Fraction
from math import prod
import json

N = 300
EDGE_N = 48
P = {0: 1, 4: 1, 5: 1, 6: 1, 7: 1, 8: 1, 12: 1}
WEIGHTS = (3, 5, 6, 8)


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def mul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            out[i+j] = out.get(i+j, 0) + x*y
    return {i: x for i, x in out.items() if x}


def add(*polys):
    out = {}
    for a in polys:
        for i, x in a.items():
            out[i] = out.get(i, 0) + x
    return {i: x for i, x in out.items() if x}


def shift(a, n):
    return {i+n: x for i, x in a.items()}


def factor(d):
    return {0: 1, d: -1}


def remainder(a, b):
    a = {i: Fraction(x) for i, x in a.items() if x}
    db, lc = max(b), b[max(b)]
    while a and max(a) >= db:
        da = max(a)
        c = a[da] / lc
        for i, x in b.items():
            a[i+da-db] = a.get(i+da-db, 0) - c*x
            if not a[i+da-db]:
                del a[i+da-db]
    return a


def series(weights, numerator, cutoff):
    out = [0] * (cutoff+1)
    out[0] = 1
    for d in weights:
        for n in range(d, cutoff+1):
            out[n] += out[n-d]
    return [sum(x*out[n-i] for i, x in numerator.items() if i <= n)
            for n in range(cutoff+1)]


# Exact polynomial certificate for the five-case simplification.
lhs = add(mul(factor(4), factor(8)),
          shift(mul(factor(6), factor(8)), 4),
          mul({6: 1, 7: 1, 8: 1}, mul(factor(3), factor(4))),
          shift(mul(factor(3), factor(8)), 10))
rhs = mul(mul(factor(4), factor(5)), P)
require(lhs == rhs, 'Five-case rational simplification failed')
D = {0: 1}
for d in WEIGHTS:
    D = mul(D, factor(d))
expected_D = {0: 1, 3: -1, 5: -1, 6: -1, 9: 1, 11: 2,
              13: 1, 16: -1, 17: -1, 19: -1, 22: 1}
require(D == expected_D, 'Denominator expansion failed')

# Explicit irreducible cyclotomic factors and their exact remainders.
cyclotomic = {
    1: {0: -1, 1: 1},
    2: {0: 1, 1: 1},
    3: {0: 1, 1: 1, 2: 1},
    4: {0: 1, 2: 1},
    5: {0: 1, 1: 1, 2: 1, 3: 1, 4: 1},
    6: {0: 1, 1: -1, 2: 1},
    8: {0: 1, 4: 1},
}
remainders = {d: remainder(P, f) for d, f in cyclotomic.items()}
require(all(remainders.values()), 'Unexpected numerator cancellation')
a, b = P, D
while b:
    a, b = b, remainder(a, b)
require(max(a) == 0, 'Numerator and denominator are not coprime')

closed = series(WEIGHTS, P, N)
regions = [[0] * (N+1) for _ in range(5)]
T = (2*N)//5
cone_vectors = 0
for a in range(T+1):
    for b in range(T//2+1):
        for c in range(T+1):
            lo = max(a+b, 2*b, a+c)
            hi = min(T, (N+2*a+2*b+c)//5)
            if lo > hi:
                continue
            if a >= b:
                region = 0 if c <= b else 1
            elif c <= a:
                region = 2
            elif c <= 2*b-a:
                region = 3
            else:
                region = 4
            for t in range(lo, hi+1):
                weight = 5*t-2*a-2*b-c
                require(0 <= weight <= N, 'Cone enumeration bound failed')
                regions[region][weight] += 1
                cone_vectors += 1
cone = [sum(v[n] for v in regions) for n in range(N+1)]
require(cone == closed, 'Direct four-variable cone coefficients disagree')
region_formulas = [
    ((3, 5, 5, 6), {0: 1}),
    ((3, 4, 5, 5), {4: 1}),
    ((5, 5, 6, 8), {8: 1}),
    ((5, 5, 6, 8), {6: 1, 7: 1}),
    ((4, 5, 5, 6), {10: 1}),
]
for i, (weights, num) in enumerate(region_formulas):
    require(regions[i] == series(weights, num, N),
            'Disjoint case %d failed' % (i+1))

# A structurally independent actual-arc computation. Form anti-isomorphism
# orbits of the 10 arcs, assign nonnegative values to the first 5 orbits,
# solve for the last singleton from its endpoint balance, and check all 5
# vertex balances. This does not use the (a,b,c,t) cone parametrization.
arcs = [(i, (i+j) % 5) for j in (1, 2) for i in range(5)]
anti = lambda e: ((-e[1]) % 5, (-e[0]) % 5)
orbits, seen = [], set()
for e in arcs:
    if e not in seen:
        orbit = sorted({e, anti(e)})
        orbits.append(orbit)
        seen.update(orbit)
require(len(orbits) == 6 and len(orbits[-1]) == 1,
        'Unexpected anti-isomorphism arc orbits')
columns = [[sum((u == v) - (w == v) for u, w in orbit)
            for v in range(5)] for orbit in orbits]
sizes = [len(orbit) for orbit in orbits]
edge = [0] * (EDGE_N+1)
edge_prefixes = 0
balance = [0] * 5
values = [0] * 6
last = columns[-1]
pivot = next(v for v in range(5) if last[v])


def visit(i, used):
    global edge_prefixes
    if i == 5:
        edge_prefixes += 1
        target = -balance[pivot]
        if target % last[pivot]:
            return
        value = target // last[pivot]
        weight = used + sizes[-1]*value
        if value < 0 or weight > EDGE_N:
            return
        if any(balance[v] + last[v]*value for v in range(5)):
            return
        values[5] = value
        literal = {e: values[j] for j, orbit in enumerate(orbits)
                   for e in orbit}
        require(all(literal[e] == literal[anti(e)] for e in arcs),
                'Literal anti-invariance failed')
        require(all(sum(literal[e] for e in arcs if e[0] == v) ==
                    sum(literal[e] for e in arcs if e[1] == v)
                    for v in range(5)), 'Literal balance failed')
        require(sum(literal.values()) == weight, 'Literal weight failed')
        edge[weight] += 1
        return
    column = columns[i]
    for value in range((EDGE_N-used)//sizes[i]+1):
        values[i] = value
        for v in range(5):
            balance[v] += value*column[v]
        visit(i+1, used+value*sizes[i])
        for v in range(5):
            balance[v] -= value*column[v]


visit(0, 0)
require(edge == closed[:EDGE_N+1], 'Actual-arc coefficients disagree')

recurrence = []
for n in range(N+1):
    recurrence.append(P.get(n, 0) - sum(c*recurrence[n-i]
        for i, c in D.items() if i and i <= n))
require(recurrence == closed, 'Exact recurrence disagrees')
require(all(P.get(i, 0) == P.get(12-i, 0) for i in range(13)),
        'Palindromic numerator failed')
A4 = Fraction(sum(P.values()), prod(WEIGHTS))
A3 = A4*sum(Fraction(d-1, 2) for d in WEIGHTS) - Fraction(
    sum(i*c for i, c in P.items()), prod(WEIGHTS))
require(A4/6 == Fraction(7, 4320), 'Leading asymptotic failed')
require(A4+A3/2 == Fraction(7, 288), 'Second asymptotic term failed')

print(json.dumps({
    'status': 'all exact checks passed',
    'cone_cutoff': N,
    'cone_vectors': cone_vectors,
    'disjoint_case_vectors': [sum(x) for x in regions],
    'actual_arc_cutoff': EDGE_N,
    'actual_arc_prefix_assignments': edge_prefixes,
    'actual_arc_fixed_flows': sum(edge),
    'coefficients_0_through_48': closed[:49],
    'coefficient_120': closed[120],
    'coefficient_240': closed[240],
    'coefficient_300': closed[300],
    'denominator_polynomial': D,
    'cyclotomic_remainders': {d: {i: str(c) for i, c in r.items()}
                              for d, r in remainders.items()},
    'leading_asymptotic': '7*n^3/4320 + 7*n^2/288 + O(n)',
    'minimal_quasiperiod': 120,
}, sort_keys=True, indent=2))
