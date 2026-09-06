#!/usr/bin/env python3
"""Exact all-scale Smith hostile to four-node metric-only interpolation.

Universe: all 923 nonempty minors of each six-by-six residual matrix for
base nodes A=(0,1,2,5), B=(0,1,3,4), independently by Bareiss and Leibniz;
full integer Smith forms at scales 2^e, 0<=e<=20; all node permutations at
e=3; optional complete translated-to-zero four-node head, ordered by height.
This does not assume that sorted distance multisets classify distance trees.
"""
import argparse
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import prod
import sympy as s
from sympy.matrices.normalforms import smith_normal_form

GATES = 0
BASES = ((0,1,2,5),(0,1,3,4))


def need(test, detail):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(detail)


def vp(x, p=2):
    need(x != 0, 'valuation requires nonzero input')
    out = 0
    while x % p == 0:
        x //= p
        out += 1
    return out


def bareiss(matrix):
    a = [list(row) for row in matrix]
    n, sign, previous = len(a), 1, 1
    for k in range(n-1):
        pivot = next((r for r in range(k,n) if a[r][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k],a[pivot] = a[pivot],a[k]
            sign = -sign
        value = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator = a[i][j]*value-a[i][k]*a[k][j]
                need(numerator % previous == 0, 'Bareiss exact division')
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = value
    return sign*a[-1][-1]


@lru_cache(None)
def signed_permutations(n):
    return tuple((q,(-1)**sum(q[i]>q[j] for i in range(n) for j in range(i+1,n)))
                 for q in permutations(range(n)))


def leibniz(matrix):
    return sum(sign*prod(matrix[i][q[i]] for i in range(len(matrix)))
               for q,sign in signed_permutations(len(matrix)))


def observer(nodes):
    size = 2*len(nodes)
    return [[x**j if r == 0 else j*x**(j-1) if j else 0
             for j in range(size)] for x in nodes for r in (0,1)]


def literal_smith(nodes):
    d = smith_normal_form(s.Matrix(observer(nodes)), domain=s.ZZ)
    factors = tuple(abs(int(d[i,i])) for i in range(len(d.tolist())))
    need(prod(factors) == prod(abs(x-y)**4 for x,y in combinations(nodes,2)),
         ('full confluent determinant',nodes))
    need(all(factors[i+1] % factors[i] == 0 for i in range(len(factors)-1)),
         ('Smith factor divisibility',nodes))
    return factors


def distances(nodes, p=2):
    return tuple(vp(nodes[j]-nodes[i],p)
                 for i,j in combinations(range(len(nodes)),2))


@lru_cache(None)
def canonical_metric(labelled):
    pairs = tuple(combinations(range(4),2))
    values = {ij:x for ij,x in zip(pairs,labelled)}
    return min(tuple(values[tuple(sorted((q[i],q[j])))] for i,j in pairs)
               for q in permutations(range(4)))


def expected_costs(which):
    return (((2,0),(1,1)),((4,0),(3,2)),((7,2),(6,6)),
            ((12,4),(11,6+which)),((17,7),),((24,12),))


def predicted(which,e):
    sums = (0,)+tuple(min(slope*e+intercept for slope,intercept in row)
                     for row in expected_costs(which))
    return (0,0)+tuple(sums[i+1]-sums[i] for i in range(6))


def minor_certificate(which):
    nodes = BASES[which]
    residual = [row[2:] for row in observer(nodes)[2:]]
    rows = []
    for size in range(1,7):
        costs = {}
        count = 0
        for rr in combinations(range(6),size):
            for cc in combinations(range(6),size):
                matrix = [[residual[i][j] for j in cc] for i in rr]
                determinant = bareiss(matrix)
                need(determinant == leibniz(matrix), ('independent determinant',which,rr,cc))
                count += 1
                if not determinant:
                    continue
                slope = sum(j+2 for j in cc)-sum(i%2 for i in rr)
                intercept = vp(determinant)
                if slope not in costs or intercept < costs[slope][0]:
                    costs[slope] = (intercept,rr,cc,determinant)
        frontier = {d:x for d,x in costs.items()
                    if not any(d2<=d and x2[0]<=x[0] and (d2,x2[0])!=(d,x[0])
                               for d2,x2 in costs.items())}
        need(set((d,x[0]) for d,x in frontier.items()) == set(expected_costs(which)[size-1]),
             ('all-scale Pareto envelope',which,size,frontier))
        rows.append((size,count,sorted(frontier.items())))
        print('MINOR_ENVELOPE',which,size,count,sorted(frontier.items()))
    need(sum(row[1] for row in rows) == 923, 'complete nonempty minor universe')
    return rows


def complete_head(height):
    seen, twins, cases = {}, [], 0
    for top in range(3,height+1):
        for a,b in combinations(range(1,top),2):
            nodes = (0,a,b,top)
            profile = tuple(vp(x) for x in literal_smith(nodes))
            key = canonical_metric(distances(nodes))
            if key in seen and seen[key][1] != profile:
                twins.append((seen[key],(nodes,profile),key))
            else:
                seen.setdefault(key,(nodes,profile))
            cases += 1
    print('COMPLETE_HEAD',height,'cases',cases,'isometry_types',len(seen),'twins',len(twins))
    if twins:
        print('FIRST_HEIGHT_ORDERED_METRIC_TWIN',twins[0])
    if height >= 40:
        need(bool(twins) and twins[0][1][0][-1] == 40, 'first-height hostile boundary')
    return twins


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--height',type=int,default=0)
    args = parser.parse_args()
    certificates = [minor_certificate(which) for which in range(2)]
    scale_rows = []
    for e in range(21):
        profiles = []
        for which,base in enumerate(BASES):
            nodes = tuple(2**e*x for x in base)
            profile = tuple(vp(x) for x in literal_smith(nodes))
            need(profile == predicted(which,e), ('all-scale full Smith replay',which,e,profile))
            profiles.append(profile)
        need(distances(tuple(2**e*x for x in BASES[0])) ==
             distances(tuple(2**e*x for x in (1,0,3,4))), 'explicit labelled isometry')
        need((profiles[0] == profiles[1]) == (e<=2), 'exact first failure depth')
        need(profiles[0][-1] == profiles[1][-1] == 7*e+5, 'largest precision is blind')
        scale_rows.append((e,profiles))
        print('SCALE',e,profiles)
    for which,base in enumerate(BASES):
        reference = literal_smith(tuple(8*x for x in base))
        for order in permutations(base):
            need(literal_smith(tuple(8*x for x in order)) == reference, 'full node permutation audit')
    # Every corresponding three-node restriction has the same p-Smith profile.
    a,b = tuple(8*x for x in BASES[0]),tuple(8*x for x in (1,0,3,4))
    for indices in combinations(range(4),3):
        pa = tuple(vp(x) for x in literal_smith(tuple(a[i] for i in indices)))
        pb = tuple(vp(x) for x in literal_smith(tuple(b[i] for i in indices)))
        need(pa == pb, ('all three-node restrictions blind',indices,pa,pb))
    # A determinant-only or largest-loss checker must miss this witness.
    need(sum(predicted(0,3)) == sum(predicted(1,3)) and predicted(0,3) != predicted(1,3),
         'determinant-only mutation is detected')
    twins = complete_head(args.height) if args.height else []
    payload = (certificates,scale_rows,twins)
    print('SEMANTIC_SHA256',sha256(json.dumps(payload).encode()).hexdigest())
    print('EXPLICIT_GATES',GATES)
    print('RESULT=PASS metric-only four-node Smith claim REFUTED; three-node theorem unaffected')


if __name__ == '__main__':
    main()
