#!/usr/bin/env python3
"""Portable exact audit. Standard library only; no producer imports or files needed."""
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement, permutations, product
from math import prod
import json

checks = 0

def need(value, label):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(label)

def determinant(a):
    """Leibniz formula, deliberately independent of elimination."""
    n = len(a)
    out = Q(0)
    for p in permutations(range(n)):
        sign = (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        out += sign * prod((a[i][p[i]] for i in range(n)), start=Q(1))
    return out

def minor(a, i, j):
    return [[v for c, v in enumerate(row) if c != j] for r, row in enumerate(a) if r != i]

def inverse(a):
    """Adjugate/cofactor formula, not producer Gaussian elimination."""
    z = determinant(a)
    need(z != 0, 'inverse is defined')
    return [[(-1) ** (i+j) * determinant(minor(a, j, i)) / z for j in range(len(a))] for i in range(len(a))]

def multiply(a, b):
    return [[sum((x*y for x, y in zip(row, col)), Q(0)) for col in zip(*b)] for row in a]

def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Q(0))

def edge(k, i, j):
    return [Q(t == i) - Q(t == j) for t in range(k-1)]

def quadratic(a, b):
    return dot(b, [dot(row, b) for row in a])

def laplacian(k, weights):
    out = [[Q(0) for _ in range(k-1)] for _ in range(k-1)]
    for (i, j), w in weights.items():
        b = edge(k, i, j)
        for a in range(k-1):
            for c in range(k-1):
                out[a][c] += w*b[a]*b[c]
    return out

def connected(k, edges):
    seen = {0}
    while True:
        extended = seen | {j for i, j in edges if i in seen} | {i for i, j in edges if j in seen}
        if extended == seen:
            return len(seen) == k
        seen = extended

def tree_count(k, weights):
    out = Q(0)
    for edges in combinations(weights, k-1):
        if connected(k, edges):
            out += prod((weights[e] for e in edges), start=Q(1))
    return out

def feasible(v, s):
    return all(quadratic(s, edge(len(v), i, j)) <= (Q(1, v[i])+Q(1, v[j]))**2 for i, j in combinations(range(len(v)), 2))

def certificate(v, weights):
    k = len(v); d = k-1
    need(all(w >= 0 for w in weights.values()), 'nonnegative weights')
    support = [e for e in weights if weights[e] > 0]
    need(connected(k, support), 'positive support connected')
    l = laplacian(k, weights)
    need(all(determinant([row[:j] for row in l[:j]]) > 0 for j in range(1, d+1)), 'L positive definite by Sylvester')
    s = inverse(l)
    need(multiply(l, s) == [[Q(i == j) for j in range(d)] for i in range(d)], 'adjugate inverse')
    need(feasible(v, s), 'every pair bound')
    resistance = {}
    budget = Q(0)
    for i, j in combinations(range(k), 2):
        z = quadratic(s, edge(k, i, j)); resistance[i, j] = z
        cap = (Q(1, v[i])+Q(1, v[j]))**2
        w = weights.get((i, j), Q(0))
        need(w == 0 or z == cap, 'weighted pair contact')
        budget += w*cap
    need(budget == d, 'trace budget equals dimension')
    need(tree_count(k, weights) == determinant(l), 'tree enumeration equals determinant')
    for h in product((-1, 0, 1), repeat=d):
        if any(h):
            need(any(dot(edge(k, *e), h) != 0 for e in support), 'nonzero translation hits a contact')
    return s, l, resistance

# Direct quotient map and interval characterization of P[-1,1]^k.
quotient_cases = 0
for v in [(1, 2, 4), (1, 2, 3, 5), (2, 2, 3), (3, 6, 9, 12)]:
    k = len(v); d = k-1; norm = sum(x*x for x in v)
    g = [[Q(v[i])*(Q(i == j)-Q(v[j]*v[j], norm)) for j in range(d)] for i in range(k)]
    gram = multiply(list(map(list, zip(*g))), g)
    need(determinant(gram) == Q(prod(x*x for x in v), norm), 'quotient volume factor')
    for z in product((-1, 1), repeat=k):
        x = [Q(z[i])-Q(v[i]*sum(a*b for a, b in zip(v, z)), norm) for i in range(k)]
        y = [x[i]/v[i]-x[-1]/v[-1] for i in range(d)]
        need([dot(row, y) for row in g] == x, 'gauge inverse on cube projections')
        need(y == [Q(z[i], v[i])-Q(z[-1], v[-1]) for i in range(d)], 'quotient ignores projection direction')
        need(all(abs((y+[Q(0)])[i]-(y+[Q(0)])[j]) <= Q(1,v[i])+Q(1,v[j]) for i,j in combinations(range(k),2)), 'projected vertices meet facets')
        quotient_cases += 1
    for short_y in product((Q(-2), Q(-1,2), Q(0), Q(1,3), Q(2)), repeat=d):
        y = short_y+(Q(0),)
        lower = max(-Q(1,v[i])-y[i] for i in range(k))
        upper = min(Q(1,v[i])-y[i] for i in range(k))
        pair_ok = all(abs(y[i]-y[j]) <= Q(1,v[i])+Q(1,v[j]) for i,j in combinations(range(k),2))
        need(pair_ok == (lower <= upper), 'pair facets iff common cube lift interval')
        if pair_ok:
            t = (lower+upper)/2
            need(all(abs(v[i]*(y[i]+t)) <= 1 for i in range(k)), 'construct cube preimage')
        quotient_cases += 1

literal = [
    ((1,2,4), {(0,2):Q(16,25), (1,2):Q(16,9)}, Q(225,256), Q(75,28)),
    ((1,2,3,5), {(0,3):Q(25,36), (1,2):Q(9,31), (1,3):Q(55,31), (2,3):Q(405,124)}, Q(124,625), Q(1488,325)),
]
literal_results = []
bound_cases = 0
for v, weights, expected_det, expected_euclidean in literal:
    k = len(v); d = k-1
    s, l, resist = certificate(v, weights)
    need(determinant(s) == expected_det, 'literal covariance determinant')
    need(determinant(s)*Q(prod(x*x for x in v), sum(x*x for x in v)) == expected_euclidean, 'literal Euclidean determinant transport')
    eps = min((Q(1,v[i])+Q(1,v[j]))**2 for i,j in combinations(range(k),2))/8
    sample_covariances = []
    for t in (Q(0), Q(1,4), Q(1,2), Q(1)):
        tmat = [[(1-t)*s[i][j]+t*eps*Q(i==j) for j in range(d)] for i in range(d)]
        need(feasible(v, tmat), 'anisotropic feasible convex sample')
        sample_covariances.append(tmat)
    networks = [weights,
                {(i,j):Q(1) for i,j in combinations(range(k),2)},
                {(i,k-1):Q(i+1,3) for i in range(d)}]
    for w in networks:
        test_l = laplacian(k, w)
        c = sum(w[e]*(Q(1,v[e[0]])+Q(1,v[e[1]]))**2 for e in w)
        for tmat in sample_covariances:
            trace = sum(test_l[i][j]*tmat[j][i] for i in range(d) for j in range(d))
            need(trace <= c, 'arbitrary-network trace bound')
            need(determinant(tmat)*determinant(test_l) <= (trace/d)**d, 'AMGM nonhomothetic sample')
            need(determinant(tmat)*determinant(test_l) <= (c/d)**d, 'general bound (4)')
            bound_cases += 1
    if k == 4:
        reduced = dict(weights); del reduced[1,2]
        need(connected(k, reduced), 'extra-edge deletion connected')
        deleted = inverse(laplacian(k, reduced))
        need(not feasible(v, deleted), 'extra-edge deletion fails feasibility')
        need(s[1][2] == Q(1,25), 'literal triangle off-diagonal')
        need(s[1][2] / ((Q(1,2)+Q(1,5))*(Q(1,3)+Q(1,5))) == Q(3,28), 'triangle correlation')
    literal_results.append({'speeds':v, 'grounded_det':str(determinant(s)), 'euclidean_det':str(expected_euclidean), 'weighted_tree_sum':str(tree_count(k,weights)), 'resistances':{str(e):str(r) for e,r in resist.items()}})

# Independently classify every tree in declared bounded sets, without matrix inverses.
def paths(k, edges, r):
    adjacency = [[] for _ in range(k)]
    for a,b in edges:
        c = (r[a]+r[b])**2
        adjacency[a].append((b,c)); adjacency[b].append((a,c))
    result = {}
    for start in range(k):
        stack = [(start,-1,Q(0))]
        while stack:
            a, parent, total = stack.pop()
            if a > start:
                result[start,a] = total
            for b, weight in adjacency[a]:
                if b != parent:
                    stack.append((b,a,total+weight))
    return result

tree_rows = []
for k, upper in [(3,6),(4,5),(5,3)]:
    complete = list(combinations(range(k),2))
    trees = [es for es in combinations(complete,k-1) if connected(k,es)]
    need(len(trees) == k**(k-2), 'labelled tree enumeration control')
    # Acyclic rank criterion and every square incidence determinant.
    for es in combinations(complete,k-1):
        incidence = [[edge(k,*e)[j] for e in es] for j in range(k-1)]
        m = determinant(incidence)
        need(m*m == Q(connected(k,es)), 'incidence minor exactly tree indicator')
    cases = valid = outside = 0
    for v in combinations_with_replacement(range(1,upper+1),k):
        r = [Q(1,x) for x in v]
        cone = v[-1]**2 >= v[-1]*(v[-2]+v[-3])+v[-2]*v[-3]
        found = False
        for es in trees:
            distance = paths(k,es,r)
            ok = all(z <= (r[i]+r[j])**2 for (i,j),z in distance.items())
            degrees = [sum(i in e for e in es) for i in range(k)]
            pivots = [i for i,d in enumerate(degrees) if d == k-1]
            condition = False
            if pivots:
                p = pivots[0]; other = sorted(v[i] for i in range(k) if i != p)
                condition = v[p]**2 >= v[p]*(other[-1]+other[-2])+other[-1]*other[-2]
            need(ok == condition, 'tree criterion iff dominant star')
            if ok:
                found = True; valid += 1
                need(pivots and v[pivots[0]] == max(v), 'certifying tree has fastest centre')
            cases += 1
        need(found == cone, 'existence of tree iff dominant cone')
        outside += not cone
    tree_rows.append({'vertices':k, 'speed_box':upper, 'trees_per_vector':len(trees), 'tree_vector_cases':cases, 'valid_tree_vector_cases':valid, 'outside_cone_vectors':outside})

# Disconnected positive support cannot be used as a grounded inverse.
need(determinant(laplacian(4,{(0,1):Q(1),(1,2):Q(2),(0,2):Q(3)})) == 0, 'isolated grounded vertex is singular')
# Repeat the complete-graph matrix-tree identity with uneven positive weights.
for k in (3,4,5):
    weights = {(i,j):Q((i+1)*(j+2),i+j+1) for i,j in combinations(range(k),2)}
    need(determinant(laplacian(k,weights)) == tree_count(k,weights), 'general weighted matrix-tree identity')

print(json.dumps({'status':'PASS', 'assertions':'explicit RuntimeError checks survive python -O', 'checks':checks, 'quotient_cases':quotient_cases, 'general_bound_samples':bound_cases, 'literal_examples':literal_results, 'tree_classification':tree_rows, 'scope':'Finite exact independent controls plus separate universal proof audit; no LRC or literature-priority claim.'}, indent=2))
