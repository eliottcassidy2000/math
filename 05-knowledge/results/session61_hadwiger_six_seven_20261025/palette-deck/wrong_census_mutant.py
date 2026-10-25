#!/usr/bin/env python3
"""Exact, assertion-free controls for palette/Fano/deck transfer boundaries."""
import itertools as it
import json
from collections import Counter

CHECKS = 0
def check(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(message)

def bits(x, n):
    return tuple((x >> i) & 1 for i in range(n))

def dot(x, y):
    return (x & y).bit_count() & 1

def edges(n, mask):
    return [e for i, e in enumerate(it.combinations(range(n), 2)) if mask >> i & 1]

def solve(n, es, a):
    adj = [[] for _ in range(n)]
    for (u, v), w in zip(es, a):
        adj[u].append((v, w)); adj[v].append((u, w))
    p = [None] * n
    comps = 0
    for root in range(n):
        if p[root] is not None:
            continue
        comps += 1
        p[root] = 0
        stack = [root]
        while stack:
            u = stack.pop()
            for v, w in adj[u]:
                if p[v] is None:
                    p[v] = p[u] ^ w; stack.append(v)
                elif p[v] != p[u] ^ w:
                    return None, comps
    return p, comps

def component_count(n, es):
    return solve(n, es, [0] * len(es))[1]

cochain_cases = residual_cases = 0
for n in range(5):
    for gm in range(1 << (n * (n-1) // 2)):
        es = edges(n, gm)
        for a in it.product(range(2), repeat=len(es)):
            direct = [x for x in range(1 << n)
                      if all(((x >> u) ^ (x >> v)) & 1 == w
                             for (u, v), w in zip(es, a))]
            p, comps = solve(n, es, a)
            check(bool(direct) == (p is not None), 'forest versus exhaustive lift')
            if p is not None:
                check(len(direct) == 2 ** comps, 'component gauge count')
            cochain_cases += 1
        # The all-negative cover is always bipartite, regardless of base.
        lifted = [(2*u+b, 2*v+(b^1)) for u, v in es for b in range(2)]
        check(all((u & 1) != (v & 1) for u,v in lifted), 'cover bipartition')
        if component_count(n, es) == 1:
            balanced = solve(n, es, [1]*len(es))[0] is not None
            check(component_count(2*n, lifted) == (2 if balanced else 1),
                  'connected base section/component criterion')
        for q in it.product(range(3), repeat=n):
            residual = [(u,v) for u,v in es if q[u] == q[v]]
            actual = any(all(2*q[u]+b[u] != 2*q[v]+b[v] for u,v in es)
                         for b in it.product(range(2), repeat=n))
            check(actual == (solve(n,residual,[1]*len(residual))[0] is not None),
                  'fixed-ternary colouring lift')
            residual_cases += 1

# Every labeled Fano kernel, independently deduplicated from its simplex words.
codes = {}
for labels in it.permutations(range(1,8)):
    code = frozenset(sum(dot(a,v) << i for i,v in enumerate(labels)) for a in range(8))
    codes.setdefault(code, labels)
check(len(codes) == 30, 'thirty Fano kernels')
triples = tuple(it.combinations(range(7),3))
crossings = ((0,1,2),(0,3,4),(0,5,6))
pencils = {}
hist = Counter()
triple_checks = 0
for code, labels in codes.items():
    # The even-representative kernel of q is the simplex row space.
    kernel = set()
    for x in range(128):
        if x.bit_count() & 1:
            continue
        syndrome = 0
        for i,v in enumerate(labels):
            if x >> i & 1:
                syndrome ^= v
        if not syndrome:
            kernel.add(x)
    check(kernel == code, 'projective quotient kernel')
    check(Counter(x.bit_count() for x in code) == {0:1,4:7}, 'simplex weights')
    for i,j in it.combinations(range(7),2):
        check({((x>>i)&1,(x>>j)&1) for x in code} == set(it.product(range(2), repeat=2)),
              'identical full pair shadows')
        check(Counter(((x>>i)&1,(x>>j)&1) for x in code) ==
              {b:2 for b in it.product(range(2),repeat=2)}, 'identical uniform pair marginals')
    fano_lines = []
    for t in triples:
        proj = {tuple((x>>i)&1 for i in t) for x in code}
        dependent = labels[t[0]] ^ labels[t[1]] ^ labels[t[2]] == 0
        check(len(proj) == (4 if dependent else 8), 'triple projections recover dependence')
        correlation = sum((-1)**sum((x>>i)&1 for i in t) for x in code)
        check(correlation == (8 if dependent else 0), 'canonical-character cubic correlation')
        complement = 127 ^ sum(1 << i for i in t)
        check((complement in code) == dependent, 'projective weight-four character gate')
        if dependent:
            check(all(sum(w)%2 == 0 for w in proj), 'dependent triple is even plane')
            fano_lines.append(t)
        triple_checks += 1
    matching = tuple(sorted(tuple(i for i in t if i != 0) for t in fano_lines if 0 in t))
    check(len(matching) == 3, 'central matching')
    pencils.setdefault(matching, []).append(code)
    hist[sum(t in fano_lines for t in crossings)] += 1
check(len(pencils) == 15 and all(len(v)==2 for v in pencils.values()), 'two completions each')
check(hist == {0:15,1:12,3:2}, 'hourglass dependency census')

def permute_word(x, perm):
    return sum(((x>>i)&1) << perm[i] for i in range(len(perm)))

for matching, pair_codes in pencils.items():
    a,b = matching[0]
    perm = list(range(7)); perm[a],perm[b] = perm[b],perm[a]
    moved = frozenset(permute_word(x,perm) for x in pair_codes[0])
    check(moved == pair_codes[1], 'one endpoint flip exchanges completions')
    # This permutation fixes the matching and every ordering of its pairs.
    check(tuple(sorted(tuple(sorted((perm[a],perm[b]))) for a,b in matching)) == matching,
          'ordered-pair data cannot select completion')

# Six labels: even subsets modulo complement form W=F2^4.
def rep6(x):
    return min(x,x^63)
W = sorted({rep6(x) for x in range(64) if x.bit_count()%2 == 0})
check(len(W) == 16, 'six-branch two-torsion vector space')
lagrangians = set()
for a,b in it.combinations(W[1:],2):
    if dot(a,b) == 0:
        lagrangians.add(frozenset((0,a,b,rep6(a^b))))
match_lagrangians = set()
for matching in pencils:
    pairs = [sum(1 << (i-1) for i in p) for p in matching]
    L = frozenset([0] + [rep6(x) for x in pairs])
    check(len(L)==4 and all(dot(a,b)==0 for a in L for b in L), 'matching Lagrangian')
    match_lagrangians.add(L)
check(len(lagrangians)==15 and match_lagrangians==lagrangians, 'all fifteen matchings/Lagrangians')

# A marked three-bit palette gauge enters the six-root deck by choosing one
# representative root in each pair. Its canonical character is parity.
for x in range(8):
    root_word = ((x&1)<<0) | (((x>>1)&1)<<2) | (((x>>2)&1)<<4)
    check(root_word.bit_count()%2 == x.bit_count()%2, 'restricted six-root character')
check(sum(x.bit_count()%2 == 0 for x in range(8)) == 4, 'lost two-bit even kernel')

# Two independent odd channels cancel after the single parity observer.
triangle = [(0,1),(1,2),(0,2)]
gains = [(1,1,0)] * 3
check(solve(3,triangle,[sum(x)%2 for x in gains])[0] is not None, 'scalar observer accepts')
check(solve(3,triangle,[x[0] for x in gains])[0] is None, 'first true channel obstructed')
check(solve(3,triangle,[x[1] for x in gains])[0] is None, 'second true channel obstructed')

# Inherited induced unit-distance three-port hostile, independently rebuilt.
points = ((0,0),(1,0),(0,1),(0,-1),(2,0),(0,2))
def norm_delta(a,b):
    x,y = a[0]-b[0], a[1]-b[1]
    return x*x+x*y+y*y # alpha=(1+i*sqrt(3))/2
unit_edges = [(i,j) for i,j in it.combinations(range(6),2) if norm_delta(points[i],points[j])==1]
check(unit_edges == [(0,1),(0,2),(0,3),(1,2),(1,4),(2,5)], 'exact induced geometry')
colourings = [c for c in it.product(range(3),repeat=6) if all(c[u]!=c[v] for u,v in unit_edges)]
boundary = {c[3:] for c in colourings}
check(len(colourings)==48 and len(boundary)==24, 'geometric three-port counts')
check(boundary == {c for c in it.product(range(3),repeat=3) if len(set(c))>1}, 'exact NAE relation')
for i,j in it.combinations(range(3),2):
    check({(c[i],c[j]) for c in boundary} == set(it.product(range(3),repeat=2)), 'NAE full pair shadows')

# Generic six-palette higher-arity hostile; NOT unit-distance realizable.
# K4 joined to three independent ports extends exactly when <=2 port colors.
port_relation = set()
for b in it.product(range(6),repeat=3):
    available = set(range(6))-set(b)
    if any(len(set(c))==4 for c in it.product(available,repeat=4)):
        port_relation.add(b)
check(len(port_relation)==96, 'six-palette abstract boundary count')
check(port_relation == {b for b in it.product(range(6),repeat=3) if len(set(b))<=2}, 'six-palette at-most-two')
for i,j in it.combinations(range(3),2):
    check({(b[i],b[j]) for b in port_relation} == set(it.product(range(6),repeat=2)), 'six-palette full pair shadows')
check(not any(b in port_relation for b in it.permutations(range(6),3)), 'gluing port triangle is non-six-colourable')

print(json.dumps({
    'status':'PASS', 'always_active_checks':CHECKS,
    'finite_universes':{'all_graphs_orders_0_through_4':76,
        'edge_cochain_cases':cochain_cases,'fixed_ternary_cases':residual_cases,
        'fano_codes':len(codes),'fano_triple_projections':triple_checks,
        'central_matchings':len(pencils),'genus_two_lagrangians':len(lagrangians)},
    'hourglass_dependency_histogram':dict(sorted(hist.items())),
    'geometric_NAE':{'vertices':6,'edges':unit_edges,'three_colourings':48,'boundary_assignments':24},
    'abstract_six_colour_relation':{'boundary_assignments':len(port_relation),'unit_realization':False},
    'negative_controls':['C3 sign cover C6 always bipartite; no section',
        'one pair endpoint flip fixes ordered matching, exchanges Fano completion',
        'all pair shadows agree across 30 kernels; triple shadows distinguish',
        'even nonzero channel holonomy (1,1,0) invisible to parity',
        'pairwise NAE shadow falsely admits all-equal ports',
        'six-colour pair shadows falsely accept triangle gluing'],
},indent=2,sort_keys=True))
