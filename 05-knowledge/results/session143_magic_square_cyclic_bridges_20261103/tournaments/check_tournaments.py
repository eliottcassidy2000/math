#!/usr/bin/env python3
"""Exact, dependency-free order-five tournament census and C5-family repair.

Run: python3 check_tournaments.py > classification.json
All validations use explicit exceptions and remain active under python3 -O.
Encoding: bit k=1 means i->j for kth pair i<j in lexicographic order.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations
import json

N = 5
EDGES = tuple(combinations(range(N), 2))
EDGE_INDEX = {edge: k for k, edge in enumerate(EDGES)}
PERMS = tuple(permutations(range(N)))
FULL = (1 << len(EDGES)) - 1
IDENTITY = tuple(range(N))
CHECKS = 0

def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)

def matrix(mask):
    a = [[0] * N for _ in range(N)]
    for k, (i, j) in enumerate(EDGES):
        a[i][j] = (mask >> k) & 1
        a[j][i] = 1 - a[i][j]
    return a

def encode(a):
    return sum(a[i][j] << k for k, (i, j) in enumerate(EDGES))

def act(mask, p):
    result = 0
    for k, (i, j) in enumerate(EDGES):
        u, v, bit = p[i], p[j], (mask >> k) & 1
        if u > v:
            u, v, bit = v, u, 1 - bit
        result |= bit << EDGE_INDEX[u, v]
    return result

def conjugate(a, p):
    b = [[0] * N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            b[p[i]][p[j]] = a[i][j]
    return b

def compose(p, q):
    return tuple(p[q[i]] for i in range(N))

def order(p):
    q = IDENTITY
    for k in range(1, 121):
        q = compose(p, q)
        if q == IDENTITY:
            return k
    raise RuntimeError('permutation order exceeded S5 bound')

def cycle_type(p):
    unseen = set(range(N))
    lengths = []
    while unseen:
        x = min(unseen)
        length = 0
        while x in unseen:
            unseen.remove(x)
            length += 1
            x = p[x]
        lengths.append(length)
    return tuple(sorted(lengths))

def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(N))
             for j in range(N)] for i in range(N)]

def identity():
    return [[int(i == j) for j in range(N)] for i in range(N)]

def char_newton(a):
    """Descending coefficients, independently from exact power traces."""
    power = identity()
    traces = [N]
    for _ in range(N):
        power = mul(power, a)
        traces.append(sum(power[i][i] for i in range(N)))
    coeff = [1]
    for k in range(1, N + 1):
        numerator = -sum(coeff[k - j] * traces[j] for j in range(1, k + 1))
        require(numerator % k == 0, 'Newton coefficient nonintegral')
        coeff.append(numerator // k)
    return coeff, traces

def char_determinant(a):
    """Leibniz expansion of det(xI-A), independent of Newton/cycle counts."""
    out = [0] * (N + 1)  # ascending
    for p in PERMS:
        inversions = sum(p[i] > p[j] for i in range(N) for j in range(i + 1, N))
        sign = (-1) ** inversions
        degree = 0
        for i, j in enumerate(p):
            if i == j:
                degree += 1
            else:
                sign *= -a[i][j]
        out[degree] += sign
    return list(reversed(out))

def cycle_counts(a):
    counts = []
    for k in (3, 4, 5):
        count = 0
        for subset in combinations(range(N), k):
            first = min(subset)
            for tail in permutations(subset[1:]):
                cycle = (first,) + tail
                count += all(a[cycle[i]][cycle[(i + 1) % k]] for i in range(k))
        counts.append(count)
    return counts

def scc_sizes_ordered(a):
    reach = [[bool(a[i][j] or i == j) for j in range(N)] for i in range(N)]
    for k in range(N):
        for i in range(N):
            for j in range(N):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    classes = []
    unused = set(range(N))
    while unused:
        i = min(unused)
        c = {j for j in unused if reach[i][j] and reach[j][i]}
        classes.append(c)
        unused -= c
    # Tournament SCC quotient is a transitive tournament; source first.
    all_classes = tuple(classes)
    classes.sort(key=lambda c: -sum(a[next(iter(c))][next(iter(d))] for d in all_classes if d != c))
    require(all(a[u][v] for i, c in enumerate(classes) for d in classes[i + 1:] for u in c for v in d),
            'SCC word is source-to-sink ordered')
    return [len(c) for c in classes]

def cayley(connection):
    return [[int((j - i) % N in connection) for j in range(N)] for i in range(N)]

def linear_combo(coefficients, matrices):
    return [[sum(c * a[i][j] for c, a in zip(coefficients, matrices))
             for j in range(N)] for i in range(N)]

canonical = {}
orbits = {}
for mask in range(FULL + 1):
    if mask in canonical:
        continue
    orbit = {act(mask, p) for p in PERMS}
    rep = min(orbit)
    require(rep == mask, 'ordered first orbit representative')
    require(not (set(canonical) & orbit), 'disjoint orbits')
    orbits[rep] = orbit
    for member in orbit:
        canonical[member] = rep
require(len(canonical) == 1024 and len(orbits) == 12, '1024 labels,12 classes')

rows = []
for rep, orbit in orbits.items():
    a = matrix(rep)
    aut = [p for p in PERMS if act(rep, p) == rep]
    anti = [p for p in PERMS if act(rep, p) == (rep ^ FULL)]
    coeff, traces = char_newton(a)
    cycles = cycle_counts(a)
    require(coeff == char_determinant(a), 'Newton versus determinant')
    require(coeff == [1, 0, 0, -cycles[0], -cycles[1], -cycles[2]], 'cycle coefficient formula')
    require(len(orbit) * len(aut) == 120, 'orbit stabilizer')
    require(all(order(p) % 2 == 1 for p in aut), 'odd automorphism orders')
    require(len(anti) in (0, len(aut)), 'anti is empty or Aut coset')
    require(all(order(p) == 2 for p in anti), 'n5 anti permutations are involutions')
    for p in PERMS:
        require(encode(conjugate(a, p)) == act(rep, p), 'two permutation-action implementations')
    for m in orbit:
        b = matrix(m)
        cc, tr = char_newton(b)
        require(cc == coeff, 'whole labeled universe characteristic invariant')
        require(tr[1] == tr[2] == 0, 'loop and two-cycle trace controls')
    extended = aut + anti
    for p in extended:
        for q in extended:
            require(compose(p, q) in extended, 'extended permutation group closure')
    if len(aut) > 1 and anti:
        g = next(p for p in aut if p != IDENTITY)
        h = anti[0]
        require(compose(h, compose(g, h)) == tuple(g.index(i) for i in range(N)), 'dihedral relation')
    rows.append(dict(mask=rep, adjacency=a, score=sorted(map(sum, a)),
                     labeled_orbit_size=len(orbit), dual_mask=canonical[rep ^ FULL],
                     self_dual=bool(anti), aut_order=len(aut), anti_count=len(anti),
                     automorphisms=[list(p) for p in aut], anti_automorphisms=[list(p) for p in anti],
                     extended_group=('D' + str(2 * len(aut))) if anti else ('C' + str(len(aut))),
                     characteristic_desc=coeff, cycles_3_4_5=cycles,
                     ordered_scc_sizes=scc_sizes_ordered(a)))

require(sum(row['self_dual'] for row in rows) == 8, '8 self-dual classes')
dual_orbits = sorted({tuple(sorted({r, canonical[r ^ FULL]})) for r in orbits})
require(len(dual_orbits) == 10, '10 converse-quotient classes')
spectra = defaultdict(list)
for row in rows:
    spectra[tuple(row['characteristic_desc'])].append(row['mask'])
require(len(spectra) == 9, 'nine characteristic polynomials')
require(len({(tuple(r['score']), tuple(r['characteristic_desc'])) for r in rows}) == 12,
        'spectrum plus score complete invariant for n5')
require(len({(tuple(r['ordered_scc_sizes']), tuple(r['characteristic_desc'])) for r in rows}) == 12,
        'spectrum plus ordered SCC-size word complete invariant for n5')
require(len({(min(tuple(r['ordered_scc_sizes']), tuple(reversed(r['ordered_scc_sizes']))), tuple(r['characteristic_desc'])) for r in rows}) == 10,
        'spectrum plus SCC word modulo reversal complete for converse quotient')
require(sorted(v for v in spectra.values() if len(v) > 1) == [[2, 5, 11], [4, 9]], 'precise collisions')
regular = [r for r in rows if r['score'] == [2] * 5]
require(len(regular) == 1 and regular[0]['mask'] == 76, 'unique regular class')

# Burnside independent count on the full labeled universe and the signed S5 action.
fixed_by_type = defaultdict(set)
twisted_by_type = defaultdict(set)
burnside = 0
signed_extra = 0
for p in PERMS:
    ordinary = sum(act(m, p) == m for m in range(1024))
    twisted = sum(act(m, p) == (m ^ FULL) for m in range(1024))
    fixed_by_type[cycle_type(p)].add(ordinary)
    twisted_by_type[cycle_type(p)].add(twisted)
    require(ordinary == {(1,1,1,1,1):1024, (1,1,3):16, (5,):4}.get(cycle_type(p), 0),
            'exact ordinary Burnside fixed-set count by cycle type')
    require(twisted == (64 if cycle_type(p) == (1,2,2) else 0),
            'exact converse-twisted Burnside fixed-set count by cycle type')
    burnside += ordinary
    signed_extra += twisted
signed_burnside = burnside + signed_extra
require(burnside == 120 * 12, 'ordinary Burnside check')
require(signed_burnside == 240 * 10, 'converse Burnside check')

connections = [(1, 2), (2, 4), (3, 4), (1, 3)]
family = [encode(cayley(s)) for s in connections]
require(len(set(family)) == 4, 'four distinct circulant orientations')
require({canonical[m] for m in family} == {76}, 'one unlabeled circulant class')
affines = {tuple((a * x + b) % 5 for x in range(5)): (a, b)
           for a in range(1, 5) for b in range(5)}
family_preservers = [p for p in PERMS if {act(m, p) for m in family} == set(family)]
require(set(family_preservers) == set(affines) and len(affines) == 20, 'full family stabilizer is AGL1(5)')
translations = {tuple((x + b) % 5 for x in range(5)) for b in range(5)}
family_kernel = {p for p in family_preservers if all(act(m, p) == m for m in family)}
require(family_kernel == translations, 'family-action kernel C5')
multiplier2 = tuple(2 * x % 5 for x in range(5))
for i, m in enumerate(family):
    require(act(m, multiplier2) == family[(i + 1) % 4], 'order4 multiplier cycles family')
    require((m ^ FULL) == family[(i + 2) % 4], 'multiplier squared is converse')
    require({p for p in PERMS if act(m, p) == m} == translations, 'individual automorphism group C5')
    ext = {p for p in PERMS if act(m, p) in {m, m ^ FULL}}
    require(ext == {p for p, (a, b) in affines.items() if a in (1, 4)}, 'individual extended group D10')

# Exact adjacency algebra and radical-spectrum identities, without numerical roots.
p = cayley((1,))
a = cayley((1, 2))
powers = [identity()]
for k in range(1, 6):
    powers.append(mul(powers[-1], a))
require(a == linear_combo([1, 1], [p, mul(p, p)]), 'A=P+P^2')
require(linear_combo([-17, -18, -2, -1, 5], powers[:5]) == linear_combo([11], [p]), '11P polynomial reconstruction')
require(linear_combo([-2, -5, -5, 0, 0, 1], powers) == [[0] * 5 for _ in range(5)], 'minimal polynomial annihilation')
require(linear_combo([1, 3, 4, 2, 1], powers[:5]) == [[11] * 5 for _ in range(5)], 'quartic spectral projector')
skew = [[a[i][j] - a[j][i] for j in range(5)] for i in range(5)]
require(char_determinant(skew) == [1, 0, 10, 0, 5, 0], 'skew characteristic x(x4+10x2+5)')
# For y=2x+1, 16*q((y-1)/2)=y4+10y2+5, whose y2 roots are -5+-2sqrt5.
# Verify transformation with exact integer binomial expansion.
from math import comb
q_asc = [1, 3, 4, 2, 1]
transformed = [sum(q_asc[k] * 2 ** (4-k) * comb(k, j) * (-1) ** (k-j)
                   for k in range(j, 5)) for j in range(5)]
require(transformed == [5, 0, 10, 0, 1], 'radical spectrum identity')

output = dict(status='FINITE-EXACT; theorem proofs in RESEARCH.md; no sporadic-group representation claimed',
              universe=dict(vertices=list(range(N)), edges=[list(e) for e in EDGES],
                            labeled_tournaments=1024, isomorphism_classes=12,
                            self_dual_classes=8, converse_quotient_classes=10,
                            characteristic_polynomials=9),
              classes=rows, converse_quotient=[list(o) for o in dual_orbits],
              cospectral_collisions=[[2, 5, 11], [4, 9]],
              regular_family=dict(connection_sets=[list(s) for s in connections], masks=family,
                                  family_group='AGL(1,5)=C5:C4, order20', kernel='C5',
                                  single_tournament_aut='C5', single_tournament_extended='D10 (order10)',
                                  characteristic_desc=[1, 0, 0, -5, -5, -2],
                                  skew_characteristic_desc=[1, 0, 10, 0, 5, 0]),
              checks_passed=CHECKS)
print(json.dumps(output, indent=2, sort_keys=True))
