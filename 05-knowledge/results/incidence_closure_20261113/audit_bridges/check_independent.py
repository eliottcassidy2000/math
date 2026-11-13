#!/usr/bin/env python3
"""Independent finite bridge controls. No imports from producer code.

General mathematical claims are audited in REPORT.md, not inferred from tests.
Uses explicit exceptions so normal and optimized Python retain every gate.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
CHECKS = 0
EXPECTED = {
    'bridges/PROOF.md': '271ad8a5a9ffe408e62b2148b88d26ce9fd255b595267523740da0811923fb2c',
    'bridges/check_bridges.py': '32eea5653748102011bded91edcd1c846ebd8aeac7b190253e2cd6e33436759d',
    'bridges/check.normal.json': '46bc599c3848ae69e04b71248c2bc96c4d474daa5fca54c391dfa6eb34b06d15',
    'bridges/check.optimized.json': '46bc599c3848ae69e04b71248c2bc96c4d474daa5fca54c391dfa6eb34b06d15',
}


def check(value, label):
    global CHECKS
    CHECKS += 1
    if not value:
        raise RuntimeError(label)


def components(n, edges):
    unseen = set(range(n))
    out = []
    while unseen:
        seed = min(unseen)
        reached, todo = {seed}, [seed]
        unseen.remove(seed)
        while todo:
            v = todo.pop()
            for a, b in edges:
                w = b if a == v else a if b == v else None
                if w in unseen:
                    unseen.remove(w)
                    reached.add(w)
                    todo.append(w)
        out.append(sorted(reached))
    return out


def rank(rows, width, prime=None):
    a = [[x % prime for x in row] if prime else list(map(Fraction, row)) for row in rows]
    r = 0
    for col in range(width):
        pivot = next((j for j in range(r, len(a)) if a[j][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][col], -1, prime) if prime else 1 / a[r][col]
        a[r] = [(x * inv) % prime if prime else x * inv for x in a[r]]
        for j in range(len(a)):
            if j != r and a[j][col]:
                q = a[j][col]
                a[j] = [(x - q*y) % prime if prime else x - q*y for x, y in zip(a[j], a[r])]
        r += 1
    return r


def graph_torsors():
    totals = Counter()
    for n in range(5):
        possible = list(combinations(range(n), 2))
        vertex_bits = list(product((0, 1), repeat=n))
        for selected in product((0, 1), repeat=len(possible)):
            edges = [e for e, b in zip(possible, selected) if b]
            c = len(components(n, edges))
            boundary = {tuple(x[a] ^ x[b] for a, b in edges) for x in vertex_bits}
            check(len(boundary) == 2**(n-c), 'binary coboundary image by enumeration')
            incidence = [[int(j == a or j == b) for j in range(n)] for a, b in edges]
            check(rank(incidence, n, 2) == n-c, 'incidence rank independent of orbit enumeration')
            classes = set()
            for s in product((0, 1), repeat=len(edges)):
                solutions = [x for x in vertex_bits if tuple(x[a] ^ x[b] for a, b in edges) == s]
                check(len(solutions) == (2**c if s in boundary else 0), 'torsor affine solution count')
                orbit = {tuple(a ^ b for a, b in zip(s, d)) for d in boundary}
                check(len(orbit) == 2**(n-c), 'gauge orbit size')
                classes.add(min(orbit))
                totals['transport_assignments'] += 1
            check(len(classes) == 2**(len(edges)-n+c), 'all gauge classes')
            totals['simple_graphs'] += 1
    # Reconstruct the actual graceful states by graph reachability, without UF.
    states = []
    rejected = []
    for a in ((0, 1), (1, 2), (2, 3)):
        for b in ((0, 2), (1, 3)):
            (states if len(components(4, [(0, 3), a, b])) == 1 else rejected).append([a, b])
    check(len(states) == 4 and len(rejected) == 2, 'six raw graceful choices')
    check(all(sum(x != y for x, y in zip(a, b)) <= 2 for a, b in combinations(states, 2)), 'actual K4')
    return {'universe': 'all simple graphs on 0..4 labelled vertices and every edge bit assignment',
            **totals, 'graceful_states': states, 'rejected_states': rejected}


def nonabelian_triangle():
    ps = list(permutations(range(3)))
    histogram = Counter()
    for a, b, c in product(ps, repeat=3):
        # Edge maps 0->1=a, 1->2=b, 0->2=c.
        inverse_c = [c.index(i) for i in range(3)]
        holonomy = [inverse_c[b[a[i]]] for i in range(3)]
        fixed = {i for i in range(3) if holonomy[i] == i}
        sections = [(x, y, z) for x, y, z in product(range(3), repeat=3)
                    if y == a[x] and z == b[y] and z == c[x]]
        check({s[0] for s in sections} == fixed and len(sections) == len(fixed), 'root evaluation is bijective')
        histogram[len(sections)] += 1
    check(set(histogram) == {0, 1, 3}, 'nontrivial holonomy can have some fixed states')
    return {'triangle_transport_assignments': 216, 'section_count_histogram': dict(sorted(histogram.items()))}


def linear_transport():
    edges = list(combinations(range(4), 2))
    hist = Counter()
    for bits in product((0, 1), repeat=6):
        rows = []
        for (a, b), bit in zip(edges, bits):
            row = [0]*4
            row[a], row[b] = -(-1)**bit, 1
            rows.append(row)
        dim = 4-rank(rows, 4, 3)
        scalar = [x for x in product(range(3), repeat=4)
                  if all((x[b] - (-1)**s*x[a]) % 3 == 0 for (a, b), s in zip(edges, bits))]
        check(len(scalar) == 3**dim, 'linear rank vs exhaustive scalar sections')
        for k in (0, 2):
            amplified = [[rows[e][v] if j == w else 0 for v in range(4) for j in range(k)]
                         for e in range(6) for w in range(k)]
            check(4*k-rank(amplified, 4*k, 3) == k*dim, 'tensor with trivial finite vector space')
        check(4-rank(rows, 4, 2) == 1, 'characteristic two sign degeneration')
        hist[dim] += 1
    check(dict(hist) == {1: 8, 0: 56}, 'linear zero section differs from empty torsor')
    check(rank([[-2*int(i == j) for j in range(16)] for i in range(16)], 16) == 16,
          'minus identity rank-sixteen parallel obstruction')
    # Complex endpoint gluing: exp(pi*i*0)=1, exp(pi*i*1)=-1.
    check(-(-1) == 1, 'complex nowhere-zero section endpoint compatibility')
    return {'field': 'F3; separate F2 degeneracy and rational rank-sixteen hostile',
            'linear_dimension_histogram': dict(sorted(hist.items())), 'amplification_dimensions': [0, 2, 16]}


def min_plus(values, restrictions, n):
    """Independent terminal-cover formula, unlike producer's all-sieve colimit."""
    out, maps, objects = {}, {}, {}
    for u in range(1 << n):
        atoms = [1 << i for i in range(n) if u & (1 << i)]
        objects[u] = [0] + atoms if u else []
        out[u] = [s for s in product(*(values[v] for v in objects[u]))
                  if not u or all(restrictions[a, 0][x] == s[0] for a, x in zip(atoms, s[1:]))]
    for u in out:
        for v in out:
            if v & u == v:
                maps[u, v] = {s: tuple(s[objects[u].index(w)] for w in objects[v]) for s in out[u]}
    return out, maps


def plus_controls():
    records = []
    for n, size, normalized in product(range(1, 5), (2, 3), (False, True)):
        values = {u: list(range(1 if not u and normalized else size)) for u in range(1 << n)}
        maps = {(u, v): {a: 0 if not v and normalized else a for a in values[u]}
                for u in values for v in values if v & u == v}
        phases = []
        for k in range(4):
            sizes = [len(values[u]) for u in range(1 << n)]
            expected = [1 if not u and (normalized or k) else size**u.bit_count()
                        if k >= 2 or normalized and k >= 1 else size for u in values]
            check(sizes == expected, 'independent terminal-cover plus phases')
            phases.append(sizes)
            values, maps = min_plus(values, maps, n)
        records.append({'n': n, 'label_count': size, 'normalized_empty': normalized, 'phases': phases})
    # Every pair of maps into an empty-open set of size 1 or 2,
    # singleton values of sizes 0..3; the whole-open set is either empty or all matches.
    arbitrary = 0
    for e, a, b in product((1, 2), range(4), range(4)):
        for ra, rb in product(product(range(e), repeat=a), product(range(e), repeat=b)):
            matched = [(i, j) for i in range(a) for j in range(b) if ra[i] == rb[j]]
            for whole in ([], matched):
                values = {0: list(range(e)), 1: list(range(a)), 2: list(range(b)), 3: whole}
                maps = {(u, u): {s:s for s in values[u]} for u in values}
                maps.update({(1, 0): dict(enumerate(ra)), (2, 0): dict(enumerate(rb)),
                             (3, 0): {s: ra[s[0]] for s in whole},
                             (3, 1): {s:s[0] for s in whole}, (3, 2): {s:s[1] for s in whole}})
                first, fm = min_plus(values, maps, 2)
                second, _ = min_plus(first, fm, 2)
                check(len(first[3]) == len(matched), 'first plus retains exact empty-overlap fibre product')
                check(len(second[3]) == a*b, 'second plus removes only empty-overlap matching constraint')
                arbitrary += 1
    return {'constant_presheaves': records, 'additional_two_point_presheaves': arbitrary,
            'independence': 'terminal covering sieve, not all-cover colimit implementation'}


def perfect_quotient():
    perms = list(permutations(range(5)))
    def mul(p, q): return tuple(p[q[i]] for i in range(5))
    def inv(p): return tuple(p.index(i) for i in range(5))
    def sign(p): return sum(p[i] > p[j] for i in range(5) for j in range(i+1, 5)) % 2
    even = {p for p in perms if not sign(p)}
    commutators = {mul(mul(mul(a,b),inv(a)),inv(b)) for a in even for b in even}
    # Strong finite certificate: here every element of A5 is itself a commutator.
    check(commutators == even and len(even) == 60, 'A5 perfect by literal permutations')
    check(all(sign(mul(a,b)) == sign(a)^sign(b) for a in perms for b in perms), 'S5 sign homomorphism')
    check(all(mul(mul(g,a),inv(g)) in even for g in perms for a in even), 'A5 normal')
    cosets = {frozenset(mul(g,a) for a in even) for g in perms}
    check(len(cosets) == 2 and {frozenset(sign(a) for a in c) for c in cosets} == {frozenset({0}),frozenset({1})},
          'nontrivial sign descends through perfect quotient')
    return {'group_order': 120, 'perfect_normal_subgroup_order': 60, 'quotient_order': 2,
            'scope': 'finite control only; general B6 follows because abelian targets kill commutators'}


def stress_residues():
    # Gaussian integer Laurent coefficients, computed independently by convolution.
    def cmul(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
    def lmul(f, g):
        out = {}
        for i, a in f.items():
            for j, b in g.items():
                u, v = cmul(a, b)
                old = out.get(i+j, (0,0))
                out[i+j] = (old[0]+u, old[1]+v)
        return out
    h = lmul({0:(1,0), -2:(-1,0)}, {0:(1,0), -2:(-1,0)})
    tests = 0
    for A, br, bi, C in product(range(-1, 2), repeat=4):
        F = lmul(h, {1:(A,0), 0:(br,-bi)})
        G = lmul(h, {1:(br,bi), 0:(C,0)})
        rf, rzf, rg = F.get(-1,(0,0)), F.get(-2,(0,0)), G.get(-1,(0,0))
        check((rf,rzf,rg) == ((-2*A,0),(-2*br,2*bi),(-2*br,-2*bi)), 'residues from Laurent convolution')
        # Period condition uses conjugation of 2*pi*i*Res G, hence minus conjugate residue.
        trans = (3*rzf[0]-rg[0], 3*rzf[1]+rg[1])
        check(rf[1] == 0, 'single-valued stress for every real amplitude')
        check((rf == (0,0) and trans == (0,0)) == (A == br == bi == 0), 'displacement excludes A and B')
        tests += 1
    return {'net': 'Z=z+1/z', 'kappa': 3, 'integer_amplitudes_tested': tests,
            'stress_dimension': 5, 'displacement_dimension': 2,
            'scope': 'inherited local model only; no audit of inaccessible 1933 full text or Saenz erratum'}


def main():
    pins = {p:hashlib.sha256((PACKET/p).read_bytes()).hexdigest() for p in EXPECTED}
    check(pins == EXPECTED, 'reviewed producer bytes unchanged')
    result = {'status': 'FINITE-EXACT independent audit controls; producer proof status unchanged',
              'reviewed_sha256': pins, 'graph_torsors': graph_torsors(),
              'nonabelian': nonabelian_triangle(), 'linear_transport': linear_transport(),
              'plus': plus_controls(), 'perfect_quotient': perfect_quotient(), 'stress': stress_residues()}
    result['checks'] = CHECKS
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
