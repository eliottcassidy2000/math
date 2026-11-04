#!/usr/bin/env python3
"""Independent, standard-library exact verifier. Imports no repository producer."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import gcd
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


MATCHINGS = [((0, 1), (2, 3), (4, 5)), ((0, 1), (2, 5), (3, 4)),
             ((0, 3), (1, 2), (4, 5)), ((0, 5), (1, 2), (3, 4)),
             ((0, 5), (1, 4), (2, 3))]


def connected(vertices, edges):
    adj = {v: [] for v in vertices}
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    reached, stack = set(), [min(vertices)]
    while stack:
        v = stack.pop()
        if v not in reached:
            reached.add(v)
            stack.extend(adj[v])
    return reached == set(vertices)


def determinant(a):
    a = [[Fraction(x) for x in row] for row in a]
    ans = Fraction(1)
    for j in range(len(a)):
        p = next((i for i in range(j, len(a)) if a[i][j]), None)
        if p is None:
            return 0
        if p != j:
            a[p], a[j] = a[j], a[p]
            ans = -ans
        q = a[j][j]
        ans *= q
        for i in range(j + 1, len(a)):
            f = a[i][j] / q
            for k in range(j, len(a)):
                a[i][k] -= f * a[j][k]
    need(ans.denominator == 1, 'nonintegral determinant')
    return int(ans)


def smith_invariants(a):
    n = len(a)
    previous, invariants = 1, []
    for k in range(1, n + 1):
        divisor = 0
        for rows in combinations(range(n), k):
            for cols in combinations(range(n), k):
                divisor = gcd(divisor, abs(determinant([[a[i][j] for j in cols] for i in rows])))
        need(divisor % previous == 0, 'determinantal divisor')
        invariants.append(divisor // previous)
        previous = divisor
    return invariants


def patch(hub, offset=6):
    # Ports 0..5; forced subdivisions offset..offset+5; rim offset+6..offset+11.
    edges = []
    for i in range(6):
        edges.extend([(i, offset + i), (offset + i, offset + 6 + i)])
    edges.extend((offset + 6 + i, offset + 6 + (i + 1) % 6) for i in range(6))
    if hub:
        edges.extend((offset + 12, offset + 6 + i) for i in range(6))
    vertices = set(range(6)) | set(range(offset, offset + 12 + int(hub)))
    return vertices, edges


def degree_subsets(vertices, edges, boundary):
    """Enumerate edge subsets, requiring degree 2 off boundary, degree <=2 on it."""
    edges = [tuple(sorted(e)) for e in edges]
    need(len(edges) == len(set(edges)), 'duplicate graph edge')
    deg, remaining = Counter(), Counter(v for edge in edges for v in edge)
    chosen = []
    exact = set(vertices) - set(boundary)
    def viable(v):
        return deg[v] <= 2 and (v not in exact or deg[v] + remaining[v] >= 2)
    def recurse(i):
        if i == len(edges):
            if all(deg[v] == 2 for v in exact):
                yield tuple(chosen)
            return
        a, b = edges[i]
        remaining[a] -= 1
        remaining[b] -= 1
        deg[a] += 1
        deg[b] += 1
        if viable(a) and viable(b):
            chosen.append((a, b))
            yield from recurse(i + 1)
            chosen.pop()
        deg[a] -= 1
        deg[b] -= 1
        if viable(a) and viable(b):
            yield from recurse(i + 1)
        remaining[a] += 1
        remaining[b] += 1
    if all(remaining[v] >= 2 for v in exact):
        yield from recurse(0)


def profile(vertices, edges):
    out = Counter()
    for selected in degree_subsets(vertices, edges, range(6)):
        adj = {v: [] for v in vertices}
        for a, b in selected:
            adj[a].append(b)
            adj[b].append(a)
        unseen = {v for v in vertices if adj[v]}
        pairs, valid = [], True
        while unseen:
            component, stack = set(), [min(unseen)]
            while stack:
                v = stack.pop()
                if v not in component:
                    component.add(v)
                    stack.extend(adj[v])
            unseen -= component
            ends = tuple(sorted(v for v in component if len(adj[v]) == 1))
            if len(ends) != 2 or any(v >= 6 for v in ends):
                valid = False
                break
            pairs.append(ends)
        if valid and pairs:
            out[(tuple(len(adj[v]) for v in range(6)), tuple(sorted(pairs)))] += 1
    return out


def cycle_witnesses(vertices, edges):
    witnesses = []
    for chosen in degree_subsets(vertices, edges, ()):
        if not connected(vertices, chosen):
            continue
        adj = {v: [] for v in vertices}
        for a, b in chosen:
            adj[a].append(b)
            adj[b].append(a)
        path = [min(vertices)]
        prev, cur = path[0], min(adj[path[0]])
        while cur != path[0]:
            path.append(cur)
            nxt = next(v for v in adj[cur] if v != prev)
            prev, cur = cur, nxt
        need(len(path) == len(vertices), 'cycle does not span')
        witnesses.append(path + path[:1])
    return sorted(witnesses)


def noncrossing_matchings(xs):
    if not xs:
        yield ()
        return
    a = xs[0]
    for j in range(1, len(xs), 2):
        for inner in noncrossing_matchings(xs[1:j]):
            for outer in noncrossing_matchings(xs[j + 1:]):
                yield tuple(sorted(((a, xs[j]),) + inner + outer))


def unrestricted_matchings(xs):
    if not xs:
        yield ()
        return
    a = xs[0]
    for j in range(1, len(xs)):
        rest = xs[1:j] + xs[j+1:]
        for matching in unrestricted_matchings(rest):
            yield tuple(sorted(((a, xs[j]),) + matching))


def states():
    for d in product(range(3), repeat=6):
        ones = tuple(i for i in range(6) if d[i] == 1)
        if ones and len(ones) % 2 == 0:
            for m in noncrossing_matchings(ones):
                yield d, m


def unit_probe(state, offset=100):
    d, matching = state
    ones = [i for i in range(6) if d[i] == 1]
    tips, raw = {}, []
    for j, u in enumerate(ones):
        chain = [u]
        v = (u + 1) % 6
        while v != ones[(j + 1) % len(ones)]:
            if d[v] == 2:
                chain.append(v)
            v = (v + 1) % 6
        raw.extend(zip(chain, chain[1:]))
        tips[u] = chain[-1]
    raw.extend((tips[u], tips[v]) for u, v in matching)
    edges = []
    for i, (u, v) in enumerate(raw):
        edges.extend([(u, offset + i), (offset + i, v)])
    return set(range(6)) | set(range(offset, offset + len(raw))), edges


def gluing(a, b):
    av, ae = a
    bv, be = b
    need(av & bv == set(range(6)), 'wrong vertex overlap')
    need(not (set(map(frozenset, ae)) & set(map(frozenset, be))), 'edge overlap')
    return cycle_witnesses(av | bv, ae + be)


def cross(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def on_segment(a, b, c):
    return cross(a, b, c) == 0 and all(min(a[i], b[i]) <= c[i] <= max(a[i], b[i]) for i in range(2))


def intersects(a, b, c, d):
    if any([on_segment(a, b, c), on_segment(a, b, d), on_segment(c, d, a), on_segment(c, d, b)]):
        return True
    return (cross(a, b, c) * cross(a, b, d) < 0 and cross(c, d, a) * cross(c, d, b) < 0)


def embedding(hub):
    rim = [(2, 0), (1, 2), (-1, 2), (-2, 0), (-1, -2), (1, -2)]
    xy = {i: (3*x, 3*y) for i, (x, y) in enumerate(rim)}
    xy.update({6+i: (2*x, 2*y) for i, (x, y) in enumerate(rim)})
    xy.update({12+i: q for i, q in enumerate(rim)})
    if hub:
        xy[18] = (0, 0)
    vertices, edges = patch(hub)
    need(len(vertices) == len(set(xy.values())), 'coordinate collision')
    for e, f in combinations(edges, 2):
        if not set(e) & set(f):
            need(not intersects(xy[e[0]], xy[e[1]], xy[f[0]], xy[f[1]]), 'drawing crossing')
    for v in vertices - set(range(6)):
        need(all(cross(xy[i], xy[(i+1)%6], xy[v]) > 0 for i in range(6)), 'interior vertex outside boundary')
    return {str(v): xy[v] for v in sorted(xy)}


def main():
    K = [[int(connected(range(6), m+n)) for n in MATCHINGS] for m in MATCHINGS]
    need(determinant(K) == -2, 'wrong connectivity determinant')
    need(smith_invariants(K) == [1, 1, 1, 1, 2], 'wrong connectivity Smith invariants')
    a, b = patch(True), patch(False)
    pa, pb = profile(*a), profile(*b)
    wanted_a = Counter({((1,)*6, m): c for m, c in zip(MATCHINGS, [3, 1, 1, 3, 1])})
    wanted_b = Counter({((1,)*6, MATCHINGS[i]): 1 for i in (0, 3)})
    need(pa == wanted_a and pb == wanted_b, 'wrong FULL patch profiles')
    physical = []
    allstates = list(states())
    need(len(allstates) == 365, 'state count')
    for state in allstates:
        exterior = unit_probe(state)
        ca, cb = gluing(a, exterior), gluing(b, exterior)
        if state[0] == (1,)*6:
            i = MATCHINGS.index(state[1])
            need(len(ca) == sum(K[i][j]*[3,1,1,3,1][j] for j in range(5)), 'wrong hub gluing')
            need(len(cb) == sum(K[i][j]*[1,0,0,1,0][j] for j in range(5)), 'wrong ring gluing')
        else:
            need(not ca and not cb, 'unexpected other-sector cycle')
        need((len(ca)-len(cb)) == 2*int(state[0] == (1,)*6), 'physical extraction identity')
        physical.append({'degrees': state[0], 'matching': state[1], 'hub_count': len(ca), 'ring_count': len(cb)})
    unrestricted = []
    crossing_cases = 0
    for matching in unrestricted_matchings(tuple(range(6))):
        response_difference = sum(int(connected(range(6), matching + m))*c
                                  for m,c in zip(MATCHINGS,[2,1,1,2,1]))
        need(response_difference in (2,4,6), 'unrestricted parity collision')
        item = {'matching':matching, 'count_difference':response_difference}
        if matching not in MATCHINGS:
            ext = unit_probe(((1,)*6,matching))
            ca, cb = gluing(a,ext), gluing(b,ext)
            need(len(ca)-len(cb)==response_difference, 'crossing physical gluing')
            item.update({'hub_count':len(ca),'ring_count':len(cb)})
            crossing_cases += 1
        unrestricted.append(item)
    need(len(unrestricted)==15 and crossing_cases==10, 'unrestricted matching count')
    extra = []
    for hub in (False, True):
        ext = patch(hub, 100)
        ca, cb = gluing(a, ext), gluing(b, ext)
        need(len(ca)-len(cb) == 2*(9 if hub else 2), 'patch-on-patch identity')
        extra.append({'exterior': 'hub' if hub else 'ring', 'hub_count':len(ca), 'ring_count':len(cb)})
    # Exact recovery from all response residues modulo four, and failure modulo two.
    mod2_fibres = {}
    for x in product(range(2), repeat=5):
        y = tuple(sum(K[i][j]*x[j] for j in range(5)) % 2 for i in range(5))
        mod2_fibres.setdefault(y, []).append(x)
    need(len(mod2_fibres) == 16 and all(len(v) == 2 for v in mod2_fibres.values()), 'parity fibre sizes')
    for fibre in mod2_fibres.values():
        need(len({sum(x) % 2 for x in fibre}) == 2, 'extra observable fails')
    residues4 = {}
    for x in product(range(4), repeat=5):
        y = tuple(sum(K[i][j]*x[j] for j in range(5)) % 4 for i in range(5))
        residues4.setdefault(y, set()).add(tuple(v % 2 for v in x))
    need(all(len(v) == 1 for v in residues4.values()), 'mod4 recovery failure')
    # Independently derive order-three invariant Walsh evaluations from labelled signs.
    # P3 is sum x01*x02 - x01*x12 + x02*x12, corresponding to root-out orientations.
    bytype = {}
    for signs in product((-1,1), repeat=3):
        scores = [0,0,0]
        for (i,j), s in zip(((0,1),(0,2),(1,2)), signs):
            scores[i if s == 1 else j] += 1
        typ = 'cyclic' if scores == [1,1,1] else 'transitive'
        x,y,z = signs
        phi = x*y-x*z+y*z
        bytype.setdefault(typ,set()).add(phi)
    need(bytype == {'transitive': {1}, 'cyclic': {-3}}, 'Walsh derivation')
    W = [[1,1],[1,-3]]
    need(smith_invariants(W) == [1,4], 'Walsh Smith invariants')
    for modulus in (4,8):
        fibres = {}
        for x in product(range(modulus), repeat=2):
            y = tuple(sum(W[i][j]*x[j] for j in range(2)) % modulus for i in range(2))
            fibres.setdefault(y,set()).add(tuple(v%2 for v in x))
        need(all(len(f)==1 for f in fibres.values()) == (modulus==8), 'Walsh precision threshold')
    report = {'status':'FINITE-EXACT; standard library; no producer imports',
              'matchings':MATCHINGS, 'K':K, 'K_determinant':determinant(K), 'K_smith':[1,1,1,1,2],
              'patches':{}, 'all_unit_gluing_cases':physical, 'additional_gluing_cases':extra,
              'physical_gluing_cases':len(physical)+len(extra)+crossing_cases,
              'unrestricted_matching_responses':unrestricted,
              'mod2_response_fibres':len(mod2_fibres), 'mod2_fibre_size':2,
              'mod4_domain_vectors_checked':4**5, 'Walsh_order3':W, 'Walsh_order3_smith':[1,4],
              'existence_hostile':{'exterior_matching':MATCHINGS[1], 'hub_cycles':gluing(a,unit_probe(((1,)*6,MATCHINGS[1]))), 'ring_cycles':[]}}
    for name, h, graph, p in [('hub',True,a,pa),('ring',False,b,pb)]:
        report['patches'][name] = {'vertices':sorted(graph[0]), 'edges':graph[1], 'coordinates':embedding(h),
            'full_profile':[{'degrees':s[0],'matching':s[1],'count':c} for s,c in sorted(p.items())]}
    print(json.dumps(report, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
