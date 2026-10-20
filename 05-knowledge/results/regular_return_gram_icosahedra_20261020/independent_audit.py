"""Independent arc-set/brute-word audit. No imports from producer checker.
Run under Python and Python -O. All validations use explicit exceptions.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def tournament(n, code):
    return frozenset((i, j) if code & (1 << k) else (j, i)
                     for k, (i, j) in enumerate(combinations(range(n), 2)))


def code_of(arcs, n):
    return sum(1 << k for k, ij in enumerate(combinations(range(n), 2)) if ij in arcs)


def path_edges(p):
    return frozenset(zip(p, p[1:]))


def directed(arcs, p):
    return all(uv in arcs for uv in zip(p, p[1:]))


def words(arcs, n):
    return [p for p in permutations(range(n)) if directed(arcs, p)]


def degree(arcs, n):
    return tuple(sum(u == v for u, _ in arcs) for v in range(n))


def flip(arcs, undirected):
    return frozenset((v, u) if frozenset((u, v)) in undirected else (u, v) for u, v in arcs)


def complement_path(arcs, p):
    keep = path_edges(p)
    require(keep <= arcs, 'K called on nondirected path')
    return frozenset((u, v) if (u, v) in keep else (v, u) for u, v in arcs)


def kernel_literal(arcs, x, y, outside):
    # Regard xy as a prescribed supervertex; choose an outside word and insert
    # the block at each cut. Do not test the artificial internal adjacency xy.
    total = 0
    for order in permutations(outside):
        for cut in range(len(order) + 1):
            prefix = order[:cut] + (x,)
            suffix = (y,) + order[cut:]
            total += directed(arcs, prefix) and directed(arcs, suffix)
    return total


def directed_triangles(arcs, n):
    for x, y, z in combinations(range(n), 3):
        if {(x, y), (y, z), (z, x)} <= arcs:
            yield x, y, z
        elif {(x, z), (z, y), (y, x)} <= arcs:
            yield x, z, y


def triangle_response(arcs, n, tri):
    x, y, z = tri
    old = frozenset(((x, y), (y, z), (z, x)))
    new = frozenset((v, u) for u, v in old)
    punctured = arcs - old
    target = punctured | new
    one = two = 0
    outside = tuple(v for v in range(n) if v not in tri)
    for a, b, c in ((x, y, z), (y, z, x), (z, x, y)):
        rest = tuple(v for v in range(n) if v not in (a, b))
        one += kernel_literal(punctured, b, a, rest) - kernel_literal(punctured, a, b, rest)
        two += kernel_literal(punctured, c, a, outside) - kernel_literal(punctured, a, c, outside)
    classes_old = Counter(len(path_edges(p) & old) for p in words(arcs, n))
    classes_new = Counter(len(path_edges(p) & new) for p in words(target, n))
    require(classes_old[0] == classes_new[0], 'zero-edge cancellation')
    require(one == classes_new[1] - classes_old[1], 'linear component')
    require(two == classes_new[2] - classes_old[2], 'quadratic component')
    require(not classes_old[3] and not classes_new[3], 'spurious cubic term')
    return one, two


def formula_audit():
    out = []
    for n in range(1, 6):
        bank = [tournament(n, code) for code in range(1 << (n * (n - 1) // 2))]
        counts = [len(words(arcs, n)) for arcs in bank]
        arc_checks = tri_checks = 0
        for code, arcs in enumerate(bank):
            for a, b in arcs:
                target = (arcs - {(a, b)}) | {(b, a)}
                rest = tuple(v for v in range(n) if v not in (a, b))
                delta = kernel_literal(arcs, b, a, rest) - kernel_literal(arcs, a, b, rest)
                require(delta == counts[code_of(target, n)] - counts[code], 'single-arc identity')
                arc_checks += 1
            for tri in directed_triangles(arcs, n):
                triangle_response(arcs, n, tri)
                tri_checks += 1
        out.append(dict(n=n, tournaments=len(bank), arc_checks=arc_checks,
                        triangle_component_checks=tri_checks))
    return out


def all_compatible_word_geometries():
    out = []
    six_witnesses = {}
    for n in range(1, 11):
        p = tuple(range(n))
        ep = path_edges(p)
        checked = same_count = six_count = rejected = 0
        support_histogram = Counter()
        # Exhaust every endpoint configuration through n7, then every
        # matching-endpoint permutation at n8..10.
        q_bank = permutations(range(n)) if n <= 7 else ((0,) + q + (n-1,) for q in permutations(range(1, n-1)))
        for q in q_bank:
            eq = path_edges(q)
            if any((v, u) in ep for u, v in eq):
                rejected += 1
                continue
            # A compatible pair determines these arcs of a possible T.
            required = ep | frozenset((v, u) for u, v in eq - ep)
            require(not any((v, u) in required for u, v in required), 'inconsistent pair')
            arcs = set(required)
            for u, v in combinations(range(n), 2):
                if (u, v) not in arcs and (v, u) not in arcs:
                    arcs.add((u, v))
            arcs = frozenset(arcs)
            middle = complement_path(arcs, p)
            require(directed(middle, q), 'constructed second path')
            target = complement_path(middle, q)
            d, dmid, dt = degree(arcs, n), degree(middle, n), degree(target, n)
            dp = tuple((v == p[0]) - (v == p[-1]) for v in range(n))
            dq = tuple((v == q[0]) - (v == q[-1]) for v in range(n))
            require(dmid == tuple(n-1-d[v]+dp[v] for v in range(n)), 'first degree law')
            require(dt == tuple(d[v]-dp[v]+dq[v] for v in range(n)), 'second degree law')
            same = (p[0], p[-1]) == (q[0], q[-1])
            require((d == dt) == same, 'same-score iff endpoints')
            checked += 1
            if not same:
                continue
            same_count += 1
            changed = arcs - target
            size = len(changed)
            support_histogram[size] += 1
            require(size == 2 * len(ep - eq), 'support successor count')
            require(size % 2 == 0 and (size == 0 or size >= 6), 'parity/minimum support')
            require(all(sum(u == v for u, _ in changed) == sum(w == v for _, w in changed)
                        for v in range(n)), 'balanced support')
            if size == 6:
                cuts = tuple(i + 1 for i, edge in enumerate(zip(p, p[1:])) if edge not in eq)
                require(len(cuts) == 3, 'three cuts')
                i, j, k = cuts
                A, B, C, D = p[:i], p[i:j], p[j:k], p[k:]
                require(A and B and C and D and q == A + C + B + D, 'block order')
                require(not (len(B) == len(C) == 1), 'singleton reversal forbidden')
                used = set(v for uv in changed for v in uv)
                kind = 'directed_six_cycle' if len(B) > 1 and len(C) > 1 else 'two_triangles_shared_vertex'
                require(len(used) == (6 if kind == 'directed_six_cycle' else 5), 'support vertex count')
                if kind == 'directed_six_cycle':
                    require(all(sum(u == v for u, _ in changed) == 1 for v in used), 'cycle indegrees')
                else:
                    triangles = list(directed_triangles(changed, n))
                    require(len(triangles) == 2 and len(set(triangles[0]) & set(triangles[1])) == 1,
                            'two shared-vertex triangles')
                six_witnesses.setdefault(kind, dict(n=n, P=p, Q=q, changed_arcs=sorted(changed)))
                six_count += 1
        out.append(dict(n=n, q_universe='all permutations' if n <= 7 else 'all matching-endpoint permutations',
                        compatible_pairs=checked, rejected_reverse_shared_edges=rejected,
                        matching_endpoint_pairs=same_count, six_arc_pairs=six_count,
                        support_histogram=dict(sorted(support_histogram.items()))))
    require(set(six_witnesses) == {'directed_six_cycle', 'two_triangles_shared_vertex'}, 'both equality shapes')
    return dict(rows=out, equality_shape_witnesses=six_witnesses)


def regular5_graph():
    n = 5
    bank = [tournament(n, c) for c in range(1024) if degree(tournament(n, c), n) == (2,)*n]
    require(len(bank) == 24, 'regular5 count')
    weight = {}
    for arcs in bank:
        code = code_of(arcs, n)
        hp = words(arcs, n)
        w = Counter()
        for p in hp:
            middle = complement_path(arcs, p)
            for q in words(middle, n):
                target = complement_path(middle, q)
                same = (p[0], p[-1]) == (q[0], q[-1])
                require((degree(target, n) == (2,)*n) == same, 'regular5 iff')
                if same:
                    w[code_of(target, n)] += 1
        require(len(hp) == w[code] == 15 and sum(w.values()) == 25, 'regular5 masses')
        require(sorted(v for t, v in w.items() if t != code) == [2]*5, 'regular5 edge weights')
        f = []
        for v in range(n):
            x, y = [b for a, b in arcs if a == v]
            f.append(x if (x, y) in arcs else y)
        require(len(set(f)) == 5, 'f injective')
        require(arcs == frozenset((v, w) for v in range(n) for w in (f[v], f[f[v]])), 'f reconstruction')
        autos = [pi for pi in permutations(range(n)) if frozenset((pi[u], pi[v]) for u, v in arcs) == arcs]
        require(len(autos) == 5, 'regular5 automorphism count')
        weight[code] = w
    graph = {c: set(w)-{c} for c, w in weight.items()}
    for c in graph:
        for t in graph:
            require((t in graph[c]) == ((c ^ t).bit_count() == 6), 'six disagreement iff adjacency')
            require(weight[c][t] == weight[t][c], 'symmetric weights')
    north, south = 947, 76
    upper = [613, 217, 710, 341, 668]
    lower = [313, 682, 355, 410, 806]
    expected = {north: set(upper), south: set(lower)}
    for i in range(5):
        expected[upper[i]] = {north, upper[(i-1)%5], upper[(i+1)%5], lower[i], lower[(i+1)%5]}
        expected[lower[i]] = {south, lower[(i-1)%5], lower[(i+1)%5], upper[i], upper[(i-1)%5]}
    require(all(graph[c] == ns for c, ns in expected.items()), 'icosahedral certificate')
    components = []
    unseen = set(graph)
    while unseen:
        seed = min(unseen)
        queue = [seed]
        dist = {seed: 0}
        for c in queue:
            for t in graph[c]:
                if t not in dist:
                    dist[t] = dist[c]+1
                    queue.append(t)
        unseen -= set(dist)
        require(len(dist) == 12, 'sheet size')
        require(len({c.bit_count()%2 for c in dist}) == 1, 'sheet parity')
        require(sorted(Counter(dist.values()).items()) == [(0,1),(1,5),(2,5),(3,1)], 'icosahedral shells')
        require(dist[seed ^ 1023] == 3, 'converse antipode')
        components.append(dict(representative=seed, size=12, distance_shells=dict(Counter(dist.values()))))
    require(len(components) == 2, 'two sheets')
    return dict(states=24, components=components, H=15, aut=5, total_presentation_mass=25,
                neighbors_per_vertex=5, weight_per_neighbor=2, nonloop_edges_total=60)


def regular7_response_controls():
    out = []
    rows_bank = [(14,28,56,112,97,67,7), (26,28,41,112,100,67,7),
                 tuple(sum(1 << ((v+d)%7) for d in (1,2,4)) for v in range(7))]
    for rows in rows_bank:
        arcs = frozenset((u, v) for u in range(7) for v in range(7) if rows[u] & (1 << v))
        require(degree(arcs, 7) == (3,)*7, 'regular7 input')
        splits = Counter(triangle_response(arcs, 7, tri) for tri in directed_triangles(arcs, 7))
        out.append(dict(rows=rows, H=len(words(arcs, 7)),
                        triangle_splits=[dict(linear=a, quadratic=b, multiplicity=c) for (a,b),c in sorted(splits.items())]))
    return out


def regular_support_witnesses():
    out = []
    cases = [
        (5, (0,2,3,4,1), (0,3,4,2,1), 'two_triangles_shared_vertex'),
        (7, (0,1,3,4,5,6,2), (0,1,5,6,3,4,2), 'directed_six_cycle'),
        (7, (0,1,3,4,5,6,2), (0,1,6,3,4,5,2), 'two_triangles_shared_vertex'),
    ]
    for n, p, q, kind in cases:
        arcs = frozenset((u, (u+d)%n) for u in range(n) for d in range(1, (n+1)//2))
        target = complement_path(complement_path(arcs, p), q)
        changed = arcs - target
        require(degree(target, n) == ((n-1)//2,)*n, 'equality witness regularity')
        require(len(changed) == 6, 'equality witness support')
        used = set(v for uv in changed for v in uv)
        require(len(used) == (6 if kind == 'directed_six_cycle' else 5), 'equality witness shape')
        out.append(dict(n=n, source_code=code_of(arcs, n), target_code=code_of(target, n),
                        P=p, Q=q, shape=kind, changed_arcs=sorted(changed)))
    return out


def main():
    result = dict(status='PASS', formulas=formula_audit(),
                  path_pair_geometry=all_compatible_word_geometries(),
                  regular5=regular5_graph(), regular7=regular7_response_controls(),
                  regular_support_witnesses=regular_support_witnesses())
    result['independent_source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
