"""Independent small checks for the endpoint/regular-return proofs.
No assert statements: normal and -O execute the same checks.
"""
from collections import Counter
from itertools import combinations, permutations
import hashlib
import json
from pathlib import Path


def need(x, message):
    if not x:
        raise RuntimeError(message)


def rows(n, code):
    a = [0] * n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if code >> k & 1 else (j, i)
        a[u] |= 1 << v
    return tuple(a)


def encode(a):
    return sum((a[i] >> j & 1) << k
               for k, (i, j) in enumerate(combinations(range(len(a)), 2)))


def literal(a):
    return sum(all(a[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(a))))


def paths(a, start=None, end=None):
    n, full = len(a), (1 << len(a)) - 1

    def go(p, mask):
        if mask == full:
            if end is None or p[-1] == end:
                yield tuple(p)
            return
        todo = a[p[-1]] & ~mask
        if end is not None and mask != (full ^ (1 << end)):
            todo &= ~(1 << end)
        while todo:
            bit = todo & -todo
            todo -= bit
            yield from go(p + [bit.bit_length() - 1], mask | bit)

    for s in range(n) if start is None else (start,):
        yield from go([s], 1 << s)


def reverse(a, edges):
    b = list(a)
    for u, v in edges:
        need(a[u] >> v & 1, 'reverse missing arc')
        b[u] ^= 1 << v
        b[v] ^= 1 << u
    return tuple(b)


def complement_path(a, p):
    n, full = len(a), (1 << len(a)) - 1
    b = [full ^ (1 << i) ^ a[i] for i in range(n)]
    for u, v in zip(p, p[1:]):
        need(a[u] >> v & 1, 'complement path must be directed')
        b[u] ^= 1 << v
        b[v] ^= 1 << u
    return tuple(b)


def endpoints(a):
    n = len(a)
    end = [[0] * n for _ in range(1 << n)]
    start = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        end[1 << i][i] = start[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            todo = a[v] & ~mask if end[mask][v] else 0
            while todo:
                bit = todo & -todo
                todo -= bit
                end[mask | bit][bit.bit_length() - 1] += end[mask][v]
            todo = sum(1 << u for u in range(n) if a[u] >> v & 1) & ~mask if start[mask][v] else 0
            while todo:
                bit = todo & -todo
                todo -= bit
                start[mask | bit][bit.bit_length() - 1] += start[mask][v]
    return end, start


def kernel(tables, x, y, w):
    end, start = tables
    total, a = 0, w
    while True:
        total += end[a | (1 << x)][x] * start[(w ^ a) | (1 << y)][y]
        if not a:
            return total
        a = (a - 1) & w


def triangles(a):
    for x, y, z in combinations(range(len(a)), 3):
        if a[x] >> y & 1 and a[y] >> z & 1 and a[z] >> x & 1:
            yield x, y, z
        elif a[y] >> x & 1 and a[x] >> z & 1 and a[z] >> y & 1:
            yield x, z, y


def triangle_delta(a, triangle):
    x, y, z = triangle
    d = list(a)
    for u, v in ((x, y), (y, z), (z, x)):
        d[u] &= ~(1 << v)
    tables = endpoints(tuple(d))
    full = (1 << len(a)) - 1
    w = full ^ (1 << x) ^ (1 << y) ^ (1 << z)
    linear = quadratic = 0
    for u, v, other in ((x, y, z), (y, z, x), (z, x, y)):
        outside = full ^ (1 << u) ^ (1 << v)
        linear += kernel(tables, v, u, outside) - kernel(tables, u, v, outside)
        quadratic += kernel(tables, other, u, w) - kernel(tables, u, other, w)
    return linear, quadratic


def graph_components(graph):
    seen, out = set(), []
    for v in graph:
        if v in seen:
            continue
        seen.add(v)
        stack, component = [v], []
        while stack:
            u = stack.pop()
            component.append(u)
            for w in graph[u] - seen:
                seen.add(w)
                stack.append(w)
        out.append(sorted(component))
    return out


def structural_f(a):
    f = []
    for v, mask in enumerate(a):
        out = [x for x in range(5) if mask >> x & 1]
        need(len(out) == 2, 'regular5 degree')
        x, y = out
        f.append(x if a[x] >> y & 1 else y)
    need(len(set(f)) == 5, 'intrinsic map permutation')
    for v in range(5):
        need(a[v] == (1 << f[v]) | (1 << f[f[v]]), 'intrinsic cycle reconstruction')
    return f


def main():
    result = {'endpoint_controls': [], 'return_graphs': [], 'regular7_triangles': []}
    for n in (3, 4, 5):
        bank = [rows(n, c) for c in range(1 << (n * (n - 1) // 2))]
        counts = [literal(a) for a in bank]
        arcs = triangle_count = 0
        changes = Counter()
        for code, a in enumerate(bank):
            tables = endpoints(a)
            need(sum(tables[0][-1]) == counts[code], 'endpoint DP vs permutations')
            full = (1 << n) - 1
            for x in range(n):
                for y in range(n):
                    if not a[x] >> y & 1:
                        continue
                    b = reverse(a, ((x, y),))
                    w = full ^ (1 << x) ^ (1 << y)
                    delta = kernel(tables, y, x, w) - kernel(tables, x, y, w)
                    need(delta == counts[encode(b)] - counts[code], 'single arc formula')
                    arcs += 1
            for triple in triangles(a):
                x, y, z = triple
                b = reverse(a, ((x, y), (y, z), (z, x)))
                linear, quadratic = triangle_delta(a, triple)
                delta = counts[encode(b)] - counts[code]
                need(delta == linear + quadratic, 'punctured triangle formula')
                need([r.bit_count() for r in a] == [r.bit_count() for r in b], 'triangle scores')
                triangle_count += 1
                changes[delta] += 1
        result['endpoint_controls'].append({'n': n, 'labelled_tournaments': len(bank),
            'single_arc_checks': arcs, 'triangle_checks': triangle_count,
            'triangle_H_increments': dict(sorted(changes.items()))})

    for n in (3, 5):
        r = (n - 1) // 2
        bank = [rows(n, c) for c in range(1 << (n * (n - 1) // 2))]
        bank = [a for a in bank if all(row.bit_count() == r for row in a)]
        graph, weights, support = {}, {}, Counter()
        for a in bank:
            code = encode(a)
            graph[code], weights[code] = set(), Counter()
            if n == 5:
                structural_f(a)
            for p in paths(a):
                b = complement_path(a, p)
                expected = [r + (v == p[0]) - (v == p[-1]) for v in range(n)]
                need([row.bit_count() for row in b] == expected, 'first degree defect')
                # Enumerate *all* Q and test iff, not only endpoint-filtered Q.
                for q in paths(b):
                    c = complement_path(b, q)
                    reg = all(row.bit_count() == r for row in c)
                    same = (p[0], p[-1]) == (q[0], q[-1])
                    need(reg == same, 'two-step return iff endpoints agree')
                    if not reg:
                        continue
                    target = encode(c)
                    weights[code][target] += 1
                    size = (code ^ target).bit_count()
                    need(size % 2 == 0, 'return parity')
                    need(size == 0 or size >= 6, 'minimum support six')
                    support[size] += 1
                    if code != target:
                        graph[code].add(target)
            need(weights[code][code] == literal(a), 'identity presentation count')
        for a in graph:
            for b in graph[a]:
                need(weights[a][b] == weights[b][a], 'symmetric edge presentation mass')
        components = graph_components(graph)
        need(sorted(map(len, components)) == ([1, 1] if n == 3 else [12, 12]), 'components')
        need(all(len({x.bit_count() % 2 for x in component}) == 1 for component in components), 'component signs')
        if n == 5:
            need(all(len(v) == 5 for v in graph.values()), 'degree five')
            need(all(weights[a][b] == 2 for a in graph for b in graph[a]), 'weight two')
            need(all(sum(w.values()) == 25 for w in weights.values()), 'mass25')
            north, south = 947, 76
            upper = [613, 217, 710, 341, 668]
            lower = [313, 682, 355, 410, 806]
            certificate = {north: set(upper), south: set(lower)}
            for i in range(5):
                certificate[upper[i]] = {north, upper[(i-1)%5], upper[(i+1)%5], lower[i], lower[(i+1)%5]}
                certificate[lower[i]] = {south, lower[(i-1)%5], lower[(i+1)%5], upper[i], upper[(i-1)%5]}
            need(all(graph[v] == neighbors for v, neighbors in certificate.items()), 'icosahedron exact adjacency')
            a = rows(5, north)
            orbit_table = []
            for p in paths(a, start=0):
                q = list(paths(complement_path(a, p), p[0], p[-1]))
                orbit_table.append({'P': p, 'Q': q})
            result['icosahedron_certificate'] = {'north': north, 'south': south,
                'upper': upper, 'lower': lower, 'edge_count': sum(map(len, certificate.values()))//2,
                'path_orbits': orbit_table}
        result['return_graphs'].append({'n': n, 'states': len(bank),
            'component_sizes': sorted(map(len, components)),
            'degrees': dict(sorted(Counter(map(len, graph.values())).items())),
            'presentation_support': dict(sorted(support.items()))})

    examples = [(14, 28, 56, 112, 97, 67, 7),
                (26, 28, 41, 112, 100, 67, 7),
                tuple(sum(1 << ((i+d)%7) for d in (1, 2, 4)) for i in range(7))]
    for a in examples:
        h = literal(a)
        changes, splits = Counter(), []
        for triple in triangles(a):
            x, y, z = triple
            b = reverse(a, ((x, y), (y, z), (z, x)))
            delta = literal(b) - h
            linear, quadratic = triangle_delta(a, triple)
            need(delta == linear + quadratic, 'regular7 triangle response')
            changes[delta] += 1
            splits.append({'triangle': triple, 'linear': linear, 'quadratic': quadratic, 'delta': delta})
        result['regular7_triangles'].append({'rows': a, 'H': h,
            'changes': dict(sorted(changes.items())), 'endpoint_splits': splits})
    need(sorted(x['H'] for x in result['regular7_triangles']) == [171, 175, 189], 'regular7 controls')
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['status'] = 'PASS'
    output = json.dumps(result, sort_keys=True, indent=2) + '\n'
    Path(__file__).with_suffix('.json').write_text(output)
    print(output, end='')


if __name__ == '__main__':
    main()
