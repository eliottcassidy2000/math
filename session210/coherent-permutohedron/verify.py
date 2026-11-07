#!/usr/bin/env python3
"""Self-contained exact checks; no assertion statements and no external data."""
import itertools as it
import json
import math
from collections import Counter
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def tournament(n, code):
    a = [[False] * n for _ in range(n)]
    for k, (i, j) in enumerate(it.combinations(range(n), 2)):
        a[i][j] = bool(code >> k & 1)
        a[j][i] = not a[i][j]
    return a


def energy(a, p):
    return sum(not a[p[i]][p[j]] for i in range(len(p)) for j in range(i + 1, len(p)))


def subset_counts(a):
    n = len(a)
    end = [[0] * n for _ in range(1 << n)]
    h = [0] * (1 << n)
    h[0] = 1
    for s in range(1, 1 << n):
        for v in range(n):
            if s >> v & 1:
                rest = s ^ (1 << v)
                end[s][v] = 1 if not rest else sum(end[rest][u] for u in range(n) if rest >> u & 1 and a[u][v])
        h[s] = sum(end[s])
    return h


def ordered_partitions(s):
    if not s:
        yield ()
        return
    b = s
    while b:
        for rest in ordered_partitions(s ^ b):
            yield (b,) + rest
        b = (b - 1) & s


def full_census():
    result = []
    first_local_nonglobal = None
    for n in range(1, 6):
        perms = list(it.permutations(range(n)))
        faces = list(ordered_partitions((1 << n) - 1))
        spectrum = set()
        tested_edges = 0
        for c in range(1 << math.comb(n, 2)):
            a = tournament(n, c)
            hp = subset_counts(a)
            es = {p: energy(a, p) for p in perms}
            minima = []
            rhs = [0] * n
            for p in perms:
                forward = sum(a[p[i]][p[i + 1]] for i in range(n - 1))
                sink = True
                for i in range(n - 1):
                    q = p[:i] + (p[i + 1], p[i]) + p[i + 2:]
                    delta = es[q] - es[p]
                    require(delta == (1 if a[p[i]][p[i + 1]] else -1), 'gradient edge')
                    sink = sink and delta > 0
                    tested_edges += 1
                if sink:
                    minima.append(p)
                for d in range(forward + 1):
                    rhs[d] += math.comb(forward, d)
            require(len(minima) == hp[-1], 'sink/endpoint count')
            require(hp[-1] % 2 == 1, 'finite odd count control')
            lhs = [0] * n
            for face in faces:
                lhs[n - len(face)] += math.prod(hp[b] for b in face)
            require(lhs == rhs, 'face incidence polynomial')
            require(sum((-1) ** d * lhs[d] for d in range(n)) == hp[-1], 'Euler sink identity')
            spectrum.add(hp[-1])
            sink_es = Counter(es[p] for p in minima)
            if first_local_nonglobal is None and len(sink_es) > 1:
                first_local_nonglobal = {'n': n, 'code': c, 'H': hp[-1], 'sink_energy_histogram': dict(sorted(sink_es.items())), 'global_energy': min(es.values()), 'paths': [{'path': p, 'energy': es[p]} for p in minima]}
        result.append({'n': n, 'labelled_tournaments': 1 << math.comb(n, 2), 'ordered_faces': len(faces), 'directed_edge_checks': tested_edges, 'spectrum': sorted(spectrum)})
    return result, first_local_nonglobal


def one_flip_a(n):
    a = [[i < j for j in range(n)] for i in range(n)]
    a[0][n - 1], a[n - 1][0] = False, True
    return a


def qbinomial_rows(m):
    # Polynomial coefficient lists, recurrence G(m,r)=G(m-1,r)+q^(m-r)G(m-1,r-1).
    rows = [[1]]
    for k in range(1, m + 1):
        new = []
        for r in range(k + 1):
            a = rows[r] if r < k else []
            b = [0] * (k - r) + rows[r - 1] if r else []
            z = [0] * max(len(a), len(b))
            for i, v in enumerate(a): z[i] += v
            for i, v in enumerate(b): z[i] += v
            new.append(z)
        rows = new
    return rows


def evalpoly(p, q):
    return sum(c * q ** i for i, c in enumerate(p))


def flip_family():
    result = []
    for n in range(3, 13):
        m = n - 2
        a = one_flip_a(n)
        sinkhist = Counter({1: 1})
        rankpolys = [Counter() for _ in range(m + 1)]
        paths = []
        for bits in range(1 << m):
            A = tuple(i + 1 for i in range(m) if bits >> i & 1)
            B = tuple(i + 1 for i in range(m) if not bits >> i & 1)
            p = A + (n - 1, 0) + B
            inv = sum(x > y for x in A for y in B)
            require(all(a[p[i]][p[i + 1]] for i in range(n - 1)), 'one-flip path')
            require(energy(a, p) == m + inv, 'one-flip energy')
            sinkhist[m + inv] += 1
            rankpolys[len(A)][inv] += 1
            paths.append(p)
        polys = qbinomial_rows(m)
        require(len(set(paths)) == 2 ** m, 'one-flip injectivity')
        for r in range(m + 1):
            require(dict(rankpolys[r]) == {i: v for i, v in enumerate(polys[r]) if v}, 'Gaussian-binomial distribution')
        require(max(sinkhist) == n * n // 4 - 1, 'sharp family barrier')
        if n <= 8:
            actual = Counter(energy(a, p) for p in it.permutations(range(n)) if all(a[p[i]][p[i + 1]] for i in range(n - 1)))
            require(actual == sinkhist, 'exhaustive family path count')
        result.append({'n': n, 'H': sum(sinkhist.values()), 'sink_energy_histogram': dict(sorted(sinkhist.items())), 'gaussian_rank_polynomials': polys, 'galois_counts_q2_q3': [sum(evalpoly(p, q) for p in polys) for q in (2, 3)]})
    return result


def rref(a, p, m):
    a = [list(v) for v in a if any(v)]
    r = 0
    for c in range(m):
        k = next((i for i in range(r, len(a)) if a[i][c]), None)
        if k is None: continue
        a[r], a[k] = a[k], a[r]
        u = pow(a[r][c], -1, p)
        a[r] = [v * u % p for v in a[r]]
        for i in range(len(a)):
            if i != r:
                u = a[i][c]
                a[i] = [(a[i][j] - u * a[r][j]) % p for j in range(m)]
        r += 1
        if r == len(a): break
    return tuple(tuple(v) for v in a[:r])


def finite_fields():
    result = []
    for q, maxm in ((2, 5), (3, 4)):
        for m in range(maxm + 1):
            images = set()
            byrank = Counter()
            for bits in range(1 << m):
                A = [i for i in range(m) if bits >> i & 1]
                B = [i for i in range(m) if not bits >> i & 1]
                slots = [(i, b) for i, a in enumerate(A) for b in B if b < a]
                for values in it.product(range(q), repeat=len(slots)):
                    rows = [[int(j == a) for j in range(m)] for a in A]
                    for (i, b), x in zip(slots, values): rows[i][b] = x
                    image = rref(rows, q, m)
                    require(image not in images, 'decorated paths injective into subspaces')
                    images.add(image)
                    byrank[len(A)] += 1
            expected = [evalpoly(poly, q) for poly in qbinomial_rows(m)]
            require([byrank[r] for r in range(m + 1)] == expected, 'finite field rank counts')
            # Independent surjectivity: close the zero subspace under adjoining every vector.
            if m <= 3:
                allspaces = {()}
                todo = [()]
                vectors = list(it.product(range(q), repeat=m))
                while todo:
                    rows = todo.pop()
                    for v in vectors:
                        image = rref(rows + (v,), q, m)
                        if image not in allspaces:
                            allspaces.add(image)
                            todo.append(image)
                require(images == allspaces, 'independent subspace closure')
            result.append({'q': q, 'm': m, 'subspaces': len(images), 'rank_counts': expected, 'independent_closure': m <= 3})
    return result


def arbitrary_sink_construction(n, k):
    perms = list(it.permutations(range(n)))
    even = [p for p in perms if sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n)) % 2 == 0]
    S = set(even[:k])
    require(len(S) == k, 'independent sink size')
    dist = {p: 0 for p in S}
    todo = list(S)
    for p in todo:
        for i in range(n - 1):
            v = p[:i] + (p[i + 1], p[i]) + p[i + 2:]
            if v not in dist:
                dist[v] = dist[p] + 1
                todo.append(v)
    sinks = []
    directions = {}
    for p in perms:
        vals = []
        for i in range(n - 1):
            v = p[:i] + (p[i + 1], p[i]) + p[i + 2:]
            vals.append((dist[v], v))
            # From higher tuple to lower: strictly decreasing lexicographic potential.
            if (dist[p], p) > (dist[v], v):
                a, b = p[i:i + 2]
                key = tuple(sorted((a, b)))
                directions.setdefault(key, set()).add(int(v.index(key[0]) < v.index(key[1])))
        if all((dist[p], p) < value for value in vals): sinks.append(p)
    require(set(sinks) == S, 'arbitrary acyclic prescribed sink set')
    require(any(len(s) == 2 for s in directions.values()), 'gap example must break parallel coherence')
    return {'n': n, 'sinks': k, 'orientation_rule': 'decrease (distance to first k even permutations, permutation lex order)', 'noncoherent_parallel_classes': [list(pair) for pair, s in sorted(directions.items()) if len(s) == 2]}


def main():
    census, hostile = full_census()
    out = {'status': 'PASS', 'code_convention': 'lexicographic pairs i<j; bit 1 means i beats j', 'complete_census': census, 'minimal_local_not_global': hostile, 'one_extreme_flip': flip_family(), 'finite_field_bijection': finite_fields(), 'ordinary_acyclic_gap_hostiles': [arbitrary_sink_construction(4, 7), arbitrary_sink_construction(5, 21)]}
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
