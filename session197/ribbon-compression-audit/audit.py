"""Independent checks for compression, arbitrary-part path merge, and sparsity."""
import itertools
import json
import math
import pathlib
import random
import sys
from functools import lru_cache

ROOT = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'source' / 'ribbon'))
import transfer


def adjacency(n, bits):
    graph = [[False] * n for i in range(n)]
    for i, j in itertools.combinations(range(n), 2):
        graph[i][j] = not (bits & 1)
        graph[j][i] = bool(bits & 1)
        bits >>= 1
    return graph


def path_counts(graph):
    n = len(graph)
    @lru_cache(None)
    def ending(vertices, last):
        rem = vertices ^ (1 << last)
        if not rem:
            return 1
        return sum(ending(rem, prev) for prev in range(n) if rem & (1 << prev) and graph[prev][last])
    def count(vertices):
        # The product convention is used only for empty parts in this helper.
        return sum(ending(vertices, last) for last in range(n) if vertices & (1 << last)) if vertices else 1
    return count


def sccs(graph):
    # Independently compute reachability using repeated depth-first searches.
    n = len(graph)
    reach = []
    for start in range(n):
        found = {start}
        stack = [start]
        while stack:
            u = stack.pop()
            for v in range(n):
                if graph[u][v] and v not in found:
                    found.add(v)
                    stack.append(v)
        reach.append(found)
    remaining = set(range(n))
    result = []
    while remaining:
        first = min(remaining)
        component = {v for v in remaining if v in reach[first] and first in reach[v]}
        result.append(component)
        remaining -= component
    return result


def triangle_at(graph, center, width):
    window = range(max(0, center - 2 * width), min(len(graph), center + 2 * width + 1))
    for a, b in itertools.permutations(window, 2):
        if graph[center][a] and graph[a][b] and graph[b][center]:
            return {center, a, b}
    raise AssertionError(('triangle missing', center, width))


def check_compression(graph, h=None):
    n = len(graph)
    width = max((j - i for i in range(n) for j in range(i + 1, n) if graph[j][i]), default=0)
    count = path_counts(graph) if h is None else None
    h = count((1 << n) - 1) if h is None else h
    parts = sccs(graph)
    keep = sorted(v for c in parts if len(c) > 1 for v in c)
    if not keep:
        assert h == 1
        return {'width': width, 'original': n, 'reduced': 0, 'h': h}
    assert h > 1 and width > 0
    reduced = [[graph[i][j] for j in keep] for i in keep]
    m = len(keep)
    assert max((j - i for i in range(m) for j in range(i + 1, m) if reduced[j][i]), default=0) <= width
    if count:
        assert count(sum(1 << v for v in keep)) == h
    triangles = [triangle_at(reduced, v, width) for v in range(m)]
    assert all(max(t) - min(t) <= 2 * width for t in triangles)
    centers = range(0, m, 4 * width + 1)
    chosen = [triangles[c] for c in centers]
    assert sum(map(len, chosen)) == len(set().union(*chosen))
    assert h >= 3 ** len(chosen)
    k, power = 0, 1
    while 3 * power <= h:
        k += 1
        power *= 3
    assert m <= (4 * width + 1) * k
    return {'width': width, 'original': n, 'reduced': m, 'h': h, 'packed': len(chosen)}


def paths(graph, vertices):
    return [p for p in itertools.permutations(vertices) if all(graph[u][v] for u, v in zip(p, p[1:]))]


def merge(graph, p, q):
    # Index cursors, separately implemented from the source's list popping.
    i = j = 0
    out = []
    while i < len(p) and j < len(q):
        if graph[p[i]][q[j]]:
            out.append(p[i])
            i += 1
        else:
            out.append(q[j])
            j += 1
    return tuple(out) + p[i:] + q[j:]


def main():
    report = {}
    tournaments = partitions = merges = 0
    for n in range(1, 7):
        for bits in range(1 << (n * (n - 1) // 2)):
            graph = adjacency(n, bits)
            check_compression(graph)
            tournaments += 1
            if n > 5:
                continue
            count = path_counts(graph)
            full = (1 << n) - 1
            # Every unordered nonempty two-part vertex partition, including
            # interleaved parts; require vertex zero in the first part.
            for left in range(1, full, 2):
                right = full ^ left
                assert count(full) >= count(left) * count(right)
                ls = [v for v in range(n) if left & (1 << v)]
                rs = [v for v in range(n) if right & (1 << v)]
                out = set()
                for p, q in itertools.product(paths(graph, ls), paths(graph, rs)):
                    z = merge(graph, p, q)
                    assert all(graph[u][v] for u, v in zip(z, z[1:]))
                    assert tuple(v for v in z if left & (1 << v)) == p
                    assert tuple(v for v in z if right & (1 << v)) == q
                    out.add(z)
                    merges += 1
                assert len(out) == count(left) * count(right)
                partitions += 1
    report['all_ordered_tournaments_through_n6_compression'] = tournaments
    report['all_arbitrary_partitions_through_n5'] = partitions
    report['distinct_merge_input_pairs_checked'] = merges
    print('exhaustive checks', report, flush=True)

    rng = random.Random(197197)
    large = []
    for width in range(1, 5):
        for n in (1, 30, 80, 200):
            word = [rng.randrange(1 << width) for i in range(n)]
            a = transfer.graph(word, width)
            h = transfer.count(word, width)
            row = check_compression(a, h)
            keep = sorted(v for c in sccs(a) if len(c) > 1 for v in c)
            if keep:
                b = [[a[i][j] for j in keep] for i in keep]
                # Encode the induced graph independently after deleting and
                # reindexing, then count with the already-audited transfer.
                recoded = [sum((1 << (width - (j - i))) for i in range(max(0, j - width), j) if b[j][i]) for j in range(len(b))]
                assert transfer.graph(recoded, width) == b
                assert transfer.count(recoded, width) == h
            large.append(row)
    report['large_random_compression_checks'] = large

    # H=1 and width=0: every tournament here is transitive, all SCCs are
    # singleton, and a nonempty canonical witness is retained separately.
    report['width0_h1_orders'] = []
    for n in (1, 2, 10, 100):
        row = check_compression(transfer.graph([0] * n, 0), 1)
        assert row['reduced'] == 0
        assert transfer.count([0], 0) == 1
        report['width0_h1_orders'].append(n)

    # Arbitrarily many transitive padding vertices show why compression is
    # essential. Three central vertices form a cycle; everything else is a
    # singleton SCC in an ordered join.
    padding = []
    for n in (3, 30, 100, 300):
        word = [0] * n
        start = (n - 3) // 2
        word[start + 1] = word[start + 2] = 1
        graph = transfer.graph(word, 1)
        assert transfer.count(word, 1) == 3
        row = check_compression(graph, 3)
        assert row['reduced'] == 3
        padding.append(row)
    report['unbounded_padding_h3_examples'] = padding

    sequence = [1, 1, 1]
    for n in range(3, 501):
        sequence.append(sum(sequence[-3:]))
        # Integer equivalent of h_n >= 3^((n-2)/3).
        assert sequence[n] ** 3 >= 3 ** (n - 2)
    bound = 100000
    monoid = {1}
    for h in sequence[3:]:
        if h > bound:
            break
        for value in list(monoid):
            while value * h <= bound:
                value *= h
                monoid.add(value)
    tests = []
    for X in (3, 10, 100, 1000, 10000, 100000):
        count = sum(x <= X for x in monoid)
        limit = math.exp(math.pi * math.sqrt(2 * math.log(X, 3)))
        assert count <= limit
        tests.append({'X': X, 'A': count, 'upper_bound': limit})
    report['width1_subpower_bound_checks'] = tests
    (ROOT / 'AUDIT.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
