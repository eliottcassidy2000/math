"""Independent literal audits of the two oriented path-count compilers.

No imports from the producer. Input universe: every loopless labelled digraph
through order four, plus the specified positive integer construction targets.
Normal and optimized Python must produce identical stdout.
"""
from collections import Counter
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


def path_total(a):
    n = len(a)
    if n == 0:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, weight in enumerate(dp[mask]):
            if weight:
                for w in range(n):
                    if a[v][w] and not (mask >> w & 1):
                        dp[mask | (1 << w)][w] += weight
    return sum(dp[-1])


def literal_paths(a):
    return [p for p in permutations(range(len(a)))
            if all(a[u][v] for u, v in zip(p, p[1:]))]


def literal_cycles(a):
    n = len(a)
    if n < 2:
        return 0
    return sum(all(a[u][v] for u, v in zip(p, p[1:] + p[:1]))
               for q in permutations(range(1, n)) for p in [(0,) + q])


def oriented(a):
    return all(not a[i][i] for i in range(len(a))) and all(
        a[i][j] + a[j][i] <= 1
        for i in range(len(a)) for j in range(i))


def strong(a):
    n = len(a)
    for start in range(n):
        reached = {start}
        todo = [start]
        while todo:
            u = todo.pop()
            for v in range(n):
                if a[u][v] and v not in reached:
                    reached.add(v)
                    todo.append(v)
        if len(reached) != n:
            return False
    return True


def split_all(a, source):
    n = len(a)
    size = 2 * n + int(source)
    out = [[0] * size for _ in range(size)]
    for v in range(n):
        out[2 * v][2 * v + 1] = 1
        if source:
            out[-1][2 * v] = 1
        for w in range(n):
            # Diagonal entries are intentionally discarded.
            if v != w and a[v][w]:
                out[2 * v + 1][2 * w] = 1
    return out


def strong_integer_graph(h):
    require(h >= 3, 'Strong compiler domain h >= 3')
    k = (h - 2).bit_length()
    a, b = k, k + 1
    out = [[0] * (k + 2) for _ in range(k + 2)]
    for i in range(k):
        for j in range(i + 1, k):
            out[i][j] = 1
        out[i][a] = 1
        if (h - 2) >> (k - 1 - i) & 1:
            out[b][i] = 1
    out[a][b] = 1
    return out


def compressed_source(g):
    k = len(g) - 2
    a, b = k, k + 1
    d = [[0] * (k + 1) for _ in range(k + 1)]
    for i in range(k):
        for j in range(k):
            d[i][j] = g[i][j]
        d[i][k] = g[i][a]
        d[k][i] = g[b][i]
    return d


def matrix_text(g):
    return [''.join(map(str, row)) for row in g]


def main():
    generic = []
    for n in range(1, 5):
        pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
        count = 0
        cycle_histogram = Counter()
        for mask in range(1 << len(pairs)):
            d = [[0] * n for _ in range(n)]
            for bit, (i, j) in enumerate(pairs):
                d[i][j] = mask >> bit & 1
            p = literal_paths(d)
            c = literal_cycles(d)
            g = split_all(d, True)
            require(oriented(g), 'Generic output oriented')
            require(path_total(g) == len(p), 'Generic path total')
            require(path_total(split_all(d, False)) == len(p) + n * c,
                    'Unforced split path-plus-cycle boundary')
            for q in p:
                lifted = (2 * n,) + tuple(z for v in q for z in (2 * v, 2 * v + 1))
                require(len(set(lifted)) == len(g) and all(
                    g[u][v] for u, v in zip(lifted, lifted[1:])),
                    'Literal lifted path')
            # Loop changes are ignored by the compiler, as they are by H.
            for v in range(n):
                d[v][v] = 1
            require(split_all(d, True) == g, 'Input loops ignored')
            cycle_histogram[c] += 1
            count += 1
        generic.append(dict(order=n, labelled_inputs=count,
                            cycle_count_histogram=dict(sorted(cycle_histogram.items()))))

    # Count both all targets in a modest interval and chosen binary boundaries.
    targets = sorted(set(range(3, 513)) | {
        2 ** k + delta for k in range(9, 13) for delta in (0, 1, 2, 3)})
    examples = []
    literal_checked = 0
    total_graphs = 0
    for h in targets:
        g = strong_integer_graph(h)
        n = len(g)
        require(oriented(g), 'Strong compiler orientation')
        require(strong(g), 'Strong compiler connectivity')
        require(path_total(g) == h, 'Strong compiler count')
        d = compressed_source(g)
        require(path_total(d) == h - 1, 'Compressed source total')
        require(literal_cycles(d) == 1 if len(d) <= 8 else True,
                'Unique source Hamilton cycle at direct-check orders')
        if n <= 8:
            paths = literal_paths(g)
            require(len(paths) == h, 'Independent target permutation count')
            a, b = n - 2, n - 1
            without = [p for p in paths if (a, b) not in list(zip(p, p[1:]))]
            require(without == [(b,) + tuple(range(a)) + (a,)],
                    'Exactly one path avoiding forced connector')
            literal_checked += 1
        if h in (3, 4, 5, 7, 21, 127, 128, 129, 130, 511, 512, 513, 514):
            examples.append(dict(target=h, order=n, rows=matrix_text(g)))
        total_graphs += 1

    require(path_total([[0, 0], [0, 0]]) == 0, 'Zero-total control')
    require(path_total([[0]]) == 1, 'One-total control')
    # Two directed triangles sharing one vertex: a strong oriented H=2 control.
    h2 = [[0] * 5 for _ in range(5)]
    for u, v in [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]:
        h2[u][v] = 1
    require(oriented(h2) and strong(h2) and path_total(h2) == 2,
            'Inherited H2 strong control')
    print(json.dumps(dict(
        status='PASS', script_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        generic_universe='Every loopless labelled digraph on n=1,2,3,4; all-loop controls',
        generic=generic, generic_inputs=sum(q['labelled_inputs'] for q in generic),
        compact_targets=targets, compact_graphs=total_graphs,
        compact_target_dp=True, compact_source_dp=True,
        compact_target_permutation_checks=literal_checked,
        compact_examples=examples,
        controls=['H0 edgeless pair', 'H1 singleton', 'H2 two triangles sharing one vertex',
                  'Unsourced split adds n times the Hamilton cycle count',
                  'Loops are discarded before compilation'],
        scope='Finite independent checks of all-order combinatorial proofs; no tournament spectrum closure'),
        indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
