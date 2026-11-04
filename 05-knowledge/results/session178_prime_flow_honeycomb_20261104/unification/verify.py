#!/usr/bin/env python3
"""Exact finite corroboration of the observer boundary-debt theorem."""
from fractions import Fraction
from itertools import combinations, product
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def rank(matrix):
    if not matrix:
        return 0
    a = [[Fraction(x) for x in row] for row in matrix]
    row = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        factor = a[row][col]
        a[row] = [x / factor for x in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                factor = a[i][col]
                a[i] = [x - factor*y for x, y in zip(a[i], a[row])]
        row += 1
        if row == len(a):
            break
    return row


def partitions(n):
    if n == 0:
        yield ()
        return
    def grow(prefix):
        if len(prefix) == n:
            yield tuple(prefix)
            return
        for label in range(max(prefix)+2):
            yield from grow(prefix + [label])
    yield from grow([0])


def boundary(n, edges, z):
    d = [0]*n
    for (u, v), amount in zip(edges, z):
        d[u] -= amount
        d[v] += amount
    return d


def observer(d, labels):
    result = [0]*(max(labels)+1)
    for value, label in zip(d, labels):
        result[label] += value
    return result


def paths_from(n, edges, root):
    paths = {root: [0]*len(edges)}
    queue = [root]
    for u in queue:
        for i, (a, b) in enumerate(edges):
            if a == u and b not in paths:
                v, sign = b, 1
            elif b == u and a not in paths:
                v, sign = a, -1
            else:
                continue
            paths[v] = paths[u].copy()
            paths[v][i] += sign
            queue.append(v)
    return paths


def check_graph(n, edges, labels):
    y = max(labels)+1
    cols = [boundary(n, edges, [int(i == j) for i in range(len(edges))])
            for j in range(len(edges))]
    b = [[col[v] for col in cols] for v in range(n)]
    qb = [[sum(col[v] for v in range(n) if labels[v] == w)
           for col in cols] for w in range(y)]
    require(rank(b) == n-1, 'original graph incidence rank')
    require(rank(qb) == y-1, 'quotient graph incidence rank')
    require((len(edges)-rank(qb))-(len(edges)-rank(b)) == n-y,
            'cycle rank difference')
    anchors = {}
    for v, label in enumerate(labels):
        anchors.setdefault(label, v)
    witnessed = 0
    for v, label in enumerate(labels):
        root = anchors[label]
        if v == root:
            continue
        z = paths_from(n, edges, root)[v]
        d = boundary(n, edges, z)
        expected = [int(i == v)-int(i == root) for i in range(n)]
        require(d == expected, 'integral surjectivity witness')
        require(not any(observer(d, labels)), 'witness is observed cycle')
        witnessed += 1
    require(witnessed == n-y, 'anchor basis cardinality')


def valuations(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0)+1
    return out


def sigma(n):
    return sum(d for d in range(1, n+1) if n % d == 0)


def main():
    graphs = cases = 0
    for n in range(1, 5):
        pairs = list(combinations(range(n), 2))
        for states in product(range(3), repeat=len(pairs)):
            edges = [(a, b) if state == 1 else (b, a)
                     for (a, b), state in zip(pairs, states) if state]
            if len(paths_from(n, edges, 0)) != n:
                continue
            graphs += 1
            for labels in partitions(n):
                check_graph(n, edges, labels)
                cases += 1
    # Parallel arrows and genuine loops, outside simple-graph enumeration.
    for edges in [[(0, 0), (0, 1), (0, 1), (1, 0)], [(0, 0)]]:
        n = 1+max(max(e) for e in edges)
        for labels in partitions(n):
            check_graph(n, edges, labels)
    # Collatz 1 -> 4 -> 2 -> 1; parity identifies vertices 4 and 2.
    edges = [(0, 1), (1, 2), (2, 0)]
    check_graph(3, edges, (0, 1, 1))
    require(boundary(3, edges, [0, 1, 0]) == [0, -1, 1], 'Collatz debt')
    require(boundary(3, edges, [1, 1, 1]) == [0, 0, 0], 'Collatz closure')
    require(Fraction(4, 1)*Fraction(2, 4)*Fraction(1, 2) == 1, 'Collatz product')
    # Multiplicative second difference has nonzero universal debt.
    require(boundary(3, [(0, 1), (1, 2)], [1, -1]) == [-1, 2, -1],
            'valuation-loss debt')
    require(Fraction(4, 2)/Fraction(8, 4) == 1, 'valuation-loss holonomy')
    # Exact component accounting for every N through 200.
    for n in range(1, 201):
        factors = valuations(n)
        predicted = {}
        product_sigma = 1
        for p, a in factors.items():
            value = sum(p**j for j in range(a+1))
            product_sigma *= value
            for q, b in valuations(value).items():
                predicted[q] = predicted.get(q, 0)+b
            predicted[p] = predicted.get(p, 0)-a
        actual = valuations(sigma(n))
        for p, a in factors.items():
            actual[p] = actual.get(p, 0)-a
        require({p:a for p,a in predicted.items() if a} ==
                {p:a for p,a in actual.items() if a}, 'sigma prime balance')
        require(product_sigma == sigma(n), 'sigma multiplicativity')
    # N=6 perfect, but its sigma reaction sum is not a state circulation.
    # Vertex order 1,2,3,4.
    require(boundary(4, [(1, 2), (2, 3), (1, 0)], [1, 1, 1]) == [1, -2, 0, 1],
            'perfect reaction debt')
    require(Fraction(3, 2)*Fraction(4, 3)*Fraction(1, 2) == 1,
            'perfect arithmetic balance')
    print(json.dumps({'status': 'PASS', 'connected_oriented_graphs': graphs,
                      'graph_partition_cases': cases, 'sigma_inputs': 200,
                      'integer_path_witnesses': True,
                      'parallel_loop_collatz_and_loss_controls': True}, sort_keys=True))


if __name__ == '__main__':
    main()
