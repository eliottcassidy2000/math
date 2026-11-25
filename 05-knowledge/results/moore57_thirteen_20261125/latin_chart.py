"""Build normalized Moore branch charts and test their reciprocal Latin rows."""
from itertools import combinations
from math import factorial, prod

from sympy import integer_nthroot

from moore_chart import hoffman_edge


def chart_from_graph(graph, root):
    heads = sorted(graph[root])
    slots = [graph[head] - {root} for head in heads]
    anchor = sorted(slots[0])
    columns = []
    for sites in slots[1:]:
        columns.append([next(iter(graph[x] & sites)) for x in anchor])
    permutations = {}
    for i, source in enumerate(columns):
        for j, target in enumerate(columns):
            if i != j:
                indices = {vertex: pos for pos, vertex in enumerate(target)}
                permutations[i, j] = tuple(indices[next(iter(graph[x] & slots[j+1]))]
                                           for x in source)
    return permutations


def passes_chart(q, permutations):
    labels = range(q)
    identity = tuple(labels)
    expected = {(i, j) for i in labels for j in labels if i != j}
    if set(permutations) != expected:
        return False
    if any(len(row) != q or set(row) != set(labels) for row in permutations.values()):
        return False

    def latin(rows):
        return all(len(set(column)) == q for column in zip(*rows))

    for j in labels:
        if not latin([identity] + [permutations[i, j] for i in labels if i != j]):
            return False
    for i, j in combinations(labels, 2):
        forward, backward = permutations[i, j], permutations[j, i]
        if any(backward[forward[x]] != x for x in labels):
            return False
        rows = [identity, forward]
        for k in labels:
            if k != i and k != j:
                rows.append(tuple(permutations[k, j][permutations[i, k][x]] for x in labels))
        if not latin(rows):
            return False
    return True


def replay():
    petersen_vertices = list(combinations(range(5), 2))
    petersen = {v: {w for w in petersen_vertices if set(v).isdisjoint(w)}
                for v in petersen_vertices}
    hoffman = {v: {w for w in range(50) if hoffman_edge(v, w)} for v in range(50)}
    for graph in (petersen, hoffman):
        for root in graph:
            q = len(graph[root])-1
            chart = chart_from_graph(graph, root)
            assert passes_chart(q, chart)
            chart[0, 1] = tuple(range(q))
            assert not passes_chart(q, chart)
        print("All rooted charts and corrupted controls passed on", len(graph), "vertices")

    q, scale = 56, 10**8
    upper = []
    for r in range(1, q+1):
        floor, exact = integer_nthroot(factorial(r)**q * scale**r, r)
        upper.append(floor + (not exact))
    assert prod(upper)**28 * 10**32412 < (scale**q * factorial(q)**q)**28
    print("Certified probability bound for order 56")


if __name__ == "__main__":
    replay()
