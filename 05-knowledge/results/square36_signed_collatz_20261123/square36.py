"""Certificates for square-color potentials, local exits and the edge (6,30)."""
from collections import Counter
from itertools import combinations
from math import isqrt


def square(value):
    return isqrt(value)**2 == value


def graph_edges(n):
    return {(a, b) for a, b in combinations(range(1, n + 1), 2) if square(a + b)}


def components(vertices, edges):
    neighbors = {v: set() for v in vertices}
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    remaining = set(neighbors)
    groups = []
    while remaining:
        root = min(remaining)
        remaining.remove(root)
        stack, group = [root], {root}
        while stack:
            fresh = neighbors[stack.pop()] & remaining
            group.update(fresh)
            stack.extend(fresh)
            remaining.difference_update(fresh)
        groups.append(frozenset(group))
    return groups


def difference_components(n):
    vertices = range(1, n)
    return components(vertices, [(a, b) for a, b in combinations(vertices, 2)
                                 if b - a >= 2 and square(a + b + 1)])


def all_factors(vertices, edges, degree):
    """Enumerate each degree-regular spanning subgraph once by include/exclude."""
    def visit(kept, free):
        while True:
            required = {v: degree for v in vertices}
            for a, b in kept:
                required[a] -= 1
                required[b] -= 1
            incident = {v: set() for v in vertices}
            for e in free:
                for x in e:
                    incident[x].add(e)
            take, drop = set(), set()
            for v in vertices:
                count = required[v]
                if count < 0 or count > len(incident[v]):
                    return
                if count == 0:
                    drop.update(incident[v])
                elif count == len(incident[v]):
                    take.update(incident[v])
            if take & drop:
                return
            if not take and not drop:
                break
            kept = kept | take
            free = free - take - drop
        if not free:
            yield frozenset(kept)
            return
        edge = min(free, key=lambda e: (min(len(incident[x]) for x in e), e))
        yield from visit(kept | {edge}, free - {edge})
        yield from visit(kept, free - {edge})
    yield from visit(set(), edges)


def unit(value, p):
    level = 0
    while value % p == 0:
        value //= p
        level += 1
    sign = 1 if pow(value, (p - 1) // 2, p) == 1 else -1
    return level, sign


def layer_exits(factor, layer):
    frequency = Counter()
    for u, v in factor:
        if (u in layer) != (v in layer):
            frequency[layer[u if u in layer else v]] += 1
    return frequency


def main():
    old = difference_components(28)
    assert len(old) == 2 and 7 in old[0] and 20 in old[1]
    assert len(difference_components(29)) == 1
    g = lambda x: (9*x - 3)//20
    for a, b in graph_edges(28):
        root = isqrt(a+b)
        assert g(a) + g(b) == root*(root-1)//2
    print("difference components at 28", [sorted(s) for s in old])
    vertices = range(1, 37)
    edges = graph_edges(36)
    matchings = list(all_factors(vertices, edges, 1))
    factors = list(all_factors(vertices, edges, 2))
    cycles = [f for f in factors if len(components(vertices, f)) == 1]
    evens = [f for f in factors if all(len(c)%2 == 0 for c in components(vertices, f))]
    assert (len(edges),len(matchings),len(factors),len(evens),len(cycles)) == (57,295,70,59,31)
    obstructed = (6,30)
    assert edges - set.union(*(set(f) for f in cycles)) == {obstructed}
    assert sum(obstructed in f for f in matchings) == 15
    assert sum(obstructed in f for f in evens) == 0
    supports = Counter(next(c for c in components(vertices, f) if 6 in c)
                       for f in factors if obstructed in f)
    assert supports == {frozenset((6,19,30)):3,
                        frozenset((2,4,6,10,15,17,19,21,23,26,30,32,34)):6}
    walks = ((6,19,30), (6,19,17,32,4,21,15,10,26,23,2,34,30),
             (6,10,26,23,2,34,15,21,4,32,17,19,30))
    for walk in walks:
        adjacent = tuple(zip(walk, walk[1:] + walk[:1]))
        block = {tuple(sorted((u,v))) for u,v in adjacent}
        assert sum(block <= f for f in factors) == 3
        assert sum((-1)**i * (u+v) for i,(u,v) in enumerate(adjacent)) == 12
    pairs = [f | s for f,s in combinations(matchings,2) if not f&s]
    assert set(pairs) == set(evens) and len(pairs) == 95
    assert sum((9,16) in f for f in cycles) == 21 and sum((7,9) in f for f in cycles) == 10
    assert all((19,30) in f for f in cycles)
    for p, population in ((3,(4,4)), (5,(3,3)), (7,(3,2))):
        layer = {x:unit(x,p)[1] for x in vertices if unit(x,p)[0] == 1}
        sizes = Counter(layer.values())
        assert (sizes[1], sizes[-1]) == population
        for f in matchings:
            exits = layer_exits(f,layer)
            if p%4 == 3:
                assert exits[1]-exits[-1] == sizes[1]-sizes[-1]
            else:
                assert all((exits[sign]-sizes[sign])%2 == 0 for sign in (-1,1))
    print("matchings / two-factors / even two-factors / Hamilton cycles",
          len(matchings),len(factors),len(evens),len(cycles))
    print("disjoint matching pairs", len(pairs), "odd components at (6,30)", list(supports.values()))


if __name__ == "__main__":
    main()
