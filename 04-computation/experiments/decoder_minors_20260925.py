"""Exact certificates for the +2/doubling graph and decoder minor boundaries.

Run: python 04-computation/experiments/decoder_minors_20260925.py
No input filters. Core certificates use only the standard library. NetworkX,
when available, supplies an independent planarity cross-check on N=1..64.
"""
from itertools import combinations


def operation_edges(n):
    return {tuple(sorted((x, x + 2))) for x in range(1, n - 1)} | {
        (x, 2 * x) for x in range(1, n // 2 + 1)
    }


def v2(n):
    return (n & -n).bit_length() - 1


ROTATION15 = {
    1: [3, 2], 2: [4, 1], 3: [1, 5, 6], 4: [6, 8, 2],
    5: [3, 10, 7], 6: [12, 8, 4, 3], 7: [5, 14, 9],
    8: [4, 6, 10], 9: [7, 11], 10: [8, 12, 5],
    11: [9, 13], 12: [14, 10, 6], 13: [11, 15],
    14: [7, 12], 15: [13],
}

PATHS16 = [
    [5, 3, 6], [5, 10], [5, 7, 14],
    [8, 6], [8, 10], [8, 16, 14],
    [12, 6], [12, 10], [12, 14],
]


def check_rotation(rotation, edges):
    darts = {(x, y) for e in edges for x, y in (e, e[::-1])}
    assert {(x, y) for x, ys in rotation.items() for y in ys} == darts
    assert all(len(ys) == len(set(ys)) for ys in rotation.values())
    visited = set()
    faces = []
    for dart in sorted(darts):
        if dart in visited:
            continue
        face = []
        current = dart
        while current not in visited:
            visited.add(current)
            face.append(current)
            x, y = current
            neighbors = rotation[y]
            current = (y, neighbors[(neighbors.index(x) + 1) % len(neighbors)])
        assert current == dart
        faces.append(face)
    assert visited == darts
    reached = {next(iter(rotation))}
    while True:
        grown = reached | {y for x in reached for y in rotation[x]}
        if grown == reached:
            break
        reached = grown
    assert reached == set(rotation)
    chi = len(rotation) - len(edges) + len(faces)
    assert chi == 2
    return faces, chi


def check_k33_paths(paths, edges):
    left, right = {5, 8, 12}, {6, 10, 14}
    branch = left | right
    assert {(p[0], p[-1]) for p in paths} == {(a, b) for a in left for b in right}
    interiors = set()
    for path in paths:
        assert len(path) == len(set(path))
        assert all(tuple(sorted(e)) in edges for e in zip(path, path[1:]))
        assert not (set(path[1:-1]) & branch)
        assert not (set(path[1:-1]) & interiors)
        interiors.update(path[1:-1])
    return sorted(interiors)


def main():
    print("STATUS: exact certificates; no Collatz convergence conclusion")
    print("GRAPH: vertices 1..N; undirected simple edges {x,x+2}, {x,2x}")
    print("UNIVERSE: certificates N=15,16; valuation n=1..256; optional cross-check N=1..64")
    edges15 = operation_edges(15)
    faces, chi = check_rotation(ROTATION15, edges15)
    print(f"G15: V=15 E={len(edges15)} F={len(faces)} Euler={chi}; planar rotation verified")
    print("G15 face boundary lengths:", sorted(map(len, faces)))
    interiors = check_k33_paths(PATHS16, operation_edges(16))
    print("G16: K3,3 branch shores [5,8,12] and [6,10,14]")
    print("G16 paths:", PATHS16)
    print("G16 pairwise-disjoint internal vertices:", interiors)
    print("ALL-N CONCLUSION: G_N planar iff N<=15 (nested family plus fixed certificates)")
    counts = {"odd_to_odd": 0, "seam_to_sea": 0, "sea_to_seam": 0}
    for n in range(1, 257):
        r, s = v2(n), v2(n + 2)
        if r == 0:
            assert s == 0
            counts["odd_to_odd"] += 1
        elif r == 1:
            assert s >= 2
            counts["seam_to_sea"] += 1
        else:
            assert s == 1
            counts["sea_to_seam"] += 1
        q = n >> r
        k = (q - 1) // 2
        assert n == (1 << r) * (2 * k + 1)
        assert 2 * n == (1 << (r + 1)) * (2 * k + 1)
        assert n + (1 << (r + 1)) == (1 << r) * (2 * (k + 1) + 1)
    print("VALUATION +2 controls:", counts)
    print("ROW-SCALED +2^(r+1): quadrant-grid coordinate identities verified n=1..256")
    # A path-contraction certificate keeps arithmetic word and intermediate guard.
    x = 4
    trace = [x]
    for letter in "DDO":
        if letter == "D":
            x *= 2
        else:
            assert x % 6 == 4
            x = (x - 1) // 3
        trace.append(x)
    assert trace == [4, 8, 16, 5]
    print("GUARDED CONTRACTED PATH: word DDO; trace", trace)
    # Underlying complete graph loses all tournament orientations.
    assert len(list(combinations(range(5), 2))) == 10
    assert len([(i, j) for i in range(3) for j in range(3, 6)]) == 9
    print("UNDERLYING TOURNAMENT: K5 subgraph from order5; K3,3 subgraph from order6")
    try:
        import networkx as nx
    except ImportError:
        print("NETWORKX: unavailable; core independent certificates already passed")
    else:
        for n in range(1, 65):
            graph = nx.Graph()
            graph.add_nodes_from(range(1, n + 1))
            graph.add_edges_from(sorted(operation_edges(n)))
            assert nx.check_planarity(graph)[0] == (n <= 15)
        print("NETWORKX", nx.__version__, "independent planarity cross-check N=1..64: PASS")
    print("PASS")


if __name__ == "__main__":
    main()
