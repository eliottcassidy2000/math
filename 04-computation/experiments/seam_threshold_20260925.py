"""Exact local threshold and triangular-ten controls; no conjecture transfer.

Reproduction: python 04-computation/experiments/seam_threshold_20260925.py
Core checks use only the standard library. Optional NetworkX independently
checks the six planarity flags at N=14,15,16.
"""
from collections import Counter
from itertools import combinations
from math import isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def graph(n, kind):
    adjacency = {v: set() for v in range(1, n + 1)}
    for u, v in combinations(adjacency, 2):
        hit = (isqrt(u + v) ** 2 == u + v) if kind == "Q" else (v == u + 2 or v == 2 * u)
        if hit:
            adjacency[u].add(v)
            adjacency[v].add(u)
    return adjacency


def edges(adjacency):
    return {frozenset((u, v)) for u in adjacency for v in adjacency[u]}


def components(adjacency):
    rest, result = set(adjacency), []
    while rest:
        seen = {rest.pop()}
        todo = list(seen)
        while todo:
            v = todo.pop()
            new = adjacency[v] & rest
            rest -= new
            seen |= new
            todo.extend(new)
        result.append(seen)
    return result


def hamilton_paths(adjacency):
    """Unpruned vertex DFS, fixing a forced leaf when one exists."""
    n = len(adjacency)
    if n == 1:
        return [[1]]
    if any(not neighbors for neighbors in adjacency.values()):
        return []
    leaves = [v for v in adjacency if len(adjacency[v]) == 1]
    require(leaves, "this finite universe uses a forced leaf")
    start, found = min(leaves), []

    def visit(path, unseen):
        if not unseen:
            found.append(path)
        else:
            for v in sorted(adjacency[path[-1]] & unseen):
                visit(path + [v], unseen - {v})

    visit([start], set(adjacency) - {start})
    return found


def main():
    print("STATUS: elementary structural proofs plus exact finite controls")
    print("UNIVERSES: Q1..Q16; local G14..G16; all1024 orientations of K5; all perfect matchings of L(K5)")
    counts = {n: len(hamilton_paths(graph(n, "Q"))) for n in range(1, 17)}
    require(counts == {n: int(n in (1, 15, 16)) for n in range(1, 17)}, "small Q counts")
    print("Q_N Hamiltonian paths up to reversal, N1..16:", counts)
    for n in (14, 15, 16):
        for kind in ("Q", "G"):
            a = graph(n, kind)
            comp = len(components(a))
            print(kind + str(n), "V", n, "E", len(edges(a)), "components", comp,
                  "cycle_rank", len(edges(a)) - n + comp,
                  "leaves", sorted(v for v in a if len(a[v]) == 1),
                  "new_neighbors", sorted(a[n]))
    q14, q15, q16 = (graph(n, "Q") for n in (14, 15, 16))
    require(len(components(q14)) == 1 and len(edges(q14)) == 13, "Q14 tree")
    arms = [[3, 1, 8], [3, 6, 10], [3, 13, 12, 4, 5, 11, 14, 2, 7, 9]]
    require({frozenset(e) for p in arms for e in zip(p, p[1:])} == edges(q14), "Q14 claw arms")
    hp15 = [8, 1, 15, 10, 6, 3, 13, 12, 4, 5, 11, 14, 2, 7, 9]
    require(hamilton_paths(q15) == [hp15], "unique Q15 path")
    require(edges(q15) - {frozenset(e) for e in zip(hp15, hp15[1:])}
            == {frozenset((1, 3))}, "unique omitted chord")
    require(hamilton_paths(q16) == [hp15 + [16]], "Q16 endpoint extension")
    print("Q14 three arms:", arms)
    print("Q15 unique path:", hp15, "; omitted edge1--3")
    print("Q ear1--15--10 replaces its third leaf; G ear8--16--14 completes the recorded K3,3")

    pairs = list(combinations(range(1, 6), 2))
    line = {i: {j for j in range(10) if j != i and set(pairs[i]) & set(pairs[j])}
            for i in range(10)}
    disjoint = {i: set(range(10)) - line[i] - {i} for i in range(10)}
    require(len(edges(line)) == 30 and len(edges(disjoint)) == 15, "L/Petersen edges")
    require(all(len(line[i]) == 6 and len(disjoint[i]) == 3 for i in range(10)), "L/Petersen regularity")
    require(not any(all((w in line[u]) == (w in line[v]) for w in range(10) if w not in (u, v))
                    for u, v in combinations(range(10), 2)), "L(K5) pair-module hostile")
    print("T4=10=edges(K5)=vertices(L(K5)); L edges30 degree6; complement edges15 degree3")
    print("L(K5) has no two-vertex module despite its perfect matchings")
    matchings = []

    def match(rem, chosen):
        if not rem:
            matchings.append(chosen)
            return
        u = min(rem)
        for v in sorted(line[u] & rem):
            match(rem - {u, v}, chosen + [(u, v)])

    match(set(range(10)), [])
    matching_profiles, orientations = Counter(), Counter()
    regular_matchings = 0
    for matching in matchings:
        owners, center_counts = {}, Counter()
        for i, j in matching:
            center, = set(pairs[i]) & set(pairs[j])
            owners[i] = owners[j] = center
            center_counts[center] += 1
        mask = sum(1 << i for i, (u, v) in enumerate(pairs) if owners[i] == u)
        orientations[mask] += 1
        scores = Counter(owners.values())
        require(all(scores[v] % 2 == 0 for v in range(1, 6)), "even outgoing scores")
        profile = tuple(sorted(center_counts.values()))
        matching_profiles[profile] += 1
        if profile == (1, 1, 1, 1, 1):
            regular_matchings += 1
            # The five connected pair branch sets give a complete ordinary minor.
            require(all(any(j in line[i] for i in A for j in B)
                        for A, B in combinations(matching, 2)), "K5 branch-set quotient")
    even_masks, regular_masks = set(), set()
    for mask in range(1 << 10):
        scores = Counter()
        for i, (u, v) in enumerate(pairs):
            scores[u if mask & (1 << i) else v] += 1
        if all(scores[v] % 2 == 0 for v in range(1, 6)):
            even_masks.add(mask)
            if all(scores[v] == 2 for v in range(1, 6)):
                regular_masks.add(mask)
                require(orientations[mask] == 1, "regular orientation unique paired decoder")
            else:
                require(sorted(scores[v] for v in range(1, 6)) == [0, 2, 2, 2, 4], "nonregular score profile")
                require(orientations[mask] == 3, "source four-arcs pairing count")
    require(set(orientations) == even_masks, "matching/orientation surjection")
    require((len(matchings), regular_matchings, len(even_masks), len(regular_masks)) == (144, 24, 64, 24), "full counts")
    print("All L(K5) perfect matchings:", len(matchings), "center profiles", sorted(matching_profiles.items()))
    print("Decoded even-score K5 orientations64 = regular24 + source-C3-sink40")
    print("Matching identity144 =24*1+40*3; all24 regular matchings give connected-pair K5 minor models")
    triangle_edges = [{frozenset(e) for e in combinations(triple, 2)}
                      for triple in combinations(range(1, 6), 3)]
    packs = sum(len(set.union(*triple)) == 9 for triple in combinations(triangle_edges, 3))
    require(packs == 0, "three edge-disjoint triangles hostile")
    print("K5 cannot partition into three edge-disjoint triangle edge-sets plus one edge:0 packs")
    try:
        import networkx as nx
    except ImportError:
        print("Optional NetworkX unavailable; prior G15/G16 certificates remain the planarity proof")
    else:
        flags = {}
        for n in (14, 15, 16):
            for kind in ("Q", "G"):
                a = graph(n, kind)
                g = nx.Graph(a)
                flags[kind + str(n)] = nx.check_planarity(g)[0]
        require(flags == {"Q14": True, "Q15": True, "Q16": True,
                          "G14": True, "G15": True, "G16": False}, "planarity flags")
        print("Independent NetworkX", nx.__version__, "planarity:", flags)
    print("PASS; no square-sum/Collatz/minor equivalence claimed")


if __name__ == "__main__":
    main()
