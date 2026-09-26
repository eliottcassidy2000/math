"""Exact controls for square-sum endpoints, local switches, and affine lifts.

Run from the repository root: python 04-computation/experiments/crossroads_family_20260926_squares.py
No third-party dependencies; checks remain active under python -O.
"""
from collections import defaultdict, deque
from itertools import combinations, combinations_with_replacement
from math import isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def square_graph(n):
    return {v: {u for u in range(1, n + 1)
                if u != v and isqrt(u + v) ** 2 == u + v}
            for v in range(1, n + 1)}


def edge_set(g):
    return {tuple(sorted((u, v))) for u in g for v in g[u]}


def connected(g, vertices=None):
    vertices = set(g) if vertices is None else set(vertices)
    if not vertices:
        return True
    reached = {next(iter(vertices))}
    todo = list(reached)
    while todo:
        for u in (g[todo.pop()] & vertices) - reached:
            reached.add(u)
            todo.append(u)
    return reached == vertices


def hamilton_paths(g):
    """Enumerate unoriented paths with only connectivity/endpoint pruning."""
    n = len(g)
    if n == 1:
        return [tuple(g)]
    if not connected(g):
        return []
    leaves = sorted(v for v in g if len(g[v]) == 1)
    if len(leaves) > 2:
        return []
    starts = leaves[:1] if leaves else sorted(g)
    answers = []

    def dfs(path, remaining):
        v = path[-1]
        if not remaining:
            if leaves or path[0] < path[-1]:
                answers.append(tuple(path))
            return
        available = remaining | {v}
        degrees = [len(g[u] & available) for u in remaining]
        if 0 in degrees or degrees.count(1) > 1:
            return
        if not connected(g, available):
            return
        for u in sorted(g[v] & remaining):
            dfs(path + [u], remaining - {u})

    for v in starts:
        dfs([v], set(g) - {v})
    return answers


def unicyclic_prediction(g):
    require(connected(g) and len(edge_set(g)) == len(g), "not unicyclic")
    core = set(g)
    degree = {v: len(g[v]) for v in g}
    queue = deque(v for v in g if degree[v] == 1)
    while queue:
        v = queue.popleft()
        core.remove(v)
        for u in g[v] & core:
            degree[u] -= 1
            if degree[u] == 1:
                queue.append(u)
    if any(len(g[v]) > 2 for v in set(g) - core):
        return 0
    if any(len(g[v]) > 3 for v in core):
        return 0
    attachments = [v for v in core if len(g[v]) == 3]
    if len(attachments) == 0:
        return len(core)
    if len(attachments) == 1:
        return 2
    if len(attachments) == 2:
        return int(attachments[1] in g[attachments[0]])
    return 0


def four_cycles(g):
    """Canonical undirected cycles, found independently from common neighbours."""
    result = set()
    for a, c in combinations(sorted(g), 2):
        for b, d in combinations(sorted(g[a] & g[c]), 2):
            word = (a, b, c, d)
            rotations = [word[k:] + word[:k] for k in range(4)]
            rev = tuple(reversed(word))
            rotations += [rev[k:] + rev[:k] for k in range(4)]
            result.add(min(rotations))
    return sorted(result)


def valid_path(g, path, cycle=False):
    if len(path) != len(g) or set(path) != set(g):
        return False
    pairs = list(zip(path, path[1:]))
    if cycle:
        pairs.append((path[-1], path[0]))
    return all(v in g[u] for u, v in pairs)


def cycle_edges(path):
    return {tuple(sorted((path[i], path[(i + 1) % len(path)])))
            for i in range(len(path))}


def kernel_dimension(g):
    """Number of bipartite connected components (isolated vertices included)."""
    unseen = set(g)
    answer = 0
    while unseen:
        root = min(unseen)
        colours = {root: 0}
        queue = [root]
        bipartite = True
        while queue:
            v = queue.pop()
            unseen.discard(v)
            for u in g[v]:
                if u not in colours:
                    colours[u] = 1 - colours[v]
                    queue.append(u)
                elif colours[u] == colours[v]:
                    bipartite = False
        answer += bipartite
    return answer


def main():
    print("S1. Independent Hamiltonian census through25, paths up to reversal")
    expected = {1: 1, 15: 1, 16: 1, 17: 1, 23: 3, 25: 10}
    saved = {}
    for n in range(1, 26):
        g = square_graph(n)
        paths = hamilton_paths(g)
        require(len(paths) == expected.get(n, 0), f"count Q{n}")
        require(all(valid_path(g, p) for p in paths), "invalid path")
        saved[n] = paths
        print(n, "paths", len(paths), "leaves", [v for v in g if len(g[v]) == 1],
              "endpoints", sorted({tuple(sorted((p[0], p[-1]))) for p in paths}))
    print("Q15 path", saved[15][0])
    hostile23 = next(p for p in saved[23] if 22 in (p[0], p[-1]))
    require(len(square_graph(23)[22]) == 2, "degree2 hostile")
    print("Q23 degree2 endpoint hostile", hostile23)

    print("S2. Unicyclic classification, every labelled simple graph on3..6 vertices")
    universe = 0
    for n in range(3, 7):
        count = 0
        for edges in combinations(list(combinations(range(1, n + 1), 2)), n):
            g = {v: set() for v in range(1, n + 1)}
            for a, b in edges:
                g[a].add(b)
                g[b].add(a)
            if connected(g):
                count += 1
                require(len(hamilton_paths(g)) == unicyclic_prediction(g),
                        f"unicyclic failure {n} {edges}")
        universe += count
        print(n, "connected unicyclic labelled graphs", count)
    print("unicyclic controls", universe)
    for n in range(15, 19):
        require(unicyclic_prediction(square_graph(n)) == len(saved[n]), "Q unicyclic")

    print("S3. Triangle and C4 birth, exact induced universe Q1..Q46")
    first_triangle = None
    first_four = None
    for n in range(1, 47):
        g = square_graph(n)
        triangles = [p for p in combinations(range(1, n + 1), 3)
                     if all(v in g[u] for u, v in combinations(p, 2))]
        c4 = four_cycles(g)
        if triangles and first_triangle is None:
            first_triangle = (n, triangles)
        if c4 and first_four is None:
            first_four = (n, c4)
    require(first_triangle == (30, [(6, 19, 30)]), "triangle birth")
    require(first_four == (46, [(1, 3, 46, 35)]), "C4 birth")
    print("first triangle", first_triangle)
    print("first C4", first_four)
    pairs = defaultdict(list)
    for pair in combinations_with_replacement([k * k for k in range(2, 10)], 2):
        pairs[sum(pair)].append(pair)
    collisions = {s: pair for s, pair in pairs.items() if len(pair) > 1}
    require(collisions == {85: [(4, 81), (36, 49)]}, "square pair bank")
    print("complete equal-pair-sum bank for squares4..81", collisions)

    print("S4. First possible Hamiltonian-cycle two-edge switch, Q46")
    c1 = (1, 3, 33, 16, 9, 40, 24, 25, 39, 42, 22, 27, 37, 12, 13, 36,
          28, 21, 43, 38, 11, 14, 35, 46, 18, 31, 5, 44, 20, 29, 7, 2,
          34, 30, 6, 19, 45, 4, 32, 17, 8, 41, 23, 26, 10, 15)
    split = c1.index(35)
    c2 = (1,) + tuple(reversed(c1[1:split + 1])) + c1[split + 1:]
    g = square_graph(46)
    require(valid_path(g, c1, True) and valid_path(g, c2, True), "switch witnesses")
    removed, added = cycle_edges(c1) - cycle_edges(c2), cycle_edges(c2) - cycle_edges(c1)
    require(removed == {(1, 3), (35, 46)} and added == {(1, 35), (3, 46)}, "switch edges")
    print("cycle1", c1)
    print("cycle2", c2)
    print("removed", sorted(removed), "added", sorted(added))
    # A C4 is necessary, not sufficient, for a cycle-preserving exchange.
    c6 = (1, 2, 3, 4, 5, 6)
    bad_edges = (cycle_edges(c6) - {(1, 2), (4, 5)}) | {(1, 5), (2, 4)}
    bad = {v: set() for v in c6}
    for a, b in bad_edges:
        bad[a].add(b)
        bad[b].add(a)
    require(all(len(bad[v]) == 2 for v in bad) and not connected(bad), "C4 hostile")
    print("hostile exchange splits6-cycle into two triangles", sorted(bad_edges))

    print("S5. Label-preserving affine lifts / signless-incidence kernel")
    print("n, bipartite-component kernel dimension",
          [(n, kernel_dimension(square_graph(n))) for n in range(1, 20)])
    require(kernel_dimension(square_graph(14)) == 1, "tree offset")
    require(kernel_dimension(square_graph(15)) == 0, "odd cycle kills offset")
    p = saved[15][0]
    offsets = {v: (-1) ** i for i, v in enumerate(p)}
    lifted = tuple(4 * v + offsets[v] for v in p)
    require(len(set(lifted)) == 15 and min(lifted) > 0, "lift positivity/injectivity")
    require(all(a + b == 4 * (p[i] + p[i + 1])
                for i, (a, b) in enumerate(zip(lifted, lifted[1:]))), "edge lift")
    require(4 * (1 + 3) + offsets[1] + offsets[3] == 14, "wholegraph hostile")
    print("path square-scale4, alternating offset1", lifted)
    print("unused edge1-3 has lifted sum14; original labelled scaled sum16")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
