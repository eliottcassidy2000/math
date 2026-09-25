"""Exact controls for collatz_bugs_20260925.md; standard library only."""
from collections import deque
from math import gcd


def T(n):
    return n // 2 if n % 2 == 0 else (3 * n + 1) // 2


def oddpart(n):
    while n % 2 == 0:
        n //= 2
    return n


def path_to_four(n):
    path = [n]
    while path[-1] != 4:
        assert len(path) < 10000 and path[-1] not in (1, 2)
        nxt = T(path[-1])
        assert nxt not in path
        path.append(nxt)
    return path


def meet_distance(px, py):
    iy = {v: i for i, v in enumerate(py)}
    for i, v in enumerate(px):
        if v in iy:
            return i, iy[v], v
    raise AssertionError("Missing meet")


def zero_one_bfs(graph, x):
    dist, queue = {x: 0}, deque([x])
    while queue:
        v = queue.popleft()
        for w, cost in graph[v]:
            candidate = dist[v] + cost
            if candidate < dist.get(w, float("inf")):
                dist[w] = candidate
                (queue.append if cost else queue.appendleft)(w)
    return dist


def main():
    for I in range(2, 5001, 2):
        legal = (I - 1) % 3 == 0
        assert legal == (I % 6 == 4)
        if legal:
            a, b = (I - 1) // 3, (4 * I - 1) // 3
            assert a % 2 == b % 2 == 1
            assert T(a) == I // 2 and T(b) == 2 * I
            assert b == 4 * a + 1
            assert oddpart(3 * a + 1) == oddpart(3 * b + 1)
            assert gcd(a, b) == 1 and (4 * I) % 6 == 4
    assert (8 - 1) % 3 != 0  # hostile: a non-bug even identity
    print("Admissibility: all 2500 even I in [2,5000]; I=8 rejected.")

    for a in range(1, 200, 2):
        ladder = [a]
        for _ in range(6):
            ladder.append(4 * ladder[-1] + 1)
        for j, b in enumerate(ladder):
            assert 3 * b + 1 == 4**j * (3 * a + 1)
            for k in range(j + 1, len(ladder)):
                assert ((4 ** (k - j) - 1) // 3) % gcd(b, ladder[k]) == 0
    assert gcd(5, 85) == 5  # hostile: nonconsecutive rungs need not be coprime
    print("Prime overlap: 100 odd seeds, 7 rungs each; nonadjacent hostile 5,85.")

    ids = list(range(4, 299, 6))
    paths = {I: path_to_four(I) for I in ids}
    graph = {}
    for path in paths.values():
        for v in path:
            graph.setdefault(v, set())
        for child, parent in zip(path, path[1:]):
            graph[child].add((parent, 1))  # against inverse-Collatz arrow
            graph[parent].add((child, 0))
    bfs = {I: zero_one_bfs(graph, I) for I in ids}
    for X in ids:
        for Y in ids:
            a, b, c = meet_distance(paths[X], paths[Y])
            assert bfs[X][Y] == a and bfs[Y][X] == b
            assert a - b == len(paths[X]) - len(paths[Y])
            assert (a == 0 or b == 0) == (X in paths[Y] or Y in paths[X])
            if a and b:
                assert c % 3 == 2
    assert meet_distance(paths[10], paths[16]) == (2, 1, 8)
    assert T(1) == 2 and T(2) == 1  # hostile to full-graph tree terminology
    print(f"Distance: {len(ids)**2} ordered pairs, {len(graph)}-vertex orbit-union tree;")
    print("  LCA formula agrees with independent 0-1 BFS; example (10,16)=(2,1), meet=8.")
    print("  Root-cycle hostile 1<->2 retained; domain excludes 1,2.")
    print("FINITE-EXACT only for the stated universes; no global Collatz test claimed.")


if __name__ == "__main__":
    main()
