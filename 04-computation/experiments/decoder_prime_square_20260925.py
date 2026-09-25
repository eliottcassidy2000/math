"""Exact square-sum halving and endpoint-aware Q24 certificate.

Pure Python; no external data. Run from repository root. Every check is an
explicit exception rather than assert, so -O does not weaken verification.
"""
from collections import Counter
from math import isqrt


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def square(n):
    return n >= 0 and isqrt(n) ** 2 == n


def graph(n, sheet=1):
    return {x: {y for y in range(1, n + 1)
                if x != y and (x + y) % sheet == 0
                and square((x + y) // sheet)}
            for x in range(1, n + 1)}


def graph_from_squares(n):
    adj = {x: set() for x in range(1, n + 1)}
    for s in range(2, isqrt(2 * n - 1) + 1):
        for x in range(max(1, s * s - n), min(n, (s * s - 1) // 2) + 1):
            y = s * s - x
            adj[x].add(y)
            adj[y].add(x)
    return adj


def v2(n):
    return (n & -n).bit_length() - 1


def bracket(n):
    if n == 1:
        return 0
    m = 1
    while (2 * m + 1) ** 2 < n:
        m += 1
    return m


def prime(n):
    return n >= 2 and all(n % p for p in range(2, isqrt(n) + 1))


def endpoint_propagation(adj, ends):
    """Only sound degree equations and a proper-cycle contradiction.

    The endpoints are INPUT, never inferred from host degree 2.
    Returns certificate events, or None if these rules cannot decide.
    """
    edges = {tuple(sorted((x, y))) for x in adj for y in adj[x]}
    chosen, excluded, events = set(), set(), []
    incident = {v: {e for e in edges if v in e} for v in adj}
    target = {v: 1 if v in ends else 2 for v in adj}
    while True:
        changed = False
        for v in sorted(adj):
            yes = incident[v] & chosen
            avail = incident[v] - excluded
            if len(yes) > target[v] or len(avail) < target[v]:
                events.append(("degree contradiction", v, len(yes),
                               len(avail), target[v]))
                return events
            if len(yes) == target[v]:
                new = avail - yes
                if new:
                    excluded |= new
                    events.append(("exclude", v, sorted(new)))
                    changed = True
            if len(avail) == target[v]:
                new = avail - chosen
                if new:
                    chosen |= new
                    events.append(("include", v, sorted(new)))
                    changed = True
        # A cycle among already mandatory edges cannot lie in a path.
        remaining = set(adj)
        while remaining:
            todo = [remaining.pop()]
            component = set(todo)
            while todo:
                v = todo.pop()
                nxt = {w for e in incident[v] & chosen for w in e}
                nxt &= remaining
                remaining -= nxt
                component |= nxt
                todo.extend(nxt)
            if len(component) >= 3 and all(
                    len(incident[v] & chosen) == 2 for v in component):
                events.append(("forced cycle", sorted(component)))
                return events
        if not changed:
            return None


def paths_from_leaf(adj):
    """Independent plain vertex-DFS; no degree-two forcing/pruning.

    Every visited prefix is a simple path. The single prune checks that no
    remaining vertex is isolated in remaining vertices plus current end.
    Starting at a fixed leaf counts undirected Hamiltonian paths once.
    """
    n = len(adj)
    leaf = min(v for v in adj if len(adj[v]) == 1)
    masks = {v: sum(1 << (u - 1) for u in adj[v]) for v in adj}
    counts, witnesses = Counter(), {}

    def visit(v, remaining, path):
        if not remaining:
            counts[v] += 1
            witnesses.setdefault(v, path)
            return
        choices = masks[v] & remaining
        while choices:
            bit = choices & -choices
            choices -= bit
            u = bit.bit_length()
            rest = remaining ^ bit
            available = rest | bit
            probe = rest
            valid = True
            while probe:
                zbit = probe & -probe
                probe -= zbit
                if not masks[zbit.bit_length()] & available:
                    valid = False
                    break
            if valid:
                visit(u, rest, path + [u])

    visit(leaf, ((1 << n) - 1) ^ (1 << (leaf - 1)), [leaf])
    return counts, witnesses


def main():
    print("decoder_prime_square_20260925: exact local controls")
    for n in range(1, 129):
        for sheet in (1, 2):
            source = graph(2 * n, sheet)
            target = graph(n, 3 - sheet)
            for x in target:
                check({y // 2 for y in source[2 * x] if y % 2 == 0}
                      == target[x], ("halving", n, sheet, x))
    print("two-sheet halving n=1..128, both sheets: PASS")
    adj = graph(256)
    check(adj == graph_from_squares(256), "independent edge generator")
    for x in adj:
        for y in adj[x]:
            if v2(x) != v2(y):
                check(min(v2(x), v2(y)) % 2 == 0,
                      ("dyadic layer condition", x, y))
    # 1+3=4 is square, but 2+6=8 is not; 4+12=16 is square.
    check(3 in graph(3)[1] and 6 not in graph(6)[2]
          and 12 in graph(12)[4], "doubling hostile / quadrupling control")
    print("dyadic layers x,y<=256; independent graph constructor: PASS")
    print("doubling hostile 1--3 -> 2--6 (8 nonsquare); x4 -> 4--12: PASS")
    eligible = [n for n in range(2, 10001) if bracket(n) == bracket(2 * n)]
    check(eligible == [2, 3, 4, 10, 11, 12], "escape set")
    print("doubling same odd-square bracket, n<=10000:", eligible)
    print("prime members:", [n for n in eligible if prime(n)])
    print("Q24 leaf / Q31 degree horizon:", sorted(graph(30)[18]),
          sorted(graph(31)[18]))
    check(graph(30)[18] == {7} and graph(31)[18] == {7, 31}, "18 leaf")
    check(min(map(len, graph(31).values())) == 2, "Q31 min degree")

    a24 = graph(24)
    check(a24 == graph_from_squares(24), "Q24 adjacency")
    check(a24[18] == {7}, "forced endpoint")
    for e in range(1, 25):
        if e == 18:
            continue
        certificate = endpoint_propagation(a24, {18, e})
        check(certificate is not None, ("Q24 certificate incomplete", e))
        print("Q24 endpoint", e, "certificate:", certificate)
    tri = {1: {2, 3}, 2: {1, 3}, 3: {1, 2}}
    check(endpoint_propagation(tri, {1, 2}) is None, "triangle hostile")
    check(endpoint_propagation(graph(23), {18, 22}) is None,
          "corrected Q23 endpoint hostile")
    print("proper endpoint controls: triangle and Q23 endpoint22 PASS")

    expected = {14: 0, 15: 1, 16: 1, 17: 1, 18: 0, 19: 0,
                20: 0, 21: 0, 22: 0, 23: 3, 24: 0, 25: 10}
    for n in range(14, 26):
        counts, witnesses = paths_from_leaf(graph_from_squares(n))
        check(sum(counts.values()) == expected[n], ("plain DFS", n, counts))
        print("independent plain DFS", n, "paths", sum(counts.values()),
              "other endpoints", sorted(counts.items()))
        if n in (23, 25):
            for path in witnesses.values():
                check(set(path) == set(range(1, n + 1))
                      and len(path) == n
                      and all(square(x + y) for x, y in zip(path, path[1:])),
                      "path witness")
            print("witness", n, witnesses[min(witnesses)])
    check(graph(25)[25] == {11, 24}, "Q25 connector")
    print("Q25 new vertex25 neighbours = [11,24]; sums36,49")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
