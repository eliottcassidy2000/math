#!/usr/bin/env python3
"""Construct an embedded black paired-quotient 4-cycle at every n>=8.

All verification is literal on labelled tournaments. No isomorphism routine is
needed: the four distinct cyclic-triangle counts distinguish merged nodes.
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def schema(n):
    return [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1)]


def tournament(mask, n):
    adj = [0] * n
    for i in range(n - 1):
        adj[i] |= 1 << (i + 1)
    for k, (x, y) in enumerate(schema(n)):
        i, j = n - x, n - y
        if mask >> k & 1:
            adj[j] |= 1 << i
        else:
            adj[i] |= 1 << j
    return adj


def edge_set(path):
    return set(zip(path, path[1:]))


def reverse_path(adj, path):
    result = adj[:]
    for u, v in zip(path, path[1:]):
        assert adj[u] >> v & 1
        result[u] ^= 1 << v
        result[v] ^= 1 << u
    return result


def triangle_count(adj):
    n = len(adj)
    return n * (n - 1) * (n - 2) // 6 - sum(
        (d := row.bit_count()) * (d - 1) // 2 for row in adj
    )


def is_blue(adj, path):
    reversal = {path[i]: path[-1 - i] for i in range(len(path))}
    return all(
        ((adj[u] >> v) & 1) != ((adj[reversal[u]] >> reversal[v]) & 1)
        for u in range(len(adj)) for v in range(u + 1, len(adj))
    )


def validate(adj, p, q):
    n = len(adj)
    assert sorted(p) == sorted(q) == list(range(n))
    assert len({p[0], p[-1], q[0], q[-1]}) == 4
    assert edge_set(p).isdisjoint(edge_set(q))
    assert all(adj[u] >> v & 1 for u, v in edge_set(p) | edge_set(q))
    for u in range(n):
        assert not (adj[u] >> u & 1)
        for v in range(u):
            assert (adj[u] >> v & 1) != (adj[v] >> u & 1)
    scores = [row.bit_count() for row in adj]
    a = scores[p[0]] - scores[p[-1]] - 1
    b = scores[q[0]] - scores[q[-1]] - 1
    assert a > 0 and b > 0 and a != b
    assert a % 2 == b % 2 == (n - 1) % 2
    tp = reverse_path(adj, p)
    tq = reverse_path(adj, q)
    tpq = reverse_path(tp, q)
    assert tpq == reverse_path(tq, p)
    c = triangle_count(adj)
    actual = [triangle_count(z) for z in [adj, tp, tpq, tq]]
    assert actual == [c, c + a, c + a + b, c + b]
    assert len(set(actual)) == 4
    assert not any([is_blue(adj, p), is_blue(tp, q),
                    is_blue(tq, p), is_blue(adj, q)])
    return {"n": n, "a": a, "b": b, "triangle_counts": actual}


def insert_four(adj, p, q):
    n = len(adj)
    endpoints = {p[0], p[-1], q[0], q[-1]}
    choice = next(
        (a, b, c, d)
        for a, b in zip(p, p[1:])
        for c, d in zip(q, q[1:])
        if len({a, b, c, d}) == 4 and not ({a, b, c, d} & endpoints)
    )
    a, b, c, d = choice
    x0, x1, x2, x3 = range(n, n + 4)
    result = adj[:] + [0] * 4

    def orient(u, v):
        result[u] |= 1 << v
        result[v] &= ~(1 << u)

    # A default orientation for cross arcs, then the explicitly prescribed arcs.
    for old in range(n):
        for new in range(n, n + 4):
            orient(old, new)
    for start in [p[0], q[0]]:
        for new in range(n, n + 4):
            orient(start, new)
    for end in [p[-1], q[-1]]:
        for new in range(n, n + 4):
            orient(new, end)
    for u, v in [(x0, x1), (x1, x2), (x2, x3),
                 (x2, x0), (x0, x3), (x3, x1),
                 (a, x0), (x3, b), (c, x2), (x1, d)]:
        orient(u, v)
    pp = p[:p.index(a) + 1] + [x0, x1, x2, x3] + p[p.index(b):]
    qq = q[:q.index(c) + 1] + [x2, x0, x3, x1] + q[q.index(d):]
    return result, pp, qq


def construct(n):
    if n < 8:
        raise ValueError("The uniform construction starts at order eight.")
    base = 8 + ((n - 8) % 4)
    seed = json.loads((HERE / f"score_square_seed_n{base}.json").read_text())
    adj, p, q = tournament(seed["t"], base), seed["P"], seed["Q"]
    while len(adj) < n:
        adj, p, q = insert_four(adj, p, q)
    return adj, p, q


def main():
    rows = []
    for n in list(range(8, 101)) + [128, 129, 130, 131, 200, 201, 202, 203]:
        adj, p, q = construct(n)
        rows.append(validate(adj, p, q))
    result = {"status": "all literal checks passed", "checked_orders": len(rows), "rows": rows}
    (HERE / "black_quadrilateral_construction_checks.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}))


if __name__ == "__main__":
    main()
