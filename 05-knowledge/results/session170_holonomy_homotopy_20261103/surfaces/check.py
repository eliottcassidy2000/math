#!/usr/bin/env python3
"""Exact finite controls for the S170 surface-capacity note."""
from collections import deque
from functools import lru_cache
import json
from pathlib import Path


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def compose(p, q):
    return tuple(p[q[x]] for x in range(len(p)))


@lru_cache(None)
def twist(n, a, x):
    if n == 0:
        return 0
    h = 1 << (n - 1)
    u, A = a % h, a // h
    y, t = x % h, x // h
    delta = lambda z: int(z != 0)
    return (twist(n - 1, u, y)
            ^ (t & (delta(u) ^ delta(y) ^ delta(u ^ y)))
            ^ (A & delta(y)) ^ (A & t))


def left(n, a):
    D = 1 << n
    return tuple(((s ^ twist(n, a, x)) * D + (a ^ x))
                 for s in range(2) for x in range(D))


def diagonal(n, f):
    D = 1 << n
    return tuple(((s ^ f(x)) * D + x)
                 for s in range(2) for x in range(D))


def orbit(gens, initial):
    seen = {initial}
    queue = deque([initial])
    while queue:
        x = queue.popleft()
        for g in gens:
            y = g[x]
            if y not in seen:
                seen.add(y)
                queue.append(y)
    return seen


def group_order(gens):
    identity = tuple(range(len(gens[0])))
    seen, queue = {identity}, deque([identity])
    while queue:
        p = queue.popleft()
        for g in gens:
            q = compose(g, p)
            if q not in seen:
                seen.add(q)
                queue.append(q)
    return len(seen)


def surface_generators(n):
    """All pairs commute, while their union generates the inherited G_n."""
    h = 1 << (n - 1)
    low = [left(n, 1 << j) for j in range(n - 1)]
    edge = [diagonal(n, lambda x, j=j: int(x % h in (0, 1 << j)))
            for j in range(n - 1)]
    top, dt = left(n, h), diagonal(n, lambda x: x // h)
    return [(top, edge[0]), (low[0], dt)] + list(zip(low[1:], edge[1:]))


def main():
    report = {"status": "PASS", "surface_generators": [],
              "prime_surface_capacity": []}
    for n in range(3, 9):
        pairs = surface_generators(n)
        gens = [g for pair in pairs for g in pair]
        require(len(pairs) == n, "Wrong handle count")
        require(all(compose(a, b) == compose(b, a) for a, b in pairs),
                f"Surface relator failure n={n}")
        require(len(orbit(gens, 0)) == 1 << (n + 1), "Not transitive")
        row = {"n": n, "handles": n, "sheet_count": 1 << (n + 1),
               "each_handle_commutes": True, "transitive": True}
        if n <= 4:
            order = group_order(gens)
            require(order == 1 << ((1 << (n - 1)) + n), "Wrong full group order")
            row["independently_enumerated_group_order"] = order
        report["surface_generators"].append(row)

    for p, r in [(3, 1), (3, 2), (5, 1), (5, 2), (7, 1), (11, 1), (13, 1)]:
        q = p ** r
        V = (q*q-1)//2
        m = q//p
        c = q*(q*q-1)//(p*(p*p-1))
        gp = 1+(p*p-1)*(p-6)//24
        b = (m-1)*V-c+1
        rank = 2*c*gp+b
        require(gp >= 0 and b >= 0, "Invalid normalization data")
        if q == 7:
            require((c, gp, b, rank, rank//2) == (1, 3, 0, 6, 3), "q=7 Klein quartic regression")
        report["prime_surface_capacity"].append(
            {"q": q, "p": p, "components": c, "component_genus": gp,
             "free_circles": b, "H1_rank": rank,
             "max_full_transport_n": rank//2 if rank >= 6 else None})

    output = json.dumps(report, indent=2, sort_keys=True)+"\n"
    Path(__file__).with_name("results.json").write_text(output)
    print(output, end="")


if __name__ == "__main__":
    main()
