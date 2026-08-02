#!/usr/bin/env python3
"""Exact controls for THM-3123.

The proof of THM-3123 is finite max-flow duality.  This companion checks the
orientation and all boundary conventions exhaustively in the first branching
partition posets, verifies the normalized filter identity on a genuine signed
alphabet bank, and records two minimal hostile controls.
"""

from __future__ import annotations

from collections import deque
from hashlib import sha256
from itertools import permutations, product


def partitions(n: int, cap: int | None = None):
    if n == 0:
        yield ()
        return
    if cap is None or cap > n:
        cap = n
    for first in range(cap, 0, -1):
        for tail in partitions(n - first, first):
            yield (first,) + tail


def hasse_edges(types):
    """Edges fine -> coarse, obtained by merging exactly two blocks."""
    index = {mu: i for i, mu in enumerate(types)}
    edges = set()
    for i, mu in enumerate(types):
        for a in range(len(mu)):
            for b in range(a + 1, len(mu)):
                nu = list(mu)
                merged = nu[a] + nu[b]
                nu = [nu[j] for j in range(len(nu)) if j not in (a, b)]
                nu.append(merged)
                nu = tuple(sorted(nu, reverse=True))
                edges.add((i, index[nu]))
    return tuple(sorted(edges))


def upper_closures(types, edges):
    q = len(types)
    out = [[] for _ in range(q)]
    for u, v in edges:
        out[u].append(v)
    ans = []
    for root in range(q):
        seen = {root}
        todo = [root]
        while todo:
            u = todo.pop()
            for v in out[u]:
                if v not in seen:
                    seen.add(v)
                    todo.append(v)
        ans.append(frozenset(seen))
    return tuple(ans)


def all_upsets(types, closures):
    q = len(types)
    answer = []
    for mask in range(1 << q):
        chosen = frozenset(i for i in range(q) if (mask >> i) & 1)
        if all(closures[i] <= chosen for i in chosen):
            answer.append(chosen)
    return tuple(answer)


def maxflow_feasible(c, edges):
    """Whether c is a nonnegative Hasse boundary."""
    assert sum(c) == 0
    q = len(c)
    source, sink = q, q + 1
    total = sum(max(0, -x) for x in c)
    if total == 0:
        return True

    cap = {}
    adj = [[] for _ in range(q + 2)]

    def add(u, v, z):
        if (u, v) not in cap:
            cap[u, v] = 0
            cap[v, u] = 0
            adj[u].append(v)
            adj[v].append(u)
        cap[u, v] += z

    for i, x in enumerate(c):
        if x < 0:
            add(source, i, -x)
        elif x > 0:
            add(i, sink, x)
    for u, v in edges:
        add(u, v, total)

    value = 0
    while True:
        parent = [-1] * (q + 2)
        parent[source] = source
        todo = deque([source])
        while todo and parent[sink] < 0:
            u = todo.popleft()
            for v in adj[u]:
                if parent[v] < 0 and cap[u, v] > 0:
                    parent[v] = u
                    todo.append(v)
        if parent[sink] < 0:
            break
        aug = total
        v = sink
        while v != source:
            u = parent[v]
            aug = min(aug, cap[u, v])
            v = u
        v = sink
        while v != source:
            u = parent[v]
            cap[u, v] -= aug
            cap[v, u] += aug
            v = u
        value += aug
    return value == total


def upset_feasible(c, upsets):
    return all(sum(c[i] for i in upset) >= 0 for upset in upsets)


def monomial_symmetric(mu, alphabet):
    if len(mu) > len(alphabet):
        return 0
    exponents = mu + (0,) * (len(alphabet) - len(mu))
    total = 0
    for alpha in set(permutations(exponents)):
        term = 1
        for x, power in zip(alphabet, alpha):
            term *= x**power
        total += term
    return total


def exhaustive_duality_controls():
    records = []
    for n, radius in ((2, 2), (3, 2), (4, 2), (5, 1)):
        types = tuple(partitions(n))
        edges = hasse_edges(types)
        closures = upper_closures(types, edges)
        upsets = all_upsets(types, closures)
        checked = passing = 0
        digest = sha256()
        for c in product(range(-radius, radius + 1), repeat=len(types)):
            if sum(c) != 0:
                continue
            checked += 1
            flow = maxflow_feasible(c, edges)
            dual = upset_feasible(c, upsets)
            assert flow == dual
            passing += int(flow)
            digest.update(bytes(x + radius for x in c))
            digest.update(bytes((int(flow),)))
        records.append(
            f"duality:N={n}:types={len(types)}:edges={len(edges)}:"
            f"upsets={len(upsets)}:vectors={checked}:passing={passing}:"
            f"digest={digest.hexdigest()}"
        )
    return records


def hostile_controls():
    types4 = tuple(partitions(4))
    edges4 = hasse_edges(types4)
    closure4 = upper_closures(types4, edges4)
    upsets4 = all_upsets(types4, closure4)
    i4 = {mu: i for i, mu in enumerate(types4)}
    c4 = [0] * len(types4)
    for mu, value in {
        (4,): 1,
        (3, 1): -1,
        (2, 2): -1,
        (2, 1, 1): 1,
    }.items():
        c4[i4[mu]] = value
    principal = [sum(c4[j] for j in closure4[i]) for i in range(len(types4))]
    bad = frozenset(i4[mu] for mu in ((4,), (3, 1), (2, 2)))
    assert min(principal) == 0
    assert sum(c4[i] for i in bad) == -1
    assert bad in upsets4
    assert not maxflow_feasible(c4, edges4)

    types3 = tuple(partitions(3))
    edges3 = hasse_edges(types3)
    c3 = [-1, 3, -2]  # order: (3),(2,1),(1,1,1)
    assert types3 == ((3,), (2, 1), (1, 1, 1))
    assert not maxflow_feasible(c3, edges3)

    return (
        "principal_hostile:N=4:vector=(1,-1,-1,1,0):"
        "principal_min=0:nonprincipal_upset=-1:PASS",
        "kernel_hostile:N=3:vector=(-1,3,-2):coarsest_upset=-1:PASS",
    )


def filter_response_control():
    n = 4
    types = tuple(partitions(n))
    edges = hasse_edges(types)
    upsets = all_upsets(types, upper_closures(types, edges))
    q_alphabet = (1, 2, 3, 4)
    signed_bank = (
        (2, (1, 2, 3)),
        (-1, (1, 1, 4)),
        (3, (2, 5)),
    )
    q = {mu: monomial_symmetric(mu, q_alphabet) for mu in types}
    phi = {
        mu: sum(eps * monomial_symmetric(mu, alphabet) for eps, alphabet in signed_bank)
        for mu in types
    }
    hq = sum(q.values())
    phih = sum(phi.values())
    g = {mu: phih * q[mu] - phi[mu] * hq for mu in types}
    assert sum(g.values()) == 0
    digest = sha256()
    for upset in upsets:
        lhs = sum(g[types[i]] for i in upset)
        fq = sum(q[types[i]] for i in upset)
        phif = sum(phi[types[i]] for i in upset)
        rhs = phih * fq - phif * hq
        assert lhs == rhs
        digest.update(f"{sorted(upset)}:{lhs};".encode())
    return (
        f"filter_identity:N=4:upsets={len(upsets)}:hq={hq}:phih={phih}:"
        f"digest={digest.hexdigest()}:PASS"
    )


def main():
    for line in exhaustive_duality_controls():
        print(line)
    for line in hostile_controls():
        print(line)
    print(filter_response_control())
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
