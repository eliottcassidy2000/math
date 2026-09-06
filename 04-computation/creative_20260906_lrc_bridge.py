"""Exact divisor-component compiler for walk cancellation; standalone controls."""
from itertools import combinations, product
from math import gcd, lcm
from fractions import Fraction
from functools import reduce
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise ArithmeticError(message)


def component(adj, allowed, root):
    if root not in allowed:
        return frozenset()
    seen = {root}
    todo = [root]
    while todo:
        for v in adj[todo.pop()]:
            if v in allowed and v not in seen:
                seen.add(v)
                todo.append(v)
    return frozenset(seen)


def subset_spectrum(values, adj, x, y):
    """Independent baseline: all connected vertex subsets with both endpoints."""
    n = len(values)
    D = gcd(values[x], values[y])
    ans = {}
    for bits in range(1 << n):
        if not ((bits >> x) & 1 and (bits >> y) & 1):
            continue
        S = frozenset(i for i in range(n) if (bits >> i) & 1)
        if component(adj, S, x) != S:
            continue
        d = reduce(gcd, (values[i] for i in S))
        J = D // d
        ans[J] = max(ans.get(J, 0), len(S))
    return ans


def divisor_spectrum(values, adj, x, y):
    """Divisor/component implementation, independent of connected-subset search."""
    D = gcd(values[x], values[y])
    ans = {}
    for h in range(1, D + 1):
        if D % h:
            continue
        S = component(adj, {i for i, w in enumerate(values) if w % h == 0}, x)
        if y in S and reduce(gcd, (values[i] for i in S)) == h:
            ans[D // h] = len(S)
    return ans


def factor_free_spectrum(values, adj, x, y):
    """At most 2**(n-2) generated gcd states; no integer factorization."""
    D = gcd(values[x], values[y])
    states = {D}
    for i, w in enumerate(values):
        if i not in (x, y):
            states |= {gcd(h, w) for h in states}
    ans = {}
    closures = {}
    for h in states:
        S = component(adj, {i for i, w in enumerate(values) if w % h == 0}, x)
        if y in S and reduce(gcd, (values[i] for i in S)) == h:
            ans[D // h] = len(S)
            closures[D // h] = S
    return ans, closures


def tree_walk(adj, S, x, y):
    parent = {x: None}
    queue = [x]
    for u in queue:
        for v in sorted(adj[u]):
            if v in S and v not in parent:
                parent[v] = u
                queue.append(v)
    children = {u: [] for u in S}
    for v, u in parent.items():
        if u is not None:
            children[u].append(v)
    path = [y]
    while path[-1] != x:
        path.append(parent[path[-1]])
    path.reverse()
    forward = dict(zip(path, path[1:]))
    walk = [x]
    def traverse(u, finish):
        for v in children[u]:
            if finish and forward.get(u) == v:
                continue
            walk.append(v)
            traverse(v, False)
            walk.append(u)
        if finish and u != y:
            v = forward[u]
            walk.append(v)
            traverse(v, True)
    traverse(x, True)
    return walk


def packet(values):
    A = B = L = 1
    ratio = Fraction(1)
    for u, v in zip(values, values[1:]):
        d = gcd(u, v)
        A *= u // d
        B *= v // d
        ratio *= Fraction(v, u)
        L = lcm(L, ratio.denominator)
    return L // (A // gcd(A, B))


def graph(n, edges):
    adj = [set() for _ in range(n)]
    for u, v in edges:
        adj[u].add(v)
        adj[v].add(u)
    return adj


def check(values, adj, x, y):
    direct = subset_spectrum(values, adj, x, y)
    by_divisor = divisor_spectrum(values, adj, x, y)
    compiled, closures = factor_free_spectrum(values, adj, x, y)
    need(direct == by_divisor == compiled, "exact spectrum and maximum support agreement")
    need(len(compiled) <= 2 ** (len(values) - 2), "height-independent state bound")
    for J, S in closures.items():
        walk = tree_walk(adj, S, x, y)
        need(walk[0] == x and walk[-1] == y and set(walk) == set(S), "witness endpoints and support")
        need(all(v in adj[u] for u, v in zip(walk, walk[1:])), "literal graph walk")
        need(len(walk) - 1 <= 2 * len(S) - 3, "bounded witness length")
        need(packet([values[i] for i in walk]) == J, "independent rational walk packet")
    return len(compiled)


def main():
    cases = 0
    edges4 = list(combinations(range(4), 2))
    for edge_bits in range(1 << len(edges4)):
        adj = graph(4, [e for i, e in enumerate(edges4) if (edge_bits >> i) & 1])
        for values in product(range(1, 5), repeat=4):
            for x, y in combinations(range(4), 2):
                check(values, adj, x, y)
                cases += 1
    families = [
        [(0, i) for i in range(1, 5)],
        [(i, i + 1) for i in range(4)],
        [(i, (i + 1) % 5) for i in range(5)],
        list(combinations(range(5), 2)),
    ]
    for edges in families:
        adj = graph(5, edges)
        for values in combinations(range(1, 9), 5):
            for x, y in combinations(range(5), 2):
                check(values, adj, x, y)
                cases += 1
    hostile_values = (6, 2, 3)
    hostile_graph = graph(3, [(0, 1), (1, 2)])
    hostile = factor_free_spectrum(hostile_values, hostile_graph, 0, 2)[0]
    need(hostile == {3: 3}, "actual atlas hostile requires cancellation three")
    star_values = (1, 2, 3, 5, 7)
    star_graph = graph(5, [(0, i) for i in range(1, 5)])
    star, closures = factor_free_spectrum(star_values, star_graph, 4, 3)
    walk = tree_walk(star_graph, closures[1], 4, 3)
    need(star == {1: 5} and len(walk) - 1 == 6, "repeated walk sees five vertices; every simple endpoint path sees three")
    scaled = tuple(10**80 * v for v in hostile_values)
    need(factor_free_spectrum(scaled, hostile_graph, 0, 2)[0] == hostile, "common scaling covariance at unbounded integer height")
    print("STATUS: PASS; exact cancellation spectra and bounded walk certificates")
    print("FINITE UNIVERSE:", cases, "labelled graph/endpoint cases; all 4-vertex graphs, labels 1..4, plus four 5-vertex graph families with distinct labels 1..8")
    print("HOSTILE:", json.dumps({"values": hostile_values, "spectrum": hostile}, sort_keys=True))
    print("REPEATED-WALK CONTROL:", json.dumps({"values": star_values, "spectrum": star, "walk_values": [star_values[i] for i in walk]}, sort_keys=True))
    print("LRC SCOPE: exact compiler for inherited conditional endpoint-walk test; universal qualifying component existence remains OPEN")
    print("ACTIVE GATES:", GATES)


if __name__ == "__main__":
    main()
