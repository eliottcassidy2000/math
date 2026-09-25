#!/usr/bin/env python3
"""Seed pre-screen: 'Berggren-antichain orbits'.  Along an accelerated orbit
x_0 -> x_1 -> ... of U(x)=(qx+s)/2^v, form the edge root pairs
P_i=(max(x_i,x_{i+1}), min(x_i,x_{i+1})) (skip degenerate equal pairs).
Question: for the plus map (q=3,s=+1) on positive odd starts, is any pair
P_i, P_j (i<j) Berggren-comparable (ancestor/descendant or equal)?
Controls: SHEET (q=3,s=-1: cycles 5-7 and 17-...), DRIFT (q=5,s=+1: cycles
{1,3}, {13,33,83}, {17,43,27}).  Orbits are followed until they enter a cycle
(cycle detected) or reach 1, with a step cap for divergent-looking 5n+1 orbits.
Ancestor test reused from procgen_incoming_20260925_audit_02 (own implementation)."""
import sys
from procgen_incoming_20260925_audit_02_berggren import is_ancestor


def U(n, q, s):
    m = q * n + s
    while m % 2 == 0:
        m //= 2
    return m


def orbit_edges(x, q, s, cap=400, bound=10 ** 30):
    xs = [x]
    seen = {x: 0}
    while True:
        y = U(xs[-1], q, s)
        if y in seen or len(xs) > cap or y > bound:
            xs.append(y)
            break
        seen[y] = len(xs)
        xs.append(y)
        if y == 1 and q == 3 and s == 1:
            break
    edges = []
    for a, b in zip(xs, xs[1:]):
        if a != b:
            edges.append((max(a, b), min(a, b)))
    return edges


def comparable(P, Q):
    if P == Q:
        return "equal"
    if P[0] < Q[0] and P[1] <= Q[1] and is_ancestor(P, Q):
        return "anc"
    if Q[0] < P[0] and Q[1] <= P[1] and is_ancestor(Q, P):
        return "desc"
    return None


def screen(q, s, xmax, maxdist):
    hits = {}
    examples = []
    checked = 0
    for x in range(3, xmax + 1, 2):
        E = orbit_edges(x, q, s)
        if not E:
            continue
        P = E[0]
        for d, Q in enumerate(E[1:1 + maxdist], start=1):
            checked += 1
            c = comparable(P, Q)
            if c:
                hits[(c, d)] = hits.get((c, d), 0) + 1
                if len(examples) < 12:
                    examples.append((x, d, P, Q, c))
    return checked, hits, examples


if __name__ == "__main__":
    xmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20001
    maxdist = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    for (q, s, name) in ((3, 1, "PLUS 3n+1"), (3, -1, "MINUS 3n-1 (SHEET)"), (5, 1, "5n+1 (DRIFT)")):
        checked, hits, ex = screen(q, s, xmax, maxdist)
        print("%s: odd starts 3..%d, first-edge vs later-edge pairs checked %d, comparable hits %s"
              % (name, xmax, checked, dict(sorted(hits.items())) if hits else 0), flush=True)
        for e in ex:
            print("    example start=%d distance=%d P=%s Q=%s %s" % e)
