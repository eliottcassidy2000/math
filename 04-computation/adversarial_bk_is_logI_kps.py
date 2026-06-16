#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adversarial_bk_is_logI_kps.py

THE CRUX. The worker asserts (docstring + derivation):
    log I(Omega,z) = sum_{k>=1} b_k z^k     with b_k = MAYER connected-cluster integral
                                            b_1=alpha_1, b_2=-|E|, b_3=P2-t.
But earlier I found Taylor coeffs of log I are c_k (c_2 = -|E|-alpha_1/2 != b_2).

So EITHER:
  (i) the graphical b_k do NOT sum to log I  (worker's labeling is wrong), or
  (ii) they DO and my Taylor c_k are something else.

I settle this by directly comparing, for a concrete graph, the brute-force Mayer
connected-cluster b_k against the Taylor c_k of log I, for k=1,2,3.

The standard statistical-mechanics fact is:
   log Xi(z) = sum_{k>=1} b_k z^k     where  b_k = (1/k!) sum over connected
   clusters of k DISTINCT particles of the product of Mayer f-functions.
The DISTINCT-particle restriction is exactly the difference. Let me check whether
the brute-force "distinct ordered tuples" b_k reproduce c_k (the true log I Taylor)
or the worker's -|E| values.
"""
import itertools, math
from fractions import Fraction


def brute_ursell_bk(adj, k):
    V = len(adj)
    pos_edges = [(i, j) for i in range(k) for j in range(i + 1, k)]
    conn_graphs = []
    for mask in range(1 << len(pos_edges)):
        edges = [pos_edges[i] for i in range(len(pos_edges)) if (mask >> i) & 1]
        if k == 1:
            conn_graphs.append(edges)
            continue
        par = list(range(k))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]; x = par[x]
            return x
        for (a, b) in edges:
            par[find(a)] = find(b)
        if len({find(i) for i in range(k)}) == 1:
            conn_graphs.append(edges)
    total = Fraction(0)
    for tup in itertools.permutations(range(V), k):
        for edges in conn_graphs:
            if all(tup[b] in adj[tup[a]] for (a, b) in edges):
                total += Fraction((-1) ** len(edges))
    return total / math.factorial(k)


def taylor_log_poly(alpha, K):
    a = [Fraction(c) for c in alpha] + [Fraction(0)] * (K + 2)
    d = [Fraction(0)] * (K + 1)
    for m in range(1, K + 1):
        s = Fraction(m) * a[m]
        for k in range(1, m):
            s -= Fraction(k) * d[k] * a[m - k]
        d[m] = s / Fraction(m)
    return d[1:K + 1]


def indep_poly(adj):
    V = len(adj)
    full = frozenset(range(V))
    memo = {}
    def poly(active):
        if not active:
            return (1,)
        if active in memo:
            return memo[active]
        v = max(active, key=lambda u: len(adj[u] & active))
        rest = active - {v}
        p1 = list(poly(frozenset(rest)))
        p2 = poly(frozenset(rest - adj[v]))
        out = p1 + [0] * (len(p2) + 1 - len(p1)) if len(p1) < len(p2) + 1 else list(p1)
        for i, c in enumerate(p2):
            out[i + 1] += c
        memo[active] = tuple(out)
        return memo[active]
    return list(poly(full))


def mkadj(V, edges):
    adj = [set() for _ in range(V)]
    for (a, b) in edges:
        adj[a].add(b); adj[b].add(a)
    return adj


def main():
    # Test graphs
    graphs = {
        "K3": (3, [(0,1),(0,2),(1,2)]),
        "P3 path": (3, [(0,1),(1,2)]),
        "single edge + isolated": (3, [(0,1)]),
        "C4": (4, [(0,1),(1,2),(2,3),(3,0)]),
        "K4": (4, [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]),
        "paw (triangle+pendant)": (4, [(0,1),(1,2),(2,0),(2,3)]),
    }
    print(f"{'graph':26s} {'k':>2s} {'brute b_k':>10s} {'Taylor c_k':>12s} {'b_k==c_k?':>10s}")
    for name, (V, E) in graphs.items():
        adj = mkadj(V, E)
        alpha = indep_poly(adj)
        c = taylor_log_poly(alpha, 4)
        for k in (1, 2, 3):
            bk = brute_ursell_bk(adj, k)
            ck = c[k - 1]
            print(f"{name:26s} {k:2d} {str(bk):>10s} {str(ck):>12s} {str(bk==ck):>10s}")
        print()
    print("VERDICT LOGIC:")
    print(" If brute b_k == Taylor c_k  => the Mayer cluster series IS log I, and the")
    print("    worker's 'b_2 = -|E|' is then a DIFFERENT object (NOT the log-I coeff).")
    print(" If brute b_k == worker's -|E| (i.e. != c_k) => the worker's b_k is the")
    print("    cluster integral but does NOT sum to log I; 'log I = sum b_k z^k' is wrong.")


if __name__ == "__main__":
    main()
