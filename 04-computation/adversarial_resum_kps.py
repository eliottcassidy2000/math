#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adversarial_resum_kps.py

Final crux check. Reconstruct exp(sum_k b_k z^k) using the WORKER's graphical b_k
(b1=alpha1, b2=-|E|, b3=P2-t, and higher brute-force b_k) and compare to I(Omega,z).
If exp(sum b_k z^k) != I, then "log I = sum b_k z^k" (the worker's stated identity)
is FALSE -- the b_k are NOT the cumulants of log I.

We compute brute-force b_k up to k=5 for a small graph (K3, C4) and reconstruct.
"""
import itertools, math
from fractions import Fraction


def brute_ursell_bk(adj, k):
    V = len(adj)
    if k > V:
        return Fraction(0)
    pos_edges = [(i, j) for i in range(k) for j in range(i + 1, k)]
    conn_graphs = []
    for mask in range(1 << len(pos_edges)):
        edges = [pos_edges[i] for i in range(len(pos_edges)) if (mask >> i) & 1]
        if k == 1:
            conn_graphs.append(edges); continue
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


def indep_poly(adj):
    V = len(adj)
    full = frozenset(range(V)); memo = {}
    def poly(active):
        if not active:
            return (1,)
        if active in memo:
            return memo[active]
        v = max(active, key=lambda u: len(adj[u] & active))
        rest = active - {v}
        p1 = list(poly(frozenset(rest)))
        p2 = poly(frozenset(rest - adj[v]))
        out = p1 + [0]*(len(p2)+1-len(p1)) if len(p1) < len(p2)+1 else list(p1)
        for i, c in enumerate(p2):
            out[i+1] += c
        memo[active] = tuple(out); return memo[active]
    return list(poly(full))


def poly_from_series(coeffs, deg):
    """exp(sum coeffs[k-1] z^k) as power series to z^deg."""
    I = [Fraction(0)] * (deg + 1); I[0] = Fraction(1)
    for m in range(1, deg + 1):
        s = Fraction(0)
        for k in range(1, m + 1):
            if k <= len(coeffs):
                s += Fraction(k) * Fraction(coeffs[k - 1]) * I[m - k]
        I[m] = s / Fraction(m)
    return I


def mkadj(V, edges):
    adj = [set() for _ in range(V)]
    for (a, b) in edges:
        adj[a].add(b); adj[b].add(a)
    return adj


def main():
    graphs = {
        "K3": (3, [(0,1),(0,2),(1,2)]),
        "C4": (4, [(0,1),(1,2),(2,3),(3,0)]),
        "K4": (4, [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]),
    }
    for name, (V, E) in graphs.items():
        adj = mkadj(V, E)
        alpha = indep_poly(adj)
        deg = len(alpha) - 1
        bk = [brute_ursell_bk(adj, k) for k in range(1, deg + 1)]
        Irec_from_b = poly_from_series(bk, deg)
        print(f"--- {name} ---")
        print(f"  I(Omega,z) true coeffs (alpha)      = {alpha}")
        print(f"  graphical Mayer b_k (k=1..{deg})     = {[str(x) for x in bk]}")
        print(f"  exp(sum b_k z^k) reconstructed I    = {[str(x) for x in Irec_from_b]}")
        match = all(Irec_from_b[k] == Fraction(alpha[k]) for k in range(deg + 1))
        print(f"  exp(sum b_k z^k) == I(Omega,z)?     = {match}")
        print()
    print("If 'match' is FALSE, then 'log I(Omega,z) = sum_k b_k z^k' (worker's")
    print("stated derivation identity) is mathematically WRONG. The graphical b_k")
    print("are a continuum-gas cluster formula that does NOT exponentiate to the")
    print("LATTICE-gas partition function I(Omega,z). Only the analytic c_k do.")


if __name__ == "__main__":
    main()
