#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adversarial_mayer_check_kps.py  (kind-pasteur ADVERSARIAL VERIFIER, 2026-06-15)

Independent re-derivation of the mayer-ocf worker's claims. I do NOT import the
worker's code. Everything below is written fresh:

  C1. OCF: H(T) = I(Omega,2) = #directed Hamiltonian paths.
  C2. b_2 = -|E(Omega)| = alpha_2 - C(alpha_1,2).
  C3. The Mayer/Ursell b_k computed by BRUTE FORCE over connected subgraphs
      (NOT by closed form) gives b_1=alpha_1, b_2=-|E|, b_3=P2-t.
  C4. exp(sum_k c_k z^k) == I(Omega,z) as a formal power series.
  C5. "log H = sum_k c_k 2^k" literal numeric sum DIVERGES (worker's caveat).
  C6. H_free = 3^alpha_1 >= H, equality iff Omega edgeless.
  C7. Forbidden: I(K3,2)=7; 7 and 21 absent from achievable H at n<=6.
"""

import itertools, math, random
from fractions import Fraction


def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    for bits in range(1 << m):
        A = [[0] * n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def odd_cycles(A, n):
    """All directed odd cycles as (canonical_seq, vertexset). Independent impl:
    enumerate by choosing vertex subset then directed Hamiltonian cycle of that subset."""
    out = []
    seen = set()
    for L in range(3, n + 1, 2):
        for verts in itertools.combinations(range(n), L):
            # directed cycles on this exact vertex set, using each vertex
            mn = verts[0]  # smallest, fix as anchor to avoid rotation dups
            others = [v for v in verts if v != mn]
            for perm in itertools.permutations(others):
                seq = (mn,) + perm
                ok = all(A[seq[t]][seq[(t + 1) % L]] == 1 for t in range(L))
                if ok:
                    # canonical: rotation already anchored at min; but reflection is a
                    # DIFFERENT directed cycle (opposite orientation) so do NOT dedupe it.
                    if seq not in seen:
                        seen.add(seq)
                        out.append((seq, frozenset(verts)))
    return out


def conflict_graph(cycles):
    V = len(cycles)
    adj = [set() for _ in range(V)]
    for a in range(V):
        for b in range(a + 1, V):
            if cycles[a][1] & cycles[b][1]:
                adj[a].add(b)
                adj[b].add(a)
    return adj


def indep_poly(adj):
    """alpha[k] = #independent sets of size k. Brute force over all subsets for
    small graphs; pivot recursion for larger. Use recursion (exact)."""
    V = len(adj)
    import sys
    sys.setrecursionlimit(100000)
    from functools import lru_cache
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


def H_hampaths(A, n):
    adj = [[j for j in range(n) if A[i][j]] for i in range(n)]
    cnt = 0
    def dfs(v, vis, depth):
        nonlocal cnt
        if depth == n:
            cnt += 1
            return
        for w in adj[v]:
            if not (vis >> w) & 1:
                dfs(w, vis | (1 << w), depth + 1)
    for s in range(n):
        dfs(s, 1 << s, 1)
    return cnt


def brute_ursell_bk(adj, k):
    """COMPLETELY INDEPENDENT brute-force Ursell b_k for the hard-core gas.
    b_k = (1/k!) * sum_{ordered k-tuples of DISTINCT vertices}
                   sum_{connected graphs G on those k labelled positions,
                        with every G-edge actually present in Omega} (-1)^{|E(G)|}.
    This is the literal Mayer connected-cluster definition with f=-1 on edges."""
    V = len(adj)
    # enumerate all connected labelled graphs on k nodes once (as edge subsets of K_k)
    pos_edges = [(i, j) for i in range(k) for j in range(i + 1, k)]
    conn_graphs = []  # list of (edge_list)
    for mask in range(1 << len(pos_edges)):
        edges = [pos_edges[i] for i in range(len(pos_edges)) if (mask >> i) & 1]
        # connectivity check on k nodes
        if k == 1:
            conn_graphs.append(edges)
            continue
        par = list(range(k))
        def find(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for (a, b) in edges:
            par[find(a)] = find(b)
        roots = {find(i) for i in range(k)}
        if len(roots) == 1:
            conn_graphs.append(edges)
    total = Fraction(0)
    for tup in itertools.permutations(range(V), k):
        for edges in conn_graphs:
            # check each graph-edge maps to an Omega edge
            ok = True
            for (a, b) in edges:
                if tup[b] not in adj[tup[a]]:
                    ok = False
                    break
            if ok:
                total += Fraction((-1) ** len(edges))
    return total / math.factorial(k)


def features(adj):
    V = len(adj)
    E = sum(len(adj[v]) for v in range(V)) // 2
    P2 = sum(len(adj[v]) * (len(adj[v]) - 1) // 2 for v in range(V))
    t = 0
    for a in range(V):
        for b in adj[a]:
            if b <= a:
                continue
            for c in adj[b] & adj[a]:
                if c > b:
                    t += 1
    return E, P2, t


def formal_cumulants(alpha, K):
    a = [Fraction(c) for c in alpha] + [Fraction(0)] * (K + 2)
    c = [Fraction(0)] * (K + 1)
    for m in range(1, K + 1):
        s = Fraction(m) * a[m]
        for k in range(1, m):
            s -= Fraction(k) * c[k] * a[m - k]
        c[m] = s / Fraction(m)
    return c[1:K + 1]


def poly_from_cumulants(c, deg):
    I = [Fraction(0)] * (deg + 1)
    I[0] = Fraction(1)
    for m in range(1, deg + 1):
        s = Fraction(0)
        for k in range(1, m + 1):
            if k <= len(c):
                s += Fraction(k) * Fraction(c[k - 1]) * I[m - k]
        I[m] = s / Fraction(m)
    return I


def main():
    print("=" * 70)
    print("ADVERSARIAL INDEPENDENT CHECK of mayer-ocf claims")
    print("=" * 70)

    fails = dict(ocf=0, b2_negE=0, b2_alpha=0, brute_b1=0, brute_b2=0,
                 brute_b3=0, exp_cum=0, free=0, free_iff=0)
    cnt = 0
    achievable_H = {}  # n -> set of H values
    for n in (3, 4, 5, 6):
        Hset = set()
        for A in all_tournaments(n):
            cnt += 1
            cyc = odd_cycles(A, n)
            adj = conflict_graph(cyc)
            alpha = indep_poly(adj)
            a1 = alpha[1] if len(alpha) > 1 else 0
            a2 = alpha[2] if len(alpha) > 2 else 0
            E, P2, t = features(adj)
            H_ocf = sum(alpha[k] * 2 ** k for k in range(len(alpha)))
            H_dir = H_hampaths(A, n)
            Hset.add(H_ocf)
            if H_ocf != H_dir:
                fails['ocf'] += 1
            if Fraction(-E) != Fraction(a2) - Fraction(a1 * (a1 - 1), 2):
                fails['b2_alpha'] += 1
            # brute-force Ursell on a SUBSET (expensive) - only do for n<=5 fully,
            # and a sample at n=6 (graphs can be large)
            do_brute = (n <= 5)
            if do_brute and len(adj) <= 14:
                b1 = brute_ursell_bk(adj, 1)
                b2 = brute_ursell_bk(adj, 2)
                if b1 != Fraction(a1):
                    fails['brute_b1'] += 1
                if b2 != Fraction(-E):
                    fails['brute_b2'] += 1
                if len(adj) <= 11:
                    b3 = brute_ursell_bk(adj, 3)
                    if b3 != Fraction(P2) - Fraction(t):
                        fails['brute_b3'] += 1
            # exp(sum c_k z^k) == I formal
            deg = len(alpha) - 1
            if deg >= 1:
                cum = formal_cumulants(alpha, deg)
                Irec = poly_from_cumulants(cum, deg)
                if any(Irec[k] != Fraction(alpha[k]) for k in range(deg + 1)):
                    fails['exp_cum'] += 1
            H_free = 3 ** a1
            if not (H_free >= H_ocf):
                fails['free'] += 1
            if (H_free == H_ocf) != (E == 0):
                fails['free_iff'] += 1
        achievable_H[n] = Hset
        print(f"n={n}: {2**(n*(n-1)//2)} tournaments, achievable H = {sorted(Hset)}")

    print()
    print("FAILURE COUNTS (want all 0 for verified claims):")
    for k, v in fails.items():
        print(f"  {k:12s}: {v}")

    print()
    print("FORBIDDEN-VALUE CHECK:")
    # I(K3,2): K3 conflict graph alpha = [1,3,0]  (3 vertices pairwise adjacent => alpha_2=0)
    IK3_2 = 1 + 3 * 2 + 0 * 4
    print(f"  I(K3, z) = 1 + 3z + 0 z^2 ; I(K3,2) = {IK3_2}  (want 7)")
    allH = set()
    for n in (3, 4, 5, 6):
        allH |= achievable_H[n]
    print(f"  7 in achievable H (n<=6)?  {7 in allH}")
    print(f"  21 in achievable H (n<=6)? {21 in allH}")
    print(f"  union achievable H (n<=6): {sorted(allH)}")

    # Paley T_7 cross-check
    print()
    print("PALEY T_7 CHECK (QR={1,2,4} mod 7):")
    n = 7
    QR = {1, 2, 4}
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and ((j - i) % 7) in QR:
                A[i][j] = 1
    cyc = odd_cycles(A, n)
    adj = conflict_graph(cyc)
    alpha = indep_poly(adj)
    # count by length
    c3 = sum(1 for s, vs in cyc if len(vs) == 3)
    c5 = sum(1 for s, vs in cyc if len(vs) == 5)
    c7 = sum(1 for s, vs in cyc if len(vs) == 7)
    H = sum(alpha[k] * 2 ** k for k in range(len(alpha)))
    Hd = H_hampaths(A, n)
    print(f"  odd cycles: c3={c3} c5={c5} c7={c7} total alpha_1={alpha[1]}")
    print(f"  alpha (first few) = {alpha[:4]}")
    print(f"  H(OCF) = {H}, H(direct paths) = {Hd}  (canon says 189)")

    # radius of convergence numeric (literal divergence at z=2)
    import numpy as np
    deg = len(alpha) - 1
    coeffs = [float(alpha[deg - k]) for k in range(deg + 1)]
    roots = np.roots(coeffs)
    rmin = min(abs(r) for r in roots)
    cum = formal_cumulants(alpha, 80)
    print(f"  smallest |root| of I = R = {rmin:.6f} ; z=2 {'OUTSIDE' if rmin<2 else 'inside'}")
    for K in (10, 40, 80):
        s = sum(float(cum[k - 1]) * 2.0 ** k for k in range(1, K + 1))
        print(f"    literal sum_(k<=){K} c_k 2^k = {s:.3e}  (log H = {math.log(H):.4f})")


if __name__ == "__main__":
    main()
