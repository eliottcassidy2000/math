#!/usr/bin/env python3
"""Auditor proof-5134-FH.  The session's construction of J(h1, h2, s) (docstring of every_base.py):
    J(h1, h2, s) = R(h1) => TT_s => R(h2)^converse, with the arc between the two apexes reversed,
    R(h) = TT_k + apex, k = bit length of 2h, apex -> c_i iff the i-th bit of 2h (most significant first) is 1.
Rebuilt here from that description (own code) and evaluated with the auditor's subset DP and, for orders <= 9, the
explicit enumeration of h1_packing_bruteforce.py.  Also: J is isomorphic to T(W) (.) B(W') of THM-5070 (2) with
W = 2^s1 (2 h1 + 1) - 1, W' = 2^s2 (2 h2 + 1) - 1 for every split s1 + s2 = s (same packing polynomial, same score
sequence, same number of Hamiltonian paths; isomorphism checked by canonical form for orders <= 8).
"""
import itertools
import sys
import os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h1_packing_bruteforce import (psi, build_T, build_B, glue, is_tournament, poly_dp, packing_poly_explicit,
                                   ham_paths, poly_add, poly_mul, evalp)


def R(h):
    u = 2 * h
    k = u.bit_length()
    adj = [[0] * (k + 1) for _ in range(k + 1)]
    for i in range(k):
        for j in range(i + 1, k):
            adj[i][j] = 1
        if (u >> (k - 1 - i)) & 1:
            adj[k][i] = 1
        else:
            adj[i][k] = 1
    return adj, k


def J(h1, h2, s):
    a1, ap1 = R(h1)
    a2, ap2 = R(h2)
    n1, n2 = len(a1), len(a2)
    n = n1 + s + n2
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = 1                      # ordered sum: earlier beats later
    for i in range(n1):
        for j in range(n1):
            if i != j:
                adj[i][j] = a1[i][j]
    off = n1 + s
    for i in range(n2):
        for j in range(n2):
            if i != j:
                adj[off + i][off + j] = a2[j][i]        # converse of R(h2)
    A, B = ap1, off + ap2
    adj[A][B] = 0
    adj[B][A] = 1
    return adj


def canon(adj):
    n = len(adj)
    best = None
    for perm in itertools.permutations(range(n)):
        key = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or key < best:
            best = key
    return best


def main():
    cnt = 0
    fail = 0
    iso_cnt = 0
    iso_fail = 0
    by_s = {}
    for h1 in range(0, 24):
        for h2 in range(0, 24):
            for s in range(0, 4):
                adj = J(h1, h2, s)
                n = len(adj)
                if n > 12:
                    continue
                assert is_tournament(adj)
                c = poly_dp(adj)
                seam = 2 ** s * (2 * h1 + 1) * (2 * h2 + 1) - (1 if s == 0 else 0)
                want = poly_add(poly_mul([1, psi(h1)], [1, psi(h2)]), [0, seam // 2])
                ok = (c == want)
                if n <= 9:
                    ce, _ = packing_poly_explicit(adj)
                    ok = ok and (ce == want)
                ok = ok and evalp(c, 2) == int(ham_paths(np.array(adj, dtype=np.int64), n))
                cnt += 1
                by_s[s] = by_s.get(s, 0) + 1
                if not ok:
                    fail += 1
                if n <= 8:
                    cj = canon(adj)
                    for s1 in range(0, s + 1):
                        W = 2 ** s1 * (2 * h1 + 1) - 1
                        Wp = 2 ** (s - s1) * (2 * h2 + 1) - 1
                        aT, u = build_T(W)
                        aB, w = build_B(Wp)
                        g = glue(aT, u, aB, w)
                        iso_cnt += 1
                        if len(g) != n or canon(g) != cj:
                            iso_fail += 1
    print(f"session construction J(h1, h2, s) = R(h1) => TT_s => R(h2)^conv, apex arc reversed; h1, h2 < 24, s <= 3, order <= 12:")
    print(f"  {cnt} joins (by s: {dict(sorted(by_s.items()))}); I(Omega(J), z) = (1 + psi1 z)(1 + psi2 z) + z (2^s y1 y2 - [s=0])/2 "
          f"and I(2) = #Hamiltonian paths: failures {fail}")
    print(f"  isomorphism J(h1,h2,s) = T(2^s1 (2h1+1) - 1) (.) B(2^s2 (2h2+1) - 1) for every split (orders <= 8, canonical forms): "
          f"{iso_cnt} comparisons, failures {iso_fail}")


if __name__ == "__main__":
    main()
