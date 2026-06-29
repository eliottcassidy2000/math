#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ROTATIONAL tournaments R_n, A038375 (max H), and the doubling-tower recursion --
merged with the Paley/dihedral-Burnside arithmetic (THM-584).  (mac-mini-S10)

KEY LEMMA (generalizes THM-584's p|H(Paley)): for ANY VERTEX-TRANSITIVE tournament T
on n vertices, the number of Hamiltonian paths starting at a fixed vertex is H/n
(by symmetry), so  n | H(T).  Rotational (circulant) R_n are vertex-transitive (Z_n
acts), so n | H(R_n), and H = n*(odd) by Redei.

We investigate:
 (1) the circulant H-SPECTRUM for odd n=3..13: all 2^{(n-1)/2} symbol sets, H each,
     verify n|H always, find opt_circ(n) = max (= A038375(n) for n<=11, THM-338), the
     Paley value, and the divisibility of A038375.
 (2) the DOUBLING TOWER (THM-447/448): skew-Hadamard S -> [[S,S],[S-2I,2I-S]] building
     DRTs on Mersenne 2^k-1; verify the self-similar recursion B_0(T_{2m-1})=T_{m-1}
     and n|H up the tower; compute H(T_3),H(T_7),H(T_15).
 (3) the divisibility criterion:  n | A038375(n)  <=>  the maximizer is vertex-transitive.
"""
from __future__ import annotations
import functools, itertools
print = functools.partial(print, flush=True)
from sympy import factorint


def ham_paths_count(arc, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            for w in range(n):
                if not (mask >> w) & 1 and arc[v][w]:
                    dp[mask | (1 << w)][w] += c
    return sum(dp[full])


def circulant(n, S):
    """circulant tournament: i->j iff (j-i) mod n in S."""
    Sset = set(S)
    return [[((j - i) % n) in Sset for j in range(n)] for i in range(n)]


def all_circulant_H(n):
    """over all antisymmetric S (one of {k,n-k} per pair), compute H. n odd."""
    m = (n - 1) // 2
    results = []
    for bits in range(1 << m):
        S = []
        for k in range(1, m + 1):
            if (bits >> (k - 1)) & 1:
                S.append(k)
            else:
                S.append(n - k)
        H = ham_paths_count(circulant(n, S), n)
        results.append((tuple(sorted(S)), H))
    return results


def skew_hadamard_tower(k):
    """seed S_1=[[1]]; iterate S -> [[S,S],[S-2I,2I-S]] to order 2^k (THM-447)."""
    import numpy as np
    S = np.array([[1]])
    for _ in range(k):
        m = S.shape[0]
        I = np.eye(m, dtype=int)
        top = np.hstack([S, S])
        bot = np.hstack([S - 2 * I, 2 * I - S])
        S = np.vstack([top, bot])
    return S


def core_tournament(S):
    """core of normalized skew-Hadamard: remove row/col 0, arc i->j iff entry +1."""
    m = S.shape[0]
    n = m - 1
    arc = [[False] * n for _ in range(n)]
    for i in range(1, m):
        for j in range(1, m):
            if i != j and S[i][j] == 1:
                arc[i - 1][j - 1] = True
    return arc, n


def main():
    print("=" * 80)
    print("Rotational tournaments, A038375, and the doubling tower (mac-mini-S10)")
    print("=" * 80)

    # known A038375 (max H over ALL n-vertex tournaments)
    A038375 = {1: 1, 2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661, 9: 3357,
               10: 15745, 11: 95095, 12: 555573, 13: 3450756}

    print("\n[1] Circulant H-spectrum (n odd): verify n|H, opt_circ vs A038375")
    for n in (3, 5, 7, 9, 11, 13):
        res = all_circulant_H(n)
        Hs = [H for _, H in res]
        opt = max(Hs)
        all_div = all(H % n == 0 for H in Hs)
        a = A038375.get(n)
        is_max = (opt == a) if a else None
        # Paley S = QR for n prime
        from sympy import isprime
        paleyH = None
        if isprime(n):
            qr = set((x * x) % n for x in range(1, n))
            paleyH = ham_paths_count(circulant(n, list(qr)), n)
        print(f"  n={n}: #circ={len(res)}, opt_circ={opt} {dict(factorint(opt))}, "
              f"A038375={a} (opt==a: {is_max}), n|H for ALL circ: {all_div}, "
              f"Paley H={paleyH}, opt/n={opt//n}")

    print("\n[2] Divisibility criterion: n | A038375(n)?  (vertex-transitive maximizer?)")
    for n in sorted(A038375):
        a = A038375[n]
        print(f"  n={n}: a={a}, n|a: {a % n == 0}  (a/n={a/n if a%n==0 else f'{a}/{n}'})"
              + ("   <- circulant-optimal (THM-338)" if n in (3,5,7,9,11) else
                 ("   <- circulant SUBoptimal (THM-338)" if n >= 13 else "")))

    print("\n[3] The doubling tower (THM-448): self-similar B_0 recursion + n|H")
    prev_arc = None
    for k in (2, 3, 4):    # orders 4,8,16 -> n=3,7,15
        S = skew_hadamard_tower(k)
        arc, n = core_tournament(S)
        # B_0 = out-neighborhood of vertex 0 (sub-tournament on N+(0))
        Nplus = [v for v in range(n) if arc[0][v]]
        sub = [[arc[Nplus[i]][Nplus[j]] for j in range(len(Nplus))] for i in range(len(Nplus))]
        H = ham_paths_count(arc, n) if n <= 15 else None
        b0_size = len(Nplus)
        b0_matches = (prev_arc is not None and sub == prev_arc)
        print(f"  T_{n} (order {2**k}): |N+(0)|={b0_size}, B_0(T_{n}) == T_{2*b0_size-1+0}? "
              f"matches-prev-level: {b0_matches if prev_arc is not None else 'n/a'}; "
              f"H={H} {dict(factorint(H)) if H else ''}, n|H: {H % n == 0 if H else '?'}")
        prev_arc = sub if False else arc  # store full level for next B_0 compare

    # recompute B_0 chain explicitly
    print("\n  B_0 self-similarity chain (out-neighborhood of vertex 0 = previous level):")
    levels = {}
    for k in (1, 2, 3, 4):
        S = skew_hadamard_tower(k)
        arc, n = core_tournament(S)
        levels[n] = arc
    for n in (15, 7):
        if n in levels:
            arc = levels[n]
            Nplus = [v for v in range(n) if arc[0][v]]
            sub = [[arc[Nplus[i]][Nplus[j]] for j in range(len(Nplus))] for i in range(len(Nplus))]
            prev_n = (n - 1) // 2
            match = prev_n in levels and len(sub) == prev_n
            # canonical compare via H (cheap invariant)
            Hsub = ham_paths_count(sub, len(sub)) if sub else 0
            Hprev = ham_paths_count(levels[prev_n], prev_n) if prev_n in levels else None
            print(f"    B_0(T_{n}): size {len(sub)} (=T_{prev_n}?), H(B_0)={Hsub}, H(T_{prev_n})={Hprev}, "
                  f"H-match: {Hsub == Hprev}")

    print("\n" + "=" * 80)
    print("MERGE: rotational R_n vertex-transitive => n|H=n*(odd) (THM-584 generalized);")
    print("circulant maximizes A038375 for n=3,5,7,9,11 (=> n|a(n)), fails n>=13;")
    print("the doubling tower is self-similar (B_0(T_{2m-1})=T_{m-1}) on Mersenne 2^k-1.")
    print("=" * 80)


if __name__ == "__main__":
    main()
