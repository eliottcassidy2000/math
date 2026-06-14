#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The first non-spectral OCF ingredient: H = 1 + 2(c3+c5) + 4D at n=6, where D =
#{vertex-disjoint triangle pairs}, and D (not c3,c5) is the non-spectral part.
kind-pasteur-2026-06-13-S5.  Completes THM-499 with the exact mechanism.

THM-499: H is not spectrally determined from n=6 (cospectral tournaments, distinct
H).  At n=6 the only odd cycles are 3- and 5-cycles; both c3=tr(A^3)/3 and
c5=tr(A^5)/5 are SPECTRAL (THM-118), so alpha_1 = c3+c5 is fixed within a cospectral
class.  The OCF H = I(Omega,2) = sum_k alpha_k 2^k.  At n=6, alpha_2 = # pairs of
VERTEX-DISJOINT odd cycles = # pairs of disjoint triangles (two triangles fill all 6
vertices; a 3-cycle and a 5-cycle can't be vertex-disjoint in 6 vertices; alpha_3
needs >=9 vertices).  So the prediction:

    H(T) = 1 + 2*(c3+c5) + 4*D,   D = #{unordered pairs of vertex-disjoint triangles},

with c3,c5 SPECTRAL and D the FIRST non-spectral ingredient (the cospectral/distinct-H
tournaments of THM-499 differ exactly in D).

VERIFY (exhaustive n=6, exact): the formula holds for every tournament; D varies
within cospectral classes; H = 25+4*alpha_2 reproduces the THM-499 witness {25,29,33}.
Also check n=5 (D=0 always, since disjoint triangles need 6 vertices) => H=1+2*c3 (c5=0
contributes? c5 at n=5 can be >0; alpha_1 = c3 + c5).  Actually at n=5 alpha_2=0 (no
room for two disjoint cycles), so H = 1 + 2*(c3+c5) exactly => H spectral at n=5
(consistent with THM-499 part 1).
"""

import sys, itertools
from collections import defaultdict
import numpy as np
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def all_tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        A = np.zeros((n, n), dtype=np.int64)
        adj = [[False] * n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i, j] = 1; adj[i][j] = True
            else:
                A[j, i] = 1; adj[j][i] = True
        yield A, adj


def H_count(A, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            cur = dp[mask][v]
            if not cur:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v, w]:
                    dp[mask | (1 << w)][w] += cur
    return sum(dp[full][v] for v in range(n))


def triangles(adj, n):
    """list of frozensets of vertices forming a directed 3-cycle."""
    tris = []
    for a, b, c in itertools.combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            tris.append((a, b, c))
    return tris


def disjoint_triangle_pairs(adj, n):
    tris = triangles(adj, n)
    D = 0
    for i in range(len(tris)):
        for j in range(i + 1, len(tris)):
            if not (set(tris[i]) & set(tris[j])):
                D += 1
    return D


def cycle_counts(A, n):
    c3 = int(np.trace(A @ A @ A)) // 3
    A2 = A @ A; A4 = A2 @ A2
    c5 = int(np.trace(A4 @ A)) // 5
    return c3, c5


def spectral_sig(A, n):
    sig = []; P = np.eye(n, dtype=np.int64)
    for k in range(1, n + 1):
        P = P @ A; sig.append(int(np.trace(P)))
    return tuple(sig)


def main():
    for n in (5, 6):
        bad = 0; total = 0
        by_sig = defaultdict(lambda: defaultdict(set))  # sig -> {} ; track (c3,c5)->H and D
        sig_D = defaultdict(set)
        for A, adj in all_tournaments(n):
            total += 1
            H = H_count(A, n)
            c3, c5 = cycle_counts(A, n)
            D = disjoint_triangle_pairs(adj, n)
            pred = 1 + 2 * (c3 + c5) + 4 * D
            if pred != H:
                bad += 1
                if bad <= 3:
                    print(f"   n={n} MISMATCH: H={H} but 1+2({c3}+{c5})+4*{D}={pred}", flush=True)
            sig = spectral_sig(A, n)
            sig_D[sig].add(D)
        # does D vary within cospectral classes? (the non-spectral part)
        split_D = sum(1 for s, ds in sig_D.items() if len(ds) > 1)
        print(f"=== n={n} ({total} tournaments) ===", flush=True)
        print(f"   H = 1 + 2(c3+c5) + 4D  holds: {'ALL' if bad==0 else str(bad)+' FAIL'} "
              f"({total-bad}/{total})", flush=True)
        print(f"   cospectral classes where D varies (D is the NON-spectral part): {split_D}", flush=True)
        if n == 6:
            print(f"   => at n=6: c3,c5 spectral (=tr/k), D NOT spectral; H resolves cospectral", flush=True)
            print(f"      tournaments exactly through D (alpha_2 = disjoint-triangle pairs).", flush=True)
        if n == 5:
            print(f"   => at n=5: D=0 always (disjoint triangles need 6 vtx), so H=1+2(c3+c5) is", flush=True)
            print(f"      fully spectral -- matching THM-499 part 1 (H spectral at n<=5).", flush=True)
        print(flush=True)


if __name__ == "__main__":
    main()
