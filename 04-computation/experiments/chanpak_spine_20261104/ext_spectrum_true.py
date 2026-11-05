#!/usr/bin/env python3
"""ext_spectrum_true.py -- the TRUE one-vertex extension spectrum E(X) = {H(X + v) : v a new vertex with any cut},
via the path-pair decomposition: a Hamiltonian path of X + v is (P1, v, P2) with P1 a path of X ending at u -> v
(or empty) and P2 a path of X starting at w with v -> w (or empty), V(P1) u V(P2) = V(X).
Checks the formula by brute force, then tests the reduction lemma: does E(X) contain every odd number in [h, 2h - 1]?"""
import sys, subprocess, itertools
import numpy as np

def tournaments(n):
    out = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True).stdout.split()
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for line in out:
        adj = [0] * n
        for (i, j), b in zip(pairs, line):
            if b == '1': adj[i] |= 1 << j
            else: adj[j] |= 1 << i
        yield adj

def end_start_tables(n, adj):
    """e[S][u] = # Hamiltonian paths of X[S] ending at u; s[S][u] = starting at u."""
    N = 1 << n
    e = np.zeros((N, n), dtype=np.int64); s = np.zeros((N, n), dtype=np.int64)
    radj = [sum(1 << u for u in range(n) if (adj[u] >> v) & 1) for v in range(n)]
    for v in range(n): e[1 << v][v] = 1; s[1 << v][v] = 1
    for S in range(1, N):
        for v in range(n):
            c = e[S][v]
            if c:
                nxt = adj[v] & ~S
                while nxt:
                    u = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
                    e[S | (1 << u)][u] += c
            c = s[S][v]
            if c:
                nxt = radj[v] & ~S
                while nxt:
                    u = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
                    s[S | (1 << u)][u] += c
    return e, s

def ext_spectrum(n, adj):
    """H(X + v) for every cut A (bitmask over V(X): v -> A, V \\ A -> v), as an array over A."""
    e, s = end_start_tables(n, adj)
    N = 1 << n; full = N - 1
    cuts = np.arange(N, dtype=np.int64)
    Amat = ((cuts[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int64)     # N x n: A membership
    EB = e @ (1 - Amat).T          # N(S) x N(cuts): sum_{u in S \ A} e[S][u]
    SA = s @ Amat.T                # N(S) x N(cuts): sum_{w in S cap A} s[S][w]
    comp = full - cuts             # complement index of S
    H = (EB[1:full] * SA[comp[1:full]]).sum(axis=0)       # S nonempty, proper
    H = H + SA[full] + EB[full]                             # P1 empty / P2 empty
    return H

def ham_paths(n, adj):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c: continue
            nxt = adj[v] & ~S
            while nxt:
                u = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
                dp[S | (1 << u)][u] += c
    return sum(dp[full])

def strong(n, adj):
    def reach_all(a):
        r = 1; fr = [0]
        while fr:
            u = fr.pop(); nx = a[u] & ~r
            while nx:
                v = (nx & -nx).bit_length() - 1; nx &= nx - 1; r |= 1 << v; fr.append(v)
        return r == (1 << n) - 1
    radj = [sum(1 << u for u in range(n) if (adj[u] >> v) & 1) for v in range(n)]
    return reach_all(adj) and reach_all(radj)

if __name__ == "__main__":
    # brute-force check of the formula on all tournaments of order <= 5, all cuts
    bad = 0; tot = 0
    for n in range(2, 6):
        for adj in tournaments(n):
            H = ext_spectrum(n, adj)
            for A in range(1 << n):
                adj2 = [adj[u] | (0 if (A >> u) & 1 else (1 << n)) for u in range(n)] + [A]
                tot += 1
                if ham_paths(n + 1, adj2) != H[A]: bad += 1
    print("formula check: %d (X, cut) pairs, orders <= 5: mismatches %d" % (tot, bad), flush=True)
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for n in range(3, NMAX + 1):
        rows = []
        for adj in tournaments(n):
            H = ext_spectrum(n, adj); h = int(H[0]); E = set(int(x) for x in H)
            assert min(E) == h
            # largest L with every odd in [h, L] in E
            L = h
            while L + 2 in E: L += 2
            rows.append((h, L, max(E), strong(n, adj)))
        # per value h: best realiser (max L/h) and whether some realiser has L >= 2h-1
        by_h = {}
        for h, L, M, st in rows:
            by_h.setdefault(h, []).append((L, M, st))
        good = sum(1 for h, lst in by_h.items() if max(L for L, M, st in lst) >= 2 * h - 1)
        bad_h = [h for h, lst in sorted(by_h.items()) if max(L for L, M, st in lst) < 2 * h - 1]
        allstrong_good = sum(1 for h, L, M, st in rows if st and L >= 2 * h - 1); nstrong = sum(1 for r in rows if r[3])
        print("n=%d: %d tournaments, %d distinct h; h with SOME realiser X having E(X) ⊇ odd[h, 2h-1]: %d of %d; failing h: %s%s; strong X satisfying it individually: %d of %d; min over X of L/h = %.3f, min of max/h = %.3f" % (
            n, len(rows), len(by_h), good, len(by_h), bad_h[:12], "..." if len(bad_h) > 12 else "", allstrong_good, nstrong,
            min(L / h for h, L, M, st in rows), min(M / h for h, L, M, st in rows)), flush=True)
