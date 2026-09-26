#!/usr/bin/env python3
"""tournament_insertion_slots_20260926.py -- the insertion-slot identity for Hamiltonian paths.

For a tournament T on n vertices, a vertex v, and a Hamiltonian path P = x_1 ... x_{n-1} of T - v,
the number of positions at which v can be inserted into P to give a Hamiltonian path of T is
   slots_v(P) = 1 + #{ i : v -> x_i and x_{i+1} -> v }
(prepend iff v -> x_1, append iff x_{n-1} -> v, insert between x_i, x_{i+1} iff x_i -> v -> x_{i+1};
with b_i = [v -> x_i] the count is #01 + b_1 + (1 - b_{n-1}) = 1 + #10 by telescoping). Hence
   H(T) = H(T - v) + sum_{P in HP(T - v)} #{ i : v -> x_i, x_{i+1} -> v }.
This script verifies the identity exhaustively for all labelled tournaments with n <= 6, on random
tournaments for n = 7..10, and tests whether the correction sum is always even (which would give
Redei's theorem, H(T) odd, by induction). Session collatz-exponent-atlas-20260926 (opus).
"""
import itertools, random, sys


def ham_paths(n, adj):
    """all Hamiltonian paths of the tournament on range(n); adj[i][j] = 1 iff i -> j."""
    paths = []
    def ext(path, used):
        if len(path) == n:
            paths.append(tuple(path)); return
        last = path[-1]
        for y in range(n):
            if not used & (1 << y) and adj[last][y]:
                path.append(y); ext(path, used | (1 << y)); path.pop()
    for s in range(n):
        ext([s], 1 << s)
    return paths


def count_ham_paths(n, adj):
    """DP count of Hamiltonian paths (number of orderings), for larger n."""
    N = 1 << n
    dp = [[0] * n for _ in range(N)]
    for s in range(n):
        dp[1 << s][s] = 1
    for mask in range(N):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if c == 0:
                continue
            for y in range(n):
                if not mask & (1 << y) and adj[last][y]:
                    dp[mask | (1 << y)][y] += c
    return sum(dp[N - 1])


def check(n, adj, v):
    """returns (H(T), H(T-v), correction sum, slot-identity ok)"""
    others = [u for u in range(n) if u != v]
    idx = {u: k for k, u in enumerate(others)}
    m = n - 1
    sub = [[adj[others[i]][others[j]] for j in range(m)] for i in range(m)]
    paths = ham_paths(m, sub)
    HT = count_ham_paths(n, adj)
    corr = 0
    ok = True
    for P in paths:
        xs = [others[k] for k in P]
        b = [adj[v][x] for x in xs]
        slots = (1 if b[0] else 0) + (1 if not b[-1] else 0) + sum(1 for i in range(m - 1) if not b[i] and b[i + 1])
        tens = sum(1 for i in range(m - 1) if b[i] and not b[i + 1])
        if slots != 1 + tens:
            ok = False
        corr += tens
    return HT, len(paths), corr, ok


def main():
    random.seed(20260926)
    for n in range(3, 7):
        pairs = list(itertools.combinations(range(n), 2))
        bad_slots = bad_identity = odd_corr = total = 0
        for bits in range(1 << len(pairs)):
            adj = [[0] * n for _ in range(n)]
            for k, (i, j) in enumerate(pairs):
                if bits >> k & 1:
                    adj[i][j] = 1
                else:
                    adj[j][i] = 1
            for v in range(n):
                HT, Hsub, corr, ok = check(n, adj, v)
                total += 1
                if not ok:
                    bad_slots += 1
                if HT != Hsub + corr:
                    bad_identity += 1
                if corr % 2 == 1:
                    odd_corr += 1
        print("n=%d: %d (tournament, v) pairs; slot formula failures %d; identity failures %d; odd correction sums %d" % (n, total, bad_slots, bad_identity, odd_corr))
    for n in range(7, 11):
        bad = odd = 0
        trials = 200 if n <= 8 else 40
        for _ in range(trials):
            adj = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    if random.random() < 0.5:
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1
            v = random.randrange(n)
            HT, Hsub, corr, ok = check(n, adj, v)
            if not ok or HT != Hsub + corr:
                bad += 1
            if corr % 2 == 1:
                odd += 1
        print("n=%d: %d random (tournament, v): failures %d; odd correction sums %d" % (n, trials, bad, odd))


if __name__ == '__main__':
    main()
