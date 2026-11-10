#!/usr/bin/env python3
"""AUDIT: Theorem C (2) as the restriction to a vertex star of the repository's digit-1 law
(THM-478 / THM-466; Grinberg-Stanley: H = 1 + 2 c_odd mod 4, c_odd = number of directed odd cycles).
An odd cycle through the new vertex v is v -> i -> ... -> j -> v with i in X, j not in X and an i-to-j path of R on
an EVEN number of vertices.  Hence, with P_ij = number of directed paths of R from i to j with an even number of
vertices,
      q_R(X) = (H(R+v_X) - H(R))/2 = sum_{i in X, j notin X} P_ij   (mod 2),
      beta_R(i,j) = mu_R(i,j)/2 = P_ij + P_ji                         (mod 2).
Checked here on every labelled tournament of orders 2..5 and on seeded random samples of orders 6, 7, 8, for every
out-set X (auditor's own code)."""
from itertools import permutations
def tournaments(m):
    prs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    for code in range(1 << len(prs)):
        adj = [[0] * m for _ in range(m)]
        for k, (i, j) in enumerate(prs):
            if code >> k & 1: adj[i][j] = 1
            else: adj[j][i] = 1
        yield adj
def ham_dp(n, adj):
    out = [sum(1 << w for w in range(n) if adj[v][w]) for v in range(n)]
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = dp[S][v]
            if c:
                nx = out[v] & ~S
                while nx:
                    w = (nx & -nx).bit_length() - 1; nx &= nx - 1; dp[S | (1 << w)][w] += c
    return sum(dp[(1 << n) - 1])
def even_paths(m, adj):
    """P[i][j] = number of simple directed paths i -> ... -> j on an even number of vertices (DFS)"""
    P = [[0] * m for _ in range(m)]
    def dfs(start, v, used, cnt):
        if cnt % 2 == 0 and v != start: P[start][v] += 1
        for w in range(m):
            if adj[v][w] and not used >> w & 1: dfs(start, w, used | 1 << w, cnt + 1)
    for s in range(m): dfs(s, s, 1 << s, 1)
    return P
import random
def sample(m, reps, rng):
    for _ in range(reps):
        adj = [[0] * m for _ in range(m)]
        for i in range(m):
            for j in range(i + 1, m):
                if rng.random() < 0.5: adj[i][j] = 1
                else: adj[j][i] = 1
        yield adj
ok = True; rng = random.Random(478)
for m in range(2, 9):
    bad = 0; cnt = 0
    for adj in (tournaments(m) if m <= 5 else sample(m, {6: 300, 7: 60, 8: 12}[m], rng)):
        cnt += 1; H0 = ham_dp(m, adj); P = even_paths(m, adj); vals = []
        for X in range(1 << m):
            a2 = [row[:] + [0] for row in adj] + [[0] * (m + 1)]
            for u in range(m):
                if X >> u & 1: a2[m][u] = 1
                else: a2[u][m] = 1
            vals.append(ham_dp(m + 1, a2))
        for X in range(1 << m):
            cut = sum(P[i][j] for i in range(m) for j in range(m) if X >> i & 1 and not X >> j & 1)
            if ((vals[X] - H0) // 2 - cut) % 2: bad += 1
        for i in range(m):
            for j in range(i + 1, m):
                mu = (vals[1 << i] - H0) + (vals[1 << j] - H0) - (vals[(1 << i) | (1 << j)] - H0)
                if (mu // 2 - P[i][j] - P[j][i]) % 2: bad += 1
    ok &= (bad == 0)
    print(f"m={m}: {cnt} labelled tournaments ({'all' if m <= 5 else 'seeded random sample'}), all out-sets and pairs: failures {bad}", flush=True)
print("ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
