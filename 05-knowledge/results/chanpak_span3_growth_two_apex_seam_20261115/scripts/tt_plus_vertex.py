#!/usr/bin/env python3
"""tt_plus_vertex.py: Hamiltonian path counts of a transitive tournament TT_m plus one vertex u whose out-neighbourhood
is a prescribed set S of ranks (1 = source ... m = sink); u loses to the other m - |S| vertices.  For m = 5..11 and every
S with |S| <= 3 (and, for m <= 9, every S) the maximum H and its argmax are printed, with the counts of the ribbon
maximisers S = {r, r+2}.  (Context for THM-5243 (2): the literal span-3 maximiser at length n is TT_{n+1} + u with
S = {n-4, n-2}.)"""
import itertools
def hpaths(adj):
    n = len(adj); full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            for w in range(n):
                if not mask >> w & 1 and adj[v][w]: dp[mask | 1 << w][w] += c
    return sum(dp[full])
def tt_plus(m, S):
    n = m + 1; adj = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(i + 1, m): adj[i][j] = 1          # rank i+1 beats rank j+1
    for r in range(1, m + 1):
        if r in S: adj[m][r - 1] = 1
        else: adj[r - 1][m] = 1
    return adj
for m in range(5, 12):
    sizes = range(0, m + 1) if m <= 9 else range(0, 4)
    best = (0, None); best3 = (0, None); table = {}
    for k in sizes:
        for S in itertools.combinations(range(1, m + 1), k):
            Hc = hpaths(tt_plus(m, set(S))); table[S] = Hc
            if Hc > best[0]: best = (Hc, S)
            if k <= 3 and Hc > best3[0]: best3 = (Hc, S)
    pairs = {r: table[(r, r + 2)] for r in range(1, m - 1)}
    print(f"m = {m}: max over |S| <= 3: H = {best3[0]} at S = {best3[1]}; " + (f"max over all S: H = {best[0]} at S = {best[1]}; " if m <= 9 else "") + f"H for S = {{r, r+2}}: {pairs}")

# The ribbon-admissible sets: in a literal span-3 word of length n (vertices v_0..v_n) each v_j is beaten by at least
# j - 3 older vertices, so if T_W is transitive the last three added vertices have ranks >= n-4, n-3, n-2 (m = n+1 ranks),
# and the new vertex beats a subset S of them.  Max of H over all S within the last six ranks, |S| <= 3:
print("maximum of H(TT_m + u) over S within the last six ranks {m-5..m}, |S| <= 3 (ribbon-admissible S are among these):")
for m in range(6, 12):
    cands = []
    for k in range(0, 4):
        for S in itertools.combinations(range(max(1, m - 5), m + 1), k):
            cands.append((hpaths(tt_plus(m, set(S))), S))
    mx = max(c[0] for c in cands); arg = [S for Hc, S in cands if Hc == mx]
    print(f"  m = {m}: max H = {mx} at S = {arg} (as ranks counted from the sink: {[tuple(m + 1 - r for r in S) for S in arg]})")
# closed form check: H(TT_m + u) = 5 * 2^(m-r-2) - 1 for S = {r, r+2}, r + 2 < m
ok = all(hpaths(tt_plus(m, {r, r + 2})) == 5 * 2 ** (m - r - 2) - 1 for m in range(4, 12) for r in range(1, m - 2))
print("closed form H = 5*2^(m-r-2) - 1 for S = {r, r+2}, r + 2 < m, checked for m <= 11:", ok)
