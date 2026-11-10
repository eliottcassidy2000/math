#!/usr/bin/env python3
"""AUDIT check (d): independent verification of printed witnesses (output of d_witness.c).
Each input line:  h n m <gentourng line> <vertices beaten by z ...>
T is rebuilt from the gentourng line (upper triangle row by row, '1' = i -> j for i < j), z is added as vertex m with
out-set X, and both Hamiltonian path counts are recomputed with a plain subset DP that knows nothing about cut
quadratics.  For m + 1 <= 9 the count is additionally confirmed by enumerating all permutations."""
import sys
from itertools import permutations
def build(m, line):
    adj = [[0] * (m + 1) for _ in range(m + 1)]; k = 0
    for i in range(m):
        for j in range(i + 1, m):
            if line[k] == "1": adj[i][j] = 1
            else: adj[j][i] = 1
            k += 1
    assert k == len(line)
    return adj
def ham(n, adj):
    out = [sum(1 << w for w in range(n) if adj[v][w]) for v in range(n)]
    dp = [dict() for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v, c in dp[S].items():
            nx = out[v] & ~S
            while nx:
                w = (nx & -nx).bit_length() - 1; nx &= nx - 1
                T = S | (1 << w); dp[T][w] = dp[T].get(w, 0) + c
    return sum(dp[(1 << n) - 1].values())
ok = True
for l in sys.stdin:
    p = l.split()
    if len(p) < 4: continue
    h, n, m, line = int(p[0]), int(p[1]), int(p[2]), p[3]; X = [int(x) for x in p[4:]]
    adj = build(m, line)
    hT = ham(m, [r[:m] for r in adj[:m]])
    for i in range(m):
        if i in X: adj[m][i] = 1
        else: adj[i][m] = 1
    hTz = ham(m + 1, adj)
    extra = ""
    if m + 1 <= 9:
        brute = sum(1 for P in permutations(range(m + 1)) if all(adj[P[t]][P[t + 1]] for t in range(m)))
        extra = f"; permutation count {brute}"; ok &= (brute == hTz)
    good = (hT == h and hTz == n); ok &= good
    print(f"h={h} n={n} (n/h = {n / h:.4f}): base order {m}, T = {line}, z beats {X}:  H(T) = {hT}, H(T+z) = {hTz}{extra}  -> {'CONFIRMED' if good else 'MISMATCH'}")
print("ALL WITNESSES CONFIRMED" if ok else "SOME WITNESS FAILED")
