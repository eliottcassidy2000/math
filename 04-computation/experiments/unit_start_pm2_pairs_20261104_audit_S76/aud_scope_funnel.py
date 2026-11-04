"""aud_scope_funnel.py - is the A/B +-2 step transportable beyond Session189's D_m wrapper?

D_m(U, r) sees U only through (h, F).  The zero-ready / FUNNELS F3 calculus sees (h, e, P(2), F, S_w(2)).
For A and B these are (29, 10, 580, 157, 34) and (31, 9, 644, 157, 31).  This script evaluates a few other literal
one-vertex and two-vertex attachments at the root and prints H for A and B, to show the difference is 2 only for the
TT_m-plus-root wrapper.
"""
import itertools


def ham_count(n, adj):
    # subset DP
    dp = [[0] * n for _ in range(1 << n)]
    for x in range(n):
        dp[1 << x][x] = 1
    for S in range(1, 1 << n):
        for x in range(n):
            c = dp[S][x]
            if not c:
                continue
            for y in range(n):
                if not (S >> y) & 1 and adj[x][y]:
                    dp[S | (1 << y)][y] += c
    return sum(dp[(1 << n) - 1])


def from_rows(R):
    n = len(R)
    return [[1 if (R[a] >> b) & 1 else 0 for b in range(n)] for a in range(n)]


def extend(adj, beats_new, new_beats):
    """add vertex x = n: new_beats = set of old vertices x beats; others beat x"""
    n = len(adj)
    out = [row[:] + [0] for row in adj] + [[0] * (n + 1)]
    for v in range(n):
        if v in new_beats:
            out[n][v] = 1
        else:
            out[v][n] = 1
    return out


A = from_rows([8, 73, 3, 4, 15, 31, 61])
B = from_rows([2, 84, 1, 7, 13, 31, 61])
rest = set(range(1, 7))
# attachments printed below: Session189 wrappers D_1, D_2; a sole exit x -> 0 only (THM-4862 (b));
# x with 0 -> x -> U - 0 (the m = 0 wrapper, H = h + F - s); x beating exactly the root's out-neighbours.


def wrapper(adj, m):
    n = len(adj)
    N = n + m + 1
    out = [[0] * N for _ in range(N)]
    for a in range(n):
        for b in range(n):
            out[a][b] = adj[a][b]
    w = n + m
    for i in range(m):
        z = n + i
        for v in range(n):
            out[z][v] = 1
        for j in range(i + 1, m):
            out[z][n + j] = 1
        out[w][z] = 1
    for v in range(1, n):
        out[w][v] = 1
    out[0][w] = 1
    return out


for name, U in (("A", A), ("B", B)):
    print(name, "H =", ham_count(7, U),
          "| D_1:", ham_count(9, wrapper(U, 1)),
          "| D_2:", ham_count(10, wrapper(U, 2)),
          "| + sole exit x -> 0:", ham_count(8, extend(U, None, {0})),
          "| + x with 0 -> x -> U-0:", ham_count(8, extend(U, None, rest)),
          "| + x beating N+(0):", ham_count(8, extend(U, None, {v for v in range(7) if U[0][v]})))
