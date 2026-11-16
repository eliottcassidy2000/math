"""Collatz on [1..N] with an ordering: arcs follow the order, except Collatz edges whose parent comes first are
reversed (child -> parent).  numeric order: only n -> n/2 survive as back arcs (doubling rows);
reverse numeric: only c -> 3c+1 survive (tripling hooks).  Also: all linear extensions of a Collatz BFS prefix."""
import subprocess, itertools, sys
sys.setrecursionlimit(10000)
def Hgen(n, back):          # vertices 0..n-1 in order; back = set of (j, i) with j > i meaning arc j -> i
    # encode as a general tournament for a small DP here (n <= 20)
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if (j, i) in back: out[j] |= 1 << i
            else: out[i] |= 1 << j
    N = 1 << n; dp = [[0] * n for _ in range(N)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, N):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if c:
                t = out[v] & ~S
                while t:
                    u = (t & -t).bit_length() - 1; t &= t - 1
                    dp[S | 1 << u][u] += c
    return sum(dp[N - 1])
def coll_parent(m):
    return m // 2 if m % 2 == 0 else 3 * m + 1
res_num, res_rev = [], []
for N in range(1, 17):
    nodes = list(range(1, N + 1))
    pos = {v: i for i, v in enumerate(nodes)}
    back = {(pos[m], pos[coll_parent(m)]) for m in nodes if m > 1 and coll_parent(m) in pos and pos[coll_parent(m)] < pos[m]}
    res_num.append(Hgen(N, back))
    nodes_r = list(range(N, 0, -1)); pos = {v: i for i, v in enumerate(nodes_r)}
    back = {(pos[m], pos[coll_parent(m)]) for m in nodes_r if m > 1 and coll_parent(m) in pos and pos[coll_parent(m)] < pos[m]}
    res_rev.append(Hgen(N, back))
print("numeric order on [1..N] (back arcs 2k -> k, the doubling rows):", res_num)
print("reverse numeric on [1..N] (back arcs c -> 3c+1, the hooks):", res_rev)
