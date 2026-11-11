#!/usr/bin/env python3
"""S93 coordinator: the counts x(h) = 2h + 1 - 2q(h) are Hamiltonian path counts of one-apex tournaments (brute force).

One-apex tournament of a word (b_1, ..., b_n): transitive spine s_1 -> s_2 -> ... -> s_n (s_i -> s_j for i < j) and
an apex u with u -> s_i iff b_i = 1.  H is computed by the subset dynamic programme.  For every n <= NMAX the multiset
{H(word) : words of length n} is compared with {g(v) : 0 <= v < 2^n} (g(0) = 1, g(2v+1) = g(v),
g(2v) = g(v) + v + (v mod 2)), and the set of values with {x(h) : h < 2^(n-1)} (x(h) = g(2h)).
"""
import sys
from collections import Counter


def q_of(h):
    q = 0
    b = 0
    while (h >> b) > 0:
        if (h >> b) & 1:
            q += ((h >> (b + 1)) + 1) >> 1
        b += 1
    return q


def g_of(v):
    return v + 1 - 2 * q_of(v) - (v & 1)


def ham(n, word):
    V = n + 1
    u = n
    adj = [[0] * V for _ in range(V)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = 1
        if (word >> i) & 1:
            adj[u][i] = 1
        else:
            adj[i][u] = 1
    dp = [[0] * V for _ in range(1 << V)]
    for x in range(V):
        dp[1 << x][x] = 1
    for S in range(1, 1 << V):
        row = dp[S]
        for x in range(V):
            c = row[x]
            if c:
                for y in range(V):
                    if adj[x][y] and not (S >> y) & 1:
                        dp[S | (1 << y)][y] += c
    return sum(dp[(1 << V) - 1])


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    ok = True
    # recursion of g against the closed form
    g = {0: 1}
    for v in range(1, 1 << 12):
        g[v] = g[v >> 1] if v & 1 else g[v >> 1] + (v >> 1) + ((v >> 1) & 1)
    ok &= all(g[v] == g_of(v) for v in range(1 << 12))
    print("g(v) = v + 1 - 2 q(v) - (v mod 2) for v < 4096:", ok)
    for n in range(1, nmax + 1):
        hs = Counter(ham(n, wd) for wd in range(1 << n))
        gs = Counter(g_of(v) for v in range(1 << n))
        xs = set(2 * h + 1 - 2 * q_of(h) for h in range(1 << (n - 1)))
        same_multiset = hs == gs
        same_set = set(hs) == xs
        ok &= same_multiset and same_set
        print(f"n = {n}: multiset of H over the 2^n words = multiset of g(v), v < 2^n: {same_multiset}; "
              f"set of values = {{x(h) : h < 2^(n-1)}}: {same_set} ({len(xs)} values, max {max(xs)})", flush=True)
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
