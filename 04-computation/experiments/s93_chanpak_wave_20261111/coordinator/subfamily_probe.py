#!/usr/bin/env python3
"""S93 coordinator probe (scratch): structured sub-families of the two-apex family G(2,0).

Spine s_0 -> ... -> s_(n-1) transitive, apexes u -> v, a_i = [u -> s_i], b_i = [v -> s_i].  H by brute-force subset DP
(n <= 11).  Families: b = f(a) letterwise or by a fixed pattern.  Printed: the set of values against simple closed
forms in the integer A = (a_0 ... a_(n-1))_2.
"""
import sys
from collections import Counter


def ham(n, a, b, uv=1):
    V = n + 2
    u, v = n, n + 1
    adj = [[0] * V for _ in range(V)]
    for i in range(n):
        for j in range(i + 1, n):
            adj[i][j] = 1
        if a[i]:
            adj[u][i] = 1
        else:
            adj[i][u] = 1
        if b[i]:
            adj[v][i] = 1
        else:
            adj[i][v] = 1
    if uv:
        adj[u][v] = 1
    else:
        adj[v][u] = 1
    full = (1 << V) - 1
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
    return sum(dp[full])


def bits(A, n):
    return [(A >> (n - 1 - i)) & 1 for i in range(n)]      # a_0 = most significant bit


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    fams = {
        'twin b=a': lambda a: a,
        'anti-twin b=1-a': lambda a: [1 - t for t in a],
        'b=0': lambda a: [0] * len(a),
        'b=1': lambda a: [1] * len(a),
        'b=shift (b_i=a_(i+1), last 0)': lambda a: a[1:] + [0],
        'b=shift (b_i=a_(i-1), first 1)': lambda a: [1] + a[:-1],
        'b=1-shift (b_i=1-a_(i+1), last 1)': lambda a: [1 - t for t in a[1:]] + [1],
        'b=1-shift (b_i=1-a_(i-1), first 0)': lambda a: [0] + [1 - t for t in a[:-1]],
    }
    for name, f in fams.items():
        print("=== family", name)
        for n in range(1, nmax + 1):
            vals = {}
            for A in range(1 << n):
                a = bits(A, n)
                vals[A] = ham(n, a, f(a))
            vs = sorted(set(vals.values()))
            odd_missing = [x for x in range(1, max(vs) + 1, 2) if x not in set(vs)]
            print(f"  n = {n}: {len(vs)} distinct of {1 << n}; max {max(vs)}; "
                  f"H(A) for A = 0..: {[vals[A] for A in range(min(1 << n, 16))]}"
                  f"{' ...' if n > 4 else ''}; missing odd below max: {len(odd_missing)}")


if __name__ == "__main__":
    main()
