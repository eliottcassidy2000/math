"""Referee side check (H): feedback arc numbers of the examples behind 'every extremal and failing example so far is
nearly transitive' (results file 11.2) and of the explicit examples found in this audit.

Exact minimum feedback arc set of a tournament by DP over vertex subsets (placing vertices from the top of the order:
f(S) = min_v f(S - v) + #{x in S - v : v -> x}).  Cross-checked against brute force over all n! orders for n <= 9.
Adjacency convention of the session scripts: bit j of adj[i] set  <=>  i -> j.
"""
import itertools, time
import numpy as np


def fas_dp(n, adj):
    f = [0] * (1 << n)
    for S in range(1, 1 << n):
        best = 10 ** 9
        T = S
        while T:
            lb = T & -T
            v = lb.bit_length() - 1
            T ^= lb
            R = S ^ lb
            c = f[R] + bin(adj[v] & R).count('1')
            if c < best:
                best = c
        f[S] = best
    return f[(1 << n) - 1]


def fas_brute(n, adj):
    best = 10 ** 9
    for p in itertools.permutations(range(n)):
        c = 0
        for i in range(n):
            for j in range(i):
                if (adj[p[i]] >> p[j]) & 1:     # later vertex beats an earlier one: backward arc
                    c += 1
        if c < best:
            best = c
    return best


def from_out_sets(sets):
    return [sum(1 << x for x in s) for s in sets]


def transitive_plus_adj(n, arcs):
    adj = [0] * n
    for a in range(n):
        for b in range(a):
            adj[a] |= 1 << b
    for lo, hi in arcs:
        adj[hi] &= ~(1 << lo)
        adj[lo] |= 1 << hi
    return adj


if __name__ == '__main__':
    t0 = time.time()
    ex = [
        ('(SS) n = 9 counterexample 1 (session; audit class 2)', 9, [502, 508, 504, 241, 480, 448, 384, 256, 8]),
        ('(SS) n = 9 counterexample 2 (session; audit class 4)', 9, [250, 508, 497, 500, 480, 448, 384, 256, 1]),
        ('D(11,5) double chain', 11, transitive_plus_adj(11, [(0, 5), (5, 10)])),
        ('W <= 0 example, n = 10, arcs (0,4),(1,9)', 10, transitive_plus_adj(10, [(0, 4), (1, 9)])),
        ('G some-fail example, n = 8, arcs (0,5),(1,6),(5,7)', 8, transitive_plus_adj(8, [(0, 5), (1, 6), (5, 7)])),
        ('A some-fail example n = 7 (labelled code 551)', 7, from_out_sets([[1, 2, 3, 6], [5], [1], [1, 2], [0, 1, 2, 3], [0, 2, 3, 4], [1, 2, 3, 4, 5]])),
        ('P some-fail example n = 7 (labelled code 161)', 7, from_out_sets([[1, 6], [3], [0, 1], [0, 2], [0, 1, 2, 3], [0, 1, 2, 3, 4], [1, 2, 3, 4, 5]])),
    ]
    for name, n, adj in ex:
        # sanity: it is a tournament
        ok = all(((adj[i] >> j) & 1) + ((adj[j] >> i) & 1) == 1 for i in range(n) for j in range(i + 1, n)) and all(not (adj[i] >> i) & 1 for i in range(n))
        f = fas_dp(n, adj)
        fb = fas_brute(n, adj) if n <= 9 else None
        print(f'{name}: tournament {ok}; feedback arc number (subset DP) = {f}' + (f'; brute force over {n}! orders = {fb}' if fb is not None else ''), flush=True)
    print(f'[{time.time() - t0:.1f} s]')
