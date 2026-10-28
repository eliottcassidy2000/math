"""Referee: THM-4836 (coherent carrier sets).  Enumerates tournament types n <= 7 from scratch,
computes Sigma(T) subset {L, C, S}, tallies classes, checks the formulas, and computes t_n
(Burnside over cycle types) and s_n (S = 1 - 1/T) through n = 9.  Run: python tournaments.py"""
import itertools
import time
from fractions import Fraction
from math import factorial, gcd
from collections import Counter

import numpy as np


def canon_factory(n):
    perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    idx = np.full((n, n), -1, dtype=np.int64)
    for k, (i, j) in enumerate(pairs):
        idx[i, j] = idx[j, i] = k
    pm = np.stack([idx[perms[:, i], perms[:, j]] for (i, j) in pairs], axis=1)
    flip = np.stack([(perms[:, i] > perms[:, j]).astype(np.int64) for (i, j) in pairs], axis=1)

    def canon(adj):
        b = np.array([adj[i][j] for (i, j) in pairs], dtype=np.int64)
        codes = ((b[None, :] ^ flip) << pm).sum(axis=1)
        return int(codes.min())
    return canon, pairs


def decode(code, n, pairs):
    adj = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if (code >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def types(n):
    if n == 1:
        return [[[0]]]
    prev = types(n - 1)
    canon, pairs = canon_factory(n)
    seen = {}
    for T in prev:
        for outset in range(2 ** (n - 1)):
            adj = [row[:] + [0] for row in T] + [[0] * n]
            for v in range(n - 1):
                if (outset >> v) & 1:
                    adj[n - 1][v] = 1
                else:
                    adj[v][n - 1] = 1
            c = canon(adj)
            if c not in seen:
                seen[c] = adj
    return list(seen.values())


def ham_path(adj, verts):
    verts = list(verts)
    m = len(verts)
    if m == 1:
        return True
    dp = [[False] * m for _ in range(1 << m)]
    for i in range(m):
        dp[1 << i][i] = True
    for mask in range(1 << m):
        for i in range(m):
            if dp[mask][i]:
                for j in range(m):
                    if not (mask >> j) & 1 and adj[verts[i]][verts[j]]:
                        dp[mask | (1 << j)][j] = True
    return any(dp[(1 << m) - 1])


def ham_cycle(adj, verts):
    verts = list(verts)
    m = len(verts)
    if m < 3:
        return False
    dp = [[False] * m for _ in range(1 << m)]
    dp[1][0] = True
    for mask in range(1 << m):
        if not mask & 1:
            continue
        for i in range(m):
            if dp[mask][i]:
                for j in range(m):
                    if not (mask >> j) & 1 and adj[verts[i]][verts[j]]:
                        dp[mask | (1 << j)][j] = True
    full = (1 << m) - 1
    return any(dp[full][i] and adj[verts[i]][verts[0]] for i in range(1, m))


def condensation_sizes(adj, n):
    reach = [[bool(adj[i][j]) or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    if reach[k][j]:
                        reach[i][j] = True
    comps = {}
    for v in range(n):
        key = frozenset(w for w in range(n) if reach[v][w] and reach[w][v])
        comps[key] = True
    comps = list(comps)
    # order: component X before Y if X -> Y
    comps.sort(key=lambda X: -sum(1 for Y in comps if Y != X and reach[next(iter(X))][next(iter(Y))]))
    return [len(X) for X in comps]


def sigma(adj, n):
    S = set()
    if ham_path(adj, range(n)):
        S.add("L")
    if any(ham_cycle(adj, [w for w in range(n) if w != v]) for v in range(n)):
        S.add("C")
    outd = [sum(r) for r in adj]
    if any(d == n - 1 for d in outd) or any(d == 0 for d in outd):
        S.add("S")
    return frozenset(S)


def burnside_t(n):
    """number of tournament types on n vertices (Davis 1954)."""
    tot = Fraction(0)
    for lam in partitions(n):
        if any(p % 2 == 0 for p in lam):
            continue
        z = 1
        for k, c in Counter(lam).items():
            z *= k ** c * factorial(c)
        orbits = sum((p - 1) // 2 for p in lam) + sum(gcd(lam[i], lam[j])
                                                      for i in range(len(lam)) for j in range(i + 1, len(lam)))
        tot += Fraction(2 ** orbits, z)
    return int(tot)


def partitions(n, maxp=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield []
        return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p):
            yield [p] + rest


if __name__ == "__main__":
    tser = [1] + [burnside_t(n) for n in range(1, 11)]
    # S = 1 - 1/T
    inv = [Fraction(0)] * 11
    inv[0] = Fraction(1)
    for k in range(1, 11):
        inv[k] = -sum(tser[i] * inv[k - i] for i in range(1, k + 1))
    sser = [0] + [int(-inv[k]) for k in range(1, 11)]
    print("t_n (n=0..10):", tser)
    print("s_n (n=1..10):", sser[1:])
    pred = {}
    for n in range(4, 11):
        S3sq_T = 0
        S3 = [0, 0, 0] + sser[3:]
        for a in range(3, n + 1):
            for b in range(3, n + 1 - a):
                S3sq_T += S3[a] * S3[b] * tser[n - a - b]
        pred[n] = (sser[n], 2 * sser[n - 1], 2 * tser[n - 1] - tser[n - 2] - 2 * sser[n - 1], S3sq_T,
                   tser[n] - sser[n] - 2 * tser[n - 1] + tser[n - 2])
        print(f"formula n={n}: {{L,C}}={pred[n][0]} {{L,C,S}}={pred[n][1]} {{L,S}}={pred[n][2]} "
              f"{{L}}=[x^n]S3^2T={pred[n][3]} (t-s-2t'+t''={pred[n][4]}) total={sum(pred[n][:4])} t_n={tser[n]}")
    for n in range(4, 8):
        t0 = time.time()
        Ts = types(n)
        tally = Counter()
        sing_ok = True
        strong_count = 0
        for adj in Ts:
            sg = sigma(adj, n)
            tally[tuple(sorted(sg))] += 1
            sizes = condensation_sizes(adj, n)
            strong = len(sizes) == 1
            strong_count += strong
            predicted_singleton = (not strong) and sizes[0] >= 3 and sizes[-1] >= 3
            if (sg == frozenset(["L"])) != predicted_singleton:
                sing_ok = False
            if "L" not in sg:
                print("NO HAMILTONIAN PATH?!")
        print(f"enumeration n={n}: types {len(Ts)} strong {strong_count} tally {dict(tally)} "
              f"singleton characterisation holds: {sing_ok} [{time.time() - t0:.1f}s]", flush=True)
