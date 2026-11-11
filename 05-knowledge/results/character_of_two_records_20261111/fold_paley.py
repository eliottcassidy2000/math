#!/usr/bin/env python3
"""The folded (quotient) tournament of a Paley tournament: Z(P_q) on the pairs {a, -a}, a a non-zero square, with
a -> b iff exactly one of a -> b, a -> -b holds in P_q.  Checks Z = P_q[QR] through a -> a^2, and compares Z with the
Paley or rotational tournament of order (q - 1)/2 by scores, H and (for small orders) isomorphism.
Usage: fold_paley.py"""
import sys
from paley_tools import field, paley_out
def ham(out, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            m = out[v] & ~S
            while m:
                w = (m & -m).bit_length() - 1; dp[S | (1 << w)][w] += c; m &= m - 1
    return sum(dp[full])
def iso(A, B, n):
    arcA = [[(A[a] >> b) & 1 for b in range(n)] for a in range(n)]; arcB = [[(B[a] >> b) & 1 for b in range(n)] for a in range(n)]
    perm = [-1] * n; used = [False] * n
    def rec(i):
        if i == n: return True
        for c in range(n):
            if used[c]: continue
            if all(arcA[i][j] == arcB[c][perm[j]] for j in range(i)):
                perm[i] = c; used[c] = True
                if rec(i + 1): return True
                perm[i] = -1; used[c] = False
        return False
    return rec(0)
for q in (7, 11, 19, 23, 27, 31):
    add, neg, mul = field(q); out = paley_out(q)
    sq = sorted({mul(a, a) for a in range(1, q)}); m = len(sq); idx = {a: i for i, a in enumerate(sq)}
    arc = lambda a, b: (out[a] >> b) & 1
    Z = [0] * m
    for a in sq:
        for b in sq:
            if a != b and (arc(a, b) ^ arc(a, neg(b))): Z[idx[a]] |= 1 << idx[b]
    assert all(((Z[i] >> j) & 1) != ((Z[j] >> i) & 1) for i in range(m) for j in range(i + 1, m))
    # induced tournament on the squares, transported by a -> a^2
    ok = all(((Z[idx[a]] >> idx[b]) & 1) == arc(mul(a, a), mul(b, b)) for a in sq for b in sq if a != b)
    D1 = [a for a in sq if arc(a, neg(a))]
    HZ = ham(Z, m) if m <= 15 else None
    line = "q = %2d: folded tournament Z on %2d pairs; Z = P_q[squares] via a -> a^2: %s; scores %s; |D_1| = %d; H(Z) = %s" % (
        q, m, ok, sorted(set(bin(z).count('1') for z in Z)), len(D1), HZ)
    if m % 4 == 3 and m in (3, 7, 11, 19, 23):
        P = paley_out(m); line += "; isomorphic to P_%d: %s (H(P_%d) = %s)" % (m, iso(Z, P, m) if m <= 11 else "?", m, ham(P, m) if m <= 15 else "?")
    print(line)
