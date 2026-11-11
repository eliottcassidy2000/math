#!/usr/bin/env python3
"""S92 (THM-5134 (G)): the join family F1 in every base b.

For a tournament T let Omega(T) be the conflict graph of its directed odd cycles (two cycles conflict when they share
a vertex) and I(Omega, z) = sum_k alpha_k z^k its independence polynomial (alpha_k = number of sets of k pairwise
disjoint odd cycles).  H(T) = I(Omega, 2) (the odd-cycle formula); H_b(T) = I(Omega, b) is the count "in base b"
(THM-4940).  When every alpha_k <= b - 1 the base-b digits of H_b(T) are alpha_0, alpha_1, alpha_2, ...: no carry.

Claims checked here on explicit tournaments (exhaustive enumeration of the odd cycles and of their packings):
  (1) one-apex R(h):  I = 1 + psi(h) z                       (all odd cycles pass through the apex);
  (2) join J(h1, h2, s) = R(h1) => TT_s => R(h2)^converse with the arc between the apexes reversed:
          I = (1 + psi1 z)(1 + psi2 z) + c3 z,     c3 = (2^s y1 y2 - [s = 0]) / 2,
      c3 being the number of odd cycles through the reversed arc; hence
          N_b(h1, h2, s) = (1 + b psi1)(1 + b psi2) + b c3,        N_2 = x1 x2 + 2^s y1 y2 - [s = 0].
  (3) balanced states 1 + b psi(h) = 2h + 1 (X_b = y): the powers of two for b = 2; h = 3, 45, 117, 309, ... for
      b = 3; none but h = 0 for b >= 4 (psi(h) > h/2 for h >= 1).
  (4) by value, ALL triples with N_b <= 2^24 on the sheet s >= 1 (until the audit of 2026-11-11 the census was
      truncated at h < 2^21 and used an unsound skip), for b = 2, 3, 4, 5, 6, 10 and the primes p = 3, 5, 7, 11, 13: N_b = 1
      modulo b always; for p not dividing b every triple with a p-divisible smaller state (p | X_b and p | y) lies
      in the class 0 (exact), and the class factors follow the first-order law (class 0: 1 + (p - 1) F, the others:
      1 - F, F the share of those triples) with a remainder below 0.2 (largest at the slow prime 7); the second
      digit (N_b - 1)/b modulo b is flat up to 8 % (slowest at the bases divisible by 3).
  (5) three-digit readings: the joins with alpha_1 <= b - 1 and alpha_2 <= b - 1 in base 10.
Exit status 1 on failure.
"""
import itertools
import math
import sys

import numpy as np
from numba import njit


def psi(h):
    if h == 0:
        return 0
    k = h >> 1
    return psi(k) + (k if h % 2 == 0 else (k >> 1) + 1)


def one_apex(h):
    u = 2 * h
    m = u.bit_length()
    adj = [[0] * (m + 1) for _ in range(m + 1)]
    for i in range(m):
        for j in range(i + 1, m):
            adj[i][j] = 1
        if (u >> (m - 1 - i)) & 1:
            adj[m][i] = 1
        else:
            adj[i][m] = 1
    return adj, m


def f1_join(h1, h2, s):
    a1, m1 = one_apex(h1)
    a2, m2 = one_apex(h2)
    n1, n2 = m1 + 1, m2 + 1
    n = n1 + s + n2
    adj = [[0] * n for _ in range(n)]
    for i in range(n1):
        for j in range(n1):
            adj[i][j] = a1[i][j]
    for i in range(n2):
        for j in range(n2):
            adj[n1 + s + i][n1 + s + j] = a2[j][i]          # converse
    for i in range(n1):
        for j in range(n1, n):
            adj[i][j] = 1
    for i in range(n1, n1 + s):
        for j in range(i + 1, n):
            adj[i][j] = 1
    ap1, ap2 = m1, n1 + s + m2
    adj[ap1][ap2] = 0
    adj[ap2][ap1] = 1
    return adj, n


def odd_cycles(adj, n):
    """list of vertex masks, one entry per directed odd cycle (cycles counted from their least vertex)."""
    out = []
    for start in range(n):
        # paths from start through vertices > start
        stack = [(start, 1 << start, 1)]
        while stack:
            v, mask, length = stack.pop()
            for w in range(start + 1, n):
                if adj[v][w] and not (mask >> w) & 1:
                    stack.append((w, mask | (1 << w), length + 1))
            if length >= 3 and length % 2 == 1 and adj[v][start]:
                out.append(mask)
    return out


def independence_poly(masks):
    """alpha_k = number of k-sets of pairwise disjoint masks."""
    alpha = [1]
    level = [(m, i) for i, m in enumerate(masks)]
    while level:
        alpha.append(len(level))
        nxt = []
        for m, i in level:
            for j in range(i + 1, len(masks)):
                if not (m & masks[j]):
                    nxt.append((m | masks[j], j))
        level = nxt
    return alpha


def ham_paths(adj, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = dp[S][v]
            if c:
                for w in range(n):
                    if adj[v][w] and not (S >> w) & 1:
                        dp[S | (1 << w)][w] += c
    return sum(dp[(1 << n) - 1])


@njit(cache=True)
def census_b(q, b, X, K, primes, cnt, div, digit2):
    """sheet s >= 1, N_b = (1 + b psi1)(1 + b psi2) + b 2^(s-1) y1 y2 <= X, ordered pairs h1 <= h2 counted twice.
    cnt[j, c]: classes modulo primes[j]; div[j, 0/1]: triples with p-divisible smaller state, and those of them in
    the class 0; digit2[c]: classes of (N_b - 1)/b modulo b (the second digit)."""
    H = q.shape[0]
    bad = 0
    open_loops = 0
    for s in range(1, K + 1):
        coef = b << (s - 1)
        if coef > X:
            break
        for h1 in range(H):
            y1 = 2 * h1 + 1
            X1 = 1 + b * (h1 - q[h1])
            # sound bounds only (audit, MISTAKE of 2026-11-11): psi(h) >= h // 2, so for h2 >= h1
            # N >= X1 (1 + b (h2 // 2)) + coef y1 y2 >= X1 (1 + b (h1 // 2)) + coef y1^2; the diagonal value
            # N(h1, h1) is NOT a lower bound (psi is not monotone).
            if X1 * (1 + b * (h1 // 2)) + coef * y1 * y1 > X:
                if (1 + b * (h1 // 2)) * (1 + b * (h1 // 2)) + coef * y1 * y1 > X:
                    break
                continue
            closed = False
            for h2 in range(h1, H):
                y2 = 2 * h2 + 1
                if (1 + b * (h2 // 2)) * X1 + coef * y1 * y2 > X:
                    closed = True
                    break
                X2 = 1 + b * (h2 - q[h2])
                N = X1 * X2 + coef * y1 * y2
                if N > X:
                    continue
                w = 1 if h2 == h1 else 2
                if N % b != 1:
                    bad += 1
                digit2[((N - 1) // b) % b] += w
                for j in range(primes.shape[0]):
                    p = primes[j]
                    cnt[j, N % p] += w
                    if h2 != h1 and X1 % p == 0 and y1 % p == 0:
                        div[j, 0] += w
                        if N % p == 0:
                            div[j, 1] += w
            if not closed:
                open_loops += 1          # the table of states is too short: the census would be truncated
    if open_loops > 0:
        return -1
    return bad


def main():
    ok = True
    # (1), (2): independence polynomials on explicit tournaments
    good1 = True
    for h in range(64):
        adj, m = one_apex(h)
        al = independence_poly(odd_cycles(adj, m + 1))
        good1 &= al == ([1, psi(h)] if psi(h) else [1])
        good1 &= sum(a << k for k, a in enumerate(al)) == 2 * psi(h) + 1
    print(f"(1) one-apex R(h), h < 64: I(Omega, z) = 1 + psi(h) z: {good1}")
    ok &= good1
    good2 = True
    nj = 0
    for h1 in range(9):
        for h2 in range(9):
            for s in range(4):
                adj, n = f1_join(h1, h2, s)
                if n > 13:
                    continue
                cyc = odd_cycles(adj, n)
                al = independence_poly(cyc)
                p1, p2 = psi(h1), psi(h2)
                c3 = (((2 * h1 + 1) * (2 * h2 + 1)) << s) - (1 if s == 0 else 0)
                assert c3 % 2 == 0
                c3 //= 2
                want = [1, p1 + p2 + c3, p1 * p2]
                while want and want[-1] == 0:
                    want.pop()
                good2 &= al == want
                if nj % 37 == 0:
                    good2 &= ham_paths(adj, n) == sum(a << k for k, a in enumerate(al))
                nj += 1
    print(f"(2) joins J(h1, h2, s) with at most 13 vertices ({nj} joins): "
          f"I(Omega, z) = (1 + psi1 z)(1 + psi2 z) + c3 z, c3 = (2^s y1 y2 - [s = 0])/2: {good2}")
    ok &= good2
    # (3) balanced states
    NB = 24
    q = np.zeros(1 << NB, dtype=np.int64)
    for k in range(NB):
        lo, hi = 1 << k, 1 << (k + 1)
        idx = np.arange(lo, hi, dtype=np.int64)
        half = idx >> 1
        q[lo:hi] = q[half] + np.where((idx & 1) == 1, (half + 1) >> 1, 0)
    h = np.arange(1 << NB, dtype=np.int64)
    ps = h - q
    print("(3) balanced states 1 + b psi(h) = 2h + 1 below 2^24:")
    for b in (2, 3, 4, 5, 6):
        sol = h[b * ps == 2 * h]
        print(f"    b = {b}: {len(sol)} states: {sol[:12].tolist()}{' ...' if len(sol) > 12 else ''}")
        if b == 2:
            ok &= all(int(v) & (int(v) - 1) == 0 for v in sol.tolist())
        if b >= 4:
            ok &= sol.tolist() == [0]
    ok &= bool(np.all(2 * ps[1:] > h[1:]))
    print(f"    psi(h) > h/2 for every 1 <= h < 2^24: {bool(np.all(2 * ps[1:] > h[1:]))}")
    # (4) census by value in several bases
    K = 24
    X = 1 << K
    primes = np.array([3, 5, 7, 11, 13], dtype=np.int64)
    q32 = q[: (1 << 22) + 2].astype(np.int64)          # every word that can occur: coef y1 y2 <= X gives h2 < X/4
    print(f"(4) sheet s >= 1, N_b <= 2^{K}: class factors p * share at the class 0 and [min, max] over the others; "
          f"F p^2 = share of p-divisible smaller states times p^2; rem = remainder of the first-order law")
    for b in (2, 3, 4, 5, 6, 10):
        cnt = np.zeros((len(primes), 16), np.int64)
        div = np.zeros((len(primes), 2), np.int64)
        digit2 = np.zeros(b, np.int64)
        bad = census_b(q32, b, X, K, primes, cnt, div, digit2)
        if bad < 0:
            print(f"    b = {b}: the table of states is too short (truncated census)")
            ok = False
            continue
        tot = cnt[0].sum()
        cells = []
        for j, p in enumerate(primes.tolist()):
            f = p * cnt[j, :p] / tot
            if b % p == 0:
                # N = 1 modulo p always
                ok &= cnt[j, 1] == tot
                cells.append(f"p={p}: N = 1 (mod p)")
            else:
                ok &= div[j, 0] == div[j, 1]                       # divisible smaller state => class 0 (exact)
                fs = div[j, 0] / tot                               # share of triples with p-divisible smaller state
                pred0, predr = 1 + (p - 1) * fs, 1 - fs            # first-order law (THM-5128 (C2) in base b)
                rem = max(abs(f[0] - pred0), abs(f[1:] - predr).max())
                ok &= rem < 0.2
                cells.append(f"p={p}: {f[0]:.3f} [{f[1:].min():.3f}, {f[1:].max():.3f}] F p^2={fs * p * p:.2f} rem={rem:.3f}")
        d2 = b * digit2 / digit2.sum()
        ok &= bad == 0 and (b == 2 or abs(d2 - 1).max() < 0.08)
        print(f"    b = {b:2d}: {int(tot):9d} triples; N = 1 (mod b): {bad == 0}; second digit flat to "
              f"{abs(d2 - 1).max():.3f}; " + "; ".join(cells))
    # (5) three-digit readings in base 10
    seen = {}
    for h1 in range(40):
        for h2 in range(h1, 40):
            for s in range(0, 6):
                p1, p2 = psi(h1), psi(h2)
                c3 = ((((2 * h1 + 1) * (2 * h2 + 1)) << s) - (1 if s == 0 else 0)) // 2
                a1, a2 = p1 + p2 + c3, p1 * p2
                if a1 <= 9 and a2 <= 9:
                    seen.setdefault((a2, a1), []).append((h1, h2, s))
    print(f"(5) base 10 without carry: the joins give the decimal strings (alpha_2)(alpha_1)1 for (alpha_2, alpha_1) in "
          f"{sorted(seen)}")
    print(f"    e.g. {[(k, v[0]) for k, v in sorted(seen.items())][:8]}")
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
