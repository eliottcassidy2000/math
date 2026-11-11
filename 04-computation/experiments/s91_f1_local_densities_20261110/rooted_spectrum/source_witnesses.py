#!/usr/bin/env python3
"""S91: explicit singleton-marked sources for the 9-grid at levels s = 7, 8, 9 (open in the Session261 packet).

A source is a rooted tournament (B, y); D = H(B), U = Yminus(B, y) = number of ordered covers (P, R) of B by
vertex-disjoint directed paths with R starting at y (P may be empty).  The tail tournament T_t(B, y), t >= 1, has
vertices B, sinks s_1 -> ... -> s_t (every vertex of B beats every s_i) and an apex z with z -> y and every other
vertex -> z.  Tail law (Session261, constructors/NINE_OFFSET_EXTENSION.md; THM-5126 (1')):  H(T_t) = D + 2^t U.

Search (numba): every rooted tournament of order m <= 9 is (class representative of order m - 1, cut of the new
vertex), root = the new vertex; the first hit of each target (D, U) is recorded in increasing order.
Targets: D = 2^s + 1, U = 9 * 2^j with j <= s, s = 2..9 (levels 2..6 are Session261's; 7, 8, 9 are new here).

Verification of every witness, independent of the search kernel (pure Python, exact integers):
  * H(B) = D by the subset DP and, for order <= 9, by direct enumeration of all permutations;
  * U from the definition (enumeration of all directed paths R from y, U = sum of H(B - V(R)));
  * the explicit tail tournaments T_t for t = 1, 2, 3: H(T_t) = D + 2^t U by the subset DP;
  * the grid values N(9, r, s) = 1 + 2^s + 9 * 2^(r + s) for r = 1, 2 are H(T_t) with t = r + s - j.
Exit status 1 on any failure.  Input: reps_le8.npz (class representatives of orders 1..8, colex pair coding of
core.py: for i < j the bit j(j-1)/2 + i is 1 iff i -> j).
"""
import itertools
import os
import sys

import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))


@njit(cache=True)
def scan(n, reps, targets, found, counts):
    """targets[k] = (D, U); found[k] = (code, cut) of the first rooted tournament of order n + 1 hitting it."""
    full = (1 << n) - 1
    f = np.zeros((1 << n, n), np.int64)
    g = np.zeros((1 << n, n), np.int64)
    hS = np.zeros(1 << n, np.int64)
    out = np.zeros(n, np.int64)
    Q = np.zeros((n, n), np.int64)
    ym = np.zeros(n, np.int64)
    nhit = 0
    for r in range(reps.shape[0]):
        code = reps[r]
        for i in range(n):
            out[i] = 0
        for j in range(1, n):
            base = j * (j - 1) // 2
            for i in range(j):
                if (code >> (base + i)) & 1:
                    out[i] |= (1 << j)
                else:
                    out[j] |= (1 << i)
        hS[0] = 1
        for S in range(1, full + 1):
            tot = 0
            for v in range(n):
                f[S, v] = 0
                g[S, v] = 0
                if ((S >> v) & 1) == 0:
                    continue
                R = S ^ (1 << v)
                if R == 0:
                    f[S, v] = 1
                    g[S, v] = 1
                    tot += 1
                    continue
                a = 0
                b = 0
                for w in range(n):
                    if (R >> w) & 1:
                        if (out[w] >> v) & 1:
                            a += f[R, w]
                        else:
                            b += g[R, w]
                f[S, v] = a
                g[S, v] = b
                tot += a
            hS[S] = tot
        for a in range(n):
            ym[a] = 0
            for b in range(n):
                Q[a, b] = 0
        for S in range(1, full + 1):
            C = full ^ S
            for a in range(n):
                if (S >> a) & 1:
                    ym[a] += g[S, a] * hS[C]
                    if f[S, a]:
                        for b in range(n):
                            if (C >> b) & 1:
                                Q[a, b] += f[S, a] * g[C, b]
        HB = hS[full]
        for cut in range(1 << n):
            H = 0
            U = HB
            for b in range(n):
                if (cut >> b) & 1:
                    H += g[full, b]
                    U += ym[b]
                else:
                    H += f[full, b]
                    for c in range(n):
                        if (cut >> c) & 1:
                            H += Q[b, c]
            for k in range(targets.shape[0]):
                if targets[k, 0] == H and targets[k, 1] == U:
                    counts[k] += 1
                    if found[k, 0] < 0:
                        found[k, 0] = code
                        found[k, 1] = cut
                        nhit += 1
    return nhit


def adjacency(code, n, cut):
    """out-masks of the tournament of order n + 1: base coded by `code`, new vertex n beats exactly `cut`."""
    out = [0] * (n + 1)
    for j in range(1, n):
        base = j * (j - 1) // 2
        for i in range(j):
            if (code >> (base + i)) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    for b in range(n):
        if (cut >> b) & 1:
            out[n] |= 1 << b
        else:
            out[b] |= 1 << n
    return out


def ham_dp(out, verts):
    """number of Hamiltonian paths of the sub-tournament on the vertex list `verts` (subset DP, exact)."""
    m = len(verts)
    if m == 0:
        return 1
    adj = [[(out[verts[i]] >> verts[j]) & 1 for j in range(m)] for i in range(m)]
    dp = [[0] * m for _ in range(1 << m)]
    for v in range(m):
        dp[1 << v][v] = 1
    for S in range(1, 1 << m):
        row = dp[S]
        for v in range(m):
            c = row[v]
            if c:
                for w in range(m):
                    if not (S >> w) & 1 and adj[v][w]:
                        dp[S | (1 << w)][w] += c
    return sum(dp[(1 << m) - 1])


def ham_perm(out, n):
    cnt = 0
    for p in itertools.permutations(range(n)):
        ok = True
        for i in range(n - 1):
            if not (out[p[i]] >> p[i + 1]) & 1:
                ok = False
                break
        cnt += ok
    return cnt


def slope_by_definition(out, n, y):
    """U = sum over directed paths R starting at y of H(B - V(R)) (R = (y) included)."""
    cache = {}

    def H(mask):
        if mask not in cache:
            cache[mask] = ham_dp(out, [v for v in range(n) if (mask >> v) & 1])
        return cache[mask]

    full = (1 << n) - 1
    total = 0
    stack = [(y, 1 << y)]
    while stack:
        v, used = stack.pop()
        total += H(full ^ used)
        m = out[v] & ~used & full
        while m:
            w = (m & -m).bit_length() - 1
            m &= m - 1
            stack.append((w, used | (1 << w)))
    return total


def tail_tournament(out, n, y, t):
    """out-masks of T_t(B, y): vertices 0..n-1 = B, n..n+t-1 = sinks in transitive order, n+t = apex z."""
    N = n + t + 1
    z = n + t
    o = [0] * N
    for v in range(n):
        o[v] = out[v]
        for i in range(t):
            o[v] |= 1 << (n + i)
        if v != y:
            o[v] |= 1 << z
    for i in range(t):
        for k in range(i + 1, t):
            o[n + i] |= 1 << (n + k)
        o[n + i] |= 1 << z
    o[z] = 1 << y
    return o


def main():
    z = np.load(os.path.join(HERE, "reps_le8.npz"))
    reps = {n: z["n%d" % n].astype(np.int64) for n in range(1, 9)}
    targets = []
    for s in range(2, 10):
        for j in range(0, s + 1):
            targets.append(((1 << s) + 1, 9 << j, s, j))
    T = np.array([[d, u] for d, u, _, _ in targets], dtype=np.int64)
    found = -np.ones((len(targets), 2), dtype=np.int64)
    order = [0] * len(targets)
    hits = {}
    for n in range(1, 9):
        before = found[:, 0].copy()
        counts = np.zeros(len(targets), dtype=np.int64)
        scan(n, reps[n], T, found, counts)
        for k in range(len(targets)):
            if counts[k]:
                hits.setdefault(k, {})[n + 1] = int(counts[k])
        for k in range(len(targets)):
            if before[k] < 0 <= found[k, 0]:
                order[k] = n + 1
    ok = True
    print("level s: least j <= s with a source (2^s + 1, 9 * 2^j) of order <= 9; least order; witness")
    best = {}
    for k, (d, u, s, j) in enumerate(targets):
        if found[k, 0] >= 0 and s not in best:
            best[s] = k
    for s in range(2, 10):
        if s not in best:
            print(f"  s = {s}: no source with j <= s at order <= 9")
            ok = False
            continue
        k = best[s]
        d, u, _, j = targets[k]
        n = order[k] - 1
        out = adjacency(int(found[k, 0]), n, int(found[k, 1]))
        m = n + 1
        y = n
        rows = "(" + ",".join(str(x) for x in out) + ")"
        print(f"  s = {s}: D = {d}, U = {u} = 9 * 2^{j}, order {m}, root y = {y}, out-masks B = {rows}")
        print(f"      rooted tournaments (class of the base, cut) with this pair, by order: {hits.get(k, {})}")
        hb = ham_dp(out, list(range(m)))
        hp = ham_perm(out, m)
        ub = slope_by_definition(out, m, y)
        good = (hb == d and hp == d and ub == u)
        tails = []
        for t in (1, 2, 3):
            ht = ham_dp(tail_tournament(out, m, y, t), list(range(m + t + 1)))
            tails.append(ht)
            good &= (ht == d + (u << t))
        grid = []
        for r in (1, 2):
            t = r + s - j
            target = 1 + (1 << s) + 9 * (1 << (r + s))
            if t <= 3:
                good &= tails[t - 1] == target
                grid.append((r, t, target, "explicit"))
            else:
                good &= d + (u << t) == target
                grid.append((r, t, target, "tail law"))
        print(f"      H(B): DP {hb}, permutations {hp}; U by definition {ub}; H(T_1), H(T_2), H(T_3) = {tails} "
              f"(D + 2^t U = {[d + (u << t) for t in (1, 2, 3)]}); grid N(9, r, {s}), r = 1, 2: {grid}  "
              f"{'ok' if good else 'FAIL'}")
        ok &= good
    # all pairs (s, j) found, as a table (least order)
    print("orders of the first witness for every (s, j), j <= s ('.' = none at order <= 9):")
    for s in range(2, 10):
        cells = []
        for k, (d, u, s2, j) in enumerate(targets):
            if s2 == s:
                cells.append(f"j={j}:{order[k] if found[k, 0] >= 0 else '.'}")
        print(f"  s = {s}: " + " ".join(cells))
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
