#!/usr/bin/env python3
"""Verifier: balanced arcs (M9) and module law / c_22 (M4, M5), own code.  usage: python v4_balanced_modules.py"""
import sys, time, random, itertools
from collections import Counter
import numpy as np
from numba import njit, prange
from vlib import *


@njit(cache=True)
def bal_one(n, m, pi, pj, res):
    adj = np.zeros(n, np.int64); radj = np.zeros(n, np.int64)
    adj_from_mask(n, m, pi, pj, adj, radj)
    f = np.zeros((1 << n, n), np.int64); g = np.zeros((1 << n, n), np.int64); Q = np.zeros((n, n), np.int64)
    h = fgQ(n, adj, f, g, Q)
    nb = 0; na_ = 0; tot = 0
    for x in range(n):
        for y in range(n):
            if x != y and (adj[x] >> y) & 1:
                mu = Q[x, y]
                tot += mu
                if 2 * mu == h - 1:
                    nb += 1
                if 2 * mu == h + 1:
                    na_ += 1
    res[0] = h; res[1] = nb; res[2] = na_; res[3] = is_strong(n, adj, radj); res[4] = tot


@njit(parallel=True, cache=True)
def bal_scan(n, pi, pj):
    N = 1 << (n * (n - 1) // 2)
    R = np.zeros((N, 5), np.int64)
    for m in prange(N):
        bal_one(n, m, pi, pj, R[m])
    return R


@njit(cache=True)
def walks(n, adj, u, v, J, K, out):
    """out[j, k] = number of walks visiting u exactly j times (1..J), v exactly k times (1..K), every other vertex once.
    (a, b, mask) lexicographic order is a topological order of the state graph."""
    full_others = ((1 << n) - 1) ^ (1 << u) ^ (1 << v)
    dp = np.zeros((J + 1, K + 1, 1 << n, n), np.int64)
    for w in range(n):
        if w == u:
            dp[1, 0, 0, w] = 1
        elif w == v:
            if K >= 1:
                dp[0, 1, 0, w] = 1
        else:
            dp[0, 0, 1 << w, w] = 1
    for a in range(J + 1):
        for b in range(K + 1):
            for mask in range(1 << n):
                if mask & ~full_others:
                    continue
                for last in range(n):
                    c = dp[a, b, mask, last]
                    if c == 0:
                        continue
                    for w in range(n):
                        if not (adj[last] >> w) & 1:
                            continue
                        if w == u:
                            if a + 1 <= J:
                                dp[a + 1, b, mask, w] += c
                        elif w == v:
                            if b + 1 <= K:
                                dp[a, b + 1, mask, w] += c
                        elif not (mask >> w) & 1:
                            dp[a, b, mask | (1 << w), w] += c
    for a in range(J + 1):
        for b in range(K + 1):
            t = 0
            for last in range(n):
                t += dp[a, b, full_others, last]
            out[a, b] = t


@njit(parallel=True, cache=True)
def c22_scan(n, pi, pj):
    N = 1 << (n * (n - 1) // 2)
    mn = np.full(N, 1 << 30, np.int64)
    st = np.zeros(N, np.int64)
    for m in prange(N):
        adj = np.zeros(n, np.int64); radj = np.zeros(n, np.int64)
        adj_from_mask(n, m, pi, pj, adj, radj)
        if not is_strong(n, adj, radj):
            continue
        st[m] = 1
        out = np.zeros((3, 3), np.int64)
        for u in range(n):
            for v in range(n):
                if u != v:
                    walks(n, adj, u, v, 2, 2, out)
                    if out[2, 2] < mn[m]:
                        mn[m] = out[2, 2]
    return mn, st


def ordered_path_partitions(k, A):
    """P[j] = number of ordered partitions of V (|V| = k, adjacency matrix A) into j directed paths, brute force."""
    P = [0] * (k + 1)
    for perm in itertools.permutations(range(k)):
        for breaks in range(1 << (k - 1)):
            ok = True
            for i in range(k - 1):
                if not (breaks >> i) & 1 and not A[perm[i]][perm[i + 1]]:
                    ok = False; break
            if ok:
                P[bin(breaks).count('1') + 1] += 1
    return P


def rand_tournament(k, rng):
    A = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def main():
    # ---- M9
    for n in range(4, 8):
        pi, pj = pair_arrays(n)
        t0 = time.time()
        R = bal_scan(n, pi, pj)
        H = R[:, 0]; nb = R[:, 1]; na_ = R[:, 2]; st = R[:, 3]
        assert np.all(R[:, 4] == (n - 1) * H)
        vals = sorted(set(H.tolist()))
        vb = sorted(set(H[nb > 0].tolist())); va = sorted(set(H[na_ > 0].tolist()))
        neither = [v for v in vals if v not in vb and v not in va]
        print(f"n = {n}: {len(H)} tournaments ({time.time()-t0:.1f}s); with a balanced arc {int((nb>0).sum())}; with an anti-"
              f"balanced arc {int((na_>0).sum())}; strong {int(st.sum())}; strong with balanced {int(((nb>0)&(st==1)).sum())}; "
              f"strong with anti-balanced {int(((na_>0)&(st==1)).sum())}")
        print(f"   H-values {len(vals)}; admitting balanced {len(vb)}; anti-balanced {len(va)}; neither {len(neither)}: "
              f"{neither}")
        print(f"   max # balanced arcs in one tournament {int(nb.max())}; anti-balanced {int(na_.max())}; "
              f"sum of arc usages = (n-1) H for every tournament: True")
        sys.stdout.flush()
    # ---- M5: min c_22 over strong labelled tournaments and ordered port pairs
    for n in range(3, 7):
        pi, pj = pair_arrays(n)
        t0 = time.time()
        mn, st = c22_scan(n, pi, pj)
        s = st == 1
        print(f"n = {n}: strong labelled tournaments {int(s.sum())}; min over all port pairs of c_22 = {int(mn[s].min())}; "
              f"tournaments attaining it {int((mn[s] == mn[s].min()).sum())} ({time.time()-t0:.1f}s)")
        sys.stdout.flush()
    # ---- M4: two-module law, random instances, everything brute force on the module side
    rng = random.Random(94)
    fails = 0; tests = 0
    dpbuf = {}
    for t in range(300):
        n0 = rng.randint(3, 6)
        A0 = rand_tournament(n0, rng)
        u, v = rng.sample(range(n0), 2)
        k1 = rng.randint(1, 3); k2 = rng.randint(1, 3)
        C1 = rand_tournament(k1, rng); C2 = rand_tournament(k2, rng)
        # build T = T0[u := C1, v := C2]
        verts = []
        for x in range(n0):
            if x == u:
                verts += [('u', i) for i in range(k1)]
            elif x == v:
                verts += [('v', i) for i in range(k2)]
            else:
                verts.append(('o', x))
        N = len(verts)

        def arc(p, q):
            (tp, ip), (tq, iq) = p, q
            if tp == 'u' and tq == 'u':
                return C1[ip][iq]
            if tp == 'v' and tq == 'v':
                return C2[ip][iq]
            xp = u if tp == 'u' else (v if tp == 'v' else ip)
            xq = u if tq == 'u' else (v if tq == 'v' else iq)
            return A0[xp][xq]
        adj = np.zeros(N, np.int64)
        for a in range(N):
            for b in range(N):
                if a != b and arc(verts[a], verts[b]):
                    adj[a] |= 1 << b
        dp = dpbuf.setdefault(N, np.zeros((1 << N, N), np.int64))
        Hdirect = int(hp(N, adj, dp))
        a0 = np.zeros(n0, np.int64)
        for a in range(n0):
            for b in range(n0):
                if A0[a][b]:
                    a0[a] |= 1 << b
        out = np.zeros((k1 + 1, k2 + 1), np.int64)
        walks(n0, a0, u, v, k1, k2, out)
        P1 = ordered_path_partitions(k1, C1); P2 = ordered_path_partitions(k2, C2)
        Hlaw = sum(int(out[j, k]) * P1[j] * P2[k] for j in range(1, k1 + 1) for k in range(1, k2 + 1))
        tests += 1
        if Hlaw != Hdirect:
            fails += 1
    print(f"module law: {tests} random instances (context order 3..6, modules of order 1..3): failures {fails}")
    # ---- M5 (ii): mixed second difference for transitive modules in strong contexts: must be > 0 for s, t >= 2
    bad = 0; cnt = 0
    for t in range(60):
        n0 = rng.randint(3, 5)
        while True:
            A0 = rand_tournament(n0, rng)
            a0 = np.zeros(n0, np.int64); r0 = np.zeros(n0, np.int64)
            for a in range(n0):
                for b in range(n0):
                    if A0[a][b]:
                        a0[a] |= 1 << b; r0[b] |= 1 << a
            if is_strong(n0, a0, r0):
                break
        u, v = rng.sample(range(n0), 2)

        def Hst(s, tt):
            N = n0 - 2 + s + tt
            verts = []
            for x in range(n0):
                if x == u:
                    verts += [('u', i) for i in range(s)]
                elif x == v:
                    verts += [('v', i) for i in range(tt)]
                else:
                    verts.append(('o', x))
            adj = np.zeros(N, np.int64)
            for a in range(N):
                for b in range(N):
                    if a == b:
                        continue
                    (tp, ip), (tq, iq) = verts[a], verts[b]
                    if tp == tq and tp in 'uv':
                        ar = ip < iq
                    else:
                        xp = u if tp == 'u' else (v if tp == 'v' else ip)
                        xq = u if tq == 'u' else (v if tq == 'v' else iq)
                        ar = A0[xp][xq]
                    if ar:
                        adj[a] |= 1 << b
            dp = dpbuf.setdefault(N, np.zeros((1 << N, N), np.int64))
            return int(hp(N, adj, dp))
        base = {(s, tt): Hst(s, tt) for s in (1, 2, 3) for tt in (1, 2, 3)}
        for s in (2, 3):
            for tt in (2, 3):
                d = base[(s, tt)] - base[(s, 1)] - base[(1, tt)] + base[(1, 1)]
                cnt += 1
                if d < (2 ** s - 2) * (2 ** tt - 2):
                    bad += 1
    print(f"transitive modules in strong contexts: {cnt} mixed differences H(s,t)-H(s,1)-H(1,t)+H(1,1), "
          f"violations of >= (2^s-2)(2^t-2): {bad}")


if __name__ == '__main__':
    main()
