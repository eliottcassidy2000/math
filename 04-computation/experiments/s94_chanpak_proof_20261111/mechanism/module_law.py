#!/usr/bin/env python3
"""S94 lane "mechanism": checks for Theorem B of PROOFS.txt (two-module law and the mixed term).

(V1) Two-module law, random instances (exact):  H(T0[u := C1, v := C2]) = sum_{j,k} c_jk (j! p_j(C1)) (k! p_k(C2)),
     c_jk = number of walks of T0 that visit u exactly j times, v exactly k times and every other vertex once;
     j! p_j(C) = number of ordered partitions of V(C) into j nonempty directed paths.
(V2) Mixed term, exhaustive: for EVERY labelled strong tournament T0 of order n = 3..6 and every pair u != v,
     c_22(T0; u, v) >= 1 (Theorem B (b): >= 1 by Camion); reports the minimum and how often the minimum is attained.
     c_22 = (number of Hamiltonian paths of the digraph obtained by duplicating u and v into non-adjacent twins) / 4.
"""
import itertools, random, sys, time
import numpy as np
from numba import njit, prange

def ham_digraph(A):
    N = len(A)
    dp = [[0] * N for _ in range(1 << N)]
    for v in range(N): dp[1 << v][v] = 1
    for mask in range(1, 1 << N):
        row = dp[mask]
        for l in range(N):
            c = row[l]
            if not c: continue
            for w in range(N):
                if A[l][w] and not (mask >> w) & 1: dp[mask | (1 << w)][w] += c
    return sum(dp[(1 << N) - 1])

def rand_tour(n, rng):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def ordered_path_partitions(C):
    """dict j -> number of ordered partitions of V(C) into j nonempty paths."""
    n = len(C)
    # paths on each subset
    hp = {}
    for S in range(1, 1 << n):
        vs = [i for i in range(n) if (S >> i) & 1]
        sub = [[C[a][b] for b in vs] for a in vs]
        hp[S] = ham_digraph(sub)
    full = (1 << n) - 1
    out = {}
    # ordered partitions: sequences of disjoint nonempty subsets covering full
    def rec(rem, j, prod):
        if rem == 0:
            out[j] = out.get(j, 0) + prod
            return
        S = rem
        while S:
            rec(rem ^ S, j + 1, prod * hp[S])
            S = (S - 1) & rem
    rec(full, 0, 1)
    return out

def walk_counts(A, u, v, J, K):
    """c[j][k] for j <= J, k <= K: walks visiting u j times, v k times, all others once."""
    n = len(A)
    others = [x for x in range(n) if x not in (u, v)]
    idx = {x: i for i, x in enumerate(others)}
    full = (1 << len(others)) - 1
    from collections import defaultdict
    cur = defaultdict(int)
    for x in range(n):
        if x == u: cur[(0, 1, 0, x)] += 1
        elif x == v: cur[(0, 0, 1, x)] += 1
        else: cur[(1 << idx[x], 0, 0, x)] += 1
    c = [[0] * (K + 1) for _ in range(J + 1)]
    L = len(others) + J + K
    for step in range(L):
        nxt = defaultdict(int)
        for (mask, a, b, last), cnt in cur.items():
            if mask == full and a >= 1 and b >= 1:
                c[a][b] += cnt
            for w in range(n):
                if not A[last][w]: continue
                if w == u:
                    if a < J: nxt[(mask, a + 1, b, w)] += cnt
                elif w == v:
                    if b < K: nxt[(mask, a, b + 1, w)] += cnt
                else:
                    if not (mask >> idx[w]) & 1: nxt[(mask | (1 << idx[w]), a, b, w)] += cnt
        cur = nxt
    return c

def substitute(A, u, v, C1, C2):
    n = len(A); n1 = len(C1); n2 = len(C2)
    blocks = []
    for x in range(n):
        if x == u: blocks.append([('u', i) for i in range(n1)])
        elif x == v: blocks.append([('v', i) for i in range(n2)])
        else: blocks.append([('o', x)])
    verts = [p for b in blocks for p in b]
    owner = {}
    for x, b in enumerate(blocks):
        for p in b: owner[p] = x
    N = len(verts)
    B = [[0] * N for _ in range(N)]
    for a in range(N):
        for b in range(N):
            if a == b: continue
            pa, pb = verts[a], verts[b]
            xa, xb = owner[pa], owner[pb]
            if xa != xb: B[a][b] = A[xa][xb]
            elif pa[0] == 'u': B[a][b] = C1[pa[1]][pb[1]]
            else: B[a][b] = C2[pa[1]][pb[1]]
    return B

def v1(trials=300, seed=94):
    rng = random.Random(seed)
    bad = 0; done = 0
    for _ in range(trials):
        n = rng.choice([3, 4, 5, 6])
        A = rand_tour(n, rng)
        u, v = rng.sample(range(n), 2)
        n1 = rng.choice([1, 2, 3]); n2 = rng.choice([1, 2, 3])
        C1 = rand_tour(n1, rng); C2 = rand_tour(n2, rng)
        c = walk_counts(A, u, v, n1, n2)
        q1 = ordered_path_partitions(C1); q2 = ordered_path_partitions(C2)
        pred = sum(c[j][k] * q1.get(j, 0) * q2.get(k, 0) for j in range(1, n1 + 1) for k in range(1, n2 + 1))
        act = ham_digraph(substitute(A, u, v, C1, C2))
        done += 1
        if pred != act: bad += 1
    print(f"(V1) two-module law: {done} random instances (context order 3..6, modules of order 1..3), failures: {bad}")

@njit(cache=True)
def ham_dup(n, m, pi, pj, u, v):
    """Hamiltonian paths of the digraph on n+2 vertices: tournament m with u, v duplicated (twins non-adjacent)."""
    N = n + 2
    adj = np.zeros(N, np.int64)
    na = n * (n - 1) // 2
    for k in range(na):
        if (m >> k) & 1:
            a = pi[k]; b = pj[k]
        else:
            a = pj[k]; b = pi[k]
        adj[a] |= 1 << b
    # twin of u is vertex n, twin of v is vertex n+1
    for a in range(n):
        if (adj[a] >> u) & 1: adj[a] |= 1 << n
        if (adj[a] >> v) & 1: adj[a] |= 1 << (n + 1)
    adj[n] = adj[u]
    adj[n + 1] = adj[v]
    # arcs between twin of u and v-copies etc. are inherited through adj[u], adj[v] (already include bits n, n+1)
    full = (1 << N) - 1
    dp = np.zeros((1 << N, N), np.int64)
    for x in range(N):
        dp[1 << x, x] = 1
    for mask in range(1, 1 << N):
        for last in range(N):
            c = dp[mask, last]
            if c == 0: continue
            nxt = adj[last] & ~mask
            for w in range(N):
                if (nxt >> w) & 1:
                    dp[mask | (1 << w), w] += c
    s = 0
    for x in range(N):
        s += dp[full, x]
    return s

@njit(cache=True)
def is_strong(n, m, pi, pj):
    sc = np.zeros(n, np.int64)
    na = n * (n - 1) // 2
    for k in range(na):
        if (m >> k) & 1: sc[pi[k]] += 1
        else: sc[pj[k]] += 1
    sc.sort()
    acc = 0
    for k in range(1, n):
        acc += sc[k - 1]
        if acc == k * (k - 1) // 2: return False
    return True

@njit(parallel=True, cache=True)
def v2_scan(n, pi, pj):
    na = n * (n - 1) // 2
    M = 1 << na
    mins = np.full(M, 1 << 40, np.int64)
    cnts = np.zeros(M, np.int64)
    for m in prange(M):
        if not is_strong(n, m, pi, pj): continue
        best = 1 << 40
        for u in range(n):
            for v in range(u + 1, n):
                h = ham_dup(n, m, pi, pj, u, v)
                c22 = h // 4
                if h % 4 != 0: c22 = -1
                if c22 < best: best = c22
                cnts[m] += 1
        mins[m] = best
    return mins, cnts

def v2():
    for n in (3, 4, 5, 6):
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        pi = np.array([p[0] for p in pairs], np.int64); pj = np.array([p[1] for p in pairs], np.int64)
        t0 = time.time()
        mins, cnts = v2_scan(n, pi, pj)
        sel = cnts > 0
        print(f"(V2) n = {n}: strong labelled tournaments {int(sel.sum())}, port pairs {int(cnts.sum())}, "
              f"min c_22 = {int(mins[sel].min())}, divisibility-by-4 failures: {int((mins[sel] < 0).sum())}, "
              f"tournaments attaining the minimum: {int((mins[sel] == mins[sel].min()).sum())}  [{time.time()-t0:.1f}s]")

if __name__ == '__main__':
    v1()
    v2()
