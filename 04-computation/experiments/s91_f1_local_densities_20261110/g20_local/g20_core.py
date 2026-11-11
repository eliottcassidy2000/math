#!/usr/bin/env python3
"""S91 helper (g20_local): direct evaluator of the two-apex grammar G(2,0).

Tournament T(w) for a word w over {0,1,2,3} of length n: chain 1 -> 2 -> ... -> n (i -> j for i < j), apexes a, b with
a -> b, a -> i iff bit 0 of w_i, b -> i iff bit 1 of w_i (otherwise i -> apex).  Same letters as Session 194.

Evaluator (own derivation; it is the path bijection of Session 191/194): a chain-only stretch of a Hamiltonian path is
an increasing sequence, so a Hamiltonian path is  S0, p, S1, q, S2  with {p, q} = {a, b} and (S0, S1, S2) an ordered
partition of the chain into three increasing runs (any of them empty), subject to
    last(S0) -> p,   p -> first(S1),   last(S1) -> q,   q -> first(S2),   and if S1 is empty: p -> q, i.e. (p, q) = (a, b).
Scanning the chain left to right, each vertex chooses its run; the state is
    (order o, s0 in {empty, last good, last bad}, s1 in {empty, last good, last bad}, s2 in {empty, nonempty}):
2 * 3 * 3 * 2 = 36 states.  "First" conditions are checked when a run is born, "last" conditions at the end.
"""
import numpy as np
from numba import njit, prange

NS = 36


def sidx(o, s0, s1, s2):
    return 18 * o + 6 * s0 + 2 * s1 + s2


def build_tables():
    t_to = np.zeros((4, NS, 3), dtype=np.int64)
    t_cnt = np.zeros((4, NS), dtype=np.int64)
    acc = np.zeros(NS, dtype=np.int64)
    for o in range(2):
        for s0 in range(3):
            for s1 in range(3):
                for s2 in range(2):
                    s = sidx(o, s0, s1, s2)
                    if s0 != 2 and (s1 == 1 or (s1 == 0 and o == 0)):
                        acc[s] = 1
                    for c in range(4):
                        abit, bbit = c & 1, (c >> 1) & 1
                        pbit, qbit = (abit, bbit) if o == 0 else (bbit, abit)
                        tg = []
                        # the vertex joins S0: always possible; it is good iff it beats p (pbit = 0)
                        tg.append(sidx(o, 2 if pbit else 1, s1, s2))
                        # joins S1: if S1 is empty it is its first vertex and needs p -> vertex
                        if s1 != 0 or pbit:
                            tg.append(sidx(o, s0, 2 if qbit else 1, s2))
                        # joins S2: if S2 is empty it needs q -> vertex
                        if s2 != 0 or qbit:
                            tg.append(sidx(o, s0, s1, 1))
                        t_cnt[c, s] = len(tg)
                        for k, t in enumerate(tg):
                            t_to[c, s, k] = t
    init = np.zeros(NS, dtype=np.int64)
    init[sidx(0, 0, 0, 0)] = 1
    init[sidx(1, 0, 0, 0)] = 1
    return t_to, t_cnt, acc, init


T_TO, T_CNT, ACC, INIT = build_tables()


def letter_matrices():
    """M[c][s][t] = number of transitions s -> t under letter c (row vector convention: v' = v M_c)."""
    M = np.zeros((4, NS, NS), dtype=np.int64)
    for c in range(4):
        for s in range(NS):
            for k in range(T_CNT[c, s]):
                M[c, s, T_TO[c, s, k]] += 1
    return M


def H_word(word):
    """Pure-python exact evaluation (Python integers)."""
    v = [int(x) for x in INIT]
    for c in word:
        out = [0] * NS
        for s in range(NS):
            x = v[s]
            if x:
                for k in range(T_CNT[c, s]):
                    out[T_TO[c, s, k]] += x
        v = out
    return sum(v[s] for s in range(NS) if ACC[s])


def state_word(word):
    v = [int(x) for x in INIT]
    for c in word:
        out = [0] * NS
        for s in range(NS):
            x = v[s]
            if x:
                for k in range(T_CNT[c, s]):
                    out[T_TO[c, s, k]] += x
        v = out
    return v


@njit(cache=True)
def nb_step(v, c, out, t_to, t_cnt):
    for j in range(36):
        out[j] = 0
    for s in range(36):
        x = v[s]
        if x != 0:
            for k in range(t_cnt[c, s]):
                out[t_to[c, s, k]] += x


@njit(cache=True)
def nb_H(word, t_to, t_cnt, acc, init):
    v = init.copy()
    out = np.zeros(36, dtype=np.int64)
    for i in range(word.shape[0]):
        nb_step(v, word[i], out, t_to, t_cnt)
        for j in range(36):
            v[j] = out[j]
    tot = 0
    for j in range(36):
        if acc[j]:
            tot += v[j]
    return tot


@njit(cache=True)
def nb_bruteforce(word):
    """Hamiltonian paths of T(word) by the vertex-subset dynamic programme (vertices: 0 = a, 1 = b, 2.. = chain)."""
    n = word.shape[0] + 2
    adj = np.zeros((n, n), dtype=np.int64)
    adj[0, 1] = 1
    for i in range(2, n):
        c = word[i - 2]
        if c & 1:
            adj[0, i] = 1
        else:
            adj[i, 0] = 1
        if (c >> 1) & 1:
            adj[1, i] = 1
        else:
            adj[i, 1] = 1
        for j in range(i + 1, n):
            adj[i, j] = 1
    full = (1 << n) - 1
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for i in range(n):
        dp[1 << i, i] = 1
    for mask in range(1, full + 1):
        for last in range(n):
            x = dp[mask, last]
            if x == 0:
                continue
            for nx in range(n):
                if (mask >> nx) & 1 == 0 and adj[last, nx] == 1:
                    dp[mask | (1 << nx), nx] += x
    tot = 0
    for i in range(n):
        tot += dp[full, i]
    return tot


@njit(cache=True)
def _dfs_hist(prefix, n, mod, t_to, t_cnt, acc, init, hist):
    """adds to hist the counts of H mod `mod` over all words of length n beginning with `prefix`; returns max H."""
    k0 = prefix.shape[0]
    st = np.zeros((n + 1, 36), dtype=np.int64)
    for j in range(36):
        st[0, j] = init[j]
    out = np.zeros(36, dtype=np.int64)
    for i in range(k0):
        nb_step(st[i], prefix[i], st[i + 1], t_to, t_cnt)
    mx = 0
    if k0 == n:
        tot = 0
        for j in range(36):
            if acc[j]:
                tot += st[n, j]
        hist[tot % mod] += 1
        return tot
    letter = np.zeros(n + 1, dtype=np.int64)
    d = k0
    letter[d] = 0
    while True:
        if letter[d] == 4:
            d -= 1
            if d < k0:
                break
            letter[d] += 1
            continue
        nb_step(st[d], letter[d], st[d + 1], t_to, t_cnt)
        if d + 1 == n:
            tot = 0
            for j in range(36):
                if acc[j]:
                    tot += st[n, j]
            hist[tot % mod] += 1
            if tot > mx:
                mx = tot
            letter[d] += 1
        else:
            d += 1
            letter[d] = 0
    return mx


@njit(parallel=True, cache=True)
def nb_full_hist(n, mod, k0, t_to, t_cnt, acc, init):
    """histogram of H mod `mod` over all 4^n words (parallel over the 4^k0 prefixes)."""
    npre = 4 ** k0
    hists = np.zeros((npre, mod), dtype=np.int64)
    mxs = np.zeros(npre, dtype=np.int64)
    for pi in prange(npre):
        prefix = np.zeros(k0, dtype=np.int64)
        for i in range(k0):
            prefix[k0 - 1 - i] = (pi // (4 ** i)) % 4
        mxs[pi] = _dfs_hist(prefix, n, mod, t_to, t_cnt, acc, init, hists[pi])
    tot = np.zeros(mod, dtype=np.int64)
    for pi in range(npre):
        for j in range(mod):
            tot[j] += hists[pi, j]
    return tot, mxs.max()


@njit(cache=True)
def _dfs_hist_reduced(prefix, n, mod, t_to, t_cnt, acc, init, hist):
    """as _dfs_hist, but only words whose LAST letter is not 3 (the caller restricts the first letter)."""
    k0 = prefix.shape[0]
    st = np.zeros((n + 1, 36), dtype=np.int64)
    for j in range(36):
        st[0, j] = init[j]
    for i in range(k0):
        nb_step(st[i], prefix[i], st[i + 1], t_to, t_cnt)
    mx = 0
    letter = np.zeros(n + 1, dtype=np.int64)
    d = k0
    letter[d] = 0
    while True:
        lim = 3 if d == n - 1 else 4
        if letter[d] == lim:
            d -= 1
            if d < k0:
                break
            letter[d] += 1
            continue
        nb_step(st[d], letter[d], st[d + 1], t_to, t_cnt)
        if d + 1 == n:
            tot = 0
            for j in range(36):
                if acc[j]:
                    tot += st[n, j]
            hist[tot % mod] += 1
            if tot > mx:
                mx = tot
            letter[d] += 1
        else:
            d += 1
            letter[d] = 0
    return mx


@njit(parallel=True, cache=True)
def nb_reduced_hist(n, mod, t_to, t_cnt, acc, init):
    """histogram of H mod `mod` over the 9 * 4^(n-2) words of length n >= 3 with first letter != 0 and last letter
    != 3 (parallel over the 12 two-letter prefixes)."""
    hists = np.zeros((12, mod), dtype=np.int64)
    mxs = np.zeros(12, dtype=np.int64)
    for pi in prange(12):
        prefix = np.zeros(2, dtype=np.int64)
        prefix[0] = 1 + pi // 4
        prefix[1] = pi % 4
        mxs[pi] = _dfs_hist_reduced(prefix, n, mod, t_to, t_cnt, acc, init, hists[pi])
    tot = np.zeros(mod, dtype=np.int64)
    for pi in range(12):
        for j in range(mod):
            tot[j] += hists[pi, j]
    return tot, mxs.max()


def psi(h):
    r = 0
    stack = []
    while h:
        stack.append(h)
        h >>= 1
    r = 0
    for v in reversed(stack):
        half = v >> 1
        if v & 1:
            r = r + (half >> 1) + 1
        else:
            r = r + half
    return r


def x_of(h):
    return 2 * psi(h) + 1


def y_of(h):
    return 2 * h + 1
