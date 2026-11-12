#!/usr/bin/env python3
"""Lane "digits" (S93 Chan-Pak wave): common exact tools.

Tournament on n vertices 0..n-1, adj[i] = bitmask of the out-neighbours of i.
H = number of Hamiltonian paths.  All counts are exact int64 (orders <= 20 keep H far below 2^62).

ham_tables(n, adj)  -> (E, S): E[mask, v] = Hamiltonian paths of T[mask] ending at v, S[mask, v] = starting at v.
hcount(n, adj)      -> H(T).
pair_stats(n, adj)  -> h, start[v], end[v], Q[e, f], hsub[mask]
                       Q[e, f] = ordered pairs (P1, P2) of disjoint NONEMPTY paths covering T, P1 ends at e, P2 starts at f.
ext_values(n, adj)  -> array of H(T + v) over all 2^n out-sets O of the new vertex (v -> O, V - O -> v).
canon(n, adj)       -> canonical key (exact: individualisation-refinement, minimum adjacency string over all leaves).
"""
import numpy as np
from numba import njit


@njit(cache=True)
def ham_tables(n, adj):
    full = 1 << n
    E = np.zeros((full, n), np.int64)
    S = np.zeros((full, n), np.int64)
    inn = np.zeros(n, np.int64)
    for i in range(n):
        for j in range(n):
            if (adj[i] >> j) & 1:
                inn[j] |= 1 << i
    for v in range(n):
        E[1 << v, v] = 1
        S[1 << v, v] = 1
    for m in range(1, full):
        for v in range(n):
            c = E[m, v]
            if c != 0:
                out = adj[v] & ~m
                while out:
                    low = out & -out
                    u = 0
                    t = low
                    while t > 1:
                        t >>= 1
                        u += 1
                    E[m | low, u] += c
                    out ^= low
            c = S[m, v]
            if c != 0:
                out = inn[v] & ~m
                while out:
                    low = out & -out
                    u = 0
                    t = low
                    while t > 1:
                        t >>= 1
                        u += 1
                    S[m | low, u] += c
                    out ^= low
    return E, S


@njit(cache=True)
def hcount(n, adj):
    """H(T) by the subset recursion (end tables only)."""
    full = 1 << n
    E = np.zeros((full, n), np.int64)
    for v in range(n):
        E[1 << v, v] = 1
    for m in range(1, full):
        for v in range(n):
            c = E[m, v]
            if c != 0:
                out = adj[v] & ~m
                while out:
                    low = out & -out
                    u = 0
                    t = low
                    while t > 1:
                        t >>= 1
                        u += 1
                    E[m | low, u] += c
                    out ^= low
    tot = 0
    for v in range(n):
        tot += E[full - 1, v]
    return tot


@njit(cache=True)
def pair_stats(n, adj):
    full = 1 << n
    E, S = ham_tables(n, adj)
    hsub = np.zeros(full, np.int64)
    hsub[0] = 1
    for m in range(1, full):
        t = 0
        for v in range(n):
            t += E[m, v]
        hsub[m] = t
    start = S[full - 1].copy()
    end = E[full - 1].copy()
    Q = np.zeros((n, n), np.int64)
    for m in range(1, full - 1):
        comp = (full - 1) ^ m
        for e in range(n):
            ce = E[m, e]
            if ce == 0:
                continue
            for f in range(n):
                cf = S[comp, f]
                if cf != 0:
                    Q[e, f] += ce * cf
    return hsub[full - 1], start, end, Q, hsub


@njit(cache=True)
def ext_values(n, adj):
    """H(T + v) for every out-set O of the new vertex: Start(O) + h - End(O) + sum_{e notin O, f in O} Q[e, f]."""
    h, start, end, Q, hsub = pair_stats(n, adj)
    full = 1 << n
    out = np.zeros(full, np.int64)
    for O in range(full):
        t = h
        for f in range(n):
            if (O >> f) & 1:
                t += start[f] - end[f]
                for e in range(n):
                    if not (O >> e) & 1:
                        t += Q[e, f]
        out[O] = t
    return out


def adj_from_rows(rows):
    """rows: list of lists of out-neighbours."""
    n = len(rows)
    a = np.zeros(n, np.int64)
    for i, r in enumerate(rows):
        for j in r:
            a[i] |= 1 << j
    return a


def extend(n, adj, O):
    """T + v with v -> O (bitmask), V - O -> v; the new vertex has index n."""
    a = np.zeros(n + 1, np.int64)
    for i in range(n):
        a[i] = int(adj[i])
        if not (O >> i) & 1:
            a[i] |= 1 << n
    a[n] = O
    return a


def transitive(n):
    a = np.zeros(n, np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            a[i] |= 1 << j
    return a


def brute_h(n, adj):
    """Independent slow count (permutation search with pruning) for cross-checks."""
    cnt = 0
    adjl = [int(x) for x in adj]

    def rec(v, used, k):
        nonlocal cnt
        if k == n:
            cnt += 1
            return
        out = adjl[v] & ~used
        while out:
            low = out & -out
            u = low.bit_length() - 1
            rec(u, used | low, k + 1)
            out ^= low

    for s in range(n):
        rec(s, 1 << s, 1)
    return cnt


# ---------------------------------------------------------------- canonical form (exact)
def _refine(n, adjl, color):
    """Equitable refinement; colors are canonical ranks (isomorphism-invariant)."""
    ncls = len(set(color))
    while True:
        sig = []
        for v in range(n):
            outs = []
            a = adjl[v]
            while a:
                low = a & -a
                outs.append(color[low.bit_length() - 1])
                a ^= low
            outs.sort()
            sig.append((color[v], tuple(outs)))
        ranks = {s: i for i, s in enumerate(sorted(set(sig)))}
        color = [ranks[s] for s in sig]
        k = len(ranks)
        if k == ncls:
            return color
        ncls = k


def _leaf_string(n, adjl, color):
    order = sorted(range(n), key=lambda v: color[v])
    pos = [0] * n
    for i, v in enumerate(order):
        pos[v] = i
    key = 0
    for i in range(n):
        v = order[i]
        row = 0
        a = adjl[v]
        while a:
            low = a & -a
            row |= 1 << pos[low.bit_length() - 1]
            a ^= low
        key = (key << n) | row
    return key


def canon(n, adj):
    adjl = [int(x) for x in adj]
    color = [bin(a).count("1") for a in adjl]
    ranks = {s: i for i, s in enumerate(sorted(set(color)))}
    color = [ranks[c] for c in color]
    color = _refine(n, adjl, color)
    best = None

    def rec(color):
        nonlocal best
        cells = {}
        for v in range(n):
            cells.setdefault(color[v], []).append(v)
        target = None
        for c in sorted(cells):
            if len(cells[c]) > 1:
                target = cells[c]
                break
        if target is None:
            s = _leaf_string(n, adjl, color)
            if best is None or s < best:
                best = s
            return
        for v in target:
            c2 = [2 * x + (1 if (x == color[v] and u != v) else 0) for u, x in enumerate(color)]
            ranks = {s: i for i, s in enumerate(sorted(set(c2)))}
            c2 = [ranks[x] for x in c2]
            rec(_refine(n, adjl, c2))

    rec(color)
    return (n, best)


def adj_from_key(key):
    n, s = key
    a = np.zeros(n, np.int64)
    mask = (1 << n) - 1
    for i in range(n - 1, -1, -1):
        a[i] = s & mask
        s >>= n
    return a


def iso_classes(nmax):
    """Isomorphism classes of tournaments of orders 1..nmax by one-vertex extension and canonical form."""
    classes = {1: {canon(1, np.zeros(1, np.int64)): np.zeros(1, np.int64)}}
    for n in range(1, nmax):
        nxt = {}
        for key, adj in classes[n].items():
            for O in range(1 << n):
                a = extend(n, adj, O)
                k = canon(n + 1, a)
                if k not in nxt:
                    nxt[k] = adj_from_key(k)
        classes[n + 1] = nxt
    return classes


if __name__ == "__main__":
    import itertools
    import random
    # self-test: subset count against permutation search; ext_values against direct counts; class counts
    rng = random.Random(7)
    ok = True
    for n in range(1, 8):
        for _ in range(30):
            a = np.zeros(n, np.int64)
            for i in range(n):
                for j in range(i + 1, n):
                    if rng.random() < 0.5:
                        a[i] |= 1 << j
                    else:
                        a[j] |= 1 << i
            h1 = hcount(n, a)
            h2 = brute_h(n, a)
            if h1 != h2:
                ok = False
            ev = ext_values(n, a)
            for O in range(1 << n):
                if hcount(n + 1, extend(n, a, O)) != ev[O]:
                    ok = False
    print("subset count = permutation search, ext_values = direct count (orders 1..7, 30 random each):", ok)
    cl = iso_classes(7)
    print("isomorphism classes of orders 1..7:", [len(cl[n]) for n in range(1, 8)], "(expected 1, 1, 2, 4, 12, 56, 456)")
