#!/usr/bin/env python3
"""auditC_lib.py -- helpers for audit C of THM-4876 / HYP-9469 (written independently of the session's scripts).

Tournaments on 0..n-1 are lists adj[u] = bitmask of out-neighbours of u.
"""
import subprocess, itertools
from functools import lru_cache


def parse_gentourng_line(n, line):
    """gentourng default output: upper triangle row by row, '1' at (i, j), i < j, means i -> j."""
    adj = [0] * n
    p = 0
    for i in range(n):
        for j in range(i + 1, n):
            if line[p] == '1':
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            p += 1
    return adj


def tournaments(n):
    out = subprocess.run(["gentourng", "-q", str(n)], capture_output=True, text=True).stdout.split()
    need = n * (n - 1) // 2
    return [parse_gentourng_line(n, line) for line in out if len(line) == need]


def labeled_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        adj = [0] * n
        for (i, j), b in zip(pairs, bits):
            if b:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
        yield adj


def reverse(n, adj):
    r = [0] * n
    for u in range(n):
        for w in range(n):
            if (adj[u] >> w) & 1:
                r[w] |= 1 << u
    return r


def dp_end(n, adj):
    """E[S][u] = number of Hamiltonian paths of T[S] ending at u (0 if u not in S)."""
    N = 1 << n
    E = [[0] * n for _ in range(N)]
    inn = reverse(n, adj)
    for S in range(1, N):
        for u in range(n):
            if not (S >> u) & 1:
                continue
            R = S & ~(1 << u)
            if R == 0:
                E[S][u] = 1
                continue
            acc = 0
            m = R & inn[u]
            while m:
                w = (m & -m).bit_length() - 1
                m &= m - 1
                acc += E[R][w]
            E[S][u] = acc
    return E


def dp_start(n, adj):
    return dp_end(n, reverse(n, adj))


def hsub_table(n, E):
    return [sum(E[S]) if S else 1 for S in range(1 << n)]


def H(n, adj):
    return sum(dp_end(n, adj)[(1 << n) - 1])


def ham_paths_list(n, adj):
    full = (1 << n) - 1
    out = []

    def rec(u, S, p):
        if S == full:
            out.append(tuple(p))
            return
        m = adj[u] & ~S
        while m:
            w = (m & -m).bit_length() - 1
            m &= m - 1
            rec(w, S | (1 << w), p + [w])
    for u in range(n):
        rec(u, 1 << u, [u])
    return out


def odd_cycles_through(n, adj, v):
    """every directed odd cycle through v, once each, as (vertex mask, length): DFS from v."""
    out = []

    def rec(u, S, length):
        if length >= 3 and length % 2 == 1 and (adj[u] >> v) & 1:
            out.append((S, length))
        m = adj[u] & ~S & ~(1 << v)
        while m:
            w = (m & -m).bit_length() - 1
            m &= m - 1
            rec(w, S | (1 << w), length + 1)
    rec(v, 1 << v, 1)
    return out


def all_odd_cycles(n, adj):
    """every directed odd cycle once, as (vertex mask, length): start at the minimum vertex, use larger vertices."""
    out = []
    for s in range(n):
        def rec(u, S, length):
            if length >= 3 and length % 2 == 1 and (adj[u] >> s) & 1:
                out.append((S, length))
            m = adj[u] & ~S
            m &= ~((1 << (s + 1)) - 1)   # only vertices > s
            while m:
                w = (m & -m).bit_length() - 1
                m &= m - 1
                rec(w, S | (1 << w), length + 1)
        rec(s, 1 << s, 1)
    return out


def I_omega_2(cycles):
    """sum over families of pairwise vertex-disjoint cycles of 2^|family|."""
    masks = [c[0] for c in cycles]
    k = len(masks)

    @lru_cache(maxsize=None)
    def f(i, used):
        if i == k:
            return 1
        r = f(i + 1, used)
        if not (masks[i] & used):
            r += 2 * f(i + 1, used | masks[i])
        return r
    return f(0, 0)


def is_transitive(n, adj, S=None):
    if S is None:
        S = (1 << n) - 1
    verts = [u for u in range(n) if (S >> u) & 1]
    for a in verts:
        for b in verts:
            for c in verts:
                if a != b and b != c and a != c and (adj[a] >> b) & 1 and (adj[b] >> c) & 1 and (adj[c] >> a) & 1:
                    return False
    return True


def is_strong(n, adj):
    full = (1 << n) - 1

    def reach(a):
        r = 1
        fr = [0]
        while fr:
            u = fr.pop()
            m = a[u] & ~r
            while m:
                w = (m & -m).bit_length() - 1
                m &= m - 1
                r |= 1 << w
                fr.append(w)
        return r == full
    return reach(adj) and reach(reverse(n, adj))


def add_vertex(n, adj, A):
    """X + v with v = n, v -> A, V \\ A -> v."""
    adj2 = [adj[u] | (0 if (A >> u) & 1 else (1 << n)) for u in range(n)] + [A]
    return n + 1, adj2


def delete_vertex(n, adj, v):
    keep = [u for u in range(n) if u != v]
    idx = {u: i for i, u in enumerate(keep)}
    adj2 = [0] * (n - 1)
    for u in keep:
        for w in keep:
            if (adj[u] >> w) & 1:
                adj2[idx[u]] |= 1 << idx[w]
    return n - 1, adj2


def induced(n, adj, S):
    keep = [u for u in range(n) if (S >> u) & 1]
    idx = {u: i for i, u in enumerate(keep)}
    adj2 = [0] * len(keep)
    for u in keep:
        for w in keep:
            if (adj[u] >> w) & 1:
                adj2[idx[u]] |= 1 << idx[w]
    return len(keep), adj2
