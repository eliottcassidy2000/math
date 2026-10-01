"""procgen_selfie_20261001_lib.py -- helpers for the selfie-tournament lane (2026-10-01).

Conventions (match procgen_selfie_20261001_arcs.c):
  vertices 0..n-1; pairs (i,j), i<j, in lexicographic order; bit b of a tournament index x
  is the orientation of pair b: 1 means i->j, 0 means j->i.
  Tiling model (THM-474, THM-022): base path (n-1) -> (n-2) -> ... -> 0, i.e. every pair
  (i,i+1) has bit 0; tiles are the pairs (i,j) with j >= i+2; tile bit 0 = 'forward'
  (larger -> smaller), tile bit 1 = 'backward'.
"""
import itertools
import math
import os
import subprocess

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
SCRATCH = os.path.join(REPO, 'scratch', 'procgen_selfie')
CSRC = os.path.join(HERE, 'procgen_selfie_20261001_arcs.c')
CSRC_EULER = os.path.join(HERE, 'procgen_selfie_20261001_euler.c')
BIN = os.path.join(SCRATCH, 'procgen_selfie_arcs')
BIN_EULER = os.path.join(SCRATCH, 'procgen_selfie_euler')
CSRC_PALEY = os.path.join(HERE, 'procgen_selfie_20261001_paley2.c')
BIN_PALEY = os.path.join(SCRATCH, 'procgen_selfie_paley2')


def build():
    os.makedirs(SCRATCH, exist_ok=True)
    subprocess.run(['cc', '-O2', '-o', BIN, CSRC], check=True)
    subprocess.run(['cc', '-O2', '-o', BIN_EULER, CSRC_EULER], check=True)
    subprocess.run(['cc', '-O2', '-o', BIN_PALEY, CSRC_PALEY], check=True)


def run_c(args, feed=None, binary=None):
    """Run the C engine (niced); feed = argv list of a generator command piped to stdin."""
    binary = binary or BIN
    if feed is None:
        r = subprocess.run(['nice', binary] + list(args), capture_output=True, text=True, check=True)
        return r.stdout
    p1 = subprocess.Popen(feed, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    r = subprocess.run(['nice', binary] + list(args), stdin=p1.stdout, capture_output=True, text=True)
    p1.stdout.close()
    p1.wait()
    return r.stdout


def gentourng(n, extra=()):
    return ['gentourng', '-q'] + list(extra) + [str(n)]


def classes(n, extra=()):
    out = subprocess.run(gentourng(n, extra), capture_output=True, text=True, check=True).stdout
    return [s for s in out.split() if s and s[0] in '01']


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def adj_from_string(s, n):
    A = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1':
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return A


def adj_from_index(x, n):
    A = [[0] * n for _ in range(n)]
    for b, (i, j) in enumerate(pairs(n)):
        if (x >> b) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1
    return A


def string_from_adj(A):
    n = len(A)
    return ''.join('1' if A[i][j] else '0' for i in range(n) for j in range(i + 1, n))


def hp_count(A, verts=None):
    """Number of directed Hamiltonian paths of the subtournament on verts (H(empty)=1)."""
    if verts is None:
        verts = list(range(len(A)))
    verts = list(verts)
    m = len(verts)
    if m <= 1:
        return 1
    f = [[0] * m for _ in range(1 << m)]
    for i in range(m):
        f[1 << i][i] = 1
    for S in range(1, 1 << m):
        row = f[S]
        for i in range(m):
            if row[i]:
                vi = verts[i]
                for j in range(m):
                    if not (S >> j) & 1 and A[vi][verts[j]]:
                        f[S | (1 << j)][j] += row[i]
    return sum(f[(1 << m) - 1])


def hp_list(A):
    """All directed Hamiltonian paths as tuples (small n only)."""
    n = len(A)
    out = []

    def rec(path, used):
        if len(path) == n:
            out.append(tuple(path))
            return
        last = path[-1]
        for w in range(n):
            if not (used >> w) & 1 and A[last][w]:
                path.append(w)
                rec(path, used | (1 << w))
                path.pop()

    for s in range(n):
        rec([s], 1 << s)
    return out


def odd_cycles(A):
    """Directed odd cycles (length >= 3) as (frozenset, canonical tuple)."""
    n = len(A)
    cyc = []
    for k in range(3, n + 1, 2):
        for sub in itertools.combinations(range(n), k):
            s0 = sub[0]
            for perm in itertools.permutations(sub[1:]):
                seq = (s0,) + perm
                if all(A[seq[i]][seq[(i + 1) % k]] for i in range(k)):
                    cyc.append((frozenset(sub), seq))
    return cyc


def indep_sets_by_cover(cycles):
    """Dict: frozenset of covered vertices -> list of sizes |S| over independent sets S."""
    res = {}

    def rec(start, used, size):
        res.setdefault(used, []).append(size)
        for idx in range(start, len(cycles)):
            vs = cycles[idx][0]
            if not (vs & used):
                rec(idx + 1, used | vs, size + 1)

    rec(0, frozenset(), 0)
    return res


def poly_mul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a:
            for j, b in enumerate(q):
                r[i + j] += a * b
    return r


def poly_add(p, q):
    r = [0] * max(len(p), len(q))
    for i, a in enumerate(p):
        r[i] += a
    for i, a in enumerate(q):
        r[i] += a
    return r


def poly_pow(p, k):
    r = [1]
    for _ in range(k):
        r = poly_mul(r, p)
    return r


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def wht(v):
    """Unnormalised Walsh-Hadamard transform along the last axis (length 2^m)."""
    v = np.array(v, dtype=np.int64, copy=True)
    h = 1
    L = v.shape[-1]
    while h < L:
        v = v.reshape(v.shape[:-1] + (L // (2 * h), 2, h))
        a = v[..., 0, :].copy()
        b = v[..., 1, :].copy()
        v[..., 0, :] = a + b
        v[..., 1, :] = a - b
        v = v.reshape(v.shape[:-3] + (L,))
        h *= 2
    return v


def arc_counts(A):
    """(H, {arc: #HPs through arc}, start[v], end[v]) by bitmask DP (small n)."""
    n = len(A)
    N = 1 << n
    f = [[0] * n for _ in range(N)]
    g = [[0] * n for _ in range(N)]
    for v in range(n):
        f[1 << v][v] = 1
        g[1 << v][v] = 1
    for S in range(1, N):
        for v in range(n):
            if f[S][v]:
                for w in range(n):
                    if not (S >> w) & 1 and A[v][w]:
                        f[S | 1 << w][w] += f[S][v]
            if g[S][v]:
                for w in range(n):
                    if not (S >> w) & 1 and A[w][v]:
                        g[S | 1 << w][w] += g[S][v]
    full = N - 1
    H = sum(f[full])
    c = {}
    for u in range(n):
        for v in range(n):
            if A[u][v]:
                c[(u, v)] = sum(f[S][u] * g[full ^ S][v] for S in range(N)
                                if (S >> u) & 1 and not (S >> v) & 1 and f[S][u])
    return H, c, [g[full][v] for v in range(n)], [f[full][v] for v in range(n)]


def strong_components_ordered(A):
    """Strong components in the order C_1 => C_2 => ... (all arcs go forward)."""
    n = len(A)
    reach = []
    for s in range(n):
        seen = {s}
        stack = [s]
        while stack:
            v = stack.pop()
            for w in range(n):
                if A[v][w] and w not in seen:
                    seen.add(w)
                    stack.append(w)
        reach.append(seen)
    comps = []
    done = set()
    for s in range(n):
        if s in done:
            continue
        comp = frozenset(t for t in range(n) if t in reach[s] and s in reach[t])
        done |= comp
        comps.append(comp)
    # order: C_i before C_j iff arcs go C_i -> C_j; sort by size of reach set (descending)
    comps.sort(key=lambda C: -len(reach[next(iter(C))]))
    return comps


def sub_adj(A, verts):
    verts = list(verts)
    return [[A[a][b] for b in verts] for a in verts]


def paley_string(p, delete=None):
    Q = {(x * x) % p for x in range(1, p)}
    V = [v for v in range(p) if v != delete]
    return ''.join('1' if ((V[b] - V[a]) % p) in Q else '0'
                   for a in range(len(V)) for b in range(a + 1, len(V)))


def circulant_strings(n):
    """All circulant tournaments on Z_n (n odd): connection sets S with S + (-S) = Z_n minus 0."""
    reps = list(range(1, (n - 1) // 2 + 1))
    out = []
    for choice in itertools.product([0, 1], repeat=len(reps)):
        S = {r if c == 0 else n - r for r, c in zip(reps, choice)}
        out.append(''.join('1' if ((b - a) % n) in S else '0' for a in range(n) for b in range(a + 1, n)))
    return out


def cayley_z3z3_strings():
    """All Cayley tournaments on Z_3 x Z_3."""
    G = [(a, b) for a in range(3) for b in range(3)]
    nz = [g for g in G if g != (0, 0)]
    reps = []
    seen = set()
    for g in nz:
        if g in seen:
            continue
        neg = ((-g[0]) % 3, (-g[1]) % 3)
        reps.append((g, neg))
        seen |= {g, neg}
    out = []
    for choice in itertools.product([0, 1], repeat=len(reps)):
        S = {pr[c] for pr, c in zip(reps, choice)}
        s = ''
        for i in range(9):
            for j in range(i + 1, 9):
                d = ((G[j][0] - G[i][0]) % 3, (G[j][1] - G[i][1]) % 3)
                s += '1' if d in S else '0'
        out.append(s)
    return out


def drt15_doubled():
    """Doubly regular tournament on 15 vertices from the skew-Hadamard doubling of QR7."""
    p = 7
    Q = {(x * x) % p for x in range(1, p)}
    S = np.zeros((p, p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i != j:
                S[i][j] = 1 if (j - i) % p in Q else -1
    S8 = np.zeros((8, 8), dtype=int)
    S8[0, 1:] = 1
    S8[1:, 0] = -1
    S8[1:, 1:] = S
    K = np.eye(8, dtype=int) + S8
    assert (K @ K.T == 8 * np.eye(8, dtype=int)).all()
    M = np.block([[K, K], [-K.T, K.T]])
    assert (M @ M.T == 16 * np.eye(16, dtype=int)).all()
    assert ((M + M.T) == 2 * np.eye(16, dtype=int)).all()
    for i in range(1, 16):
        if M[0, i] == -1:
            M[i, :] *= -1
            M[:, i] *= -1
    T = (M[1:, 1:] == 1).astype(int)
    np.fill_diagonal(T, 0)
    co = T @ T.T
    assert all(co[i][j] == 3 for i in range(15) for j in range(15) if i != j)
    return T.tolist()


def all_sub_H(A):
    """H(T[S]) for every vertex subset S (bitmask index), H(empty) = 1, by one bitmask DP."""
    n = len(A)
    N = 1 << n
    f = [[0] * n for _ in range(N)]
    for v in range(n):
        f[1 << v][v] = 1
    out = [0] * N
    out[0] = 1
    for S in range(1, N):
        row = f[S]
        tot = 0
        for v in range(n):
            x = row[v]
            if x:
                tot += x
                for w in range(n):
                    if not (S >> w) & 1 and A[v][w]:
                        f[S | (1 << w)][w] += x
        out[S] = tot
    return out
