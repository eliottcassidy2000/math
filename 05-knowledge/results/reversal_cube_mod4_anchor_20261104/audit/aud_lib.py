"""Referee library for the THM-4864 audit (written from scratch; no packet code imported).

Conventions (as in the packet): ports x->y fixed, a_v=[v->x], c_v=[v->y].
Tournaments are lists `out` of out-neighbour bitmasks.  Queries are pairs (amask, cmask).
"""
import itertools
import os
import random
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.join(HERE, 'aud_engine.exe')


# ----------------------------------------------------------------------------- engine
def run_jobs(jobs):
    """jobs: list of (mode, n, out, wmask, queries). Returns list of result lists."""
    if not jobs:
        return []
    lines = []
    for mode, n, out, wmask, qs in jobs:
        assert len(qs) >= 1
        parts = [mode, str(n)] + [str(o) for o in out]
        if mode == 'C':
            parts.append(str(wmask))
        parts.append(str(len(qs)))
        for a, c in qs:
            parts += [str(a), str(c)]
        lines.append(' '.join(parts))
    res = subprocess.run([ENGINE], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError('engine failed: ' + res.stderr)
    outl = res.stdout.rstrip('\n').split('\n')
    if len(outl) != len(jobs):
        raise RuntimeError('engine output length mismatch')
    return [[int(t) for t in l.split()] for l in outl]


def hamF(out, queries):
    return run_jobs([('H', len(out), out, None, queries)])[0]


def gcyc(out, w, queries):
    return run_jobs([('C', len(out), out, w, queries)])[0]


# ----------------------------------------------------------------------------- tournaments
def rand_tour(n, rng):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
    return out


def all_tours(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        out = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
        yield out


def converse(out):
    n = len(out)
    new = [0] * n
    for i in range(n):
        for j in range(n):
            if out[i] >> j & 1:
                new[j] |= 1 << i
    return new


def arc(out, i, j):
    return (out[i] >> j) & 1


def transitive_from_order(order, n_local):
    """order: list of local indices b_1..b_k; b_i -> b_j for i<j."""
    out = [0] * n_local
    for p in range(len(order)):
        for q in range(p + 1, len(order)):
            out[order[p]] |= 1 << order[q]
    return out


def is_tournament(out):
    n = len(out)
    for i in range(n):
        if out[i] >> i & 1:
            return False
        for j in range(i + 1, n):
            if ((out[i] >> j) & 1) + ((out[j] >> i) & 1) != 1:
                return False
    return True


def compose(Q, interiors):
    """W = Q[A_1..A_m]. interiors[i] = out-bitmasks of child i on local indices.
    Returns (out, blocks) with blocks[i] = list of global names of child i (consecutive)."""
    blocks, start = [], 0
    for A in interiors:
        blocks.append(list(range(start, start + len(A))))
        start += len(A)
    N = start
    out = [0] * N
    m = len(Q)
    for i in range(m):
        for loc, p in enumerate(blocks[i]):
            for loc2 in range(len(interiors[i])):
                if interiors[i][loc] >> loc2 & 1:
                    out[p] |= 1 << blocks[i][loc2]
            for j in range(m):
                if j != i and arc(Q, i, j):
                    for q in blocks[j]:
                        out[p] |= 1 << q
    assert is_tournament(out)
    return out, blocks


# ----------------------------------------------------------------------------- Hamilton data on small tournaments
def ham_start_end(out):
    """Return (S, E): S[v] = # Hamilton paths starting at v, E[v] = # ending at v (exact)."""
    n = len(out)
    full = (1 << n) - 1
    # end[S][v]: paths with vertex set S ending at v; start-specific counts via per-start DP
    S_cnt = [0] * n
    E_cnt = [0] * n
    for s in range(n):
        dp = {(1 << s, s): 1}
        for size in range(1, n):
            ndp = {}
            for (S, v), val in dp.items():
                nb = out[v] & ~S & full
                while nb:
                    w = (nb & -nb).bit_length() - 1
                    nb &= nb - 1
                    key = (S | (1 << w), w)
                    ndp[key] = ndp.get(key, 0) + val
            dp = ndp
        for (S, v), val in dp.items():
            S_cnt[s] += val
            E_cnt[v] += val
    return S_cnt, E_cnt


def all_paths(out, allowed=None):
    """All directed simple paths (tuples), singletons included, inside `allowed` (bitmask)."""
    n = len(out)
    if allowed is None:
        allowed = (1 << n) - 1
    res = []
    stack = [((v,), 1 << v) for v in range(n) if allowed >> v & 1]
    while stack:
        p, S = stack.pop()
        res.append(p)
        nb = out[p[-1]] & ~S & allowed
        while nb:
            w = (nb & -nb).bit_length() - 1
            nb &= nb - 1
            stack.append((p + (w,), S | (1 << w)))
    return res


def wt(path, w):
    return sum((w >> v) & 1 for v in path) & 1


# ----------------------------------------------------------------------------- Boolean coefficients
def subq(A, C):
    """all sub-queries of the monomial with a-support A, c-support C (bitmasks)."""
    res = []
    a_bits = [1 << i for i in range(32) if A >> i & 1]
    c_bits = [1 << i for i in range(32) if C >> i & 1]
    for ra in range(len(a_bits) + 1):
        for sa in itertools.combinations(a_bits, ra):
            for rc in range(len(c_bits) + 1):
                for sc in itertools.combinations(c_bits, rc):
                    res.append((sum(sa), sum(sc)))
    return res


def coef(f, A, C):
    """Moebius coefficient of the multilinear F2 polynomial of f (dict query->bit) at monomial (A,C)."""
    s = 0
    for q in subq(A, C):
        s ^= f[q] & 1
    return s


def Lfun_queries(n, v):
    qs = set()
    for t in range(n):
        if t == v:
            continue
        qs.update(subq((1 << t) | (1 << v), 1 << v))
    return qs


def Lfun(f, Z, v):
    """L_v = sum_{t!=v} [a_t c_v] f + sum_{v->t in Z} [a_t a_v c_v] f   (mod 2)."""
    n = len(Z)
    s = 0
    for t in range(n):
        if t == v:
            continue
        s ^= coef(f, 1 << t, 1 << v)
        if arc(Z, v, t):
            s ^= coef(f, (1 << t) | (1 << v), 1 << v)
    return s


def battery_queries(n, vs):
    qs = set()
    for v in vs:
        qs |= Lfun_queries(n, v)
    return qs


def expand_query(q, kernel_blocks):
    """kernel query (A,C) on kernel indices -> actual masks via kernel_blocks[i] = list of names."""
    A, C = q
    a = c = 0
    for i, names in enumerate(kernel_blocks):
        if A >> i & 1:
            for p in names:
                a |= 1 << p
        if C >> i & 1:
            for p in names:
                c |= 1 << p
    return a, c


def half_contrast(FW, FWp):
    d = FW - FWp
    if d % 2:
        raise RuntimeError('odd contrast: Redei violated?')
    return (d // 2) % 2
