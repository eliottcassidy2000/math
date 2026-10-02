#!/usr/bin/env python3
"""procgen_tcpc_20261001_lib.py -- TCPC lane library (collatz-procgen-20260922, 2026-10-01).

Clock digraphs Clk(G, C): vertices = a finite abelian group G (Z/m by default),
arc multiset {x -> x + c : c in C} (C a multiset).  Loops (0 in C), doubled arcs
(c and -c in C), missing pairs (neither) and tournaments (C disjoint union -C = G - 0)
are all allowed.

Two independent Hamiltonian-path engines:
  * hp_c(A): the C program procgen_tcpc_20261001_hp.c (128-bit exact / bitset parity);
  * hp_py(A): a pure-Python subset DP (written separately, used as the check);
  * hp_brute(A): permutation brute force (n <= 8).
"""
import itertools
import math
import os
import subprocess
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
SCRATCH = os.path.join(REPO, 'scratch', 'procgen_tcpc')
HP_SRC = os.path.join(HERE, 'procgen_tcpc_20261001_hp.c')
HP_BIN = os.path.join(SCRATCH, 'procgen_tcpc_20261001_hp')


def ensure_hp_binary():
    os.makedirs(SCRATCH, exist_ok=True)
    if (not os.path.exists(HP_BIN)) or os.path.getmtime(HP_BIN) < os.path.getmtime(HP_SRC):
        subprocess.run(['clang', '-O2', '-o', HP_BIN, HP_SRC], check=True)
    return HP_BIN


# ----------------------------------------------------------------------------
# clocks
# ----------------------------------------------------------------------------

def clock_matrix(m, C):
    """Adjacency multiplicity matrix of Clk(Z/m, C); C is an iterable (multiset)."""
    A = [[0] * m for _ in range(m)]
    for x in range(m):
        for c in C:
            A[x][(x + c) % m] += 1
    return A


def cayley_matrix(elems, add, C):
    """Clk(G, C) for a finite abelian group given by its element list and addition."""
    idx = {g: i for i, g in enumerate(elems)}
    n = len(elems)
    A = [[0] * n for _ in range(n)]
    for g in elems:
        for c in C:
            A[idx[g]][idx[add(g, c)]] += 1
    return A


def zm_product(m1, m2):
    elems = [(a, b) for a in range(m1) for b in range(m2)]
    add = lambda x, y: ((x[0] + y[0]) % m1, (x[1] + y[1]) % m2)
    return elems, add


def classify_pairs(m, C):
    """For Clk(Z/m, C) return dict with loops, and for each unordered pair class {c,-c}
    (c != 0) its type: 'missing', 'forward' (c in C, -c not), 'backward', 'doubled'
    (multiplicities summed)."""
    mult = {}
    for c in C:
        mult[c % m] = mult.get(c % m, 0) + 1
    out = {'loops': mult.get(0, 0)}
    seen = set()
    for c in range(1, m):
        if c in seen:
            continue
        d = (-c) % m
        seen.add(c); seen.add(d)
        a, b = mult.get(c, 0), mult.get(d, 0)
        if c == d:
            t = 'involution' if a else 'missing'
        elif a and b:
            t = 'doubled'
        elif a:
            t = 'forward'
        elif b:
            t = 'backward'
        else:
            t = 'missing'
        out[(c, d)] = (t, a, b)
    return out


def is_tournament(A):
    n = len(A)
    for i in range(n):
        if A[i][i]:
            return False
        for j in range(i + 1, n):
            if A[i][j] + A[j][i] != 1:
                return False
    return True


def lex_product(T, S):
    """Lexicographic product T[S] (S the same digraph at every vertex of T).
    Vertex (u, x) -> index u*|S| + x."""
    n, m = len(T), len(S)
    A = [[0] * (n * m) for _ in range(n * m)]
    for u in range(n):
        for x in range(m):
            for v in range(n):
                for y in range(m):
                    if u == v:
                        A[u * m + x][v * m + y] = S[x][y]
                    else:
                        A[u * m + x][v * m + y] = T[u][v]
    return A


def lex_product_multi(T, Ss):
    """Composition T[S_0, ..., S_{n-1}] with possibly different S_v."""
    n = len(T)
    sizes = [len(S) for S in Ss]
    off = [0]
    for s in sizes:
        off.append(off[-1] + s)
    N = off[-1]
    A = [[0] * N for _ in range(N)]
    for u in range(n):
        for x in range(sizes[u]):
            for v in range(n):
                for y in range(sizes[v]):
                    if u == v:
                        A[off[u] + x][off[v] + y] = Ss[u][x][y]
                    else:
                        A[off[u] + x][off[v] + y] = T[u][v]
    return A


def complement(A):
    n = len(A)
    return [[(0 if i == j else (1 - (1 if A[i][j] else 0))) for j in range(n)] for i in range(n)]


def transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


# ----------------------------------------------------------------------------
# Hamiltonian paths: C engine
# ----------------------------------------------------------------------------

def hp_c(A, mode='exact'):
    ensure_hp_binary()
    n = len(A)
    inp = '%d\n' % n + '\n'.join(' '.join(str(int(x)) for x in row) for row in A) + '\n'
    r = subprocess.run([HP_BIN, mode], input=inp, capture_output=True, text=True, check=True)
    res = {'start': [0] * n, 'end': [0] * n, 'arc': {}}
    for line in r.stdout.split('\n'):
        if not line:
            continue
        p = line.split()
        if p[0] == 'H':
            res['H'] = int(p[1])
        elif p[0] == 'HC':
            res['HC'] = int(p[1])
        elif p[0] == 'S':
            res['start'][int(p[1])] = int(p[2])
        elif p[0] == 'E':
            res['end'][int(p[1])] = int(p[2])
        elif p[0] == 'C':
            res['arc'][(int(p[1]), int(p[2]))] = int(p[3])
    return res


# ----------------------------------------------------------------------------
# Hamiltonian paths: independent pure-Python engine
# ----------------------------------------------------------------------------

def hp_py(A, want_arcs=True):
    """Pure-Python subset DP.  Returns dict(H, HC, start, end, arc)."""
    n = len(A)
    N = 1 << n
    # forward: f[S] = dict v -> count of paths on S ending at v
    f = [None] * N
    g = [None] * N  # backward: paths on S starting at v
    for v in range(n):
        f[1 << v] = {v: 1}
        g[1 << v] = {v: 1}
    order = sorted(range(1, N), key=lambda s: bin(s).count('1'))
    for S in order:
        fs = f[S]
        gs = g[S]
        if fs is None and gs is None:
            continue
        for w in range(n):
            if S >> w & 1:
                continue
            T = S | (1 << w)
            if fs:
                tot = 0
                for v, c in fs.items():
                    if A[v][w]:
                        tot += c * A[v][w]
                if tot:
                    if f[T] is None:
                        f[T] = {}
                    f[T][w] = f[T].get(w, 0) + tot
            if gs:
                tot = 0
                for v, c in gs.items():
                    if A[w][v]:
                        tot += c * A[w][v]
                if tot:
                    if g[T] is None:
                        g[T] = {}
                    g[T][w] = g[T].get(w, 0) + tot
    FULL = N - 1
    fF = f[FULL] or {}
    gF = g[FULL] or {}
    res = {'H': sum(fF.values()),
           'end': [fF.get(v, 0) for v in range(n)],
           'start': [gF.get(v, 0) for v in range(n)]}
    # Hamiltonian cycles through vertex 0 (paths from 0)
    if n >= 2:
        h = {1: {0: 1}}
        for S in order:
            if not (S & 1) or S not in h:
                continue
            for v, c in h[S].items():
                for w in range(n):
                    if S >> w & 1 or not A[v][w]:
                        continue
                    T = S | (1 << w)
                    h.setdefault(T, {})
                    h[T][w] = h[T].get(w, 0) + c * A[v][w]
        hc = 0
        for v, c in h.get(FULL, {}).items():
            if v != 0 and A[v][0]:
                hc += c * A[v][0]
        res['HC'] = hc
    else:
        res['HC'] = 0
    if want_arcs:
        arc = {}
        for S in range(1, FULL):
            fs = f[S]
            if not fs:
                continue
            R = FULL ^ S
            gs = g[R]
            if not gs:
                continue
            for u, cu in fs.items():
                for v, cv in gs.items():
                    if A[u][v]:
                        arc[(u, v)] = arc.get((u, v), 0) + cu * A[u][v] * cv
        for u in range(n):
            for v in range(n):
                if u != v and A[u][v]:
                    arc.setdefault((u, v), 0)
        res['arc'] = arc
    return res


def hp_brute(A):
    n = len(A)
    H = 0
    for p in itertools.permutations(range(n)):
        w = 1
        for i in range(n - 1):
            w *= A[p[i]][p[i + 1]]
            if not w:
                break
        H += w
    return H


def path_cover_numbers(A):
    """pc[k] = number of (unordered) covers of V by k vertex-disjoint directed paths
    (a single vertex is a path; arc multiplicities weight the count)."""
    n = len(A)
    N = 1 << n
    # P[S] = number of directed paths with vertex set exactly S
    f = [dict() for _ in range(N)]
    for v in range(n):
        f[1 << v][v] = 1
    for S in sorted(range(1, N), key=lambda s: bin(s).count('1')):
        for v, c in list(f[S].items()):
            for w in range(n):
                if S >> w & 1 or not A[v][w]:
                    continue
                T = S | (1 << w)
                f[T][w] = f[T].get(w, 0) + c * A[v][w]
    P = [sum(f[S].values()) for S in range(N)]
    # Q[S][k] = number of covers of S by k paths: block containing lowest vertex
    Q = [None] * N
    Q[0] = {0: 1}
    for S in range(1, N):
        low = S & (-S)
        rest = S ^ low
        acc = {}
        sub = rest
        while True:
            blk = sub | low
            pb = P[blk]
            if pb:
                for k, c in Q[S ^ blk].items():
                    acc[k + 1] = acc.get(k + 1, 0) + pb * c
            if sub == 0:
                break
            sub = (sub - 1) & rest
        Q[S] = acc
    return Q[N - 1]


# ----------------------------------------------------------------------------
# arithmetic helpers
# ----------------------------------------------------------------------------

def is_prime(p):
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    d = 3
    while d * d <= p:
        if p % d == 0:
            return False
        d += 2
    return True


def primitive_root(p):
    phi = p - 1
    fac = [q for q in range(2, phi + 1) if phi % q == 0 and is_prime(q)]
    for g in range(2, p):
        if all(pow(g, phi // q, p) != 1 for q in fac):
            return g
    raise ValueError


def qr_set(p):
    return sorted({(x * x) % p for x in range(1, p)})


def dlog_table(g, m):
    """discrete log base g on the units mod m (g a generator)."""
    t = {}
    x = 1
    k = 0
    while True:
        if x in t:
            break
        t[x] = k
        x = (x * g) % m
        k += 1
    return t


def v2(n):
    n = abs(n)
    k = 0
    while n % 2 == 0:
        n //= 2
        k += 1
    return k


def syracuse(a):
    """Syracuse map on odd a: (3a+1)/2^v, returns (S(a), v)."""
    t = 3 * a + 1
    v = v2(t)
    return t >> v, v


# ----------------------------------------------------------------------------
# nauty automorphism group order for a simple digraph (0/1 matrix, loops ignored)
# ----------------------------------------------------------------------------

def aut_order_digraph(A):
    n = len(A)
    lines = ['n=%d $=0 d g' % n]
    for i in range(n):
        nb = [str(j) for j in range(n) if j != i and A[i][j]]
        lines.append(' '.join(nb) + (';' if i < n - 1 else '.'))
    lines.append('x')
    lines.append('q')
    r = subprocess.run(['dreadnaut'], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    out = r.stdout
    import re
    m = re.search(r'grpsize=([0-9.e+]+)', out)
    if not m:
        raise RuntimeError('dreadnaut failed: ' + out[:500])
    s = m.group(1)
    if 'e' in s:
        mant, ex = s.split('e')
        return round(float(mant) * 10 ** int(ex))
    return int(round(float(s)))


def canon_digraph(A):
    """Canonical certificate of a simple digraph via nauty dreadnaut (c + b)."""
    n = len(A)
    lines = ['n=%d $=0 d g' % n]
    for i in range(n):
        nb = [str(j) for j in range(n) if j != i and A[i][j]]
        lines.append(' '.join(nb) + (';' if i < n - 1 else '.'))
    lines += ['c', 'x', 'b', 'q']
    r = subprocess.run(['dreadnaut'], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    out = r.stdout
    # the canonically labelled graph is printed after the line containing 'canupdates' or the grpsize line;
    # take everything after the last line that starts with a digit-colon pattern of the labelling
    import re
    m = re.search(r'grpsize=([0-9.e+]+)', out)
    grp = m.group(1) if m else '?'
    # canonical adjacency lines look like '  0 : 3 5 7;'
    adj = re.findall(r'^\s*(\d+)\s*:\s*([0-9 ]*);', out, re.M)
    cert = tuple((int(a), tuple(int(x) for x in b.split())) for a, b in adj[-n:])
    return cert, grp


def twisted_clock(N, C0, C1):
    """Two-sheet clock on Z/N (N even): x -> x + d iff d in C_{x mod 2}."""
    A = [[0] * N for _ in range(N)]
    for x in range(N):
        for d in (C0 if x % 2 == 0 else C1):
            A[x][(x + d) % N] += 1
    return A
