"""procgen_sumgraph_20260926_theory.py -- finite checks for the reflection-orbit theory of sum graphs.

Session collatz-procgen-20260922, lane "sumgraph" (2026-09-26).
Pure Python, deterministic.  Uses only procgen_sumgraph_20260926_solver (own code).
"""
import sys
import os
from itertools import combinations
from math import gcd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_sumgraph_20260926_solver as SV   # noqa: E402


# ----------------------------------------------------------------------------------
# matchings and small unions (Task A)
# ----------------------------------------------------------------------------------
def msize(n, t):
    """|M_t| on [n] by Lemma R1."""
    if t < 3 or t > 2 * n - 1:
        return 0
    return (n - abs(t - n - 1)) // 2


def msize_direct(n, t):
    return sum(1 for x in range(1, n + 1) if x < t - x <= n)


def union_edges(n, T):
    E = []
    for t in T:
        for x in range(max(1, t - n), (t - 1) // 2 + 1):
            y = t - x
            if x < y <= n:
                E.append((x, y))
    return E


def is_path_graph(n, E):
    """The graph ([n], E) is itself a Hamiltonian path."""
    if len(E) != n - 1:
        return False
    deg = [0] * (n + 1)
    adj = [[] for _ in range(n + 1)]
    for x, y in E:
        deg[x] += 1
        deg[y] += 1
        adj[x].append(y)
        adj[y].append(x)
    if n > 1 and max(deg[1:]) > 2:
        return False
    seen = {1}
    st = [1]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                st.append(w)
    return len(seen) == n


def has_cycle(n, E):
    parent = list(range(n + 1))

    def f(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for x, y in E:
        rx, ry = f(x), f(y)
        if rx == ry:
            return True
        parent[rx] = ry
    return False


def a2_prediction(n):
    out = set()
    for s in (n, n + 1):
        out.add((s, s + 1))
    for s in (n - 1, n, n + 1):
        if s % 2 == 1 and s >= 3:
            out.add((s, s + 2))
    return {p for p in out if 3 <= p[0] and p[1] <= 2 * n - 1}


def a3_cond_i(n, a, b, c):
    lo, hi = max(1, c - n), min(n, a - 1)
    return all(2 * x in (a, b, c) for x in range(lo, hi + 1))


def a3_cond_ii(n, a, b, c):
    return msize(n, a) + msize(n, b) + msize(n, c) == n - 1


def a3_cond_iii(a, b, c):
    g = gcd(b - a, c - b)
    return g == 1 or (g == 2 and a % 2 == 1)


def a3_generic_table(n, a, b, c):
    """Corollary A3': explicit form of (ii) in the generic case c - a >= n."""
    m = c - a
    ub = b - n - 1
    if n % 2 == 1:
        return abs(ub) <= 1 and (m == n or (m == n + 1 and a % 2 == n % 2))
    if ub == 0:
        return m == n + 1 or (m == n and a % 2 == 0) or (m == n + 2 and a % 2 == 1)
    return 1 <= abs(ub) <= 2 and m == n and a % 2 == 1


def maxdeg_union(n, T):
    d = [0] * (n + 1)
    for x, y in union_edges(n, T):
        d[x] += 1
        d[y] += 1
    return max(d[1:]) if n >= 1 else 0


def zigzag(k):
    n = 4 * k - 1
    z = [0] * (n + 1)
    for j in range(k):
        z[4 * j + 1] = 2 * k - 2 * j
        z[4 * j + 2] = 2 * j + 1
        z[4 * j + 3] = 4 * k - 1 - 2 * j
        if 4 * j + 4 <= n:
            z[4 * j + 4] = 2 * k + 2 + 2 * j
    return z[1:]


def rotation_orbit_check(n, a, b, c):
    """Generic tight triple: along the union path, the vertices carrying a b-edge satisfy
    F(x) = x + (c-b) mod (c-a) at every return (b-edge then a/c-edge)."""
    E = union_edges(n, [a, b, c])
    adj = {v: [] for v in range(1, n + 1)}
    for x, y in E:
        adj[x].append(y)
        adj[y].append(x)
    ok = True
    for x in range(1, n + 1):
        y = b - x
        if not (1 <= y <= n and y != x):
            continue
        for z in adj[y]:
            if z == x:
                continue
            e = y + z
            if e not in (a, c):
                return False
            ok &= (z - x - (c - b)) % (c - a) == 0
    return ok


# ----------------------------------------------------------------------------------
# C_n structure (Task B)
# ----------------------------------------------------------------------------------
def pow2_le(n):
    p = 1
    while 2 * p <= n:
        p *= 2
    return p


def pow3_le(n):
    p = 1
    while 3 * p <= n:
        p *= 3
    return p


def top_layer_check(n):
    """Lemma T: every x in (M, n], M = max(P, 3^b), has neighbours exactly {2P-x} u {3^(b+1)-x if in [1,n]}."""
    P, Q = pow2_le(n), pow3_le(n)
    M = max(P, Q)
    Ts = set(SV.targets_for(n, 'pow23'))
    for x in range(M + 1, n + 1):
        nb = sorted(t - x for t in Ts if 1 <= t - x <= n and t - x != x)
        pred = [2 * P - x]
        if 1 <= 3 * Q - x <= n and 3 * Q - x != x:
            pred.append(3 * Q - x)
        if nb != sorted(pred):
            return False
    return True


def nbrs(n, family='pow23'):
    T = SV.targets_for(n, family)
    nb = [[] for _ in range(n + 1)]
    for x, y, t in SV.edges_of(n, T):
        nb[x].append(y)
        nb[y].append(x)
    return nb


def first_order_conflicts(n, family='pow23'):
    """Choke Lemma (first order).  Returns (leaves, conflicts); each conflict is (kind, centre, resolving set)
    where the resolving set is the set of possible second ends that could dissolve it (empty = never).
      over  : v has more degree-2 (non-leaf) neighbours than its need      -> the free end is one of them;
      choke : d (not a leaf) has at most one neighbour that is not saturated by degree-2 neighbours
              -> the free end is d or a saturating vertex (choke-1), or impossible (choke-0)."""
    nb = nbrs(n, family)
    deg = [len(a) for a in nb]
    L = [v for v in range(1, n + 1) if deg[v] <= 1]
    Ls = set(L)
    F = [[] for _ in range(n + 1)]
    for w in range(1, n + 1):
        if deg[w] == 2 and w not in Ls:
            for v in nb[w]:
                F[v].append(w)
    conf = []
    for v in range(1, n + 1):
        need = 1 if v in Ls else 2
        if len(F[v]) > need:
            conf.append(('over', v, frozenset(F[v]) if (len(F[v]) == need + 1 and v not in Ls) else frozenset()))
    for d in range(1, n + 1):
        if d in Ls:
            continue
        U, W = [], set()
        for v in nb[d]:
            Wv = [w for w in F[v] if w != d]
            if len(Wv) >= (1 if v in Ls else 2):
                U.append(v)
                W.update(Wv)
        avail = [v for v in nb[d] if v not in U]
        if len(avail) == 1:
            conf.append(('choke', d, frozenset(W) | {d}))
        elif len(avail) == 0:
            conf.append(('choke0', d, frozenset()))
    return L, conf


def first_order_verdict(n, family='pow23'):
    """'NONE' if the first-order conflicts already refute a Hamiltonian path, else 'OPEN'."""
    L, conf = first_order_conflicts(n, family)
    if len(L) > 2:
        return 'LOCAL', conf
    if not conf:
        return 'OPEN', conf
    if len(L) == 2:
        return 'NONE', conf
    common = None
    for kind, c, R in conf:
        common = set(R) if common is None else common & R
    return ('NONE' if not common else 'OPEN'), conf


def runs_of(ns):
    runs = []
    for n in sorted(ns):
        if runs and runs[-1][1] == n - 1:
            runs[-1][1] = n
        else:
            runs.append([n, n])
    return runs


def fmt_runs(runs):
    return ', '.join(f'{a}-{b}' if a != b else f'{a}' for a, b in runs)


def shortest_word(n, target_v, sources, maxlen=3):
    """Shortest word in the reflections r_t (t in targets of C_n with t <= n, i.e. full matchings) carrying a
    source to target_v, all intermediate points in [1, n].  Returns (length, path) or None."""
    T = [t for t in SV.targets_for(n, 'pow23') if t <= n]
    prev = {s: None for s in sources if 1 <= s <= n}
    if target_v in prev:
        return 0, [target_v]
    frontier = list(prev)
    for step in range(1, maxlen + 1):
        nf = []
        for x in frontier:
            for t in T:
                y = t - x
                if 1 <= y <= n and y != x and y not in prev:
                    prev[y] = (x, t)
                    nf.append(y)
                    if y == target_v:
                        path = [y]
                        while prev[path[-1]] is not None:
                            path.append(prev[path[-1]][0])
                        return step, path[::-1]
        frontier = nf
    return None


def boundary_word(n):
    """n* = T - v with T a top target (T > n) and v = w(s), s a power of 2 or 3, shortest w."""
    T = SV.targets_for(n, 'pow23')
    top = [t for t in T if t > n]
    sources = [2 ** j for j in range(0, 40) if 2 ** j <= n] + [3 ** j for j in range(1, 30) if 3 ** j <= n]
    best = None
    for t in top:
        r = shortest_word(n, t - n, sources)
        if r is not None and (best is None or r[0] < best[0]):
            best = (r[0], t, r[1])
    return best


# ----------------------------------------------------------------------------------
# parametric first-order families (Propositions W1-W3): endpoints are linear forms in P and B
# ----------------------------------------------------------------------------------
def pow23_set(lim):
    return set(SV.pow23(lim))


def nb_direct(n, x, Ts):
    """Neighbours of x in C_n without building the graph."""
    return {t - x for t in Ts if x < t <= x + n and t != 2 * x}


def level(a):
    B = 3 ** (a - 1)
    P = 1
    while P <= B:
        P *= 2
    return B, P   # B < P < 2B


def w_ranges(a):
    """The three families' n-ranges [lo, hi) (possibly empty), as linear forms in P, B."""
    B, P = level(a)
    q = P // 8
    r1 = (max(10 * q, 3 * B - 6 * q), min(12 * q, 3 * B - 4 * q))
    r2 = (max(2 * B + 2 * q, 2 * P - B + 2 * q, 10 * q, 3 * B - 6 * q), min(14 * q, 3 * B - 2 * q))
    r3 = (max(14 * q, 10 * q, 3 * B - 6 * q, 2 * P - B + 2 * q, 2 * B + 2 * q, 13 * q, 3 * B - 3 * q, 9 * q, 3 * B - 7 * q,
              2 * P - B + q, 2 * B + q), min(3 * B - 2 * q, 15 * q, 3 * B - q))
    return {'W1': r1, 'W2': r2, 'W3': r3}


def w_hypotheses(a, fam, n):
    """Check the hypotheses of Proposition fam at n directly from neighbour sets.  Returns True/False."""
    B, P = level(a)
    q = P // 8
    Ts = pow23_set(2 * n)

    def N(x):
        return nb_direct(n, x, Ts)

    def rigid(x):
        return len(N(x)) == 2

    def sat(v, ws):
        nv = N(v)
        return all(w in nv and rigid(w) for w in ws)
    if not (P <= n < 2 * P and B <= n < 3 * B):
        return False
    if N(P) != {3 * B - P}:
        return False
    if fam == 'W1':
        return (N(4 * q) == {B - 4 * q} and N(2 * q) == {6 * q, B - 2 * q} and sat(6 * q, [10 * q, 3 * B - 6 * q]))
    if fam == 'W2':
        return (N(2 * q) == {6 * q, B - 2 * q} and sat(6 * q, [10 * q, 3 * B - 6 * q])
                and sat(B - 2 * q, [2 * P - B + 2 * q, 2 * B + 2 * q]))
    if fam == 'W3':
        R4 = {2 * q, 10 * q, 3 * B - 6 * q, 2 * P - B + 2 * q, 2 * B + 2 * q}
        R8 = {q, 13 * q, 3 * B - 3 * q, 9 * q, 3 * B - 7 * q, 2 * P - B + q, 2 * B + q}
        return (N(2 * q) == {6 * q, B - 2 * q, 14 * q} and rigid(14 * q)
                and sat(6 * q, [10 * q, 3 * B - 6 * q]) and sat(B - 2 * q, [2 * P - B + 2 * q, 2 * B + 2 * q])
                and N(q) == {3 * q, 7 * q, B // 3 - q, B - q}
                and sat(3 * q, [13 * q, 3 * B - 3 * q]) and sat(7 * q, [9 * q, 3 * B - 7 * q])
                and sat(B - q, [2 * P - B + q, 2 * B + q]) and not (R4 & R8))
    return False
