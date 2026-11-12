"""tt.py - auditor's own tournament toolkit (independent of scripts/chain_semigroup.py).

A tournament on n vertices is a list adj of n bitmasks: bit w of adj[v] is set iff v -> w.

tables(n, adj) returns (hsub, st, en):
    hsub[S]   = number of Hamiltonian paths of the sub-tournament induced on the vertex set S  (hsub[0] = 1)
    st[S, v]  = number of those paths that START at v
    en[S, v]  = number of those paths that END at v
From these:
    H(T)        = hsub[full]
    W+_T(u)     = sum over S containing u of st[S, u] * hsub[full ^ S]        (paths P starting at u, H(T - P))
    W-_T(u)     = sum over S containing u of en[S, u] * hsub[full ^ S]
    tau(T)      = sum over all S (S = 0 included) of hsub[S] * hsub[full ^ S]
    o, e        = W+ split by the parity of |S|
"""
import numpy as np
from numba import njit


@njit(cache=False)
def _tables(n, adj):
    full = (1 << n) - 1
    en = np.zeros((1 << n, n), dtype=np.int64)
    st = np.zeros((1 << n, n), dtype=np.int64)
    hsub = np.zeros(1 << n, dtype=np.int64)
    hsub[0] = 1
    for v in range(n):
        en[1 << v, v] = 1
        st[1 << v, v] = 1
    for S in range(1, full + 1):
        tot = 0
        for last in range(n):
            c = en[S, last]
            if c != 0:
                tot += c
                nxt = adj[last] & ~S
                for w in range(n):
                    if (nxt >> w) & 1:
                        en[S | (1 << w), w] += c          # append w after last
        hsub[S] = tot
        for first in range(n):
            c = st[S, first]
            if c != 0:
                for w in range(n):
                    if ((S >> w) & 1) == 0 and ((adj[w] >> first) & 1):
                        st[S | (1 << w), w] += c          # prepend w before first
    return hsub, st, en


def tables(n, adj):
    return _tables(n, np.array(adj, dtype=np.int64))


def H(n, adj):
    hsub, st, en = tables(n, adj)
    return int(hsub[(1 << n) - 1])


def popcount(x):
    return bin(x).count("1")


def invariants(n, adj, u, tabs=None):
    """returns dict with H, Wp (=W+), Wm (=W-), o, e (W+ split by parity of the number of vertices of the path),
    Om, Em (W- split), tau_minus_root (= tau of T - u)"""
    hsub, st, en = tabs if tabs is not None else tables(n, adj)
    full = (1 << n) - 1
    Wp = Wm = o = e = Om = Em = 0
    for S in range(1, full + 1):
        if (S >> u) & 1:
            comp = int(hsub[full ^ S])
            a = int(st[S, u]) * comp
            b = int(en[S, u]) * comp
            Wp += a
            Wm += b
            if popcount(S) & 1:
                o += a
                Om += b
            else:
                e += a
                Em += b
    rest = full ^ (1 << u)
    tau = 0
    S = rest
    while True:                      # all subsets S of rest, including 0
        tau += int(hsub[S]) * int(hsub[rest ^ S])
        if S == 0:
            break
        S = (S - 1) & rest
    return dict(H=int(hsub[full]), Wp=Wp, Wm=Wm, o=o, e=e, Om=Om, Em=Em, tau_minus_root=tau)


def ordinal_sum(blocks):
    """blocks: list of (n_i, adj_i).  Returns (n, adj, offsets): every vertex of an earlier block beats every
    vertex of a later block."""
    offs = []
    n = 0
    for (m, _) in blocks:
        offs.append(n)
        n += m
    adj = [0] * n
    for bi, (m, a) in enumerate(blocks):
        later = 0
        for bj in range(bi + 1, len(blocks)):
            later |= ((1 << blocks[bj][0]) - 1) << offs[bj]
        for v in range(m):
            adj[offs[bi] + v] = (a[v] << offs[bi]) | later
    return n, adj, offs


def reverse_arc(adj, a, b):
    """make the arc between a and b point b -> a (it must currently be a -> b)"""
    assert (adj[a] >> b) & 1 and not (adj[b] >> a) & 1
    adj[a] &= ~(1 << b)
    adj[b] |= 1 << a


def transitive(m):
    """transitive tournament 0 -> 1 -> ... -> m-1 (i < j implies i -> j)"""
    return m, [((1 << m) - 1) & ~((1 << (v + 1)) - 1) for v in range(m)]


def all_tournaments(n):
    """all labeled tournaments on n vertices"""
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        adj = [0] * n
        for t, (i, j) in enumerate(pairs):
            if (bits >> t) & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
        yield adj


def one_apex(word):
    """one-apex tournament: vertex 0 = apex, vertices 1..len(word) = spine (transitive, earlier beats later);
    word[j] = 1 iff the apex beats spine vertex j (the first spine vertex is word[0])."""
    n = len(word) + 1
    adj = [0] * n
    for j in range(len(word)):
        for k in range(j + 1, len(word)):
            adj[1 + j] |= 1 << (1 + k)
        if word[j]:
            adj[0] |= 1 << (1 + j)
        else:
            adj[1 + j] |= 1
    return n, adj


def g_fun(v):
    """the one-apex function: g(0) = 1, g(2v+1) = g(v), g(2v) = g(v) + v + (v mod 2)"""
    if v == 0:
        return 1
    if v & 1:
        return g_fun(v >> 1)
    h = v >> 1
    return g_fun(h) + h + (h & 1)


def word_PQ(word):
    """P = value of the letters 1 (first spine vertex = top bit); Q = value of the letters 0 (first spine vertex =
    unit bit)"""
    P = 0
    for a in word:
        P = 2 * P + a
    Q = sum(1 << j for j, a in enumerate(word) if not a)
    return P, Q
