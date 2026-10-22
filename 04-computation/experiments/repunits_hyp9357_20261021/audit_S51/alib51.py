"""Independent audit helpers for the opus-2026-10-21-S51 packet (THM-4748, THM-4749, HYP-9359).
Written from scratch by the S51 auditor; imports nothing from the session's scripts and nothing from the S50 audit.

Convention: a tournament on n vertices is an int64 array `adj` of out-neighbour bitmasks (bit w of adj[u] set iff u -> w).
"""
import itertools
import numpy as np
import numba as nb


# ------------------------------------------------------------------ basic constructions
def transitive(m):
    """i -> j for i < j"""
    adj = np.zeros(m, dtype=np.int64)
    for i in range(m):
        for j in range(i + 1, m):
            adj[i] |= 1 << j
    return adj


def reverse_arcs(m, backs):
    """transitive order on m vertices with the arcs (y -> x) for (y, x) in backs, y > x"""
    adj = transitive(m)
    for (y, x) in backs:
        assert y > x
        adj[x] &= ~(1 << y)
        adj[y] |= 1 << x
    return adj


def one_apex(m, A):
    """vertex 0 plus a transitive chain 1 < ... < m-1; 0 loses exactly to the chain positions in A"""
    adj = np.zeros(m, dtype=np.int64)
    for i in range(1, m):
        for j in range(i + 1, m):
            adj[i] |= 1 << j
    for j in range(1, m):
        if j in A:
            adj[j] |= 1
        else:
            adj[0] |= 1 << j
    return adj


def converse(adj):
    n = len(adj)
    full = (1 << n) - 1
    return np.array([full ^ int(adj[v]) ^ (1 << v) for v in range(n)], dtype=np.int64)


def is_tournament(adj):
    n = len(adj)
    for u in range(n):
        if (int(adj[u]) >> u) & 1:
            return False
        for w in range(u + 1, n):
            a = (int(adj[u]) >> w) & 1
            b = (int(adj[w]) >> u) & 1
            if a + b != 1:
                return False
    return True


# ------------------------------------------------------------------ Hamiltonian paths
@nb.njit(cache=True)
def ham(n, adj):
    """number of Hamiltonian paths; subset DP over (vertex set, last vertex)"""
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            nb_ = adj[v] & ~S
            for u in range(n):
                if (nb_ >> u) & 1:
                    dp[S | (1 << u), u] += c
    t = 0
    for v in range(n):
        t += dp[N - 1, v]
    return t


def ham_brute(adj):
    n = len(adj)
    return sum(1 for p in itertools.permutations(range(n))
               if all((int(adj[p[i]]) >> p[i + 1]) & 1 for i in range(n - 1)))


@nb.njit(cache=True)
def ham_subsets(n, adj):
    """Hs[S] = number of Hamiltonian paths of the subtournament induced on S, for every S"""
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    Hs = np.zeros(N, dtype=np.int64)
    for S in range(1, N):
        t = 0
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            t += c
            nb_ = adj[v] & ~S
            for u in range(n):
                if (nb_ >> u) & 1:
                    dp[S | (1 << u), u] += c
        Hs[S] = t
    return Hs


@nb.njit(cache=True)
def is_strong_mask(n, adj, S):
    """is the subtournament induced on the vertex set S strongly connected?"""
    if S == 0:
        return False
    low = S & -S
    # forward closure
    r = low
    while True:
        nr = r
        for u in range(n):
            if (r >> u) & 1:
                nr |= adj[u] & S
        if nr == r:
            break
        r = nr
    if r != S:
        return False
    # backward closure
    r = low
    while True:
        nr = r
        for u in range(n):
            if (S >> u) & 1 and not (r >> u) & 1 and (adj[u] & r):
                nr |= 1 << u
        if nr == r:
            break
        r = nr
    return r == S


@nb.njit(cache=True)
def fas(n, adj):
    """minimum feedback arc number: f[S] = min over the last vertex v of f[S - v] + #(arcs v -> S - v)"""
    N = 1 << n
    f = np.empty(N, dtype=np.int64)
    f[0] = 0
    for S in range(1, N):
        best = 1 << 40
        for v in range(n):
            if (S >> v) & 1:
                T = S ^ (1 << v)
                x = adj[v] & T
                c = 0
                while x:
                    x &= x - 1
                    c += 1
                if f[T] + c < best:
                    best = f[T] + c
        f[S] = best
    return f[N - 1]


# ------------------------------------------------------------------ canonical form (individualisation-refinement, full tree)
@nb.njit(cache=True)
def _lexcmp(sig, a, b, L):
    for c in range(L):
        if sig[a, c] < sig[b, c]:
            return -1
        if sig[a, c] > sig[b, c]:
            return 1
    return 0


@nb.njit(cache=True)
def _refine(n, adj, col):
    """equitable refinement in place; colours are canonical dense ranks; returns the number of cells"""
    k = 0
    for v in range(n):
        if col[v] + 1 > k:
            k = col[v] + 1
    sig = np.zeros((n, n + 1), dtype=np.int64)
    order = np.empty(n, dtype=np.int64)
    newcol = np.empty(n, dtype=np.int64)
    while True:
        for v in range(n):
            for c in range(k + 1):
                sig[v, c] = 0
            sig[v, 0] = col[v]
            x = adj[v]
            for u in range(n):
                if (x >> u) & 1:
                    sig[v, 1 + col[u]] += 1
        for v in range(n):
            order[v] = v
        for a in range(1, n):
            x = order[a]
            b = a - 1
            while b >= 0 and _lexcmp(sig, x, order[b], k + 1) < 0:
                order[b + 1] = order[b]
                b -= 1
            order[b + 1] = x
        c = 0
        newcol[order[0]] = 0
        for a in range(1, n):
            if _lexcmp(sig, order[a], order[a - 1], k + 1) != 0:
                c += 1
            newcol[order[a]] = c
        knew = c + 1
        for v in range(n):
            col[v] = newcol[v]
        if knew == k:
            return k
        k = knew


@nb.njit(cache=True)
def _leaf_code(n, adj, col):
    order = np.empty(n, dtype=np.int64)
    for v in range(n):
        order[col[v]] = v
    w0 = 0
    w1 = 0
    e = 0
    for i in range(n):
        ai = adj[order[i]]
        for j in range(i + 1, n):
            if (ai >> order[j]) & 1:
                if e < 60:
                    w0 |= 1 << e
                else:
                    w1 |= 1 << (e - 60)
            e += 1
    return w1, w0


@nb.njit(cache=True)
def canon(n, adj):
    """canonical code (w1, w0) = minimum over all leaves of the full search tree, and |Aut| = number of minimal leaves"""
    col = np.empty(n, dtype=np.int64)
    # initial colouring: score ranks
    sc = np.empty(n, dtype=np.int64)
    for v in range(n):
        x = adj[v]
        c = 0
        while x:
            x &= x - 1
            c += 1
        sc[v] = c
    for v in range(n):
        r = 0
        seen = np.zeros(n + 1, dtype=np.int64)
        for u in range(n):
            if sc[u] < sc[v] and seen[sc[u]] == 0:
                seen[sc[u]] = 1
                r += 1
        col[v] = r
    k = _refine(n, adj, col)
    best1 = -1
    best0 = -1
    cnt = 0
    if k == n:
        w1, w0 = _leaf_code(n, adj, col)
        return w1, w0, 1
    cols = np.empty((n + 1, n), dtype=np.int64)
    cells = np.empty((n + 1, n), dtype=np.int64)
    csz = np.zeros(n + 1, dtype=np.int64)
    ptr = np.zeros(n + 1, dtype=np.int64)
    depth = 0
    for v in range(n):
        cols[0, v] = col[v]
    # target cell at depth 0
    _target(n, cols, 0, cells, csz)
    ptr[0] = 0
    newcol = np.empty(n, dtype=np.int64)
    while depth >= 0:
        if ptr[depth] >= csz[depth]:
            depth -= 1
            continue
        v = cells[depth, ptr[depth]]
        ptr[depth] += 1
        c = cols[depth, v]
        for u in range(n):
            cu = cols[depth, u]
            if cu > c:
                newcol[u] = cu + 1
            elif cu == c and u != v:
                newcol[u] = c + 1
            else:
                newcol[u] = cu
        kk = _refine(n, adj, newcol)
        if kk == n:
            w1, w0 = _leaf_code(n, adj, newcol)
            if best1 < 0 or w1 < best1 or (w1 == best1 and w0 < best0):
                best1 = w1
                best0 = w0
                cnt = 1
            elif w1 == best1 and w0 == best0:
                cnt += 1
        else:
            depth += 1
            for u in range(n):
                cols[depth, u] = newcol[u]
            _target(n, cols, depth, cells, csz)
            ptr[depth] = 0
    return best1, best0, cnt


@nb.njit(cache=True)
def _target(n, cols, d, cells, csz):
    """smallest non-singleton cell (ties: smallest colour)"""
    size = np.zeros(n, dtype=np.int64)
    for u in range(n):
        size[cols[d, u]] += 1
    bc = -1
    bs = n + 1
    for c in range(n):
        if size[c] > 1 and size[c] < bs:
            bs = size[c]
            bc = c
    m = 0
    for u in range(n):
        if cols[d, u] == bc:
            cells[d, m] = u
            m += 1
    csz[d] = m


@nb.njit(cache=True)
def decode(n, w1, w0):
    """adjacency (out-masks) of the canonical representative of a code"""
    adj = np.zeros(n, dtype=np.int64)
    e = 0
    for i in range(n):
        for j in range(i + 1, n):
            if e < 60:
                b = (w0 >> e) & 1
            else:
                b = (w1 >> (e - 60)) & 1
            if b:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            e += 1
    return adj


@nb.njit(cache=True)
def extend(n, adj, x):
    """one-vertex extension: new vertex n beats exactly the old vertices in the bitmask x"""
    new = np.zeros(n + 1, dtype=np.int64)
    for v in range(n):
        a = adj[v]
        if not (x >> v) & 1:
            a |= 1 << n
        new[v] = a
    new[n] = x
    return new


# ------------------------------------------------------------------ ear values: H of every one-vertex extension
@nb.njit(cache=True)
def path_tables(n, adj):
    """F[S, u] = paths with vertex set S ending at u; G[S, u] = paths with vertex set S starting at u"""
    N = 1 << n
    F = np.zeros((N, n), dtype=np.int64)
    G = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        F[1 << v, v] = 1
        G[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = F[S, v]
            if c != 0:
                x = adj[v] & ~S
                for u in range(n):
                    if (x >> u) & 1:
                        F[S | (1 << u), u] += c
            g = G[S, v]
            if g != 0:
                for u in range(n):
                    if not (S >> u) & 1 and (adj[u] >> v) & 1:
                        G[S | (1 << u), u] += g
    return F, G


@nb.njit(cache=True)
def ear_coeffs(n, adj):
    """H(x) = c0 + sum_u b[u] x_u - sum_{u<u'} Q[u,u'] x_u x_u' for the extension by a vertex w beating the set x.
    Derivation: w first (w -> start u in x), w last (end u not in x -> w), or P1 w P2 with end(P1) = u not in x and
    start(P2) = u' in x:  sum_S F[S,u] G[S^c,u'] (1 - x_u) x_u'."""
    N = 1 << n
    full = N - 1
    F, G = path_tables(n, adj)
    M = np.zeros((n, n), dtype=np.int64)
    for S in range(1, full):
        Sc = full ^ S
        for u in range(n):
            f = F[S, u]
            if f == 0:
                continue
            for u2 in range(n):
                g = G[Sc, u2]
                if g != 0:
                    M[u, u2] += f * g
    c0 = 0
    for u in range(n):
        c0 += F[full, u]
    b = np.zeros(n, dtype=np.int64)
    for u2 in range(n):
        a = 0
        for u in range(n):
            a += M[u, u2]
        b[u2] = G[full, u2] - F[full, u2] + a
    Q = np.zeros((n, n), dtype=np.int64)
    for u in range(n):
        for u2 in range(n):
            if u != u2:
                Q[min(u, u2), max(u, u2)] += M[u, u2]
    return c0, b, Q


@nb.njit(cache=True)
def ear_all(n, adj):
    c0, b, Q = ear_coeffs(n, adj)
    N = 1 << n
    out = np.zeros(N, dtype=np.int64)
    for x in range(N):
        h = c0
        for u in range(n):
            if (x >> u) & 1:
                h += b[u]
                for u2 in range(u + 1, n):
                    if (x >> u2) & 1:
                        h -= Q[u, u2]
        out[x] = h
    return out


# ------------------------------------------------------------------ closed form for one-apex tournaments (path decomposition)
def one_apex_paths(m, A):
    """H for vertex 0 + transitive chain 1..m-1, 0 losing exactly to A (any A subset of 1..m-1).
    A Hamiltonian path is P1 0 P2 with P1, P2 increasing:  [1 not in A] (0 first) + [m-1 in A] (0 last)
    + sum_{t in A, i not in A, 1 <= i < t} 2^(t-i-1) + #{t in A : t+1 <= m-1, t+1 not in A}."""
    A = set(A)
    h = (1 if 1 not in A else 0) + (1 if (m - 1) in A else 0)
    for t in A:
        for i in range(1, t):
            if i not in A:
                h += 2 ** (t - i - 1)
        if t + 1 <= m - 1 and (t + 1) not in A:
            h += 1
    return h


def one_apex_pairwise(A):
    """THM-4748 (4): H - 1 = sum_{t in A} 2^(t-1) - sum_{s<t in A} 2^max(t-s-1, 1), A subset of {2, ..., m-1}"""
    A = sorted(A)
    return 1 + sum(2 ** (t - 1) for t in A) - sum(2 ** max(t - s - 1, 1) for s, t in itertools.combinations(A, 2))
