"""Independent audit helpers (audit of opus-2026-10-21-S50).  Written from scratch; imports nothing from the session.

Conventions: a tournament on n vertices is a list/array `adj` of out-neighbour bitmasks (bit w of adj[u] set iff u -> w).
"""
import itertools
import numpy as np
import numba as nb


# ------------------------------------------------------------------ Hamiltonian paths: two independent algorithms
@nb.njit(cache=True)
def ham_dp(n, adj):
    """Subset DP over (covered set, last vertex), iterating over out-neighbour bits."""
    N = 1 << n
    dp = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = dp[S, v]
            if c == 0:
                continue
            free = adj[v] & ~S
            while free:
                low = free & (-free)
                u = 0
                t = low
                while t > 1:
                    t >>= 1
                    u += 1
                dp[S | low, u] += c
                free ^= low
    tot = 0
    for v in range(n):
        tot += dp[N - 1, v]
    return tot


@nb.njit(cache=True)
def ham_ie(n, adj):
    """Inclusion-exclusion: #Ham paths = sum_S (-1)^(n-|S|) * #(walks on n vertices inside S)."""
    tot = 0
    vec = np.zeros(n, dtype=np.int64)
    nxt = np.zeros(n, dtype=np.int64)
    for S in range(1, 1 << n):
        k = 0
        for v in range(n):
            if (S >> v) & 1:
                vec[v] = 1
                k += 1
            else:
                vec[v] = 0
        for step in range(n - 1):
            for u in range(n):
                nxt[u] = 0
            for v in range(n):
                if vec[v] != 0:
                    o = adj[v] & S
                    for u in range(n):
                        if (o >> u) & 1:
                            nxt[u] += vec[v]
            for u in range(n):
                vec[u] = nxt[u]
        w = 0
        for u in range(n):
            w += vec[u]
        if (n - k) % 2 == 0:
            tot += w
        else:
            tot -= w
    return tot


def ham_brute(adj):
    n = len(adj)
    return sum(1 for p in itertools.permutations(range(n))
               if all((adj[p[i]] >> p[i + 1]) & 1 for i in range(n - 1)))


def H(adj):
    return int(ham_dp(len(adj), np.array(adj, dtype=np.int64)))


# ------------------------------------------------------------------ ear H values by a quadratic form
@nb.njit(cache=True)
def end_table(n, adj):
    """F[S, v] = number of paths with vertex set S ending at v."""
    N = 1 << n
    F = np.zeros((N, n), dtype=np.int64)
    for v in range(n):
        F[1 << v, v] = 1
    for S in range(1, N):
        for v in range(n):
            c = F[S, v]
            if c == 0:
                continue
            for u in range(n):
                if not (S >> u) & 1 and (adj[v] >> u) & 1:
                    F[S | (1 << u), u] += c
    return F


@nb.njit(cache=True)
def ear_values(n, adj, radj):
    """For a tournament on n vertices (adj) with converse radj, return Hs[x] = H of the one-vertex extension in which the
    new vertex w beats exactly the vertices of the bitmask x (and loses to the rest), for every x in [0, 2^n).
    H(x) = |F_full| + (G_full - F_full + a) . x - x^T M x,  a = sum_S |F_S| G_{S^c},  M = sum_S F_S G_{S^c}^T."""
    N = 1 << n
    full = N - 1
    F = end_table(n, adj)
    G = end_table(n, radj)        # paths of the converse ending at u = paths of T starting at u
    a = np.zeros(n, dtype=np.int64)
    M = np.zeros((n, n), dtype=np.int64)
    for S in range(1, full):
        Sc = full ^ S
        fs = 0
        for u in range(n):
            fs += F[S, u]
        for u2 in range(n):
            g = G[Sc, u2]
            if g != 0:
                a[u2] += fs * g
                for u in range(n):
                    if F[S, u] != 0:
                        M[u, u2] += F[S, u] * g
    c0 = 0
    b = np.zeros(n, dtype=np.int64)
    for u in range(n):
        c0 += F[full, u]
        b[u] = G[full, u] - F[full, u] + a[u]
    Q = np.zeros((n, n), dtype=np.int64)
    for u in range(n):
        for u2 in range(u + 1, n):
            Q[u, u2] = M[u, u2] + M[u2, u]
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


def converse(adj):
    n = len(adj)
    full = (1 << n) - 1
    return [full ^ adj[v] ^ (1 << v) for v in range(n)]


def extend(adj, x):
    """Child of adj by a new last vertex w = n that beats exactly the bitmask x of old vertices."""
    n = len(adj)
    new = []
    for v in range(n):
        a = adj[v]
        if not (x >> v) & 1:
            a |= 1 << n
        new.append(a)
    new.append(x)
    return new


# ------------------------------------------------------------------ canonical form (refinement + individualisation)
def _refine(outl, col):
    n = len(col)
    while True:
        sig = [(col[u], tuple(sorted(col[w] for w in outl[u]))) for u in range(n)]
        keys = sorted(set(sig))
        rank = {k: i for i, k in enumerate(keys)}
        new = [rank[s] for s in sig]
        if len(keys) == len(set(col)):
            return new
        col = new


def canon(adj, want_aut=False):
    """Exact canonical code: minimum, over all leaves of the full individualisation-refinement tree, of the upper-
    triangular arc code in the leaf order.  The number of leaves attaining the minimum is |Aut| (free action)."""
    n = len(adj)
    outl = [[w for w in range(n) if (adj[u] >> w) & 1] for u in range(n)]
    best = [None, 0]

    def code_of(order):
        c = 0
        e = 0
        for i in range(n):
            ai = adj[order[i]]
            for j in range(i + 1, n):
                if (ai >> order[j]) & 1:
                    c |= 1 << e
                e += 1
        return c

    def search(col):
        col = _refine(outl, col)
        if len(set(col)) == n:
            order = sorted(range(n), key=lambda v: col[v])
            c = code_of(order)
            if best[0] is None or c < best[0]:
                best[0] = c
                best[1] = 1
            elif c == best[0]:
                best[1] += 1
            return
        cells = {}
        for v in range(n):
            cells.setdefault(col[v], []).append(v)
        target = min((len(vs), k) for k, vs in cells.items() if len(vs) > 1)[1]
        for v in cells[target]:
            new = [2 * c + (1 if (c == target and u != v) else 0) for u, c in enumerate(col)]
            search(new)

    search([len(outl[u]) for u in range(n)])
    if want_aut:
        return best[0], best[1]
    return best[0]


def decode(n, code):
    adj = [0] * n
    e = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> e) & 1:
                adj[i] |= 1 << j
            else:
                adj[j] |= 1 << i
            e += 1
    return adj


def is_strong(adj):
    n = len(adj)
    # reachability from 0 forward and backward
    def reach(nbr):
        seen = 1
        stack = [0]
        while stack:
            u = stack.pop()
            new = nbr(u) & ~seen
            seen |= new
            while new:
                low = new & -new
                stack.append(low.bit_length() - 1)
                new ^= low
        return seen
    full = (1 << n) - 1
    inn = [sum(1 << u for u in range(n) if (adj[u] >> v) & 1) for v in range(n)]
    return reach(lambda u: adj[u]) == full and reach(lambda u: inn[u]) == full


def c3count(adj):
    n = len(adj)
    s = [bin(a).count('1') for a in adj]
    return n * (n - 1) * (n - 2) // 6 - sum(x * (x - 1) // 2 for x in s)


# ------------------------------------------------------------------ interval chains
def chain(parts):
    """IC(p_1..p_k): spine 0..N, all arcs forward except one back arc l_(i+1) -> l_i per block."""
    N = sum(parts)
    n = N + 1
    adj = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            adj[i] |= 1 << j
    l = 0
    for p in parts:
        r = l + p
        adj[l] &= ~(1 << r)
        adj[r] |= 1 << l
        l = r
    return adj


def compositions(total, minpart=2, maxpart=None):
    if total == 0:
        yield ()
        return
    hi = total if maxpart is None else min(total, maxpart)
    for p in range(minpart, hi + 1):
        for rest in compositions(total - p, minpart, maxpart):
            yield (p,) + rest


def busch(m):
    q, r = divmod(m, 3)
    return {0: 3 * 5 ** (q - 1), 1: 5 ** q, 2: 9 * 5 ** (q - 1)}[r]


def busch_min(m):
    """min 3^x 5^y over 2x + 3y = m - 1 (computed by search, independent of the closed form)."""
    M = m - 1
    return min(3 ** x * 5 ** y for y in range(M // 3 + 1) for x in range(M // 2 + 1) if 2 * x + 3 * y == M)


def prod_law(parts):
    out = 1
    for p in parts:
        out *= 2 ** (p - 1) + 1
    return out
