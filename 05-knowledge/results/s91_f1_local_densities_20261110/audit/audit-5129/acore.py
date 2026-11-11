# audit-5129: core counting routines, written from the Setting section of THM-5129 only.
#
# Conventions.  A tournament on n vertices is an array adj of n out-masks (bit w of adj[v] set iff v -> w).
# A labelled tournament is also coded by an integer: pairs (i, j), i < j, in lexicographic order, bit k set iff i -> j.
#
#   D = H(B)                       number of Hamiltonian paths
#   E = E_B(y)                     Hamiltonian paths ending at y
#   U = U_B(y)                     ordered pairs (P, R) of vertex-disjoint directed paths, V(P) + V(R) = V(B),
#                                  R starting at y (so R is nonempty), P possibly empty
#   U_even, U_odd                  the parts with |V(R)| even, odd
#   Pi_2(C)                        ordered pairs (P, R) of vertex-disjoint paths covering C, each possibly empty
import numpy as np
from numba import njit


@njit(cache=True)
def popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@njit(cache=True)
def adj_from_code(code, n, adj):
    for i in range(n):
        adj[i] = 0
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> k) & 1:
                adj[i] |= (1 << j)
            else:
                adj[j] |= (1 << i)
            k += 1


@njit(cache=True)
def code_from_adj(adj, n):
    code = 0
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (adj[i] >> j) & 1:
                code |= (1 << k)
            k += 1
    return code


@njit(cache=True)
def tables(adj, n, hp, sp, Hs):
    """Subset dynamic programme.
    hp[m, v] = number of Hamiltonian paths of T[m] ending at v,
    sp[m, v] = number of Hamiltonian paths of T[m] starting at v,
    Hs[m]    = H(T[m]), Hs[0] = 1."""
    full = 1 << n
    fm = full - 1
    for m in range(full):
        for v in range(n):
            hp[m, v] = 0
            sp[m, v] = 0
    inm = np.zeros(n, np.int64)
    for v in range(n):
        for w in range(n):
            if (adj[v] >> w) & 1:
                inm[w] |= (1 << v)
    for v in range(n):
        hp[1 << v, v] = 1
        sp[1 << v, v] = 1
    Hs[0] = 1
    for m in range(1, full):
        tot = 0
        for v in range(n):
            if (m >> v) & 1:
                c = hp[m, v]
                if c != 0:
                    tot += c
                    ext = adj[v] & (~m) & fm
                    w = 0
                    while ext != 0:
                        if ext & 1:
                            hp[m | (1 << w), w] += c
                        ext >>= 1
                        w += 1
                c2 = sp[m, v]
                if c2 != 0:
                    ext = inm[v] & (~m) & fm
                    w = 0
                    while ext != 0:
                        if ext & 1:
                            sp[m | (1 << w), w] += c2
                        ext >>= 1
                        w += 1
        Hs[m] = tot


@njit(cache=True)
def root_stats(n, hp, sp, Hs, y, out):
    """out = (D, E, U_even, U_odd) for the root y."""
    fm = (1 << n) - 1
    ue = 0
    uo = 0
    for m in range(1, fm + 1):
        if (m >> y) & 1:
            c = sp[m, y]
            if c != 0:
                val = c * Hs[fm ^ m]
                if popcount(m) & 1:
                    uo += val
                else:
                    ue += val
    out[0] = Hs[fm]
    out[1] = hp[fm, y]
    out[2] = ue
    out[3] = uo


@njit(cache=True)
def pi2_from_tables(n, Hs):
    fm = (1 << n) - 1
    s = 0
    for m in range(fm + 1):
        s += Hs[m] * Hs[fm ^ m]
    return s


@njit(cache=True)
def ham_count(adj, n, dp):
    """H(T) by the forward subset programme only (dp is an (2^n, n) int64 scratch array)."""
    full = 1 << n
    fm = full - 1
    for m in range(full):
        for v in range(n):
            dp[m, v] = 0
    for v in range(n):
        dp[1 << v, v] = 1
    for m in range(1, full):
        for v in range(n):
            if (m >> v) & 1:
                c = dp[m, v]
                if c != 0:
                    ext = adj[v] & (~m) & fm
                    w = 0
                    while ext != 0:
                        if ext & 1:
                            dp[m | (1 << w), w] += c
                        ext >>= 1
                        w += 1
    tot = 0
    for v in range(n):
        tot += dp[fm, v]
    return tot


@njit(cache=True)
def is_strong(adj, n):
    """Strongly connected: every vertex reachable from 0 and 0 reachable from every vertex."""
    fm = (1 << n) - 1
    reach = 1
    frontier = 1
    while frontier:
        nxt = 0
        for v in range(n):
            if (frontier >> v) & 1:
                nxt |= adj[v]
        nxt &= ~reach
        reach |= nxt
        frontier = nxt
    if reach != fm:
        return False
    # reverse reachability
    reach = 1
    changed = True
    while changed:
        changed = False
        for v in range(n):
            if not ((reach >> v) & 1):
                if adj[v] & reach:
                    reach |= (1 << v)
                    changed = True
    return reach == fm


# ---------------------------------------------------------------- literal counts (no subset programme)

@njit(cache=True)
def literal_stats(adj, n, perms, E, Ue, Uo):
    """Literal enumeration of all n! vertex orders.
    A permutation with no bad consecutive pair is a Hamiltonian path; cut before position k it gives the pair
    (P, R) = (perm[:k], perm[k:]) for every k.  A permutation with exactly one bad pair, at positions (b, b+1),
    gives the single pair (perm[:b+1], perm[b+1:]).  Every pair (P, R) with R nonempty arises exactly once.
    Returns H; fills E[y], Ue[y], Uo[y]."""
    for v in range(n):
        E[v] = 0
        Ue[v] = 0
        Uo[v] = 0
    H = 0
    for p in range(perms.shape[0]):
        nbad = 0
        bad = -1
        for i in range(n - 1):
            if not ((adj[perms[p, i]] >> perms[p, i + 1]) & 1):
                nbad += 1
                bad = i
                if nbad > 1:
                    break
        if nbad == 0:
            H += 1
            E[perms[p, n - 1]] += 1
            for k in range(n):
                y = perms[p, k]
                if (n - k) & 1:
                    Uo[y] += 1
                else:
                    Ue[y] += 1
        elif nbad == 1:
            k = bad + 1
            y = perms[p, k]
            if (n - k) & 1:
                Uo[y] += 1
            else:
                Ue[y] += 1
    return H


@njit(cache=True)
def literal_H(adj, n, perms):
    H = 0
    for p in range(perms.shape[0]):
        ok = True
        for i in range(n - 1):
            if not ((adj[perms[p, i]] >> perms[p, i + 1]) & 1):
                ok = False
                break
        if ok:
            H += 1
    return H


@njit(cache=True)
def literal_pi2(adj, n, perms):
    """Ordered pairs (P, R) covering the vertex set, each possibly empty: permutation + cut position k = 0..n."""
    if n == 0:
        return 1
    tot = 0
    for p in range(perms.shape[0]):
        nbad = 0
        for i in range(n - 1):
            if not ((adj[perms[p, i]] >> perms[p, i + 1]) & 1):
                nbad += 1
                if nbad > 1:
                    break
        if nbad == 0:
            tot += n + 1
        elif nbad == 1:
            tot += 1
    return tot


def perm_array(n):
    import itertools
    return np.array(list(itertools.permutations(range(n))), dtype=np.int64).reshape(-1, max(n, 1))


# ---------------------------------------------------------------- closings

@njit(cache=True)
def head_adj(adj, n, y, out):
    """B + z: z = n, z -> y, every other vertex of B -> z."""
    for v in range(n):
        out[v] = adj[v]
        if v != y:
            out[v] |= (1 << n)
    out[n] = (1 << y)


@njit(cache=True)
def tail_adj(adj, n, y, t, out):
    """B => TT_t => z, z -> y, every other vertex of B -> z.  TT vertices n..n+t-1, z = n+t."""
    z = n + t
    ttmask = 0
    for i in range(t):
        ttmask |= (1 << (n + i))
    for v in range(n):
        out[v] = adj[v] | ttmask
        if v != y:
            out[v] |= (1 << z)
    for i in range(t):
        m = 1 << z
        for j in range(i + 1, t):
            m |= (1 << (n + j))
        out[n + i] = m
    out[z] = (1 << y)


# ---------------------------------------------------------------- strong components

@njit(cache=True)
def components(adj, n, comp):
    """comp[v] = index of the strong component of v in the condensation order (0 = initial component, which beats
    all the others).  Returns the number of components.  Uses scores of the condensation: in a tournament the
    components are linearly ordered, C_i => C_j for i < j."""
    fm = (1 << n) - 1
    # forward reachability sets
    reach = np.zeros(n, np.int64)
    for v in range(n):
        r = 1 << v
        fr = r
        while fr:
            nx = 0
            for u in range(n):
                if (fr >> u) & 1:
                    nx |= adj[u]
            nx &= ~r
            r |= nx
            fr = nx
        reach[v] = r
    # v, w in the same component iff they reach each other; order components by decreasing size of reach set
    for v in range(n):
        comp[v] = -1
    ncomp = 0
    # repeatedly take the vertex with the largest reach set among unassigned
    while True:
        best = -1
        bs = -1
        for v in range(n):
            if comp[v] < 0:
                s = popcount(reach[v])
                if s > bs:
                    bs = s
                    best = v
        if best < 0:
            break
        for w in range(n):
            if comp[w] < 0 and ((reach[best] >> w) & 1) and ((reach[w] >> best) & 1):
                comp[w] = ncomp
        ncomp += 1
    return ncomp
