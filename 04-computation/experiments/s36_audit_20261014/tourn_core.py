#!/usr/bin/env python3
"""Independent referee code (S36 audit): Hamiltonian-path statistics of tournaments.

Written from scratch for the audit of HYP-9280 (session opus-2026-10-14-S36); it does not import the session's scripts.
Methods deliberately differ from the session's subset DP:
  * Hamiltonian paths are ENUMERATED explicitly by an iterative depth-first search (and, for n <= 6, cross-checked
    against a brute-force scan of all n! vertex orders);
  * strong connectivity of T uses Landau's score criterion; strong connectivity of the forward digraph F uses a
    Warshall transitive closure; transitivity of M uses the Kendall-Babington Smith 3-cycle count.
Tournament encoding: adj[v] = bitmask of out-neighbours of v.
"""
import numpy as np
import numba as nb


@nb.njit(cache=True)
def popcount(x):
    c = 0
    while x:
        x &= x - 1
        c += 1
    return c


@nb.njit(cache=True)
def decode(n, code, pi, pj, adj):
    """code bit e = 1 means pi[e] -> pj[e], else pj[e] -> pi[e]"""
    for v in range(n):
        adj[v] = 0
    for e in range(pi.shape[0]):
        if (code >> e) & 1:
            adj[pi[e]] |= 1 << pj[e]
        else:
            adj[pj[e]] |= 1 << pi[e]


@nb.njit(cache=True)
def landau_strong(n, adj):
    s = np.empty(n, dtype=np.int64)
    for v in range(n):
        s[v] = popcount(adj[v])
    s.sort()
    acc = 0
    for k in range(1, n):
        acc += s[k - 1]
        if acc <= k * (k - 1) // 2:
            return False
    return True


@nb.njit(cache=True)
def closure_strong(n, g):
    """g = bitmask digraph; strongly connected iff the reflexive transitive closure is complete (Warshall)."""
    R = np.empty(n, dtype=np.int64)
    for v in range(n):
        R[v] = g[v] | (1 << v)
    for k in range(n):
        for i in range(n):
            if (R[i] >> k) & 1:
                R[i] |= R[k]
    full = (1 << n) - 1
    for i in range(n):
        if R[i] != full:
            return False
    return True


@nb.njit(cache=True)
def three_cycles_from_scores(n, g):
    """number of cyclic triples of a tournament g (Kendall-Babington Smith)"""
    t = n * (n - 1) * (n - 2) // 6
    for v in range(n):
        s = popcount(g[v])
        t -= s * (s - 1) // 2
    return t


@nb.njit(cache=True)
def enum_paths(n, adj, C, first, Nd, inmask):
    """Enumerate all Hamiltonian paths by iterative DFS; accumulate
    C[x, y] = #paths with x before y, first[v] = #paths starting at v,
    Nd[w, k] = #paths in which exactly k in-neighbours of w precede w.  Returns H."""
    for i in range(n):
        first[i] = 0
        for j in range(n):
            C[i, j] = 0
            Nd[i, j] = 0
    H = 0
    path = np.empty(n, dtype=np.int64)
    nxt = np.zeros(n + 1, dtype=np.int64)
    for s0 in range(n):
        path[0] = s0
        visited = 1 << s0
        depth = 1
        nxt[1] = 0
        while depth >= 1:
            if depth == n:
                H += 1
                first[path[0]] += 1
                for i in range(n):
                    a = path[i]
                    for j in range(i + 1, n):
                        C[a, path[j]] += 1
                pm = 0
                for i in range(n):
                    w = path[i]
                    Nd[w, popcount(pm & inmask[w])] += 1
                    pm |= 1 << w
                depth -= 1
                visited &= ~(1 << path[depth])
                continue
            v = path[depth - 1]
            c = nxt[depth]
            while c < n and (((visited >> c) & 1) == 1 or ((adj[v] >> c) & 1) == 0):
                c += 1
            if c < n:
                nxt[depth] = c + 1
                path[depth] = c
                visited |= 1 << c
                depth += 1
                if depth < n:
                    nxt[depth] = 0
            else:
                depth -= 1
                if depth >= 1:
                    visited &= ~(1 << path[depth])
    return H


@nb.njit(cache=True)
def perm_paths(n, adj, perms, C, first, Nd, inmask):
    """Same statistics by brute force over all vertex orders (rows of perms)."""
    for i in range(n):
        first[i] = 0
        for j in range(n):
            C[i, j] = 0
            Nd[i, j] = 0
    H = 0
    for r in range(perms.shape[0]):
        ok = True
        for i in range(n - 1):
            if ((adj[perms[r, i]] >> perms[r, i + 1]) & 1) == 0:
                ok = False
                break
        if not ok:
            continue
        H += 1
        first[perms[r, 0]] += 1
        for i in range(n):
            for j in range(i + 1, n):
                C[perms[r, i], perms[r, j]] += 1
        pm = 0
        for i in range(n):
            w = perms[r, i]
            Nd[w, popcount(pm & inmask[w])] += 1
            pm |= 1 << w
    return H


@nb.njit(cache=True)
def is_Tn(n, adj):
    """Is the tournament isomorphic to T_n = transitive with the extreme arc reversed?
    Test: some arc x -> y whose reversal gives a transitive tournament with y the source and x the sink."""
    g = np.empty(n, dtype=np.int64)
    for x in range(n):
        for y in range(n):
            if x == y or ((adj[x] >> y) & 1) == 0:
                continue
            for v in range(n):
                g[v] = adj[v]
            g[x] &= ~(1 << y)
            g[y] |= 1 << x
            if three_cycles_from_scores(n, g) == 0 and popcount(g[y]) == n - 1 and popcount(g[x]) == 0:
                return True
    return False


# statistics vector layout (per strong tournament)
#  0 H, 1 M transitive, 2 F strong, 3 #inverted arcs, 4 forward 3-cycle exists, 5 best 3-cycle cert numerator,
#  6 vertex-lemma numerator (min_w max_{u->w} C[u,w]), 7 min_w D_in numerator (sum 2C-H), 8 min_w D_out numerator,
#  9 max_w first[w], 10 #SS violations (vertices), 11 #violations of max_u C[u,w] + first[w] >= H (vertices),
#  12 delta numerator (max over pairs min(C, H-C)), 13 #3-cycles of M, 14 #3-cycles of T, 15 consistency errors,
#  16 #vertices with D_in <= 0, 17 #vertices failing the vertex lemma, 18 squeezed-triangle failures among (1,1) types,
#  19 #inverted arcs with p < 1/3
NSTAT = 20


@nb.njit(cache=True)
def analyse(n, adj, C, first, Nd, inmask, H, out, typehist):
    """Fill out[0:NSTAT] from the path statistics; add inverted-arc degree types to typehist[outdeg y * 16 + indeg x]."""
    for k in range(NSTAT):
        out[k] = 0
    out[0] = H
    err = 0
    if H % 2 == 0:
        err += 1
    # majority tournament and forward digraph
    madj = np.zeros(n, dtype=np.int64)
    fadj = np.zeros(n, dtype=np.int64)
    ninv = 0
    nthird = 0
    for x in range(n):
        for y in range(n):
            if x == y:
                continue
            if C[x, y] + C[y, x] != H:
                err += 1
            if 2 * C[x, y] > H:
                madj[x] |= 1 << y
                if (adj[x] >> y) & 1:
                    fadj[x] |= 1 << y
            if ((adj[x] >> y) & 1) and 2 * C[x, y] < H:
                ninv += 1
                if 3 * C[x, y] < H:
                    nthird += 1
                od = popcount(adj[y])
                idg = popcount(inmask[x])
                typehist[od * 16 + idg] += 1
                if od == 1 and idg == 1:
                    z = 0
                    while ((adj[y] >> z) & 1) == 0:
                        z += 1
                    if not ((inmask[x] >> z) & 1):
                        out[18] += 1
    m3 = three_cycles_from_scores(n, madj)
    out[1] = 1 if m3 == 0 else 0
    out[2] = 1 if closure_strong(n, fadj) else 0
    out[3] = ninv
    out[13] = m3
    out[14] = three_cycles_from_scores(n, adj)
    out[19] = nthird
    # 3-cycles of T: forward 3-cycle and best certificate
    best = -1
    fwd = 0
    for x in range(n):
        for y in range(n):
            if y == x or ((adj[x] >> y) & 1) == 0:
                continue
            for z in range(n):
                if z == x or z == y:
                    continue
                if ((adj[y] >> z) & 1) and ((adj[z] >> x) & 1):
                    m = C[x, y]
                    if C[y, z] < m:
                        m = C[y, z]
                    if C[z, x] < m:
                        m = C[z, x]
                    if m > best:
                        best = m
                    if 2 * C[x, y] > H and 2 * C[y, z] > H and 2 * C[z, x] > H:
                        fwd = 1
    out[4] = fwd
    out[5] = best
    # vertex statistics
    vl = 1 << 40
    dimin = 1 << 40
    domin = 1 << 40
    fmax = -1
    ssv = 0
    bfv = 0
    dneg = 0
    vlfail = 0
    for w in range(n):
        bi = -1
        din = 0
        dout = 0
        d = 0
        sumC = 0
        for u in range(n):
            if u == w:
                continue
            if (inmask[w] >> u) & 1:
                d += 1
                sumC += C[u, w]
                din += 2 * C[u, w] - H
                if C[u, w] > bi:
                    bi = C[u, w]
            else:
                dout += 2 * C[w, u] - H
        if bi < vl:
            vl = bi
        if 2 * bi < H:
            vlfail += 1
        if din <= 0:
            dneg += 1
        if din < dimin:
            dimin = din
        if dout < domin:
            domin = dout
        if first[w] > fmax:
            fmax = first[w]
        if bi + first[w] < H:
            bfv += 1
        # consistency: sum_k Nd = H, sum_k k Nd = sum_u C[u,w], Nd[w,0] = first[w]
        tot = 0
        mom = 0
        for k in range(n):
            tot += Nd[w, k]
            mom += k * Nd[w, k]
        if tot != H or mom != sumC or Nd[w, 0] != first[w]:
            err += 1
        # stochastic symmetry (SS): Pr[N >= d - k] >= Pr[N <= k] for all k
        le = 0
        for k in range(d + 1):
            le += Nd[w, k]
            ge = 0
            for j in range(d - k, d + 1):
                ge += Nd[w, j]
            if ge < le:
                ssv += 1
                break
    out[6] = vl
    out[7] = dimin
    out[8] = domin
    out[9] = fmax
    out[10] = ssv
    out[11] = bfv
    out[16] = dneg
    out[17] = vlfail
    # best balance delta
    dl = 0
    for x in range(n):
        for y in range(x + 1, n):
            m = C[x, y]
            if H - C[x, y] < m:
                m = H - C[x, y]
            if m > dl:
                dl = m
    out[12] = dl
    out[15] = err


@nb.njit(cache=True)
def inmasks(n, adj, inmask):
    for w in range(n):
        inmask[w] = 0
    for u in range(n):
        for w in range(n):
            if (adj[u] >> w) & 1:
                inmask[w] |= 1 << u
