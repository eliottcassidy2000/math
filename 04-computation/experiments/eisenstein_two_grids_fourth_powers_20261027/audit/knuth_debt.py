"""Independent audit solver (not derived from the producer scripts).

Window-restricted catalytic merge debt g_W on an integer lattice, computed with Knuth's generalisation of Dijkstra's
algorithm (Knuth 1977, superior functions): every rule x <- (y, z) has cost g(y) + g(z) + max(beta, 0), which is
>= max(g(y), g(z)) and nondecreasing in both arguments, so nodes can be finalised in increasing order of g and each
pair is combined exactly once.  This is a different algorithm from the producers' all-pairs Bellman rounds.

Lattice points are integer vectors P (rows); Q(x) = x^T G x / s with integer G, s; beta(y, z) = 2 y^T G z / s.
All nodes of every tree are restricted to the window Q <= W (the point list itself is the window).  Unit charges are
the minimal nonzero vectors of the window (callers assert that this is the true minimal norm).

Returns (g, Q) as numpy arrays aligned with P; g = INF where unreachable inside the window.
"""
import numpy as np

INF = np.int64(1) << 60


def knuth_debt(P, G, s, units=None, record_parent=False):
    P = np.ascontiguousarray(np.asarray(P, dtype=np.int64))
    N, d = P.shape
    G = np.asarray(G, dtype=np.int64)
    PG = P @ G
    Q2 = np.einsum('ij,ij->i', PG, P)
    assert np.all(Q2 % s == 0), "form not integral on the given points"
    Q = Q2 // s
    R = int(np.abs(P).max())
    B = 4 * R + 3
    assert B ** d < (1 << 62), "code overflow"
    pw = np.array([B ** k for k in range(d)], dtype=np.int64)
    codes = P @ pw
    order = np.argsort(codes, kind='stable')
    sc = codes[order]
    assert np.all(np.diff(sc) > 0), "duplicate points"
    if units is None:
        m = Q[Q > 0].min()
        units = np.where(Q == m)[0]
    T = np.full(N, INF, dtype=np.int64)
    T[units] = 0
    Tm = T.copy()
    fin = np.zeros(N, dtype=bool)
    finlist = np.empty(N, dtype=np.int64)
    nf = 0
    parent = {} if record_parent else None
    SENT = INF + 1
    while True:
        i = int(np.argmin(Tm))
        if Tm[i] >= INF:
            break
        fin[i] = True
        Tm[i] = SENT
        finlist[nf] = i
        nf += 1
        J = finlist[:nf]
        sums = codes[i] + codes[J]
        pos = np.searchsorted(sc, sums)
        pos[pos >= N] = N - 1
        ok = sc[pos] == sums
        if not ok.any():
            continue
        Jv = J[ok]
        K = order[pos[ok]]
        b2 = (P[Jv] @ PG[i]) * 2
        assert np.all(b2 % s == 0)
        beta = b2 // s
        # consistency: beta = Q(x) - Q(y) - Q(z)
        assert np.all(beta == Q[K] - Q[i] - Q[Jv])
        cost = T[i] + T[Jv] + np.maximum(beta, 0)
        nfk = ~fin[K]
        K2 = K[nfk]
        c2 = cost[nfk]
        better = c2 < T[K2]
        if better.any():
            Kb = K2[better]
            T[Kb] = c2[better]
            Tm[Kb] = c2[better]
            if record_parent:
                for kk, jj in zip(Kb.tolist(), Jv[nfk][better].tolist()):
                    parent[kk] = (i, jj)
    return T, Q, parent


def window_points_box(gen, Qfun, W):
    """Utility: filter generator of integer tuples by Qfun <= W."""
    return [p for p in gen if Qfun(p) <= W]
