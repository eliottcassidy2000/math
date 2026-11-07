"""How much of the deterministic (DC) count of THM-4890 (D) is forced by the converse involution?

For v in X_n write v = y 1^k with y in Y_(n-k) (or v = 1^n).  Its layer-converse partner is y* 1^k, with the same g.
With the split v = h 2^m + l (m = floor(n/2) - 2), a pair is of type 4 when h != h' and l != l'.
Prints, for n = 12..20: DC = all type-4 collisions; the ordered partner pairs (v, partner(v)) of type 4; the rest;
and the random-model type-4 expectation from aud_randmodel.out for comparison.
"""
import numpy as np


def gtable(N):
    g = np.zeros(1 << N, dtype=np.int64)
    g[0] = 1
    for v in range(1, 1 << N):
        u = v >> 1
        g[v] = g[u] if v & 1 else g[u] + u + (u & 1)
    return g


def rc(x, w):
    r = 0
    for i in range(w):
        r |= (1 - ((x >> i) & 1)) << (w - 1 - i)
    return r


N = 20
G = gtable(N)
for n in range(12, N + 1):
    m = n // 2 - 2
    C = 1 << m
    lo = 1 << (n - 1)
    vs = np.arange(lo, 1 << n, dtype=np.int64)
    gv = G[vs]
    order = np.argsort(gv, kind='stable')
    sv, sg = vs[order], gv[order]
    st = np.flatnonzero(np.r_[True, sg[1:] != sg[:-1]])
    en = np.r_[st[1:], len(sg)]
    dc = 0
    for i0, i1 in zip(st, en):
        if i1 - i0 < 2:
            continue
        grp = sv[i0:i1]
        hs = grp >> m; ls = grp & (C - 1)
        for a in range(len(grp)):
            dc += int(np.sum((hs != hs[a]) & (ls != ls[a])))
    # partner pairs
    part4 = 0
    for v in range(lo, 1 << n):
        if v == (1 << n) - 1:
            continue
        k = 0
        while (v >> k) & 1:
            k += 1
        y = v >> k
        j = n - k
        ys = rc(y, j)
        w = (ys << k) | ((1 << k) - 1)
        if w == v:
            continue
        assert G[w] == G[v]
        if (w >> m) != (v >> m) and (w & (C - 1)) != (v & (C - 1)):
            part4 += 1
    X = 1 << (n - 1)
    print('n=%d: DC (type-4 collisions) = %.4f |X|; of which converse-partner pairs %.4f |X|; other type-4 %.4f |X|'
          % (n, dc / X, part4 / X, (dc - part4) / X))
