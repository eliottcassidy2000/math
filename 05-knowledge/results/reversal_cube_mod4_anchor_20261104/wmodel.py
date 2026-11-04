"""Fast weighted odd-cycle parity model (same output as wmodel.gamma)."""
import numpy as np

def path_tables(Z, w):
    n = len(Z); N = 1 << n
    out = [sum(1 << j for j in range(n) if Z[i][j]) for i in range(n)]
    T = np.zeros((N, n, n), dtype=np.uint8)
    for v in range(n): T[1 << v, v, v] = 1
    for S in range(1, N):
        blk = T[S]
        for e in range(n):
            if not (S >> e) & 1: continue
            col = blk[:, e]
            if not col.any(): continue
            nx = out[e] & ~S
            while nx:
                f = (nx & -nx).bit_length() - 1; nx &= nx - 1
                T[S | (1 << f), :, f] ^= col
    W = np.array([bin(S & sum(1 << i for i in range(n) if w[i])).count('1') & 1 for S in range(N)], dtype=np.uint8)
    return T, W

def zeta_cols(F, n):
    G = F.copy()
    idx = np.arange(G.shape[0])
    for i in range(n):
        bit = 1 << i
        sel = (idx & bit) != 0
        G[sel] ^= G[idx[sel] ^ bit]
    return G

def gamma(Z, w):
    n = len(Z); N = 1 << n
    T, W = path_tables(Z, w)
    cnt = {}
    def add(A, C):
        k = (frozenset(A), frozenset(C)); cnt[k] = cnt.get(k, 0) ^ 1
    P0 = (T[W == 0].sum(axis=0) % 2).astype(np.uint8); P1 = (T[W == 1].sum(axis=0) % 2).astype(np.uint8)
    for s in range(n):
        for e in range(n):
            if P0[s, e]:
                add([e], []); add([s, e], []); add([], [e]); add([], [s, e])
            if P1[s, e]:
                add([e], []); add([e], [s])
    flat = T.reshape(N, n * n)
    comp = (N - 1) ^ np.arange(N)
    for par2 in (0, 1):
        A = flat * (W == 1 - par2)[:, None]
        Bz = zeta_cols(flat * (W == par2)[:, None], n)[comp]
        M = (A.astype(np.int64).T @ Bz.astype(np.int64)) % 2
        for i1, i2 in zip(*np.nonzero(M)):
            s, e = divmod(int(i1), n); s2, e2 = divmod(int(i2), n)
            for A0 in ([], [s]):
                for C0 in ([], [s2]):
                    add(A0 + [e2], C0 + [e])
    return {k: 1 for k, v in cnt.items() if v}

def contrast(p, q):
    return {k: 1 for k in set(p) | set(q) if p.get(k, 0) ^ q.get(k, 0)}

def Lfun(d, arc, v, verts):
    s = 0
    for u in verts:
        if u == v: continue
        s ^= d.get((frozenset([u]), frozenset([v])), 0)
        if arc(v, u): s ^= d.get((frozenset([u, v]), frozenset([v])), 0)
    return s
