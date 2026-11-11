# audit-5129, step 10: the section (B') added to THM-5129 after the audit was launched (offsets 129 and 257), and
# more generally the complete slope lists for EVERY offset D <= 257.
#
# Needed: every strong tournament with H <= 257.  Orders <= 8: class list (step 2).  Order 9: class list (step 4).
# Order 10: class list (step 7).  Order 11: a strong tournament T of order n >= 4 has a vertex v with T - v strong
# and H(T - v) <= H(T) - 2 (THM-5129 (A2)), so the strong tournaments of order 11 with H <= 257 are one-vertex
# extensions of the strong classes of order 10 with H <= 255; same for order 12.
import time
import numpy as np
from numba import njit, prange
import acore as A
from s2_classes import canon

t0 = time.time()
HMAX = 257


@njit(cache=True)
def H_only(adj, n, sp, LOG):
    full = 1 << n
    fm = full - 1
    inm = np.zeros(n, np.int64)
    for v in range(n):
        for w in range(n):
            if (adj[v] >> w) & 1:
                inm[w] |= (1 << v)
    for m in range(full):
        for v in range(n):
            sp[m, v] = 0
    for v in range(n):
        sp[1 << v, v] = 1
    for m in range(1, full):
        mm = m
        while mm:
            lb = mm & (-mm)
            v = LOG[lb]
            mm ^= lb
            c = sp[m, v]
            if c:
                ext = inm[v] & (~m) & fm
                while ext:
                    lb2 = ext & (-ext)
                    sp[m | lb2, LOG[lb2]] += c
                    ext ^= lb2
    tot = 0
    for v in range(n):
        tot += sp[fm, v]
    return tot


@njit(parallel=True, cache=True)
def extend(codes, n_old, hmax, LOG, outH, outStrong, outCanon):
    """All 2^n_old one-vertex extensions of each parent: H (0 if > hmax is not needed: store H), strong, canon."""
    n = n_old + 1
    R = codes.shape[0]
    P = 1 << n_old
    for r in prange(R):
        adj0 = np.zeros(n_old, np.int64)
        adj = np.zeros(n, np.int64)
        sp = np.zeros((1 << n, n), np.int64)
        A.adj_from_code(codes[r], n_old, adj0)
        for S in range(P):
            for v in range(n_old):
                adj[v] = adj0[v]
                if not ((S >> v) & 1):
                    adj[v] |= (1 << n_old)
            adj[n_old] = S
            h = H_only(adj, n, sp, LOG)
            outH[r, S] = h
            if h <= hmax:
                st = A.is_strong(adj, n)
                outStrong[r, S] = 1 if st else 0
                if st:
                    outCanon[r, S] = canon(adj, n, -1)


@njit(cache=True)
def stats(codes, n, Hout, Pi2, Uarr):
    adj = np.zeros(n, np.int64)
    hp = np.zeros((1 << n, n), np.int64)
    sp = np.zeros((1 << n, n), np.int64)
    Hs = np.zeros(1 << n, np.int64)
    st = np.zeros(4, np.int64)
    for r in range(codes.shape[0]):
        A.adj_from_code(codes[r], n, adj)
        A.tables(adj, n, hp, sp, Hs)
        Hout[r] = Hs[(1 << n) - 1]
        Pi2[r] = A.pi2_from_tables(n, Hs)
        for y in range(n):
            A.root_stats(n, hp, sp, Hs, y, st)
            Uarr[r, y] = st[2] + st[3]


if __name__ == "__main__":
    LOG = np.zeros(1 << 12, np.int64)
    for i in range(12):
        LOG[1 << i] = i
    rows = []          # (order, H, Pi2, [U...])
    Z = np.load("classes_le8.npz")
    for n in range(3, 9):
        H = Z["H%d" % n]; sg = Z["strong%d" % n]; P = Z["Pi2%d" % n]; U = Z["U%d" % n]
        for r in range(H.size):
            if sg[r] and H[r] <= HMAX:
                rows.append((n, int(H[r]), int(P[r]), [int(x) for x in U[r]]))
    n8 = len(rows)
    # order 9
    Z9 = np.load("classes9.npz")
    sel = (Z9["strong"] == 1) & (Z9["H"] <= HMAX)
    c9 = Z9["codes"][sel]
    H9 = np.zeros(c9.size, np.int64); P9 = np.zeros(c9.size, np.int64); U9 = np.zeros((c9.size, 9), np.int64)
    stats(c9, 9, H9, P9, U9)
    assert np.array_equal(H9, Z9["H"][sel])
    for r in range(c9.size):
        rows.append((9, int(H9[r]), int(P9[r]), [int(x) for x in U9[r]]))
    print("strong classes with H <= %d: orders <= 8: %d ; order 9: %d" % (HMAX, n8, c9.size))
    # order 10
    codes10 = np.load("o10_codes_0.npy"); st10 = np.load("o10_strong_0.npy"); Hh10 = np.load("o10_H_0.npy")
    u10, first = np.unique(codes10, return_index=True)
    sel = (st10[first] == 1) & (Hh10[first] <= HMAX)
    c10 = u10[sel]
    H10 = np.zeros(c10.size, np.int64); P10 = np.zeros(c10.size, np.int64); U10 = np.zeros((c10.size, 10), np.int64)
    stats(c10, 10, H10, P10, U10)
    assert np.array_equal(H10, Hh10[first][sel])
    for r in range(c10.size):
        rows.append((10, int(H10[r]), int(P10[r]), [int(x) for x in U10[r]]))
    print("order 10: strong classes with H <= %d: %d (least H %d); with H <= 255: %d  [%.1fs]" %
          (HMAX, c10.size, int(H10.min()), int(np.sum(H10 <= 255)), time.time() - t0), flush=True)
    # order 11
    par = c10[H10 <= HMAX - 2]
    oH = np.zeros((par.size, 1 << 10), np.int64); oS = np.zeros((par.size, 1 << 10), np.int8)
    oC = np.zeros((par.size, 1 << 10), np.int64)
    extend(par, 10, HMAX, LOG, oH, oS, oC)
    m11 = (oH <= HMAX) & (oS == 1)
    c11 = np.unique(oC[m11])
    H11 = np.zeros(c11.size, np.int64); P11 = np.zeros(c11.size, np.int64); U11 = np.zeros((c11.size, 11), np.int64)
    stats(c11, 11, H11, P11, U11)
    for r in range(c11.size):
        rows.append((11, int(H11[r]), int(P11[r]), [int(x) for x in U11[r]]))
    print("order 11: %d parents, %d labelled extensions with H <= %d, %d of them strong, %d strong classes; "
          "least H of a strong one: %d (f(11)); strong with H = 129: %d, with H = 257: %d  [%.1fs]" %
          (par.size, int(np.sum(oH <= HMAX)), HMAX, int(m11.sum()), c11.size, int(H11.min()),
           int(np.sum(H11 == 129)), int(np.sum(H11 == 257)), time.time() - t0), flush=True)
    # order 12
    par = c11[H11 <= HMAX - 2]
    oH = np.zeros((par.size, 1 << 11), np.int64); oS = np.zeros((par.size, 1 << 11), np.int8)
    oC = np.zeros((par.size, 1 << 11), np.int64)
    extend(par, 11, HMAX, LOG, oH, oS, oC)
    m12 = (oH <= HMAX) & (oS == 1)
    print("order 12: %d parents (strong classes of order 11 with H <= 255), %d labelled extensions with H <= %d, "
          "strong among them: %d ; least H over all extensions with a non-trivial cut: %d  [%.1fs]" %
          (par.size, int(np.sum(oH <= HMAX)), HMAX, int(m12.sum()),
           int(oH[:, 1:-1].min()) if par.size else -1, time.time() - t0), flush=True)
    for h in (3, 43, 129, 257):
        byo = {}
        us = set(); ps = set()
        for (o, hh, p2, Ul) in rows:
            if hh == h:
                byo[o] = byo.get(o, 0) + 1
                us.update(Ul); ps.add(p2)
        print("strong classes with H = %d by order: %s ; %d values of U, %d values of Pi_2" % (h, byo, len(us), len(ps)))
    import pickle
    with open("strong_le257.pkl", "wb") as f:
        pickle.dump(rows, f)
    print("saved strong_le257.pkl (%d strong classes with H <= %d)  [%.1fs]" % (len(rows), HMAX, time.time() - t0))
