# audit-5129, step 4: every tournament of order 9 as a one-vertex extension of a class representative of order 8
# (6880 * 256 labelled tournaments, covering all isomorphism classes).  For each: H, strong, Pi_2, canonical code,
# and (D, U) at every root.  Gives f(9), the class counts (known: 191536 tournaments, 178133 strong), the rooted
# census of order 9 and the strong tournaments of order 9 with few paths.
import time
import numpy as np
from numba import njit, prange
import acore as A
from s2_classes import canon

t0 = time.time()


@njit(parallel=True, cache=True)
def run(codes8, Hout, strong, Pi2, canon9, keys):
    n = 9
    R = codes8.shape[0]
    for r in prange(R):
        adj0 = np.zeros(8, np.int64)
        adj = np.zeros(n, np.int64)
        hp = np.zeros((1 << n, n), np.int64)
        sp = np.zeros((1 << n, n), np.int64)
        Hs = np.zeros(1 << n, np.int64)
        st = np.zeros(4, np.int64)
        A.adj_from_code(codes8[r], 8, adj0)
        for S in range(256):
            for v in range(8):
                adj[v] = adj0[v]
                if not ((S >> v) & 1):
                    adj[v] |= (1 << 8)
            adj[8] = S
            A.tables(adj, n, hp, sp, Hs)
            Hout[r, S] = Hs[(1 << n) - 1]
            strong[r, S] = 1 if A.is_strong(adj, n) else 0
            Pi2[r, S] = A.pi2_from_tables(n, Hs)
            canon9[r, S] = canon(adj, n, -1)
            for y in range(n):
                A.root_stats(n, hp, sp, Hs, y, st)
                keys[r, S, y] = (st[0] << 32) | (st[2] + st[3])


if __name__ == "__main__":
    Z = np.load("classes_le8.npz")
    codes8 = Z["codes8"]
    R = codes8.size
    Hout = np.zeros((R, 256), np.int64)
    strong = np.zeros((R, 256), np.int8)
    Pi2 = np.zeros((R, 256), np.int64)
    canon9 = np.zeros((R, 256), np.int64)
    keys = np.zeros((R, 256, 9), np.int64)
    run(codes8, Hout, strong, Pi2, canon9, keys)
    print("order-9 extension run finished [%.1fs]" % (time.time() - t0), flush=True)
    ucl, first = np.unique(canon9.ravel(), return_index=True)
    print("labelled extensions: %d ; isomorphism classes: %d (known 191536)" % (canon9.size, ucl.size))
    sflat = strong.ravel()[first] == 1
    print("strong classes: %d (known 178133)" % int(sflat.sum()))
    Hcl = Hout.ravel()[first]
    f9 = int(Hcl[sflat].min())
    print("f(9) = min H over strong tournaments of order 9 = %d ; classes attaining it: %d" %
          (f9, int(np.sum(Hcl[sflat] == f9))))
    print("least H over ALL strong labelled extensions (no canonical form used): %d" %
          int(Hout[strong == 1].min()))
    print("strong classes of order 9 with H <= 73: %d ; with H <= 123: %d ; max H at order 9: %d" %
          (int(np.sum(Hcl[sflat] <= 73)), int(np.sum(Hcl[sflat] <= 123)), int(Hcl.max())))
    # consistency: a strong class of order 9 and its H
    hs = np.sort(Hcl[sflat])
    print("ten least H of strong classes of order 9:", hs[:10].tolist())
    uk = np.unique(keys.ravel())
    np.save("pairs9.npy", uk)
    print("distinct pairs (D, U) at order 9: %d ; with D <= 73: %d" %
          (uk.size, int(np.sum((uk >> 32) <= 73))))
    # strong order-9 data with small H, one row per class
    sel = np.nonzero(sflat & (Hcl <= 225))[0]
    rr, ss = np.unravel_index(first[sel], canon9.shape)
    rows = np.zeros((sel.size, 11), np.int64)
    rows[:, 0] = Hout[rr, ss]
    rows[:, 1] = Pi2[rr, ss]
    rows[:, 2:] = keys[rr, ss, :] & 0xFFFFFFFF
    np.save("strong9_small.npy", rows)
    print("saved strong9_small.npy (%d strong classes of order 9 with H <= 225)" % sel.size)
    # all classes of order 9 (codes, H, strong) for later use
    np.savez_compressed("classes9.npz", codes=ucl, H=Hcl, strong=sflat.astype(np.int8),
                        rep=first)
    # level 7 and 9 sources: rooted classes with the given pair at order 9
    for (D, U) in ((129, 144), (513, 576), (129, 288), (257, 288), (65, 72), (2049, 2304)):
        key = (D << 32) | U
        hit = np.argwhere(keys == key)
        rc = set()
        adj0 = np.zeros(8, np.int64); adj = np.zeros(9, np.int64)
        for (r, S, y) in hit:
            A.adj_from_code(codes8[r], 8, adj0)
            for v in range(8):
                adj[v] = adj0[v]
                if not ((S >> v) & 1):
                    adj[v] |= (1 << 8)
            adj[8] = S
            rc.add(int(canon(adj, 9, int(y))))
        print("order 9, pair (%d, %d): %d labelled hits, %d rooted isomorphism classes" % (D, U, len(hit), len(rc)))
    print("elapsed %.1fs" % (time.time() - t0))
