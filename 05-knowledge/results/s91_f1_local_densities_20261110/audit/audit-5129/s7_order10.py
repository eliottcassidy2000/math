# audit-5129, step 7: an independent rooted census of order 10 for the line counts of (F).
#
# Every tournament of order 10 is a one-vertex extension P + v of a class representative P of order 9 (191,536
# representatives from step 4, 512 out-sets each).  A candidate is kept only if the new vertex v has the least
# value of an isomorphism-invariant vertex key; every class is still produced (take a vertex u of least key in a
# tournament T: T - u is isomorphic to a representative P and T to P + v with v of least key).  The kept
# candidates are counted up to isomorphism by a canonical form (known: 9,733,056 tournaments, 9,355,949 strong).
# For every kept candidate and every root: D = H(B), U = U_B(y); the lattice points (D, U 2^j) of the lines
# D + 2U = N < 2^15 are marked (reading (a): sinks appended).
import sys, time
import numpy as np
from numba import njit, prange
import acore as A
from s2_classes import canon

NMAX = 1 << 15
NTH = 3


@njit(cache=True)
def vkeys(adj, n, key):
    fm = (1 << n) - 1
    score = np.zeros(n, np.int64)
    for v in range(n):
        score[v] = A.popcount(adj[v])
    for v in range(n):
        s1 = 0
        s2 = 0
        c3 = 0
        inv = fm & (~adj[v]) & (~(1 << v))
        for w in range(n):
            if (adj[v] >> w) & 1:
                s1 += score[w]
                s2 += score[w] * score[w]
                c3 += A.popcount(adj[w] & inv)
        key[v] = ((score[v] * 128 + s1) * 1024 + s2) * 64 + c3


@njit(cache=True)
def build(adj0, S, adj):
    for v in range(9):
        adj[v] = adj0[v]
        if not ((S >> v) & 1):
            adj[v] |= (1 << 9)
    adj[9] = S


@njit(cache=True)
def accepted(adj, key):
    vkeys(adj, 10, key)
    k = key[9]
    for v in range(9):
        if key[v] < k:
            return False
    return True


@njit(parallel=True, cache=True)
def pass1(codes9, cnt):
    R = codes9.shape[0]
    for r in prange(R):
        adj0 = np.zeros(9, np.int64)
        adj = np.zeros(10, np.int64)
        key = np.zeros(10, np.int64)
        A.adj_from_code(codes9[r], 9, adj0)
        c = 0
        for S in range(512):
            build(adj0, S, adj)
            if accepted(adj, key):
                c += 1
        cnt[r] = c


@njit(cache=True)
def census_one(adj, n, sp, Hs, LOG, Uout):
    """sp[m, v] = number of Hamiltonian paths of T[m] starting at v; Hs[m] = H(T[m]); Uout[y] = U(y)."""
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
    Hs[0] = 1
    for m in range(1, full):
        tot = 0
        mm = m
        while mm:
            lb = mm & (-mm)
            v = LOG[lb]
            mm ^= lb
            c = sp[m, v]
            if c:
                tot += c
                ext = inm[v] & (~m) & fm
                while ext:
                    lb2 = ext & (-ext)
                    sp[m | lb2, LOG[lb2]] += c
                    ext ^= lb2
        Hs[m] = tot
    for y in range(n):
        Uout[y] = 0
    for m in range(1, full):
        h2 = Hs[fm ^ m]
        mm = m
        while mm:
            lb = mm & (-mm)
            v = LOG[lb]
            mm ^= lb
            Uout[v] += sp[m, v] * h2
    return Hs[fm]


@njit(parallel=True, cache=True)
def pass2(codes9, off, r_lo, r_hi, ccodes, cstrong, cH, mark, LOG, maxDU):
    for r in prange(r_lo, r_hi):
        adj0 = np.zeros(9, np.int64)
        adj = np.zeros(10, np.int64)
        key = np.zeros(10, np.int64)
        sp = np.zeros((1024, 10), np.int64)
        Hs = np.zeros(1024, np.int64)
        Uout = np.zeros(10, np.int64)
        A.adj_from_code(codes9[r], 9, adj0)
        pos = off[r] - off[r_lo]
        for S in range(512):
            build(adj0, S, adj)
            if not accepted(adj, key):
                continue
            ccodes[pos] = canon(adj, 10, -1)
            cstrong[pos] = 1 if A.is_strong(adj, 10) else 0
            D = census_one(adj, 10, sp, Hs, LOG, Uout)
            cH[pos] = D
            pos += 1
            if 3 * D < NMAX:
                for y in range(10):
                    U = Uout[y]
                    N = D + 2 * U
                    while N < NMAX:
                        mark[(N - 1) // 2, (D - 1) // 2] = True
                        U *= 2
                        N = D + 2 * U
        # (per-thread maxima are not tracked; not needed)


if __name__ == "__main__":
    t0 = time.time()
    Z9 = np.load("classes9.npz")
    codes9 = Z9["codes"]
    R = codes9.size
    cnt = np.zeros(R, np.int64)
    pass1(codes9, cnt)
    off = np.zeros(R + 1, np.int64)
    off[1:] = np.cumsum(cnt)
    tot = int(off[-1])
    print("parents: %d ; candidates: %d ; kept (new vertex of least key): %d [%.1fs]" %
          (R, R * 512, tot, time.time() - t0), flush=True)
    nchunk = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    which = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    # chunk boundaries with about equal numbers of kept candidates
    bounds = [0]
    for c in range(1, nchunk):
        bounds.append(int(np.searchsorted(off, tot * c // nchunk)))
    bounds.append(R)
    r_lo, r_hi = bounds[which], bounds[which + 1]
    k = int(off[r_hi] - off[r_lo])
    LOG = np.zeros(1024, np.int64)
    for i in range(10):
        LOG[1 << i] = i
    ccodes = np.zeros(k, np.int64)
    cstrong = np.zeros(k, np.int8)
    cH = np.zeros(k, np.int64)
    mark = np.zeros((NMAX // 2, NMAX // 6 + 2), dtype=np.bool_)
    maxDU = np.zeros(2, np.int64)
    pass2(codes9, off, r_lo, r_hi, ccodes, cstrong, cH, mark, LOG, maxDU)
    print("chunk %d of %d: parents %d..%d, %d tournaments counted [%.1fs]" %
          (which, nchunk, r_lo, r_hi, k, time.time() - t0), flush=True)
    np.save("o10_codes_%d.npy" % which, ccodes)
    np.save("o10_strong_%d.npy" % which, cstrong)
    np.save("o10_H_%d.npy" % which, cH)
    np.save("o10_mark_%d.npy" % which, np.packbits(mark, axis=1))
    print("saved chunk %d [%.1fs]" % (which, time.time() - t0))
