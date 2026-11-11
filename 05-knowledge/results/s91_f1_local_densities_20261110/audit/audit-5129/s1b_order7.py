# audit-5129, step 1b: ALL 2^21 labelled tournaments of order 7, every root (14,680,064 rooted tournaments).
#   - (D, E, U_even, U_odd) by the subset programme (validated against literal counts in step 1 for orders <= 6
#     and here on a sample of order 7)
#   - E = U_odd - U_even and U >= D on every root
#   - product law (A1) on every non-strong tournament (component data by the subset programme on the induced
#     sub-tournaments)
#   - closing laws on a sample (explicit closed tournaments of orders 8, 9, 10)
#   - the set of pairs (D, U) of order 7
import time
import numpy as np
from numba import njit, prange
import acore as A

t0 = time.time()
N7 = 1 << 21


@njit(parallel=True, cache=True)
def run(keys, fails, counts):
    n = 7
    nthr = 12
    f_par = np.zeros(nthr, np.int64); f_ge = np.zeros(nthr, np.int64); f_prod = np.zeros(nthr, np.int64)
    f_close = np.zeros(nthr, np.int64)
    c_ns = np.zeros(nthr, np.int64); c_triv = np.zeros(nthr, np.int64); c_close = np.zeros(nthr, np.int64)
    c_strong = np.zeros(nthr, np.int64)
    chunk = N7 // nthr + 1
    for th in prange(nthr):
        adj = np.zeros(n, np.int64)
        hp = np.zeros((1 << n, n), np.int64); sp = np.zeros((1 << n, n), np.int64); Hs = np.zeros(1 << n, np.int64)
        st = np.zeros(4, np.int64)
        comp = np.zeros(n, np.int64)
        sub = np.zeros(n, np.int64)
        shp = np.zeros((1 << n, n), np.int64); ssp = np.zeros((1 << n, n), np.int64); sHs = np.zeros(1 << n, np.int64)
        cH = np.zeros(n, np.int64); cPi = np.zeros(n, np.int64); cU = np.zeros(n, np.int64); csz = np.zeros(n, np.int64)
        vs = np.zeros(n, np.int64)
        cadj = np.zeros(12, np.int64)
        dp = np.zeros((1 << 10, 10), np.int64)
        Uall = np.zeros(n, np.int64)
        for code in range(th * chunk, min(N7, (th + 1) * chunk)):
            A.adj_from_code(code, n, adj)
            A.tables(adj, n, hp, sp, Hs)
            for y in range(n):
                A.root_stats(n, hp, sp, Hs, y, st)
                U = st[2] + st[3]
                Uall[y] = U
                keys[code, y] = (st[0] << 32) | U
                if st[1] != st[3] - st[2]:
                    f_par[th] += 1
                if U < st[0]:
                    f_ge[th] += 1
            D = Hs[(1 << n) - 1]
            ncomp = A.components(adj, n, comp)
            if ncomp == 1:
                c_strong[th] += 1
            else:
                c_ns[th] += 1
                Dp = 1
                for c in range(ncomp):
                    k = 0
                    for v in range(n):
                        if comp[v] == c:
                            vs[k] = v
                            k += 1
                    csz[c] = k
                    for a in range(k):
                        sub[a] = 0
                        for b in range(k):
                            if (adj[vs[a]] >> vs[b]) & 1:
                                sub[a] |= (1 << b)
                    A.tables(sub, k, shp, ssp, sHs)
                    cH[c] = sHs[(1 << k) - 1]
                    cPi[c] = A.pi2_from_tables(k, sHs)
                    Dp *= cH[c]
                    for a in range(k):
                        A.root_stats(k, shp, ssp, sHs, a, st)
                        cU[vs[a]] = st[2] + st[3]
                if Dp != D:
                    f_prod[th] += 1
                for y in range(n):
                    k = comp[y]
                    pred = cU[y]
                    if csz[k] == 1:
                        c_triv[th] += 1
                    for i in range(ncomp):
                        if i < k:
                            pred *= cH[i]
                        elif i > k:
                            pred *= cPi[i]
                    if pred != Uall[y]:
                        f_prod[th] += 1
            if code % 61 == 0:
                for y in range(n):
                    A.root_stats(n, hp, sp, Hs, y, st)
                    U = st[2] + st[3]
                    A.head_adj(adj, n, y, cadj)
                    if A.ham_count(cadj, n + 1, dp) != st[0] - st[1] + U:
                        f_close[th] += 1
                    for t in range(1, 3):
                        A.tail_adj(adj, n, y, t, cadj)
                        if A.ham_count(cadj, n + t + 1, dp) != st[0] + (1 << t) * U:
                            f_close[th] += 1
                    c_close[th] += 1
    fails[0] = f_par.sum(); fails[1] = f_ge.sum(); fails[2] = f_prod.sum(); fails[3] = f_close.sum()
    counts[0] = c_ns.sum(); counts[1] = c_triv.sum(); counts[2] = c_close.sum(); counts[3] = c_strong.sum()


if __name__ == "__main__":
    keys = np.zeros((N7, 7), np.int64)
    fails = np.zeros(4, np.int64); counts = np.zeros(4, np.int64)
    run(keys, fails, counts)
    print("order 7: %d labelled tournaments, %d rooted" % (N7, N7 * 7))
    print("strong labelled: %d ; non-strong labelled: %d ; roots in a trivial component (non-strong): %d" %
          (counts[3], counts[0], counts[1]))
    print("failures: E = U_odd - U_even: %d ; U >= D: %d ; product law (D or U): %d ; closing laws: %d "
          "(on %d sampled rooted tournaments, head and t = 1, 2)" % (fails[0], fails[1], fails[2], fails[3], counts[2]))
    uk = np.unique(keys.ravel())
    np.save("pairs7.npy", uk)
    print("distinct pairs (D, U) at order 7: %d ; offsets: %d ; max D %d ; max U %d" %
          (uk.size, np.unique(uk >> 32).size, int((uk >> 32).max()), int((uk & 0xFFFFFFFF).max())))
    # literal sample
    perms = A.perm_array(7)
    E = np.zeros(7, np.int64); Ue = np.zeros(7, np.int64); Uo = np.zeros(7, np.int64)
    adj = np.zeros(7, np.int64)
    bad = 0
    ns = 0
    for code in range(0, N7, 1009):
        A.adj_from_code(code, 7, adj)
        H = A.literal_stats(adj, 7, perms, E, Ue, Uo)
        ns += 1
        for y in range(7):
            if keys[code, y] != ((H << 32) | int(Ue[y] + Uo[y])):
                bad += 1
    print("literal recount of (D, U) on %d sampled labelled tournaments of order 7 (all roots): %d mismatches" % (ns, bad))
    print("elapsed %.1fs" % (time.time() - t0))
