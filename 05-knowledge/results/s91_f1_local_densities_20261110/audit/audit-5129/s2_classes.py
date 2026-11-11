# audit-5129, step 2: isomorphism classes of tournaments of orders 1..8 by one-vertex extension and a canonical
# form (minimum code over all relabellings that sort the vertices by an isomorphism-invariant key), then the
# statistics (H, strong, Pi_2, U(y), E(y), U_even(y)) of every class and every root.
# Known counts used as a check of the generator (OEIS A000568, A051337):
#   tournaments 1, 1, 2, 4, 12, 56, 456, 6880;  strong 1, 0, 1, 1, 6, 35, 353, 6008.
import time
import numpy as np
from numba import njit
import acore as A

t0 = time.time()


@njit(cache=True)
def next_perm(a, lo, hi):
    i = hi - 2
    while i >= lo and a[i] >= a[i + 1]:
        i -= 1
    if i < lo:
        l = lo
        r = hi - 1
        while l < r:
            tmp = a[l]; a[l] = a[r]; a[r] = tmp
            l += 1
            r -= 1
        return False
    j = hi - 1
    while a[j] <= a[i]:
        j -= 1
    tmp = a[i]; a[i] = a[j]; a[j] = tmp
    l = i + 1
    r = hi - 1
    while l < r:
        tmp = a[l]; a[l] = a[r]; a[r] = tmp
        l += 1
        r -= 1
    return True


@njit(cache=True)
def canon(adj, n, root):
    """Canonical code of the tournament (rooted at `root` if root >= 0)."""
    fm = (1 << n) - 1
    score = np.zeros(n, np.int64)
    for v in range(n):
        score[v] = A.popcount(adj[v])
    key = np.zeros(n, np.int64)
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
        if v == root:
            key[v] += (1 << 40)
    order = np.arange(n)
    # insertion sort by (key, index)
    for i in range(1, n):
        x = order[i]
        j = i - 1
        while j >= 0 and key[order[j]] > key[x]:
            order[j + 1] = order[j]
            j -= 1
        order[j + 1] = x
    gs = np.zeros(n + 1, np.int64)
    ng = 0
    gs[0] = 0
    for i in range(1, n):
        if key[order[i]] != key[order[i - 1]]:
            ng += 1
            gs[ng] = i
    ng += 1
    gs[ng] = n
    best = np.int64(1) << 62
    while True:
        code = 0
        k = 0
        for i in range(n):
            ai = adj[order[i]]
            for j in range(i + 1, n):
                if (ai >> order[j]) & 1:
                    code |= (np.int64(1) << k)
                k += 1
        if code < best:
            best = code
        g = 0
        while g < ng:
            if next_perm(order, gs[g], gs[g + 1]):
                break
            g += 1
        if g == ng:
            break
    return best


@njit(cache=True)
def extend_all(rep_codes, n_old):
    n = n_old + 1
    out = np.empty(rep_codes.shape[0] << n_old, np.int64)
    adj = np.zeros(n, np.int64)
    adj0 = np.zeros(max(n_old, 1), np.int64)
    idx = 0
    for r in range(rep_codes.shape[0]):
        A.adj_from_code(rep_codes[r], n_old, adj0)
        for S in range(1 << n_old):
            for v in range(n_old):
                adj[v] = adj0[v]
                if not ((S >> v) & 1):
                    adj[v] |= (1 << n_old)
            adj[n_old] = S
            out[idx] = canon(adj, n, -1)
            idx += 1
    return out


@njit(cache=True)
def class_stats(codes, n, Hout, strong, Pi2, Uarr, Earr, Uearr):
    adj = np.zeros(n, np.int64)
    hp = np.zeros((1 << n, n), np.int64)
    sp = np.zeros((1 << n, n), np.int64)
    Hs = np.zeros(1 << n, np.int64)
    st = np.zeros(4, np.int64)
    for r in range(codes.shape[0]):
        A.adj_from_code(codes[r], n, adj)
        A.tables(adj, n, hp, sp, Hs)
        Hout[r] = Hs[(1 << n) - 1]
        strong[r] = 1 if A.is_strong(adj, n) else 0
        Pi2[r] = A.pi2_from_tables(n, Hs)
        for y in range(n):
            A.root_stats(n, hp, sp, Hs, y, st)
            Earr[r, y] = st[1]
            Uearr[r, y] = st[2]
            Uarr[r, y] = st[2] + st[3]


if __name__ == "__main__":
    known_all = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
    known_strong = {1: 1, 2: 0, 3: 1, 4: 1, 5: 6, 6: 35, 7: 353, 8: 6008}
    reps = {1: np.array([0], dtype=np.int64)}
    for n in range(2, 9):
        c = extend_all(reps[n - 1], n - 1)
        reps[n] = np.unique(c)
        print("order %d: %d candidates -> %d classes (known %d)  [%.1fs]"
              % (n, c.size, reps[n].size, known_all[n], time.time() - t0), flush=True)
        assert reps[n].size == known_all[n]
    save = {}
    nroots = 0
    viol = 0
    for n in range(1, 9):
        codes = reps[n]
        k = codes.size
        H = np.zeros(k, np.int64); strong = np.zeros(k, np.int64); Pi2 = np.zeros(k, np.int64)
        U = np.zeros((k, n), np.int64); E = np.zeros((k, n), np.int64); Ue = np.zeros((k, n), np.int64)
        class_stats(codes, n, H, strong, Pi2, U, E, Ue)
        nroots += k * n
        viol += int(np.sum(E != (U - Ue) - Ue)) + int(np.sum(U < H[:, None]))
        assert np.all(Pi2 == U.sum(axis=1) + H)
        ns = int(strong.sum())
        s = strong == 1
        hu = set()
        for r in np.nonzero(s)[0]:
            for y in range(n):
                hu.add((int(H[r]), int(U[r, y])))
        hp2 = {(int(a), int(b)) for a, b in zip(H[s], Pi2[s])}
        fmin = int(H[s].min()) if ns else None
        print("order %d: classes %d, strong %d (known %d), f(n) = %s, #strong with H = f: %d, distinct (H,U) %d, "
              "distinct (H,Pi_2) %d, max H %d" % (n, k, ns, known_strong[n], fmin,
                                                  int(np.sum(H[s] == fmin)) if ns else 0, len(hu), len(hp2),
                                                  int(H.max())), flush=True)
        assert ns == known_strong[n]
        save["codes%d" % n] = codes; save["H%d" % n] = H; save["strong%d" % n] = strong; save["Pi2%d" % n] = Pi2
        save["U%d" % n] = U; save["E%d" % n] = E; save["Ue%d" % n] = Ue
    print("rooted (class, vertex) pairs of orders 1..8:", nroots)
    print("violations of E = U_odd - U_even or U >= D over these roots:", viol)
    np.savez_compressed("classes_le8.npz", **save)
    print("saved classes_le8.npz  [%.1fs]" % (time.time() - t0))
