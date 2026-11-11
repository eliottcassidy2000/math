#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9666 (A). The fusion states of Session264 (fusion/PROOF.md):
#   pointed T with root z:  H = H(T), h = H(T - z) (H(empty) = 1),
#   U = sum over paths R starting at z of H(T - V(R)),  V = same for paths ending at z,
#   d = sum over R starting at z of (-1)^(|R|+1) H(T - V(R)),  e = same for paths ending at z   (|R| = number of vertices).
#   fusion law (3):  H(A * B) = H h' + h H' - h h' + (U - h)(V' - h') + (d - h)(e' - h').
# Checks (brute force, own code):
#   (i)   the law (3) on all pairs of pointed tournaments of orders <= 4 x <= 4 (A * B: roots identified, all arcs
#         from A - z to B - z);
#   (ii)  parities: H, h odd; U = d and V = e (mod 2)  -> the states (h, H-h, U-h, d-h) are confined modulo 2;
#   (iii) the distribution of the states (h, H-h, U-h, d-h) modulo 3 and 5 over ALL labelled pointed tournaments of
#         order n (root = vertex 0), n <= 7: how many of the p^4 residue vectors occur, the largest and smallest
#         p^4 * share, and any linear relation modulo p.
import sys, itertools
import numpy as np
from numba import njit


@njit(cache=True)
def invariants(adj, radj, n, Hs, out):
    # adj[v]: out-neighbour mask; radj[v]: in-neighbour mask. root = vertex 0.
    full = 1 << n
    # H of all induced subtournaments
    dp = np.zeros((full, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, full):
        for v in range(n):
            c = dp[S, v]
            if c:
                o = adj[v] & ~S
                w = 0
                while o:
                    if o & 1:
                        dp[S | (1 << w), w] += c
                    o >>= 1
                    w += 1
    Hs[0] = 1
    for S in range(1, full):
        t = 0
        for v in range(n):
            t += dp[S, v]
        Hs[S] = t
    # paths starting at 0 (gp) and paths ending at 0 (gm, = paths starting at 0 in the converse)
    gp = np.zeros((full, n), dtype=np.int64)
    gm = np.zeros((full, n), dtype=np.int64)
    gp[1, 0] = 1
    gm[1, 0] = 1
    U = 0
    V = 0
    d = 0
    e = 0
    for S in range(1, full):
        if not (S & 1):
            continue
        fp = 0
        fm = 0
        for v in range(n):
            c = gp[S, v]
            if c:
                fp += c
                o = adj[v] & ~S
                w = 0
                while o:
                    if o & 1:
                        gp[S | (1 << w), w] += c
                    o >>= 1
                    w += 1
            c = gm[S, v]
            if c:
                fm += c
                o = radj[v] & ~S
                w = 0
                while o:
                    if o & 1:
                        gm[S | (1 << w), w] += c
                    o >>= 1
                    w += 1
        comp = Hs[(full - 1) ^ S]
        # popcount parity of S
        k = 0
        t = S
        while t:
            k += t & 1
            t >>= 1
        sg = 1 if k % 2 == 1 else -1
        U += fp * comp
        V += fm * comp
        d += sg * fp * comp
        e += sg * fm * comp
    out[0] = Hs[full - 1]
    out[1] = Hs[(full - 1) ^ 1]
    out[2] = U
    out[3] = V
    out[4] = d
    out[5] = e


@njit(cache=True)
def all_states(n, res):
    m = n * (n - 1) // 2
    full = 1 << n
    Hs = np.zeros(full, dtype=np.int64)
    adj = np.zeros(n, dtype=np.int64)
    radj = np.zeros(n, dtype=np.int64)
    out = np.zeros(6, dtype=np.int64)
    for code in range(1 << m):
        for v in range(n):
            adj[v] = 0
            radj[v] = 0
        k = 0
        for i in range(n):
            for j in range(i + 1, n):
                if (code >> k) & 1:
                    adj[i] |= 1 << j
                    radj[j] |= 1 << i
                else:
                    adj[j] |= 1 << i
                    radj[i] |= 1 << j
                k += 1
        invariants(adj, radj, n, Hs, out)
        for t in range(6):
            res[code, t] = out[t]


def tour_from_code(n, code):
    adj = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> k) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj


def fuse(A, B):
    # roots (vertex 0 of each) identified; all arcs from A - 0 to B - 0
    na, nb = len(A), len(B)
    N = na + nb - 1
    adj = [[0] * N for _ in range(N)]
    for i in range(na):
        for j in range(na):
            adj[i][j] = A[i][j]
    def idx(j):
        return 0 if j == 0 else na + j - 1
    for i in range(nb):
        for j in range(nb):
            if i != j:
                adj[idx(i)][idx(j)] = B[i][j]
    for i in range(1, na):
        for j in range(1, nb):
            adj[i][na + j - 1] = 1
            adj[na + j - 1][i] = 0
    return adj


def ham(adj):
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = dp[S][v]
            if c:
                for w in range(n):
                    if adj[v][w] and not (S >> w) & 1:
                        dp[S | (1 << w)][w] += c
    return sum(dp[(1 << n) - 1])


def main():
    print("# fusion_states.py (auditor hypotheses; own code)")
    states = {}
    for n in range(1, 8):
        m = n * (n - 1) // 2
        res = np.zeros((1 << m, 6), dtype=np.int64)
        all_states(n, res)
        states[n] = res
    # sanity: singleton and triangle (Session264: (H,h,Q,U,V,d,e) = (3,1,4,3,3,1,1) for the directed triangle)
    print("singleton: (H, h, U, V, d, e) = %s" % states[1][0].tolist())
    tri = [states[3][c].tolist() for c in range(8)]
    print("order 3, the 8 labelled tournaments rooted at 0: %s" % sorted(set(tuple(t) for t in tri)))
    # (i) fusion law
    tot = 0
    bad = 0
    for na in range(1, 5):
        for nb in range(1, 5):
            for ca in range(1 << (na * (na - 1) // 2)):
                A = tour_from_code(na, ca)
                Ha, ha, Ua, Va, da, ea = states[na][ca].tolist()
                for cb in range(1 << (nb * (nb - 1) // 2)):
                    B = tour_from_code(nb, cb)
                    Hb, hb, Ub, Vb, db, eb = states[nb][cb].tolist()
                    pred = Ha * hb + ha * Hb - ha * hb + (Ua - ha) * (Vb - hb) + (da - ha) * (eb - hb)
                    tot += 1
                    if ham(fuse(A, B)) != pred:
                        bad += 1
    print("(i) fusion law (3) on %d pairs of pointed tournaments of orders <= 4: %d violations" % (tot, bad))
    # (ii) parities
    for n in range(1, 8):
        r = states[n]
        H, h, U, V, d, e = (r[:, i] for i in range(6))
        print("(ii) n = %d: H odd %s, h odd %s, U = d mod 2 %s, V = e mod 2 %s; (U - h) mod 2 takes values %s; ((U-h) mod 4, (d-h) mod 4) pairs: %s" % (
            n, bool(np.all(H % 2 == 1)), bool(np.all(h % 2 == 1)), bool(np.all((U - d) % 2 == 0)), bool(np.all((V - e) % 2 == 0)),
            sorted(set(((U - h) % 2).tolist())), sorted(set(zip(((U - h) % 4).tolist(), ((d - h) % 4).tolist())))))
    # (iii) distribution modulo 3, 5, 7
    for p in (3, 5, 7):
        for n in range(4, 8):
            r = states[n]
            H, h, U, V, d, e = (r[:, i] for i in range(6))
            for name, vec in (("left  (h, H-h, U-h, d-h)", (h, H - h, U - h, d - h)), ("right (h, H-h, V-h, e-h)", (h, H - h, V - h, e - h))):
                key = ((vec[0] % p) * p + (vec[1] % p)) * p * p + (vec[2] % p) * p + (vec[3] % p)
                cnt = np.bincount(key, minlength=p ** 4)
                occ = int((cnt > 0).sum())
                share = cnt / cnt.sum()
                # linear relations mod p: brute force over nonzero coefficient vectors (c0..c3, const)
                rels = []
                if occ < p ** 4:
                    vals = [v % p for v in vec]
                    for c in itertools.product(range(p), repeat=4):
                        if not any(c):
                            continue
                        # normalise: first nonzero coefficient 1
                        fz = next(x for x in c if x)
                        if fz != 1:
                            continue
                        lin = (c[0] * vals[0] + c[1] * vals[1] + c[2] * vals[2] + c[3] * vals[3]) % p
                        u = np.unique(lin)
                        if len(u) == 1:
                            rels.append((c, int(u[0])))
                print("(iii) p = %d, n = %d, %s: %d of %d residue vectors occur; p^4 * share in [%.3f, %.3f]; zero vector: %.3f; affine relations mod p: %s" % (
                    p, n, name, occ, p ** 4, p ** 4 * share[cnt > 0].min(), p ** 4 * share.max(), p ** 4 * share[0], rels if rels else "none"))


if __name__ == "__main__":
    main()
