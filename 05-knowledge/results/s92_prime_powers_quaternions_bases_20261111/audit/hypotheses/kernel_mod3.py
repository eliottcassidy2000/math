#!/usr/bin/env python3
# Auditor "hypotheses": HYP-9666 (B), (D). Generalisation of (D) suggested by the audit:
# for EVERY two-port kernel K (any number of visits), modulo 3
#     H(K[A, B]) = c_K(|A| mod 2, |B| mod 2) * H(A) * H(B)       (mod 3),
# because P_2 = ((-1)^n + 1) H and P_m = 0 (mod 3) for m >= 3 (THM-4756) in the expansion
# H(K[A, B]) = sum N_ij P_i(A) P_j(B) of THM-5133 (A).  So module kernels of every rank collapse to rank 1 at 3.
# Brute force: every labelled K of order 3, 4 and 5, every ordered port pair, A and B over one representative of
# every isomorphism class of orders 1..4 (8 tournaments).  For each (K, ports, parity of |A|, parity of |B|) the
# residue H(K[A,B]) must be a fixed multiple of H(A) H(B) modulo 3.
import itertools, sys
import numpy as np
from numba import njit


@njit(cache=True)
def ham(adjm, n):
    full = 1 << n
    dp = np.zeros((full, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v, v] = 1
    for S in range(1, full):
        for v in range(n):
            c = dp[S, v]
            if c:
                o = adjm[v] & ~S
                w = 0
                while o:
                    if o & 1:
                        dp[S | (1 << w), w] += c
                    o >>= 1
                    w += 1
    t = 0
    for v in range(n):
        t += dp[full - 1, v]
    return t


def tour(n, code):
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


def masks(adj):
    n = len(adj)
    return np.array([sum(1 << j for j in range(n) if adj[i][j]) for i in range(n)], dtype=np.int64)


def substitute(K, x, y, A, B):
    m = len(K)
    verts = []
    for v in range(m):
        if v == x:
            verts += [(v, "A", i) for i in range(len(A))]
        elif v == y:
            verts += [(v, "B", i) for i in range(len(B))]
        else:
            verts.append((v, "v", 0))
    N = len(verts)
    adj = [[0] * N for _ in range(N)]
    for a in range(N):
        for b in range(N):
            if a == b:
                continue
            va, ka, ia = verts[a]
            vb, kb, ib = verts[b]
            if va == vb:
                adj[a][b] = (A if ka == "A" else B)[ia][ib]
            else:
                adj[a][b] = K[va][vb]
    return adj


def canon_classes(n):
    # one representative per isomorphism class (by brute-force canonical form), n <= 4
    seen = {}
    for code in range(1 << (n * (n - 1) // 2)):
        T = tour(n, code)
        best = None
        for perm in itertools.permutations(range(n)):
            key = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
            if best is None or key < best:
                best = key
        seen.setdefault(best, T)
    return list(seen.values())


def main():
    print("# kernel_mod3.py (auditor hypotheses; own code)")
    mods = []
    for n in (1, 2, 3, 4):
        for T in canon_classes(n):
            mods.append((n, T, int(ham(masks(T), n))))
    print("modules: %s" % [(n, h) for n, _, h in mods])
    for m in (3, 4, 5):
        nk = 0
        bad = 0
        three = 0
        consts = {}
        for code in range(1 << (m * (m - 1) // 2)):
            K = tour(m, code)
            for x in range(m):
                for y in range(m):
                    if x == y:
                        continue
                    nk += 1
                    c = {}
                    ok = True
                    for (na, A, hA) in mods:
                        for (nb, B, hB) in mods:
                            if m - 2 + na + nb > 10:
                                continue
                            T = substitute(K, x, y, A, B)
                            hT = int(ham(masks(T), len(T)))
                            key = (na % 2, nb % 2)
                            prod = (hA * hB) % 3
                            if prod == 0:
                                if hT % 3 != 0:
                                    ok = False
                            else:
                                # hT = c * prod mod 3  ->  c = hT * prod^-1 = hT * prod (mod 3)
                                cc = (hT * prod) % 3
                                if key in c and c[key] != cc:
                                    ok = False
                                c.setdefault(key, cc)
                    if not ok:
                        bad += 1
                    consts[tuple(sorted(c.items()))] = consts.get(tuple(sorted(c.items())), 0) + 1
        print("order %d: %d (labelled kernel, ordered ports); violations of 'H(K[A,B]) = c(parities) H(A) H(B) mod 3': %d; distinct constant tables: %d" % (
            m, nk, bad, len(consts)))


if __name__ == "__main__":
    main()
