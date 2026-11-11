#!/usr/bin/env python3
"""Auditor proof-5135, script a3: THM-5135 (3)(b), the E8 form.

Independent recomputation:
  [C] the Cartan matrix of E8 from the Bourbaki simple roots in the model D8+ (Gram matrix of explicit vectors),
      determinant 1; Q(x) = x^T C x / 2 = sum x_i^2 - sum_{edges} x_i x_j.
  [M] class counts of Q modulo 2, 4, 8 by direct enumeration of all residue vectors (numpy, chunked) and
      modulo 16 by a numba loop over all 16^8 vectors; independently by an exact dynamic programme on the
      Dynkin tree for the moduli 2, 4, ..., 128 and 3, 9, 27, 5, 25, 7, 49, 11, 13, and composite 6, 12, 15, 36;
      compared with the rank-4 seam (4-fold additive convolution of the multiplication table).
      For contrast: the forms H^3 + A (A = x^2 + xy + y^2), D8 and Z^8 at 2.
  [T] shells: #{v in E8 : v.v = 2n} from the model D8+ by generating functions, n <= 60, against 240 sigma_3(n);
      explicit enumeration of the shells n = 1, 2, 3 with conversion to root coordinates (so that the count is
      a count of integer solutions of x^T C x / 2 = n).
  [I] the identity (8 pi^4 / 3) n^3 * sigma_3(n) / (zeta(4) n^3) = 240 sigma_3(n) and the real density
      d/dn vol{x : x^T C x / 2 <= n} = (8 pi^4 / 3) n^3 (Monte Carlo free: closed form of the 8-ball).
"""
import sys
import itertools
from fractions import Fraction
from math import pi, gamma

import numpy as np

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False


def table(m):
    cnt = [0] * m
    for a in range(m):
        for b in range(m):
            cnt[(a * b) % m] += 1
    return cnt


def cyc_conv(a, b):
    m = len(a)
    out = [0] * m
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            out[(i + j) % m] += ai * bj
    return out


def seam_counts(m, r):
    t = table(m)
    cur = t[:]
    for _ in range(r - 1):
        cur = cyc_conv(cur, t)
    return cur


# ---------------------------------------------------------------- the lattice
def simple_roots():
    """Bourbaki simple roots of E8 in R^8 (times 2, to stay in integers)."""
    e = np.eye(8, dtype=np.int64)
    a = []
    a.append(e[0] + e[7] - (e[1] + e[2] + e[3] + e[4] + e[5] + e[6]))  # 2 * alpha_1
    a.append(2 * (e[0] + e[1]))                                          # 2 * alpha_2
    for i in range(1, 7):
        a.append(2 * (e[i] - e[i - 1]))                                  # 2 * alpha_3 .. alpha_8
    return np.array(a, dtype=np.int64)  # rows = 2 alpha_i


def tree_dp_counts(edges, nvert, m, diag=None):
    """Exact class counts modulo m of Q(x) = sum_i diag_i x_i^2 - sum_{edges} x_i x_j on a tree."""
    if diag is None:
        diag = [1] * nvert
    adj = {i: [] for i in range(nvert)}
    for (i, j) in edges:
        adj[i].append(j)
        adj[j].append(i)

    def conv_rows(A, B):
        # cyclic convolution over the class index, row by row (same x)
        out = np.zeros_like(A)
        for s in range(m):
            out += A * 0  # no-op, keeps dtype
        res = np.zeros((m, m), dtype=object)
        for x in range(m):
            a = A[x]
            b = B[x]
            row = [0] * m
            for i in range(m):
                ai = int(a[i])
                if ai == 0:
                    continue
                for j in range(m):
                    bj = int(b[j])
                    if bj:
                        row[(i + j) % m] += ai * bj
            res[x] = row
        return res

    def F(v, parent):
        # F[x][c] = number of assignments of the proper descendants of v with partial form value c, x_v = x
        cur = np.zeros((m, m), dtype=object)
        for x in range(m):
            cur[x, (diag[v] * x * x) % m] = 1
        for u in adj[v]:
            if u == parent:
                continue
            Fu = F(u, v)
            G = np.zeros((m, m), dtype=object)
            for x in range(m):
                for xu in range(m):
                    sh = (x * xu) % m
                    # total = c' - x xu  => c' = c + x xu ; G[x][c] += Fu[xu][c + sh]
                    G[x] += np.roll(Fu[xu], -sh)
            cur = conv_rows(cur, G)
        return cur

    root = 0
    Fr = F(root, -1)
    return [int(sum(Fr[x][c] for x in range(m))) for c in range(m)]


def brute_counts_np(C, m):
    """Direct enumeration of x^T C x / 2 over (Z/m)^8, chunked over the first two coordinates."""
    n = C.shape[0]
    rest = np.array(list(itertools.product(range(m), repeat=n - 2)), dtype=np.int64)
    cnt = np.zeros(m, dtype=np.int64)
    C = C.astype(np.int64)
    # Q = sum_i x_i^2 + sum_{i<j} C_ij x_i x_j
    for x0 in range(m):
        for x1 in range(m):
            head = np.array([x0, x1], dtype=np.int64)
            X = np.concatenate([np.broadcast_to(head, (rest.shape[0], 2)), rest], axis=1)
            q = (X * X).sum(axis=1)
            for i in range(n):
                for j in range(i + 1, n):
                    if C[i, j]:
                        q += C[i, j] * X[:, i] * X[:, j]
            cnt += np.bincount(q % m, minlength=m)
    return [int(x) for x in cnt]


if HAVE_NUMBA:
    @njit(parallel=True, cache=False)
    def brute16(C, m):
        cnt = np.zeros((m, m), dtype=np.int64)
        for x0 in prange(m):
            loc = np.zeros(m, dtype=np.int64)
            x = np.zeros(8, dtype=np.int64)
            x[0] = x0
            for x1 in range(m):
                x[1] = x1
                for x2 in range(m):
                    x[2] = x2
                    for x3 in range(m):
                        x[3] = x3
                        for x4 in range(m):
                            x[4] = x4
                            for x5 in range(m):
                                x[5] = x5
                                for x6 in range(m):
                                    x[6] = x6
                                    # partial value without x7
                                    q = 0
                                    for i in range(7):
                                        q += x[i] * x[i]
                                        for j in range(i + 1, 7):
                                            q += C[i, j] * x[i] * x[j]
                                    lin = 0
                                    for i in range(7):
                                        lin += C[i, 7] * x[i]
                                    for x7 in range(m):
                                        v = (q + x7 * x7 + lin * x7) % m
                                        if v < 0:
                                            v += m
                                        loc[v] += 1
            for c in range(m):
                cnt[x0, c] = loc[c]
        return cnt.sum(axis=0)


def main():
    out = []
    P = out.append

    # ---------------------------------------------------------------- [C]
    A2 = simple_roots()                 # rows 2*alpha_i
    G4 = A2 @ A2.T                      # 4 * Gram
    assert (G4 % 4 == 0).all()
    C = G4 // 4
    P("[C] Gram matrix of the Bourbaki simple roots (model D8+):")
    for row in C:
        P("    " + " ".join("%2d" % v for v in row))
    detC = int(round(np.linalg.det(C.astype(float))))
    P("    det = %d; diagonal all 2: %s; off-diagonal in {0,-1}: %s"
      % (detC, bool((np.diag(C) == 2).all()), bool(((C - np.diag(np.diag(C))) <= 0).all())))
    edges = [(i, j) for i in range(8) for j in range(i + 1, 8) if C[i, j] == -1]
    P("    Dynkin edges (0-based): %s (a tree with %d edges; branch node of degree 3: %s)"
      % (edges, len(edges), [v for v in range(8) if sum(1 for e in edges if v in e) == 3]))
    # exact determinant by fraction elimination
    M = [[Fraction(int(v)) for v in row] for row in C]
    det = Fraction(1)
    for i in range(8):
        piv = next(r for r in range(i, 8) if M[r][i] != 0)
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            det = -det
        det *= M[i][i]
        for r in range(i + 1, 8):
            f = M[r][i] / M[i][i]
            for c in range(i, 8):
                M[r][c] -= f * M[i][c]
    P("    exact determinant: %s" % det)

    # ---------------------------------------------------------------- [M]
    P("")
    P("[M] class counts of Q = x^T C x / 2 against the rank-4 seam")
    for m in (2, 4, 8):
        b = brute_counts_np(C, m)
        d = tree_dp_counts(edges, 8, m)
        s = seam_counts(m, 4)
        P("  modulo %d: direct enumeration of all %d^8 vectors: %s" % (m, m, b))
        P("            tree DP agrees with enumeration: %s; equals rank-4 seam counts: %s" % (b == d, b == s))
    if HAVE_NUMBA:
        b16 = [int(x) for x in brute16(C.astype(np.int64), 16)]
        d16 = tree_dp_counts(edges, 8, 16)
        s16 = seam_counts(16, 4)
        P("  modulo 16: numba enumeration of all 16^8 = %d vectors; total %d" % (16 ** 8, sum(b16)))
        P("            tree DP agrees: %s; equals rank-4 seam counts: %s" % (b16 == d16, b16 == s16))
        P("            counts by valuation v = 0..4: %s" % [b16[(2 ** v) % 16] for v in range(5)])
    for m in (32, 64, 128):
        d = tree_dp_counts(edges, 8, m)
        s = seam_counts(m, 4)
        P("  modulo %d (tree DP, exact): equals rank-4 seam counts: %s" % (m, d == s))
    oddres = []
    for m in (3, 9, 27, 5, 25, 7, 49, 11, 13, 6, 12, 15, 36):
        d = tree_dp_counts(edges, 8, m)
        s = seam_counts(m, 4)
        oddres.append((m, d == s))
    P("  odd prime powers and composite moduli (tree DP): %s" % oddres)

    # contrast forms at 2
    P("  contrast at 2 (tree DP / convolution):")
    # H^3 + A: A = x^2 + x y + y^2 ; as a 'tree' form: diag (1,1), edge with coefficient +1 -> use x -> -x: x^2 - xy + y^2
    for m in (2, 4, 8):
        tA = tree_dp_counts([(0, 1)], 2, m)             # x^2 - xy + y^2, equivalent to x^2 + xy + y^2
        h3 = seam_counts(m, 3)
        h3a = cyc_conv(h3, tA)
        s = seam_counts(m, 4)
        # D8 root lattice form: tree D8
        d8_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (5, 7)]
        d8 = tree_dp_counts(d8_edges, 8, m)
        sq = [0] * m
        for a in range(m):
            sq[(a * a) % m] += 1
        z8 = sq[:]
        for _ in range(7):
            z8 = cyc_conv(z8, sq)
        P("    modulo %d: H^3 + A equals seam: %s; D8 equals seam: %s; Z^8 (eight squares) equals seam: %s"
          % (m, h3a == s, d8 == s, z8 == s))
        if m == 2:
            P("      zeros modulo 2: E8/seam %d, H^3 + A %d, D8 %d, Z^8 %d" % (s[0], h3a[0], d8[0], z8[0]))

    # ---------------------------------------------------------------- [T]
    P("")
    P("[T] shells of E8")
    NMAX = 60
    # integer vectors: number with sum of squares = 2n (coordinate sum is then even automatically)
    L = 2 * NMAX
    one = np.zeros(L + 1, dtype=object)
    x = 0
    while x * x <= L:
        one[x * x] += 1 if x == 0 else 2
        x += 1
    cur = one.copy()
    for _ in range(7):
        nxt = np.zeros(L + 1, dtype=object)
        for i in range(L + 1):
            if cur[i] == 0:
                continue
            for j in range(L + 1 - i):
                if one[j]:
                    nxt[i + j] += cur[i] * one[j]
        cur = nxt
    r8_even = [int(cur[2 * n]) for n in range(NMAX + 1)]
    # half-integer vectors: o_i odd, sum o_i^2 = 8n, sum o_i = 0 (mod 4)
    L2 = 8 * NMAX
    odd_terms = []
    o = 1
    while o * o <= L2:
        odd_terms.append((o * o, o % 4))        # +o
        odd_terms.append((o * o, (-o) % 4))     # -o
        o += 2
    state = {(0, 0): 1}
    for _ in range(8):
        nxt = {}
        for (s, t), c in state.items():
            for (q, tt) in odd_terms:
                if s + q <= L2:
                    key = (s + q, (t + tt) % 4)
                    nxt[key] = nxt.get(key, 0) + c
        state = nxt
    half = [state.get((8 * n, 0), 0) for n in range(NMAX + 1)]
    sig3 = lambda n: sum(d ** 3 for d in range(1, n + 1) if n % d == 0)
    bad = [n for n in range(1, NMAX + 1) if r8_even[n] + half[n] != 240 * sig3(n)]
    P("  #{v in D8+ : v.v = 2n} = 240 sigma_3(n) for n = 1..%d: mismatches %s" % (NMAX, bad))
    P("  n = 1, 2, 3: %s (240 sigma_3 = %s)" % ([r8_even[n] + half[n] for n in (1, 2, 3)], [240 * sig3(n) for n in (1, 2, 3)]))

    # explicit enumeration of the shells n <= 3 and conversion to root coordinates
    vecs = []
    for v in itertools.product(range(-2, 3), repeat=8):       # integer vectors, norm <= 6
        nn = sum(t * t for t in v)
        if nn in (2, 4, 6) and sum(v) % 2 == 0:
            vecs.append(tuple(2 * t for t in v))
    for v in itertools.product((-3, -1, 1, 3), repeat=8):     # half-integers (times 2), norm <= 6 -> sum o^2 <= 24
        nn = sum(t * t for t in v)
        if nn in (8, 16, 24) and sum(v) % 4 == 0:
            vecs.append(v)
    V2 = np.array(vecs, dtype=np.int64)                        # rows = 2 v
    # v = sum x_i alpha_i  <=>  2 v = x^T (2 alpha) ; solve with exact inverse via C: x = C^-1 (alpha_i . v)
    Cinv = np.linalg.inv(C.astype(float))
    dots4 = V2 @ A2.T                                           # 4 * (v . alpha_i)
    X = (dots4 / 4.0) @ Cinv.T
    Xi = np.rint(X).astype(np.int64)
    integral = bool(np.abs(X - Xi).max() < 1e-9)
    Qv = ((Xi @ C) * Xi).sum(axis=1) // 2
    shell = {n: int((Qv == n).sum()) for n in (1, 2, 3)}
    back = bool(((Xi @ A2) == V2).all())
    P("  explicit vectors of norm 2, 4, 6 in D8+: %d; root coordinates integral: %s; reconstruct exactly: %s"
      % (len(vecs), integral, back))
    P("  integer solutions of x^T C x / 2 = n, n = 1, 2, 3: %s; max |x_i| used: %d" % (shell, int(np.abs(Xi).max())))

    # ---------------------------------------------------------------- [I]
    P("")
    P("[I] the identity and the real density")
    vol8 = pi ** 4 / gamma(5)                                  # volume of the unit 8-ball = pi^4 / 24
    # vol{x^T C x <= 2n} = vol8 (2n)^4 / sqrt(det C) = (2 pi^4 / 3) n^4 ; derivative (8 pi^4 / 3) n^3
    P("  unit 8-ball volume pi^4/24 = %.9f; vol{Q <= n} = (2 pi^4/3) n^4 -> coefficient %.9f vs %.9f"
      % (vol8, vol8 * 16, 2 * pi ** 4 / 3))
    P("  real density coefficient: 4 * (2 pi^4/3) = %.9f = 8 pi^4/3 = %.9f" % (4 * vol8 * 16, 8 * pi ** 4 / 3))
    z4 = pi ** 4 / 90
    P("  (8 pi^4/3) / zeta(4) = %.12f (exactly 8 * 90 / 3 = 240)" % ((8 * pi ** 4 / 3) / z4))
    # numerical zeta(4)
    zs = sum(1.0 / k ** 4 for k in range(1, 200001))
    P("  zeta(4) by summation %.12f, pi^4/90 = %.12f" % (zs, z4))
    # Z^8: real density pi^4 n^3 / 6, odd-prime factors (1/zeta(4)) / (1 - 2^-4) sigma_3(n)/n^3 -> 16 sigma_3(n) (n odd)
    P("  Z^8, n odd: (pi^4/6) * (90/pi^4) * (16/15) = %.12f (= 16, so the 2-adic density is 1 at odd n)" % ((pi ** 4 / 6) * (90 / pi ** 4) * (16 / 15)))
    P("  four squares, n odd: pi^2 * (6/pi^2) * (4/3) = %.12f (= 8)" % (pi ** 2 * (6 / pi ** 2) * (4 / 3)))

    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
