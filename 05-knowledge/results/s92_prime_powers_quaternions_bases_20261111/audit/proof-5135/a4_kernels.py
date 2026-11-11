#!/usr/bin/env python3
"""Auditor proof-5135, script a4: THM-5135 (4), kernels.

  [K] class law of u^T G v for the 2 x 2 kernels named in the statement, modulo 9, 25, 49 (and 3, 5, 7, 27):
        [[5,2],[2,1]]   (Gaussian kernel of THM-5133 (B): (al,be,ga,de) = (1,0,0,1), det 1),
        [[11,6],[6,3]]  (hyperbolic kernel of THM-5133 (B), det -3),
        diag(1,2)       (det 2),
      against the rank-2 law; Session261's kernel [[149,72],[72,37]] (det 329 = 7 * 47) modulo 7, 47 against the
      rank-1 law (multiplication table), and modulo 49 for the record.
      Smith-form reading: the law of a degenerate kernel is the convolution of the scaled tables.
  [D] the two rank-4 Gram matrices: determinants, and their uniform class law modulo 3, 5, 9 (direct enumeration)
      against the rank-4 law; the value at the class 0 modulo p: 1 + p^-3 - p^-4.
  [F] THM-4916's factor 1 + (3p-2)(p-1)/(2p^3) against: the uniform rank-4 value 1 + p^-3 - p^-4 (two independent
      uniform states), the value 1 of a linear functional of ONE uniform state (THM-4916's statistic is F(s) =
      K(t, s) over words s, one free state), and the uniform rank-1, rank-2 values.
"""
import sys
import itertools
from fractions import Fraction

import numpy as np


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


def scaled_table(m, s):
    """distribution of s * a * b modulo m."""
    cnt = [0] * m
    for a in range(m):
        for b in range(m):
            cnt[(s * a * b) % m] += 1
    return cnt


def brute(G, m):
    r = len(G)
    Gm = np.array(G, dtype=np.int64) % m
    grid = np.array(list(itertools.product(range(m), repeat=r)), dtype=np.int64)
    Gv = (grid @ Gm.T) % m
    cnt = np.zeros(m, dtype=np.int64)
    for u in grid:
        cnt += np.bincount((Gv @ u) % m, minlength=m)
    return [int(x) for x in cnt]


def det_exact(M):
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    d = Fraction(1)
    for i in range(n):
        piv = None
        for r in range(i, n):
            if A[r][i] != 0:
                piv = r
                break
        if piv is None:
            return Fraction(0)
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            d = -d
        d *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            for c in range(i, n):
                A[r][c] -= f * A[i][c]
    return d


def main():
    out = []
    P = out.append

    P("[K] two-port kernels")
    kernels = [("[[5,2],[2,1]]", [[5, 2], [2, 1]]), ("[[11,6],[6,3]]", [[11, 6], [6, 3]]), ("diag(1,2)", [[1, 0], [0, 2]])]
    for name, G in kernels:
        d = det_exact(G)
        res = []
        for m in (3, 5, 7, 9, 25, 27, 49):
            same = brute(G, m) == seam_counts(m, 2)
            res.append("%d:%s" % (m, "law(1)" if same else "DIFFERENT"))
        P("  G = %s, det = %s: %s" % (name, d, "  ".join(res)))
    # the degenerate case [[11,6],[6,3]] modulo 3 and 9: Smith form diag(1, 3) over Z_3
    G = [[11, 6], [6, 3]]
    for m in (3, 9, 27):
        b = brute(G, m)
        sm = cyc_conv(scaled_table(m, 1), scaled_table(m, 3))
        P("  [[11,6],[6,3]] modulo %d: counts by valuation %s; convolution of the tables a1 b1 and 3 a2 b2: %s; rank-2 law: %s"
          % (m, [b[(3 ** v) % m] for v in range(len(str(m)) + 2) if 3 ** v <= m], b == sm,
             [seam_counts(m, 2)[(3 ** v) % m] for v in range(len(str(m)) + 2) if 3 ** v <= m]))
    # Session261
    G = [[149, 72], [72, 37]]
    P("  Session261 kernel [[149,72],[72,37]]: det = %s = 7 * 47 = %d" % (det_exact(G), 7 * 47))
    for p in (7, 47):
        b = brute(G, p)
        # the law is a distribution: 4 variables against the 2 variables of the table, so the counts carry a factor p^2
        P("    modulo %d: equals the rank-1 law (p^2 times the multiplication table): %s; equals the rank-2 law: %s; G mod p nonzero: %s; G mod p = %s"
          % (p, b == [p * p * x for x in table(p)], b == seam_counts(p, 2), any(x % p for row in G for x in row),
             [[x % p for x in row] for row in G]))
    b = brute(G, 49)
    sm = cyc_conv(scaled_table(49, 1), scaled_table(49, 7))
    P("    modulo 49: equals rank-1 law: %s; rank-2 law: %s; convolution of a1 b1 and 7 a2 b2 (Smith form diag(1,7)): %s"
      % (b == [49 * 49 * x for x in table(49)], b == seam_counts(49, 2), b == sm))
    P("    modulo 49, p^k * share by valuation v = 0, 1, 2: %s (rank-1 law: %s; rank-2 law: %s)"
      % ([str(Fraction(b[(7 ** v) % 49] * 49, 49 ** 4)) for v in range(3)],
         [str(Fraction(table(49)[(7 ** v) % 49] * 49, 49 ** 2)) for v in range(3)],
         [str(Fraction(seam_counts(49, 2)[(7 ** v) % 49] * 49, 49 ** 4)) for v in range(3)]))
    for p in (3, 5, 11, 13):
        P("    modulo %d (p prime to 329): rank-2 law: %s" % (p, brute(G, p) == seam_counts(p, 2)))

    P("")
    P("[D] the two rank-4 Gram matrices")
    # K-form of THM-4916 (1d) in the coordinates (U, E, d, a): U U' + E E' - d d' + a d' + a' d
    GK = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -1, 1], [0, 0, 1, 0]]
    # the same form in the state coordinates (D, E, U, d), a = D - E - U + d : G_state = T^T GK T, (U,E,d,a) = T (D,E,U,d)
    T = np.array([[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, -1, -1, 1]], dtype=np.int64)
    GKs = T.T @ np.array(GK, dtype=np.int64) @ T
    # fusion corner of Session264 eq. (3) in (h, H-h, U-h, d-h) x (h', H'-h', V'-h', e'-h')
    GF = [[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    # the same in raw coordinates (h, H, U, d) x (h', H', V', e'):  H h' + h H' - h h' + (U-h)(V'-h') + (d-h)(e'-h')
    S = np.array([[1, 0, 0, 0], [-1, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]], dtype=np.int64)
    GFr = S.T @ np.array(GF, dtype=np.int64) @ S
    P("  K-form, coordinates (U,E,d,a): det = %s; in state coordinates (D,E,U,d): matrix %s, det = %s"
      % (det_exact(GK), GKs.tolist(), det_exact(GKs.tolist())))
    P("  fusion corner, shifted coordinates: det = %s; in raw coordinates (h,H,U,d) x (h',H',V',e'): matrix %s, det = %s"
      % (det_exact(GF), GFr.tolist(), det_exact(GFr.tolist())))
    # check the raw-coordinate fusion matrix against eq. (3) on random integers
    rng = np.random.default_rng(264)
    ok = True
    for _ in range(2000):
        h, H, U, d, h2, H2, V2, e2 = [int(x) for x in rng.integers(-50, 50, size=8)]
        lhs = H * h2 + h * H2 - h * h2 + (U - h) * (V2 - h2) + (d - h) * (e2 - h2)
        rhs = int(np.array([h, H, U, d]) @ GFr @ np.array([h2, H2, V2, e2]))
        thm = h * h2 + (H - h) * h2 + h * (H2 - h2) + (U - h) * (V2 - h2) + (d - h) * (e2 - h2)
        ok = ok and lhs == rhs == thm
    P("  Session264 eq. (3) = the form displayed in THM-5135 (4) = u^T G v with the raw matrix, 2000 random integer points: %s" % ok)
    # K-form check against (1d): K = U U' + E E' - d d' + a d' + a' d with a = D - E - U + d
    ok = True
    for _ in range(2000):
        D, E, U, d, D2, E2, U2, d2 = [int(x) for x in rng.integers(-50, 50, size=8)]
        a = D - E - U + d
        a2 = D2 - E2 - U2 + d2
        lhs = U * U2 + E * E2 - d * d2 + a * d2 + a2 * d
        rhs = int(np.array([D, E, U, d]) @ GKs @ np.array([D2, E2, U2, d2]))
        ok = ok and lhs == rhs
    P("  THM-4916 (1d) K-form = u^T G v with the state matrix, 2000 random integer points: %s" % ok)
    for name, G in (("K-form", GK), ("fusion corner", GF)):
        res = []
        for m in (2, 3, 4, 5, 9):
            b = brute(G, m)
            res.append("%d:%s" % (m, "law(1)" if b == seam_counts(m, 4) else "DIFFERENT"))
        P("  %s, uniform states, direct enumeration: %s" % (name, "  ".join(res)))
    for p in (3, 5):
        b = brute(GK, p)
        P("  K-form modulo %d: p * share of the class 0 = %s; 1 + p^-3 - p^-4 = %s"
          % (p, Fraction(b[0] * p, p ** 8), 1 + Fraction(1, p ** 3) - Fraction(1, p ** 4)))

    P("")
    P("[F] THM-4916 against the uniform values at the class 0 modulo p")
    P("   p | THM-4916 K=2: 1+(3p-2)(p-1)/(2p^3) | rank 4 uniform pair: 1+p^-3-p^-4 | rank 2: 1+1/p-1/p^2 | rank 1: 2-1/p | one uniform state, linear functional: 1")
    for p in (3, 5, 7, 11, 13):
        t = 1 + Fraction((3 * p - 2) * (p - 1), 2 * p ** 3)
        r4 = 1 + Fraction(1, p ** 3) - Fraction(1, p ** 4)
        r2 = 1 + Fraction(1, p) - Fraction(1, p * p)
        r1 = 2 - Fraction(1, p)
        P("  %2d | %-10s = %.5f | %-12s = %.5f | %-8s = %.5f | %-6s = %.5f | 1" % (p, t, float(t), r4, float(r4), r2, float(r2), r1, float(r1)))
    P("  THM-4916's quoted examples: p = 3: %s (34/27), p = 5: %s (151/125)"
      % (1 + Fraction(7 * 2, 54), 1 + Fraction(13 * 4, 250)))

    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
