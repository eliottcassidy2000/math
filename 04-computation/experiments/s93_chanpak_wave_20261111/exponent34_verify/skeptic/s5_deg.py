#!/usr/bin/env python3
"""Skeptic check 5: the degree bounds (D1), (D2), (D3) and the CRT form of (B2), by brute force.

Half split: w = s1 + s2, s1 = ceil(w/2), s2 = floor(w/2), ds = s1 - s2.
  P_u = (a, b), |a|, |b| < 2^s1, b odd;  P_v = (a', b'), |a'|, |b'| < 2^s2, b' odd;
  X = a' 2^s1 + b,  Xi = a 2^s2 + b',  M = max(|X|, |Xi|)/gcd(X, Xi).
  Deg_K(P_u)  = #{P_v : M <= K},  Deg'_K(P_v) = #{P_u : M <= K}.

Part 1 (numba, exhaustive over all P_u, P_v in the boxes, all integer K <= 2^s2):
  (D1)  Deg_K <= 4 K^2, Deg'_K <= 8 K^2
  (D2)  Deg_K, Deg'_K <= 3 * 2^s1 (3 + 2 ln K)
  (D3)  implied constant  max Deg_K / [K(1+lnK)^2 + K(1+lnK)|a|^(1/3) + K^(5/2)|a|^(-1/2)]   (a != 0)
        and the same for Deg' with |a'|.
  also  inequality (1) with the true r_s1, r_s2 (computed directly), and N(K) from the pair definition.
Part 2 (pure Python, exact rationals, small w): the CRT form.
"""
import sys, math, time
from fractions import Fraction
import numpy as np
from numba import njit, prange, set_num_threads

set_num_threads(4)

@njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

@njit(parallel=True, cache=True)
def deg_tables(s1, s2, Kmax):
    """returns cumulative tables: DegU[(a index, b index), K], DegV[(a' index, b' index), K] for K = 0..Kmax."""
    A1 = (1 << (s1 + 1)) - 1          # a in (-2^s1, 2^s1): index a + 2^s1 - 1
    B1 = 1 << s1                      # b odd in (-2^s1, 2^s1): b = 2 j - 2^s1 + 1, j = 0..2^s1-1
    A2 = (1 << (s2 + 1)) - 1
    B2 = 1 << s2
    DegU = np.zeros((A1 * B1, Kmax + 1), dtype=np.int32)
    for ia in prange(A1):
        a = ia - (1 << s1) + 1
        for jb in range(B1):
            b = 2 * jb - (1 << s1) + 1
            row = DegU[ia * B1 + jb]
            for iap in range(A2):
                ap = iap - (1 << s2) + 1
                X = ap * (1 << s1) + b
                aX = abs(X)
                for jbp in range(B2):
                    bp = 2 * jbp - (1 << s2) + 1
                    Xi = a * (1 << s2) + bp
                    aXi = abs(Xi)
                    g = gcd(aX, aXi)
                    M = (aX if aX > aXi else aXi) // g
                    if M <= Kmax:
                        row[M] += 1
            for K in range(1, Kmax + 1):
                row[K] += row[K - 1]
    DegV = np.zeros((A2 * B2, Kmax + 1), dtype=np.int32)
    for iap in prange(A2):
        ap = iap - (1 << s2) + 1
        for jbp in range(B2):
            bp = 2 * jbp - (1 << s2) + 1
            row = DegV[iap * B2 + jbp]
            for ia in range(A1):
                a = ia - (1 << s1) + 1
                Xi = a * (1 << s2) + bp
                aXi = abs(Xi)
                for jb in range(B1):
                    b = 2 * jb - (1 << s1) + 1
                    X = ap * (1 << s1) + b
                    aX = abs(X)
                    g = gcd(aX, aXi)
                    M = (aX if aX > aXi else aXi) // g
                    if M <= Kmax:
                        row[M] += 1
            for K in range(1, Kmax + 1):
                row[K] += row[K - 1]
    return DegU, DegV

def r_table(s):
    """r_s(a, b) for all (a, b) by direct enumeration of pairs; returned as dict-free dense array [a idx, b idx(all b)]"""
    n = 1 << s
    u = np.arange(n, dtype=np.int64)
    ru = np.zeros(n, dtype=np.int64)
    for i in range(s):
        ru |= ((u >> i) & 1) << (s - 1 - i)
    a = (u[:, None] - u[None, :]).ravel() + n - 1
    b = (ru[:, None] - ru[None, :]).ravel() + n - 1
    r = np.zeros((2 * n - 1, 2 * n - 1), dtype=np.int64)
    np.add.at(r, (a, b), 1)
    return r

def main():
    out = []
    def P(*a):
        s = " ".join(str(t) for t in a)
        print(s, flush=True)
        out.append(s)
    wlist = [int(t) for t in sys.argv[1].split(",")]
    outname = sys.argv[2]
    for w in wlist:
        t0 = time.time()
        s1 = (w + 1) // 2; s2 = w // 2; ds = s1 - s2
        Kmax = 1 << s2
        DegU, DegV = deg_tables(s1, s2, Kmax)
        A1 = (1 << (s1 + 1)) - 1; B1 = 1 << s1; A2 = (1 << (s2 + 1)) - 1; B2 = 1 << s2
        K = np.arange(1, Kmax + 1, dtype=np.float64)
        lnK = np.log(K)
        # (D1), (D2)
        mU = DegU[:, 1:].max(axis=0).astype(np.float64)   # max over P_u for each K
        mV = DegV[:, 1:].max(axis=0).astype(np.float64)
        d1u = float(np.max(mU / (4 * K * K))); d1v = float(np.max(mV / (8 * K * K)))
        d2u = float(np.max(mU / (3 * 2 ** s1 * (3 + 2 * lnK)))); d2v = float(np.max(mV / (3 * 2 ** s1 * (3 + 2 * lnK))))
        P("w=%d (s1=%d, s2=%d, ds=%d), all P_u (%d) and P_v (%d), all integer K <= 2^s2 = %d:" % (w, s1, s2, ds, A1 * B1, A2 * B2, Kmax))
        P("   (D1) max Deg_K/(4K^2) = %.4f ; max Deg'_K/(8K^2) = %.4f   (must be <= 1); Deg_1 max = %d, Deg'_1 max = %d"
          % (d1u, d1v, int(DegU[:, 1].max()), int(DegV[:, 1].max())))
        P("   (D2) max Deg_K/(3*2^s1(3+2lnK)) = %.4f ; max Deg'_K/(3*2^s1(3+2lnK)) = %.4f   (must be <= 1)" % (d2u, d2v))
        # (D3) implied constants
        aU = np.repeat(np.arange(A1) - (1 << s1) + 1, B1).astype(np.float64)
        aV = np.repeat(np.arange(A2) - (1 << s2) + 1, B2).astype(np.float64)
        nzU = aU != 0; nzV = aV != 0
        absU = np.abs(aU[nzU])[:, None]; absV = np.abs(aV[nzV])[:, None]
        Kr = K[None, :]; lr = lnK[None, :]
        denU = Kr * (1 + lr) ** 2 + Kr * (1 + lr) * absU ** (1 / 3) + Kr ** 2.5 * absU ** (-0.5)
        ratU = DegU[nzU][:, 1:] / denU
        denV = Kr * (1 + lr) ** 2 + Kr * (1 + lr) * absV ** (1 / 3) + Kr ** 2.5 * absV ** (-0.5)
        ratV = DegV[nzV][:, 1:] / denV
        iu = np.unravel_index(np.argmax(ratU), ratU.shape); iv = np.unravel_index(np.argmax(ratV), ratV.shape)
        P("   (D3) implied constant: max over P_u (a != 0), K of Deg_K/[K(1+lnK)^2 + K(1+lnK)|a|^(1/3) + K^(5/2)|a|^(-1/2)] = %.4f "
          "(at |a|=%d, K=%d) ; for Deg' with |a'|: %.4f (at |a'|=%d, K=%d)"
          % (float(ratU.max()), int(absU[iu[0], 0]), int(K[iu[1]]), float(ratV.max()), int(absV[iv[0], 0]), int(K[iv[1]])))
        # sharper diagnostics: per dyadic K, the max over P_u of Deg_K/K and Deg_K / (K^2), to see the true growth
        line = "   max_Pu Deg_K at K=2^i: " + " ".join("%d" % int(DegU[:, 1 << i].max()) for i in range(0, s2 + 1))
        P(line)
        line = "   max_Pv Deg'_K at K=2^i: " + " ".join("%d" % int(DegV[:, 1 << i].max()) for i in range(0, s2 + 1))
        P(line)
        # inequality (1) with true multiplicities
        r1 = r_table(s1); r2 = r_table(s2)
        # P_u index (ia, jb): a = ia - 2^s1 + 1 -> r index a + 2^s1 - 1 = ia ; b = 2 jb - 2^s1 + 1 -> index b + 2^s1 - 1 = 2 jb
        rU = r1[:, 0::2].reshape(-1).astype(np.float64)        # shape (A1, B1) flattened, b odd
        rV = r2[:, 0::2].reshape(-1).astype(np.float64)
        assert rU.size == A1 * B1 and rV.size == A2 * B2
        sumU = (rU ** 2) @ DegU[:, 1:].astype(np.float64)
        sumV = (rV ** 2) @ DegV[:, 1:].astype(np.float64)
        rhs1 = 0.5 * sumU + 0.5 * sumV
        # N(K) from the saved exact height counts of s3 (definition), if present
        try:
            cnt = np.load("cnt_w%d.npy" % w)
            NK = np.cumsum(cnt)[1:Kmax + 1].astype(np.float64)
            ok = bool(np.all(NK <= rhs1 + 1e-6))
            P("   (1): N(K) <= (1/2) sum r^2 Deg_K + (1/2) sum r'^2 Deg'_K for all K <= 2^s2: %s ; max N(K)/rhs = %.4f ; "
              "sum_{b odd} r_s1^2 = %d (6^s1 = %d), sum_{b' odd} r_s2^2 = %d (6^s2 = %d)"
              % (ok, float(np.max(NK / rhs1)), int(np.sum(rU ** 2)), 6 ** s1, int(np.sum(rV ** 2)), 6 ** s2))
            # weighted (D3)-shape: [sum r^2 Deg_K] / (6^s1 [K(1+lnK)^2 + K(1+lnK)2^(s1/3) + K^(5/2) 2^(-s2/2)])
            shape = 6 ** s1 * (K * (1 + lnK) ** 2 + K * (1 + lnK) * 2 ** (s1 / 3) + K ** 2.5 * 2 ** (-s2 / 2))
            P("   'The sum': max_K N(K)/(6^s1[K(1+lnK)^2 + K(1+lnK)2^(s1/3) + K^(5/2)2^(-s2/2)]) = %.5f ; "
              "max_K rhs(1)/(same) = %.5f" % (float(np.max(NK / shape)), float(np.max(rhs1 / shape))))
        except FileNotFoundError:
            P("   (cnt file for w=%d missing: (1) not compared)" % w)
        P("   time %.1fs" % (time.time() - t0))
    with open(outname, "w") as f:
        f.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
