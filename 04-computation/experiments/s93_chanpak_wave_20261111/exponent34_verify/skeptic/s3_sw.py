#!/usr/bin/env python3
"""Skeptic check 3: the gcd sum S_w and the height counts N(K), exactly.

 Method 1 (definition): all 4^w ordered pairs (x, x') of w-bit words, X = x - x', Xi = rev x - rev x',
          odd pairs only, M = max(|X|,|Xi|)/gcd; cnt[M] += 1.                      (w <= 13)
 Method 2 (uses (A)2, verified separately): all eps in {-1,0,1}^w with eps_0, eps_(w-1) != 0,
          weight 2^z(eps); cnt[M] += 2^z.                                          (w <= WMAX2)
 Output: S_w as an exact rational sum (printed as float), S_w / 6^(w/2), S_w / 2^(3w/2),
         N(K) for K = 2^i, the height-1 mass, and the explicit inequalities (2), (B1), (B2-shape).
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

@njit(cache=True)
def rev_bits(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r

@njit(parallel=True, cache=True)
def direct_counts(w):
    n = 1 << w
    rv = np.zeros(n, dtype=np.int64)
    for x in range(n):
        rv[x] = rev_bits(x, w)
    nthr = 4
    cnt = np.zeros((nthr, n + 1), dtype=np.int64)
    nodd = np.zeros(nthr, dtype=np.int64)
    for th in prange(nthr):
        for x in range(th, n, nthr):
            for xp in range(n):
                X = x - xp
                Xi = rv[x] - rv[xp]
                if (X & 1) and (Xi & 1):
                    aX = abs(X); aXi = abs(Xi)
                    g = gcd(aX, aXi)
                    M = max(aX, aXi) // g
                    cnt[th, M] += 1
                    nodd[th] += 1
    tot = np.zeros(n + 1, dtype=np.int64)
    for th in range(nthr):
        tot += cnt[th]
    return tot, nodd.sum()

@njit(cache=True)
def half_table(m, shiftD, w):
    """all strings in {-1,0,1}^m placed at positions shiftD .. shiftD+m-1 of a w-string."""
    n = 3 ** m
    D = np.zeros(n, dtype=np.int64)
    Ds = np.zeros(n, dtype=np.int64)
    z = np.zeros(n, dtype=np.int64)
    for idx in range(n):
        t = idx
        d = 0; ds = 0; zz = 0
        for i in range(m):
            e = t % 3 - 1
            t //= 3
            pos = shiftD + i
            d += e * (1 << pos)
            ds += e * (1 << (w - 1 - pos))
            if e == 0:
                zz += 1
        D[idx] = d; Ds[idx] = ds; z[idx] = zz
    return D, Ds, z

@njit(parallel=True, cache=True)
def eps_counts(w):
    # eps_0 = +1 fixed (the case eps_0 = -1 is the negative: same M, same weight) -> factor 2 at the end
    # eps_(w-1) = +-1 ; middle digits positions 1 .. w-2 split in two halves
    mid = w - 2
    m1 = mid // 2
    m2 = mid - m1
    D1, Ds1, z1 = half_table(m1, 1, w)
    D2, Ds2, z2 = half_table(m2, 1 + m1, w)
    n1 = D1.size; n2 = D2.size
    n = 1 << w
    nthr = 4
    cnt = np.zeros((nthr, n + 1), dtype=np.int64)
    for th in prange(nthr):
        for i2 in range(th, n2, nthr):
            for top in (-1, 1):
                baseD = 1 + D2[i2] + top * (1 << (w - 1))
                baseDs = (1 << (w - 1)) + Ds2[i2] + top
                bz = z2[i2]
                for i1 in range(n1):
                    X = baseD + D1[i1]
                    Xi = baseDs + Ds1[i1]
                    aX = abs(X); aXi = abs(Xi)
                    g = gcd(aX, aXi)
                    M = max(aX, aXi) // g
                    cnt[th, M] += 1 << (bz + z1[i1])
    tot = np.zeros(n + 1, dtype=np.int64)
    for th in range(nthr):
        tot += cnt[th]
    return 2 * tot

def report(w, cnt, P, label):
    n = 1 << w
    nodd = int(cnt.sum())
    # exact S_w
    S = Fraction(0)
    # group to keep it fast: exact rational via integer arithmetic on floats is unsafe; use Fraction on nonzero entries
    nzidx = np.nonzero(cnt)[0]
    Sf = math.fsum(int(cnt[M]) / int(M) for M in nzidx)
    if w <= 12:
        for M in nzidx:
            S += Fraction(int(cnt[M]), int(M))
        Sf_exact = float(S)
    else:
        Sf_exact = Sf
    s1 = (w + 1) // 2; s2 = w // 2
    cum = np.cumsum(cnt)
    P("%s w=%d: odd pairs=%d (4^(w-1)=%d)  S_w=%.6f  S_w/6^(w/2)=%.4f  S_w/2^(3w/2)=%.6f  S_w/(w 2^(3w/2))=%.6f  "
      "height-1 mass=%d (/6^(w/2)=%.4f)"
      % (label, w, nodd, 4 ** (w - 1) if w >= 2 else 2, Sf_exact, Sf_exact / 6 ** (w / 2), Sf_exact / 2 ** (1.5 * w),
         Sf_exact / (w * 2 ** (1.5 * w)), int(cnt[1]), int(cnt[1]) / 6 ** (w / 2)))
    # N(K) for dyadic K and the inequalities
    rows = []
    rhs2 = 0.0
    for i in range(0, w + 1):
        K = 1 << (i + 1)
        NK = int(cum[min(K, n)])
        rhs2 += NK / (1 << i)
        rows.append((i, K, NK))
    # N(K) at K = 2^i, i = 0..w
    line = "    N(2^i), i=0..%d: " % w + " ".join(str(int(cum[min(1 << i, n)])) for i in range(0, w + 1))
    P(line)
    P("    (2): S_w=%.4f <= sum_i 2^-i N(2^(i+1)) = %.4f : %s" % (Sf_exact, rhs2, Sf_exact <= rhs2 + 1e-9))
    if w >= 4:
        # (B1)-type explicit: N(K) <= 6^s1 min(8K^2, 3*2^s1 (3+2 ln K)) for 1 <= K <= 2^s2  (all integer K)
        worst1 = 0.0; worst2 = 0.0; worstB2 = 0.0; argB2 = 0
        for K in range(1, (1 << s2) + 1):
            NK = int(cum[K])
            worst1 = max(worst1, NK / (6 ** s1 * 8 * K * K))
            worst2 = max(worst2, NK / (6 ** s1 * 3 * 2 ** s1 * (3 + 2 * math.log(K))))
            r = NK / (K * 8 ** s1)
            if r > worstB2:
                worstB2 = r; argB2 = K
        P("    for all integer 1<=K<=2^s2: max N(K)/(8K^2 6^s1)=%.4f  max N(K)/(3*6^s1*2^s1(3+2lnK))=%.4f  "
          "max N(K)/(K 8^s1)=%.5f at K=%d   [the (B2) bound says the last is <= C]"
          % (worst1, worst2, worstB2, argB2))
        b1 = (6 * math.sqrt(2)) ** s1 * (82 + 12 * s2) + 2 ** (1.5 * w + 2)
        P("    (B1) explicit: S_w=%.2f <= (6 sqrt2)^s1 (82+12 s2) + 2^(3w/2+2) = %.2f : %s ; tail part 2*4^w*2^-s2=%.2f ; "
          "part of S_w from M>2^s2: %.2f ; from M<=2^s2: %.2f"
          % (Sf_exact, b1, Sf_exact <= b1, 2 * 4 ** w / 2 ** s2,
             math.fsum(int(cnt[M]) / int(M) for M in nzidx if M > (1 << s2)),
             math.fsum(int(cnt[M]) / int(M) for M in nzidx if M <= (1 << s2))))
    return Sf_exact

def main():
    out = []
    def P(*a):
        s = " ".join(str(t) for t in a)
        print(s, flush=True)
        out.append(s)
    wmax1 = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    wmax2 = int(sys.argv[2]) if len(sys.argv) > 2 else 18
    outname = sys.argv[3] if len(sys.argv) > 3 else "s3_sw.out"
    Sd = {}
    cnts = {}
    for w in range(1, wmax1 + 1):
        t0 = time.time()
        cnt, nodd = direct_counts(w)
        cnts[w] = cnt
        Sd[w] = report(w, cnt, P, "DIRECT")
        P("    time %.1fs" % (time.time() - t0))
    for w in range(3, wmax2 + 1):
        t0 = time.time()
        cnt = eps_counts(w)
        if w in cnts:
            P("EPS w=%d: counts identical to the direct enumeration: %s" % (w, np.array_equal(cnt, cnts[w])))
        else:
            report(w, cnt, P, "EPS")
        P("    time %.1fs" % (time.time() - t0))
        np.save("cnt_w%d.npy" % w, cnt)
    with open(outname, "w") as f:
        f.write("\n".join(out) + "\n")

if __name__ == "__main__":
    main()
