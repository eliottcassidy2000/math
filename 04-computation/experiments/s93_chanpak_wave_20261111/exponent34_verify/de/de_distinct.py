#!/usr/bin/env python3
"""Independent verifier (part D/E) - distinct values of x on n-bit words and the lower-bound chain of (E).

For each n: W = n-bit words (top bit 1) with no run of more than R = ceil(log2 n) + 1 equal bits, t = ceil(n/4),
prefix classes W_p (p = first t bits), d_p = distinct values on W_p, E_p = ordered colliding pairs inside W_p.
Chain of (E):
   D_all >= D_W >= (sum_p d_p) / 2^(t-1) >= |W|^2 / (2^(t-1) (|W| + E)) >= 2^(2n-4) / (2^(t-1) (2^(n-2) + E)),
with the first Cauchy-Schwarz |W_p|^2 <= d_p (|W_p| + E_p) checked class by class (empty classes included) and the
second one sum_p d_p >= |W|^2/(|W| + E).  Also: |W| >= 2^(n-2); the count of words with a long run against
(n - R) 2^(n-1-R); "a value is taken in at most 2^(t-1) classes"; the same chain for the elementary variant
t = ceil(n (beta0 - 1)/(2 beta0 - 1)); the total energy; and the optimal t for the chain.
Exact integer arithmetic (Python ints for the chain).
"""
import sys, math, time
from fractions import Fraction
import numpy as np
from numba import njit


def q_table32(nmax):
    N = 1 << nmax
    q = np.zeros(N, dtype=np.int32)
    for k in range(nmax):
        lo, hi = 1 << k, 1 << (k + 1)
        h = np.arange(lo, hi, dtype=np.int64)
        half = h >> 1
        q[lo:hi] = (q[half].astype(np.int64) + (h & 1) * ((half + 1) >> 1)).astype(np.int32)
    return q


@njit(cache=True)
def max_runs(n):
    lo = 1 << (n - 1)
    m = lo
    out = np.empty(m, dtype=np.int8)
    for i in range(m):
        h = lo + i
        best = 1
        cur = 1
        prev = h & 1
        for k in range(1, n):
            b = (h >> k) & 1
            if b == prev:
                cur += 1
                if cur > best:
                    best = cur
            else:
                cur = 1
                prev = b
        out[i] = best
    return out


def chain(n, xs, inW, t):
    lo = 1 << (n - 1)
    E = 0
    sum_dp = 0
    cs1_ok = True
    empty = 0
    nclasses = 1 << (t - 1)
    width = 1 << (n - t)
    for p in range(1 << (t - 1), 1 << t):
        a = (p << (n - t)) - lo
        v = xs[a:a + width][inW[a:a + width]]
        if v.size == 0:
            empty += 1
            continue
        u, c = np.unique(v, return_counts=True)
        c = c.astype(np.int64)
        d_p = int(u.size)
        E_p = int((c * (c - 1)).sum())
        if int(v.size) ** 2 > d_p * (int(v.size) + E_p):
            cs1_ok = False
        E += E_p
        sum_dp += d_p
    return E, sum_dp, cs1_ok, empty, nclasses


def main():
    n0 = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    n1 = int(sys.argv[2]) if len(sys.argv) > 2 else 26
    beta0 = math.log2(6) / 2 + 0.25
    q = q_table32(n1)
    print("beta0 = %.9f, (beta0-1)/(2 beta0-1) = %.9f, beta0/(2 beta0-1) = %.9f" %
          (beta0, (beta0 - 1) / (2 * beta0 - 1), beta0 / (2 * beta0 - 1)))
    bad = 0
    for n in range(n0, n1 + 1):
        t0 = time.time()
        lo = 1 << (n - 1)
        h = np.arange(lo, 2 * lo, dtype=np.int64)
        xs = (2 * h + 1 - 2 * q[lo:2 * lo].astype(np.int64)).astype(np.int64)
        del h
        assert int(xs.max()) <= 2 ** (n + 1) - 1 and int(xs.min()) >= lo + 2
        R = math.ceil(math.log2(n)) + 1
        mr = max_runs(n)
        inW = mr <= R
        Wsize = int(inW.sum())
        longrun = lo - Wsize
        ub_longrun = (n - R) * 2 ** (n - 1 - R)
        u_all, c_all = np.unique(xs, return_counts=True)
        D_all = int(u_all.size)
        E_total = int((c_all.astype(np.int64) * (c_all.astype(np.int64) - 1)).sum())
        maxmult = int(c_all.max())
        del u_all, c_all
        D_W = int(np.unique(xs[inW]).size)
        out = []
        for label, t in (('t=ceil(n/4)', -(-n // 4)),
                         ('t=ceil(n(b0-1)/(2b0-1))', math.ceil(n * (beta0 - 1) / (2 * beta0 - 1)))):
            E, sum_dp, cs1_ok, empty, ncl = chain(n, xs, inW, t)
            LB1 = Fraction(sum_dp, 2 ** (t - 1))
            LB2 = Fraction(Wsize ** 2, 2 ** (t - 1) * (Wsize + E))
            LB3 = Fraction(2 ** (2 * n - 4), 2 ** (t - 1) * (2 ** (n - 2) + E))
            ok = cs1_ok and (D_all >= D_W) and (D_W >= LB1) and (LB1 >= LB2) and (LB2 >= LB3) \
                and (sum_dp * (Wsize + E) >= Wsize ** 2) and (Wsize >= 2 ** (n - 2)) and (longrun <= ub_longrun) \
                and (ub_longrun <= 2 ** (n - 2)) and (4 * t >= n)
            if not ok:
                bad += 1
            out.append((label, t, E, sum_dp, LB1, LB2, LB3, ok, empty, ncl))
        # all-words chain (no run condition) for every t, to see the best t
        best = None
        allW = np.ones(lo, dtype=bool)
        for tt in range(1, n):
            Et, sdp, _, _, _ = chain(n, xs, allW, tt)
            lb = Fraction(lo * lo, 2 ** (tt - 1) * (lo + Et))
            if best is None or lb > best[1]:
                best = (tt, lb, Et)
        print("n=%2d R=%d |W|=%d (2^(n-2)=%d, |W|/2^(n-1)=%.4f) long-run words=%d <= (n-R)2^(n-1-R)=%d <= 2^(n-2)  "
              "D_all=%d (D_all/2^n=%.4f, D_all/2^(3n/4)=%.2f) D_W=%d  E_total(all words, ordered)=%d (=%.4f*2^(n-1)) max mult=%d"
              % (n, R, Wsize, 2 ** (n - 2), Wsize / lo, longrun, ub_longrun, D_all, D_all / 2 ** n, D_all / 2 ** (0.75 * n),
                 D_W, E_total, E_total / lo, maxmult))
        for (label, t, E, sum_dp, LB1, LB2, LB3, ok, empty, ncl) in out:
            print("     %-24s t=%d classes=%d (empty %d)  E=%d (E/2^n=%.5f)  sum d_p=%d  chain: D_W=%d >= sum d_p/2^(t-1)=%.1f >= "
                  "|W|^2/(2^(t-1)(|W|+E))=%.1f >= 2^(2n-4)/(2^(t-1)(2^(n-2)+E))=%.1f ; /2^(n-t)=%.4f ; /2^(3n/4)=%.4f  %s"
                  % (label, t, ncl, empty, E, E / 2 ** n, sum_dp, D_W, float(LB1), float(LB2), float(LB3),
                     float(LB3) / 2 ** (n - t), float(LB3) / 2 ** (0.75 * n), 'OK' if ok else 'VIOLATION'))
        print("     all words, best t for |X|^2/(2^(t-1)(|X|+E_t)): t=%d bound=%.1f (E_t=%d) ; true D_all=%d   %.1fs"
              % (best[0], float(best[1]), best[2], D_all, time.time() - t0))
        sys.stdout.flush()
    print("TOTAL VIOLATIONS:", bad)


if __name__ == '__main__':
    main()
