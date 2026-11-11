# a_constants.py -- auditor proof-5134-AE.  The constants A_(m,d)^sigma of THM-5134 (C):
#   A_(m,d)^sigma = lim_H sum_{h < H} chi_d(h) w^sigma(h),  chi_d(h) = [gcd(g(h), m) = d] - c_d,
#   w^0(h) = (1/x) ln(2 (x + y)/(x + 2 y)),  w^1(h) = (1/x) ln(1 + x/(2 y)).
# Partial sums at H = 2^k (floating point, block sums), compared with the census quantity
#   (D^(d) - c_d #T)/(2 X)   at X = 2^22, 2^24, 2^26, 2^28  (my census_L.json).
# By the proof of THM-5128 (D2) the census value is the partial sum over y1 < (X/2^s)^(1/2) less a boundary term.
# Also: the predicted limit correction of (D):  m share - L_m = (2/kappa)(sum_{d | g0} d A_(m,d)) / ln X.
import math, json, sys
import numpy as np
import numba as nb
from fractions import Fraction as Fr

HERE = "C:/tmp/wt-s92/05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/proof-5134-AE/"
K0 = 0.1615389173
K1 = 0.2125139092


def primes_of(n):
    ps = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            ps.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        ps.append(n)
    return ps


def c_d(m, d):
    v = Fr(1, d * d)
    for p in primes_of(m // d):
        v *= (1 - Fr(1, p * p))
    return v


@nb.njit(cache=True)
def gcd(a, b):
    while b:
        a, b = b, a % b
    return a


@nb.njit(cache=True)
def partial_sums(KMAX, m, divs, cd):
    # returns S[k, t, sigma] = sum_{h < 2^k} chi_{d_t}(h) w^sigma(h), and W[k, sigma] = sum_{h<2^k} w^sigma(h)
    H = 1 << KMAX
    psi = np.zeros(H, np.int64)
    nd = divs.shape[0]
    S = np.zeros((KMAX + 1, nd, 2))
    W = np.zeros((KMAX + 1, 2))
    acc = np.zeros((nd, 2))
    accw = np.zeros(2)
    comp = np.zeros((nd, 2))
    k = 0
    for h in range(0, H):
        if h >= 1:
            kk = h >> 1
            if (h & 1) == 0:
                psi[h] = psi[kk] + kk
            else:
                psi[h] = psi[kk] + (kk >> 1) + 1
        if h == (1 << k):
            for t in range(nd):
                S[k, t, 0] = acc[t, 0]
                S[k, t, 1] = acc[t, 1]
            W[k, 0] = accw[0]
            W[k, 1] = accw[1]
            k += 1
        x = 2 * psi[h] + 1
        y = 2 * h + 1
        g = gcd(gcd(x, y), m)
        w0 = math.log(2.0 * (x + y) / (x + 2.0 * y)) / x
        w1 = math.log(1.0 + x / (2.0 * y)) / x
        accw[0] += w0
        accw[1] += w1
        for t in range(nd):
            ch = -cd[t]
            if g == divs[t]:
                ch += 1.0
            # Kahan summation
            v0 = ch * w0 - comp[t, 0]
            s0 = acc[t, 0] + v0
            comp[t, 0] = (s0 - acc[t, 0]) - v0
            acc[t, 0] = s0
            v1 = ch * w1 - comp[t, 1]
            s1 = acc[t, 1] + v1
            comp[t, 1] = (s1 - acc[t, 1]) - v1
            acc[t, 1] = s1
    for t in range(nd):
        S[KMAX, t, 0] = acc[t, 0]
        S[KMAX, t, 1] = acc[t, 1]
    W[KMAX, 0] = accw[0]
    W[KMAX, 1] = accw[1]
    return S, W


def main():
    KMAX = 27
    cens = {L: json.load(open(HERE + "census_%d.json" % L)) for L in (22, 24, 26, 28)}
    for m in (9, 27, 25, 49, 45, 15):
        divs = [d for d in range(1, m + 1) if m % d == 0]
        cd = np.array([float(c_d(m, d)) for d in divs])
        S, W = partial_sums(KMAX, m, np.array(divs, np.int64), cd)
        print("=" * 110)
        print("modulus %d  divisors %s" % (m, divs))
        print("  check of the mass: sum_{h < 2^k} w^sigma(h) / (k ln 2) -> kappa_sigma:  k = %d: %.5f %.5f  (kappa %.5f %.5f; increments per digit %.6f %.6f against ln2*kappa = %.6f %.6f)"
              % (KMAX, W[KMAX, 0] / (KMAX * math.log(2)), W[KMAX, 1] / (KMAX * math.log(2)), K0, K1,
                 W[KMAX, 0] - W[KMAX - 1, 0], W[KMAX, 1] - W[KMAX - 1, 1], math.log(2) * K0, math.log(2) * K1))
        for sh in range(2):
            print("  sheet %d: partial sums of A_(m,d) at H = 2^k" % sh)
            for k in (8, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24, 26, 27):
                print("     k = %2d: " % k + "  ".join("d=%d: %+.5f" % (d, S[k, t, sh]) for t, d in enumerate(divs)))
            for L in (22, 24, 26, 28):
                dat = cens[L]["data"][str(m)]
                off = dat["off"][sh]
                dia = dat["diag"][sh]
                T = 2 * sum(sum(r) for r in off) + sum(sum(r) for r in dia)
                X = 1 << L
                vals = []
                for t, d in enumerate(divs):
                    D = 2 * sum(off[t])
                    vals.append((D - float(c_d(m, d)) * T) / (2 * X))
                print("     census X = 2^%d: (D^(d) - c_d #T)/(2X) = " % L + "  ".join("d=%d: %+.5f" % (d, v) for d, v in zip(divs, vals)))
            kap = K0 if sh == 0 else K1
            A = {d: S[KMAX, t, sh] for t, d in enumerate(divs)}
            print("     (D) with A at k = %d: correction (2/kappa) sum_{d | g0} d A_d = " % KMAX + "  ".join("g0=%d: %+.3f" % (g0, 2 / kap * sum(d * A[d] for d in divs if g0 % d == 0)) for g0 in divs)
                  + "   [to be divided by ln X; ln 2^28 = %.2f]" % (28 * math.log(2)))
            L = 28
            dat = cens[L]["data"][str(m)]
            off = dat["off"][sh]
            dia = dat["diag"][sh]
            CM = [2 * sum(off[t][r] for t in range(len(divs))) + sum(dia[t][r] for t in range(len(divs))) for r in range(m)]
            T = sum(CM)
            line = []
            for g0 in divs:
                cl = [r for r in range(m) if (math.gcd(m, r) if r else m) == g0]
                fac = m * sum(CM[r] for r in cl) / (len(cl) * T)
                lim = float(sum(d * c_d(m, d) for d in divs if g0 % d == 0))
                pred = lim + 2 / kap * sum(d * A[d] for d in divs if g0 % d == 0) / (28 * math.log(2))
                line.append("g0=%d: measured %.4f  limit %.4f  limit + 1/ln X term %.4f" % (g0, fac, lim, pred))
            print("     class factors at X = 2^28: " + ";  ".join(line))


if __name__ == "__main__":
    main()
