# v_gcd.py - independent brute-force checks for Prop. 7 (B2'): S_w, N(K) for ALL integer K, extended degree bounds
# for ALL integer K >= 1 and all points; Prop. 6 family.  Exact integer arithmetic (floats only in the bounds).
import sys, math, time
import numpy as np
from numba import njit

@njit(cache=True)
def gcd(a, b):
    a = abs(a); b = abs(b)
    while b:
        a, b = b, a % b
    return a

@njit(cache=True)
def rev_bits(h, n):
    r = 0
    for i in range(n):
        if (h >> i) & 1:
            r |= 1 << (n - 1 - i)
    return r

@njit(cache=True)
def height_hist(w):
    # histogram of M = max(|X|,|Xi|)/gcd over ordered odd pairs (x, x') of w-bit words; and S_w numerator pieces
    N = 1 << w
    rv = np.empty(N, dtype=np.int64)
    for x in range(N):
        rv[x] = rev_bits(x, w)
    hist = np.zeros(N + 1, dtype=np.int64)
    for x in range(N):
        for xp in range(N):
            X = x - xp
            Xi = rv[x] - rv[xp]
            if (X & 1) and (Xi & 1):
                g = gcd(X, Xi)
                M = max(abs(X), abs(Xi)) // g
                hist[M] += 1
    return hist

@njit(cache=True)
def deg_check(w, bound_u, bound_v):
    # for every p_u = (a, b) (|a|,|b| < 2^s1, b odd): cumulative count of (a', b') (|a'|,|b'| < 2^s2, b' odd) with M <= K
    # compared with bound_u[K] for every integer K >= 1; same for Deg' with bound_v.  Returns max ratios and #violations.
    s1 = (w + 1) // 2; s2 = w // 2
    A1 = 1 << s1; A2 = 1 << s2
    KM = 1 << w
    viol = 0
    mr_u = 0.0; mr_v = 0.0
    hist = np.zeros(KM + 1, dtype=np.int64)
    for a in range(-A1 + 1, A1):
        for b in range(-A1 + 1, A1):
            if (b & 1) == 0:
                continue
            hist[:] = 0
            for ap in range(-A2 + 1, A2):
                for bp in range(-A2 + 1, A2):
                    if (bp & 1) == 0:
                        continue
                    X = ap * A1 + b
                    Xi = a * A2 + bp
                    g = gcd(X, Xi)
                    M = max(abs(X), abs(Xi)) // g
                    hist[M] += 1
            cum = 0
            for K in range(1, KM + 1):
                cum += hist[K]
                r = cum / bound_u[K]
                if r > mr_u: mr_u = r
                if cum > bound_u[K]: viol += 1
    for ap in range(-A2 + 1, A2):
        for bp in range(-A2 + 1, A2):
            if (bp & 1) == 0:
                continue
            hist[:] = 0
            for a in range(-A1 + 1, A1):
                for b in range(-A1 + 1, A1):
                    if (b & 1) == 0:
                        continue
                    X = ap * A1 + b
                    Xi = a * A2 + bp
                    g = gcd(X, Xi)
                    M = max(abs(X), abs(Xi)) // g
                    hist[M] += 1
            cum = 0
            for K in range(1, KM + 1):
                cum += hist[K]
                r = cum / bound_v[K]
                if r > mr_v: mr_v = r
                if cum > bound_v[K]: viol += 1
    return viol, mr_u, mr_v

def main():
    wmax_hist = int(sys.argv[1]); wmax_deg = int(sys.argv[2])
    beta1 = math.log2(6) / 2 + 0.2
    print("beta1 = %.6f ; 6*2^(2/5) = %.5f ; 2*3^(1/4) = %.5f (log2 = %.5f)" % (beta1, 6 * 2 ** 0.4, 2 * 3 ** 0.25, math.log2(2 * 3 ** 0.25)))
    for w in range(4, wmax_hist + 1):
        t0 = time.time()
        s1 = (w + 1) // 2; s2 = w // 2
        hist = height_hist(w)
        Mmax = int(np.max(np.nonzero(hist)[0]))
        total = int(hist.sum())
        S = sum(int(hist[M]) / M for M in range(1, len(hist)) if hist[M])
        cum = np.cumsum(hist)
        # extended N(K) bound for every integer K >= 2^s2
        okN = True; mrN = 0.0
        for K in range(1 << s2, (1 << w) + 1):
            bd = 4 * 6 ** s1 * K * (1 + math.log(K)) * (4 * K / 2 ** s1 + 1)
            r = int(cum[K]) / bd
            mrN = max(mrN, r)
            if int(cum[K]) > bd: okN = False
        # (2): S_w <= sum_i 2^-i N(2^(i+1))
        dy = sum(2.0 ** (-i) * int(cum[min(1 << (i + 1), len(cum) - 1)]) for i in range(w))
        # range (a) display: N(K)/K vs (1+s1)^2 (6*2^0.4)^s1 for K <= 2^s2
        mra = max(int(cum[K]) / K / ((1 + s1) ** 2 * (6 * 2 ** 0.4) ** s1) for K in range(1, (1 << s2) + 1))
        # N(K)/K over all K against 6^(w/2)
        mrall = max(int(cum[K]) / K / 6 ** (w / 2) for K in range(1, (1 << w) + 1))
        print("w=%2d odd pairs=%d (4^(w-1)=%d) max height=%d (<2^w: %s)  S_w=%.4f  S_w/6^(w/2)=%.4f  S_w/(w^3 2^(beta1 w))=%.3e  dyadic sum (2)=%.1f >= S_w: %s" % (
            w, total, 4 ** (w - 1), Mmax, Mmax < (1 << w), S, S / 6 ** (w / 2), S / (w ** 3 * 2 ** (beta1 * w)), dy, dy >= S))
        print("      extended N(K) bound, every integer K >= 2^s2: %s (max ratio %.4f) ; max_K<=2^s2 N(K)/K / ((1+s1)^2 (6 2^.4)^s1) = %.4f ; max_K N(K)/K / 6^(w/2) = %.3f  (%.1fs)" % (
            okN, mrN, mra, mrall, time.time() - t0))
        sys.stdout.flush()
    for w in range(4, wmax_deg + 1):
        t0 = time.time()
        s1 = (w + 1) // 2; s2 = w // 2
        KM = 1 << w
        bu = np.zeros(KM + 1); bv = np.zeros(KM + 1)
        bu[0] = 1.0; bv[0] = 1.0
        for K in range(1, KM + 1):
            bu[K] = (2 * K / 2 ** s1 + 1) * (2 ** (s2 + 1) * (1 + math.log(K)) + K)
            bv[K] = (2 * K / 2 ** s2 + 1) * (2 ** (s1 + 1) * (1 + math.log(K)) + K)
        viol, mu, mv = deg_check(w, bu, bv)
        print("w=%2d extended degree bounds for ALL points and ALL integer K in [1, 2^w]: violations=%d  max Deg/bound=%.4f  max Deg'/bound=%.4f  (%.1fs)" % (w, viol, mu, mv, time.time() - t0))
        sys.stdout.flush()
    # Prop. 6 family
    okall = True; cnt = 0
    for w in range(8, 15):
        for tp in range(2, 9):
            s = tp - 1
            P = 1 << (tp - 1)
            for mid in range(1 << (w - 8)):
                for midp in range(1 << (w - 8)):
                    x = (0b1001 << (w - 4)) | (mid << 4) | 0b1111
                    xp = (midp << 4)
                    X = x - xp
                    Xi = int(rev_bits(x, w)) - int(rev_bits(xp, w))
                    if not (X & 1 and Xi & 1): okall = False
                    # nearest integer to P Xi / X  (X > 0)
                    c = (2 * P * Xi + X) // (2 * X)
                    g = c * X - P * Xi
                    if not ((1 << s) < c <= (1 << (s + 1)) and abs(g) <= (1 << (w - 1))):
                        okall = False
                    cnt += 1
    print("Prop. 6 family, w = 8..14, t' = 2..8 (%d instances): %s" % (cnt, okall))

if __name__ == "__main__":
    main()
