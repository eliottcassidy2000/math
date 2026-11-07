"""S84 EXPLORATORY: criterion (DC) of THM-4890 (D) measured directly.

Split m = floor(n/2) - 2, C = 2^m, P = 2^(n-m-1). With c = C_m(l), g(l) = g(c) (converse), and the parity term
par(h, c) = [h odd][c odd], the off-diagonal energy is
    E_n = |X_n| + T3 + T4,
    T3 = #{h != h', c : g(h) + h c + par(h,c) = g(h') + h' c + par(h',c)}            (deterministic)
    T4 = #{h != h', c != c' : g(c') - g(c) = t},   t = h c - h' c' - (g(h') - g(h)) + par(h,c) - par(h',c').
Model: the micro values replaced by i.i.d. draws with the empirical law of g on [0, C):
    T4_model = sum over the same (h, h', c, c') of rho(t),  rho(t) = sum_N p(N) p(N + t).
Reports R_n = T4 / T4_model (DC holds with constant R if R stays bounded) and checks E_n.
"""
import sys, time
import numpy as np
from numba import njit


@njit(cache=True)
def gtable(N):
    g = np.empty(1 << N, dtype=np.int64)
    g[0] = 1
    size = 1
    for _ in range(N):
        for v in range(size - 1, -1, -1):
            gv = g[v]
            g[2 * v + 1] = gv
            g[2 * v] = gv + v + (v & 1)
        size *= 2
    g[0] = 1
    return g


@njit(cache=True)
def rcw(x, w):
    r = 0
    for i in range(w):
        r |= (1 - ((x >> i) & 1)) << (w - 1 - i)
    return r


@njit(cache=True)
def counts(g, n, m, rho, off):
    C = 1 << m
    P = 1 << (n - m - 1)
    T3 = 0
    T4 = 0
    T4c = 0
    T4m = 0.0
    for h in range(P, 2 * P):
        gh = g[h]
        ho = h & 1
        for hp in range(P, 2 * P):
            if hp == h:
                continue
            D = g[hp] - gh
            hpo = hp & 1
            for c in range(C):
                co = c & 1
                # T3: c' = c
                if gh + h * c + ho * co == g[hp] + hp * c + hpo * co:
                    T3 += 1
                X = h * c - D + ho * co
                lo = (X - C - 1) // hp
                if lo < 0:
                    lo = 0
                hi = (X + C + 1) // hp + 1
                if hi > C - 1:
                    hi = C - 1
                for cp in range(lo, hi + 1):
                    if cp == c:
                        continue
                    t = X - hp * cp - hpo * (cp & 1)
                    if t <= -C or t >= C:
                        continue
                    if g[cp] - g[c] == t:
                        T4 += 1
                        v = h * C + rcw(c, m)
                        vp = hp * C + rcw(cp, m)
                        if rcw(v, n) == vp:
                            T4c += 1
                    T4m += rho[t + off]
    return T3, T4, T4c, T4m


def main(ns):
    N = max(ns)
    g = gtable(N)
    for n in ns:
        t0 = time.time()
        m = n // 2 - 2
        C = 1 << m
        vals = g[:C]
        p = np.bincount(vals, minlength=2 * C + 2).astype(np.float64) / C
        # rho(t) = sum_N p(N) p(N+t), t in (-C, C)
        full = np.correlate(p, p, mode='full')   # index k <-> shift k - (len-1); full[k] = sum_N p(N+ (k-(L-1))) p(N)
        L = len(p)
        off_full = L - 1
        rho = np.zeros(2 * C + 1)
        for t in range(-C, C + 1):
            rho[t + C] = full[off_full + t] if 0 <= off_full + t < len(full) else 0.0
        T3, T4, T4c, T4m = counts(g, n, m, rho, C)
        X = 1 << (n - 1)
        vv = g[X:2 * X]
        _, cnt = np.unique(vv, return_counts=True)
        En = int((cnt.astype(np.int64) ** 2).sum())
        print('n=%d m=%d C=%d P=%d: E_n/|X| = %.4f ; check |X|+T3+T4 = %d vs E_n = %d ; T3/|X| = %.4f ; '
              'T4/|X| = %.4f (converse pairs %.4f) ; T4_model/|X| = %.4f ; R_n = T4/T4_model = %.4f ; '
              'R*_n = (T4 - converse)/T4_model = %.4f  [%.1fs]'
              % (n, m, C, 1 << (n - m - 1), En / X, X + T3 + T4, En, T3 / X, T4 / X, T4c / X, T4m / X, T4 / T4m,
                 (T4 - T4c) / T4m, time.time() - t0))
        sys.stdout.flush()


if __name__ == '__main__':
    main([int(x) for x in sys.argv[1:]])
