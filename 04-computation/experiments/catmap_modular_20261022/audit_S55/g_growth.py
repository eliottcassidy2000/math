"""Audit G (cont.): c_eff (q -> 1), cusp-wise growth exponents, actual coefficient growth, q -> 1 ratios.
Auditor for opus-2026-10-22-S55.

For f = prod ((+-q^a; q^m)_inf)^pw, at q = exp(2 pi i h/k) e^(-t), t -> 0:
   log|(q^a;q^m)|  ~ -(pi^2/t) B2({h k' a/k})/(m k'^2),   k' = k/gcd(k,m)
   log|(-q^a;q^m)| ~ -(pi^2/t) B2({h k' a/k + k'/2})/(m k'^2),   B2(x) = x^2 - x + 1/6
(from log(1-z) = -sum z^j/j and Re Li2(e(x)) = pi^2 B2({x})).  Writing log|f| ~ (pi^2/6t) c_{h/k}, the coefficients grow like
exp(2 pi sqrt(c* n/6)) with c* = max_{h/k} c_{h/k} (circle method); c_{0/1} is the usual c_eff.
Usage: python g_growth.py [NMAX]"""
import sys, math
from fractions import Fraction
from math import gcd

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6000


def B2(x):
    x = x % 1
    return x * x - x + Fraction(1, 6)


def c_cusp(factors, h, k):
    tot = Fraction(0)
    for pw, sign, a, m in factors:            # sign -1: (q^a;q^m), +1: (-q^a;q^m)
        kp = k // gcd(k, m)
        x = Fraction(h * kp * a, k) + (Fraction(kp, 2) if sign == 1 else 0)
        tot += pw * (-B2(x) / (m * kp * kp))
    return 6 * tot                            # c_{h/k}


def tri(a, m):
    return [(1, -1, a, m), (1, -1, m - a, m), (1, -1, m, m)]


fam = {
    'RR G': tri(2, 5) + [(-1, -1, 1, 1)], 'RR H': tri(1, 5) + [(-1, -1, 1, 1)],
    'RS A': tri(3, 7) + [(-1, -1, 2, 2)], 'RS B': tri(2, 7) + [(-1, -1, 2, 2)], 'RS C': tri(1, 7) + [(-1, -1, 2, 2)],
    'R14 (61)': tri(6, 14) + [(-1, -1, 1, 1)], 'R14 (60)': tri(4, 14) + [(-1, -1, 1, 1)], 'R14 (59)': tri(2, 14) + [(-1, -1, 1, 1)],
    'GG 1': [(-1, -1, 1, 8), (-1, -1, 4, 8), (-1, -1, 7, 8)], 'GG 2': [(-1, -1, 3, 8), (-1, -1, 4, 8), (-1, -1, 5, 8)],
    'B9 A': tri(4, 9) + [(-1, -1, 3, 3)], 'B9 B': tri(2, 9) + [(-1, -1, 3, 3)], 'B9 C': tri(1, 9) + [(-1, -1, 3, 3)],
    'D27 A': tri(12, 27) + [(-1, -1, 1, 1)], 'D27 B': [(1, -1, 9, 9), (-1, -1, 1, 1)], 'D27 C': tri(6, 27) + [(-1, -1, 1, 1)],
    'D27 D': tri(3, 27) + [(-1, -1, 1, 1)],
    'JS': [(1, -1, 1, 8), (1, -1, 7, 8), (1, -1, 8, 8), (1, -1, 6, 16), (1, -1, 10, 16), (-1, -1, 1, 1)],
}

# ---------------- exact coefficients via partitions and sparse theta series
def partitions(n):
    p = [0] * (n + 1); p[0] = 1
    pent = []
    k = 1
    while True:
        g1, g2 = k * (3 * k - 1) // 2, k * (3 * k + 1) // 2
        if g1 > n:
            break
        sg = 1 if k % 2 else -1
        pent.append((g1, sg))
        if g2 <= n:
            pent.append((g2, sg))
        k += 1
    for i in range(1, n + 1):
        s = 0
        for g, sg in pent:
            if g > i:
                break
            s += sg * p[i - g]
        p[i] = s
    return p


P = partitions(NMAX)


def inv_qq(m):                                 # 1/(q^m;q^m)
    out = [0] * (NMAX + 1)
    for i in range(0, NMAX // m + 1):
        out[i * m] = P[i]
    return out


def theta_sparse(a, m):                        # (q^a, q^(m-a), q^m; q^m)_inf = sum_j (-1)^j q^(m j(j-1)/2 + a j)
    out = {}
    j = 0
    for j in range(-400, 401):
        e = m * j * (j - 1) // 2 + a * j
        if 0 <= e <= NMAX:
            out[e] = out.get(e, 0) + (-1) ** (j % 2)
    return out


def sparse_mul(sp, dense):
    out = [0] * (NMAX + 1)
    for e, c in sp.items():
        for i in range(0, NMAX + 1 - e):
            if dense[i]:
                out[i + e] += c * dense[i]
    return out


coef = {
    'RS A': sparse_mul(theta_sparse(3, 7), inv_qq(2)), 'RS C': sparse_mul(theta_sparse(1, 7), inv_qq(2)),
    'B9 A': sparse_mul(theta_sparse(4, 9), inv_qq(3)), 'B9 C': sparse_mul(theta_sparse(1, 9), inv_qq(3)),
    'RR G': sparse_mul(theta_sparse(2, 5), inv_qq(1)), 'D27 A': sparse_mul(theta_sparse(12, 27), inv_qq(1)),
}

print(f'{"family":9s} {"c_eff (q->1)":>12s} {"6/c_eff":>8s}   dominant cusp h/k (k <= 36): c*      | growth from exact coefficients n <= {NMAX}')
for name, fs in fam.items():
    c0 = c_cusp(fs, 0, 1)
    best = max(((c_cusp(fs, h, k), Fraction(h, k)) for k in range(1, 37) for h in range(k) if gcd(h, k) == 1), key=lambda x: x[0])
    est = ''
    if name in coef:
        a = coef[name]
        # estimate c from max |a_n| on two windows: log max ~ 2 pi sqrt(c n/6) + beta log n
        def lm(n1, n2):
            return max(math.log(abs(x)) for x in a[n1:n2] if x)
        n1, n2, n3 = NMAX // 4, NMAX // 2, NMAX
        L1, L2, L3 = lm(n1 - 60, n1), lm(n2 - 60, n2), lm(n3 - 60, n3)
        # solve L = A sqrt(n) + B log n + C with three points
        import numpy as np
        X = np.array([[math.sqrt(n), math.log(n), 1] for n in (n1, n2, n3)])
        A_, B_, C_ = np.linalg.solve(X, np.array([L1, L2, L3]))
        cfit = 6 * (A_ / (2 * math.pi)) ** 2
        signs = ''.join('+' if x > 0 else '-' if x < 0 else '0' for x in a[NMAX - 12:NMAX])
        est = f'fit c = {cfit:.4f}; log|a_n| at n = {NMAX}: {L3:.1f}; signs of last 12: {signs}'
    print(f'{name:9s} {str(c0):>12s} {str(6 / c0) if c0 else "inf":>8s}   {str(best[1]):>6s}: {str(best[0]):8s} {float(best[0]):.4f} | {est}')

# ---------------- q -> 1 ratios with Richardson extrapolation in t
import mpmath as mp
mp.mp.dps = 30


def logprod(t, fs):
    import numpy as np
    t = float(t); tot = 0.0
    for pw, sign, a, m in fs:                # here sign = +1 means (q^a;q^m), -1 means (-q^a;q^m) (after fix())
        K = int(80 / (t * m)) + 10
        n = a + m * np.arange(K, dtype=float)
        tot += pw * float(np.sum(np.log1p(-sign * np.exp(-t * n))))
    return mp.mpf(tot)


def ratio(f1, f2, t):
    return mp.exp(logprod(t, f1) - logprod(t, f2))


fix = lambda fs: [(pw, (1 if s == -1 else -1), a, m) for pw, s, a, m in fs]     # mp.qp(a, q) = (a; q)_inf: (q^a;q^m) -> a = +q^a
print('\nq -> 1 ratios (Richardson extrapolation from t = 0.004, 0.002, 0.001):')
phi = (1 + mp.sqrt(5)) / 2
targets = [('RR H / RR G', 'RR H', 'RR G', 1 / phi), ('RS B / RS A', 'RS B', 'RS A', mp.sinpi(mp.mpf(2) / 7) / mp.sinpi(mp.mpf(3) / 7)),
           ('RS C / RS A', 'RS C', 'RS A', mp.sinpi(mp.mpf(1) / 7) / mp.sinpi(mp.mpf(3) / 7)),
           ('R14 (60)/(61)', 'R14 (60)', 'R14 (61)', mp.sinpi(mp.mpf(2) / 7) / mp.sinpi(mp.mpf(3) / 7)),
           ('GG 2 / GG 1', 'GG 2', 'GG 1', mp.tan(mp.pi / 8)),
           ('B9 B / B9 A', 'B9 B', 'B9 A', mp.sinpi(mp.mpf(2) / 9) / mp.sinpi(mp.mpf(4) / 9)),
           ('B9 C / B9 A', 'B9 C', 'B9 A', mp.sinpi(mp.mpf(1) / 9) / mp.sinpi(mp.mpf(4) / 9)),
           ('D27 B / D27 A', 'D27 B', 'D27 A', mp.sinpi(mp.mpf(1) / 3) / mp.sinpi(mp.mpf(4) / 9))]
for lab, n1, n2, tgt in targets:
    r = [ratio(fix(fam[n1]), fix(fam[n2]), mp.mpf(t)) for t in ('0.004', '0.002', '0.001')]
    r1 = 2 * r[1] - r[0]; r2 = 2 * r[2] - r[1]; rr = (4 * r2 - r1) / 3
    print(f'  {lab:15s}: extrapolated {mp.nstr(rr, 10)}  target {mp.nstr(tgt, 10)}  diff {mp.nstr(rr - tgt, 2)}')
