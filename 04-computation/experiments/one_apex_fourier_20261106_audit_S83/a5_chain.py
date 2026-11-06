# Audit S83, step 5: every inequality in the proof chain of THM-4886 (A), and of the audit's two-block bound, on exact data.
# Three-block chain (per split (a, L, l), boundary bits (p, eps), theta):
#   (c0) S_n = sum of the four S_(p,eps), each computed directly from H(uvw) = g(A(x));
#   (c1) |S_(p,eps)|^2 <= K_u K_w Q;
#   (c2) |Q - diag| <= sum_{v != v'} |G_u(dC)| |G_w(dA)|          (triangle inequality on the expansion);
#   (c3) ... <= sum_{v != v'} min(K_u, 1/(2||2t dC||)) min(K_w, 1/(2||2t dA||));
#   (c4) ... <= 2^L sqrt(T_u T_w);
#   (c5) T_K <= 2 (2^L/q + 1)(3K^2 + 2qK).
# Two-block chain (split (a, m), bit p):
#   (d1) |S_p|^2 <= K sum_u |F(u)|^2;  (d2) sum_u |F(u)|^2 <= K 2^m + 2^m sum_{0<|d|<2^m} min(K, 1/(2||2t d||));
#   (d3) that sum <= 2 (2^m/q + 1)(3K + q (1 + ln q)).
import math, cmath
import numpy as np
from fractions import Fraction
from oa_engine import g_auto, IRR

def gv(A):
    return np.array([g_auto(int(x)) for x in np.atleast_1d(A)], dtype=np.int64)

def rev(x, m):
    r = 0
    for _ in range(m):
        r = (r << 1) | (x & 1); x >>= 1
    return r

def Cval(A, m):  # C(w) = sum_j (1 - w_j) 2^(j-1), w_1 = most significant letter
    return (1 << m) - 1 - rev(A, m)

def E(x):
    return np.exp(2j * np.pi * x)

def nrm(x):
    return np.abs(x - np.round(x))

def minterm(K, beta, d):
    dd = nrm(beta * d)
    with np.errstate(divide='ignore'):
        return np.minimum(K, 1.0 / (2 * dd))

thetas = [Fraction(3, 37), Fraction(5, 101), Fraction(1, 7), Fraction(11, 257), Fraction(1234, 65537),
          ('phi', IRR['phi']), ('pi-3', IRR['pi-3'])]

def convergents(x, qmax):
    """continued-fraction convergents r/q of the Fraction x with q <= qmax"""
    out = []; h0, h1, k0, k1 = 0, 1, 1, 0; y = x
    while True:
        a = math.floor(y)
        h0, h1 = h1, a * h1 + h0; k0, k1 = k1, a * k1 + k0
        if k1 > qmax:
            break
        out.append(Fraction(h1, k1))
        if y == a:
            break
        y = 1 / (y - a)
    return out

def q_for(th, n):
    if isinstance(th, Fraction):
        return (2 * th).denominator, float(th)
    b = 2 * th[1]
    fr = convergents(b, 1 << (n // 2))[-1]          # largest convergent denominator <= 2^(n/2)
    assert abs(b - fr) <= Fraction(1, fr.denominator ** 2)
    return fr.denominator, float(th[1])

worst = {k: 0.0 for k in ('c1', 'c2', 'c3', 'c4', 'c5', 'd1', 'd2', 'd3')}
viol = {k: 0 for k in worst}; c0bad = 0; ncases = 0; ncases2 = 0
for n in (10, 12, 14):
    allA = np.arange(1 << (n - 1), 1 << n, dtype=np.int64)
    Hall = gv(allA)
    for th in thetas:
        q, t = q_for(th, n)
        Sn = E(t * Hall.astype(np.float64)).sum()
        for L in range(1, n - 2):
            for a in range(2, n - L):
                l = n - a - L
                Ku, Kw = 1 << (a - 2), 1 << (l - 1)
                Av = np.arange(1 << L, dtype=np.int64); Hv = gv(Av)
                Cv = np.array([Cval(int(x), L) for x in Av], dtype=np.int64)
                epsv = (Av < (1 << (L - 1))).astype(np.int64); pv = Av & 1
                tot = 0j
                for p in (0, 1):
                    Au = (1 << (a - 1)) + p + 2 * np.arange(Ku, dtype=np.int64)
                    for eps in (0, 1):
                        Aw = (1 - eps) * (1 << (l - 1)) + np.arange(Kw, dtype=np.int64)
                        Cw = np.array([Cval(int(x), l) for x in Aw], dtype=np.int64)
                        # direct S_(p,eps)
                        Ax = (Au[:, None, None] << (L + l)) + (Av[None, :, None] << l) + Aw[None, None, :]
                        Sd = E(t * gv(Ax.ravel()).astype(np.float64)).sum()
                        tot += Sd
                        Hhat = Hv + p * epsv + pv * eps
                        M1 = E(t * np.outer(Au, Cv).astype(np.float64))
                        M2 = E(t * np.outer(Av, Cw).astype(np.float64))
                        F = (M1 * E(t * Hhat.astype(np.float64))[None, :]) @ M2
                        Q = float((np.abs(F) ** 2).sum()); diag = (1 << L) * Ku * Kw
                        ncases += 1
                        r1 = abs(Sd) ** 2 / (Ku * Kw * Q); worst['c1'] = max(worst['c1'], r1); viol['c1'] += r1 > 1 + 1e-9
                        dC = Cv[:, None] - Cv[None, :]; dA = Av[:, None] - Av[None, :]
                        off = ~np.eye(1 << L, dtype=bool)
                        Gu = np.abs(E(t * np.outer(dC[off].astype(np.float64), Au.astype(np.float64))).sum(axis=1)) if (1 << L) > 1 else np.zeros(0)
                        Gw = np.abs(E(t * np.outer(dA[off].astype(np.float64), Cw.astype(np.float64))).sum(axis=1)) if (1 << L) > 1 else np.zeros(0)
                        offabs = float((Gu * Gw).sum())
                        r2 = abs(Q - diag) / offabs if offabs > 0 else 0.0
                        worst['c2'] = max(worst['c2'], r2); viol['c2'] += abs(Q - diag) > offabs + 1e-6
                        offmin = float((minterm(Ku, 2 * t, dC[off]) * minterm(Kw, 2 * t, dA[off])).sum())
                        r3 = offabs / offmin if offmin > 0 else 0.0
                        worst['c3'] = max(worst['c3'], r3); viol['c3'] += offabs > offmin * (1 + 1e-9)
                        ds = np.concatenate([np.arange(1, 1 << L), -np.arange(1, 1 << L)])
                        Tu = float((minterm(Ku, 2 * t, ds) ** 2).sum()); Tw = float((minterm(Kw, 2 * t, ds) ** 2).sum())
                        r4 = offmin / ((1 << L) * math.sqrt(Tu * Tw)) if Tu * Tw > 0 else 0.0
                        worst['c4'] = max(worst['c4'], r4); viol['c4'] += r4 > 1 + 1e-9
                        for K, T in ((Ku, Tu), (Kw, Tw)):
                            r5 = T / (2 * ((1 << L) / q + 1) * (3 * K * K + 2 * q * K))
                            worst['c5'] = max(worst['c5'], r5); viol['c5'] += r5 > 1 + 1e-9
                c0bad += abs(tot - Sn) > 1e-6 * (1 << n)
        # two-block chain
        for a in range(2, n):
            m = n - a; K = 1 << (a - 2)
            Ay = np.arange(1 << m, dtype=np.int64); Hy = gv(Ay)
            Cy = np.array([Cval(int(x), m) for x in Ay], dtype=np.int64)
            epsy = (Ay < (1 << (m - 1))).astype(np.int64)
            for p in (0, 1):
                Au = (1 << (a - 1)) + p + 2 * np.arange(K, dtype=np.int64)
                Ax = (Au[:, None] << m) + Ay[None, :]
                Sp = E(t * gv(Ax.ravel()).astype(np.float64)).sum()
                b = E(t * (Hy - 1 + p * epsy).astype(np.float64))
                Fu = (E(t * np.outer(Au, Cy).astype(np.float64)) * b[None, :]).sum(axis=1)
                s2 = float((np.abs(Fu) ** 2).sum())
                ncases2 += 1
                r = abs(Sp) ** 2 / (K * s2); worst['d1'] = max(worst['d1'], r); viol['d1'] += r > 1 + 1e-9
                ds = np.concatenate([np.arange(1, 1 << m), -np.arange(1, 1 << m)])
                msum = float(minterm(K, 2 * t, ds).sum()) if m >= 1 else 0.0
                r = s2 / (K * (1 << m) + (1 << m) * msum); worst['d2'] = max(worst['d2'], r); viol['d2'] += r > 1 + 1e-9
                r = msum / (2 * ((1 << m) / q + 1) * (3 * K + q * (1 + math.log(q)))); worst['d3'] = max(worst['d3'], r); viol['d3'] += r > 1 + 1e-9
print(f'three-block chain: {ncases} (n, theta, split, p, eps) cases; S_n = sum of the four pieces failures: {c0bad}')
for k in ('c1', 'c2', 'c3', 'c4', 'c5'):
    print(f'  {k}: violations {viol[k]}, max lhs/rhs {worst[k]:.4f}')
print(f'two-block chain: {ncases2} (n, theta, a, p) cases')
for k in ('d1', 'd2', 'd3'):
    print(f'  {k}: violations {viol[k]}, max lhs/rhs {worst[k]:.4f}')
