# Audit S83, step 8: (i) HYP-9479 P3 (Q / diagonal at the splits (9,5,9), (10,5,10), (9,7,9), theta in {3/37, 5/101,
# 11/257, phi}, all four (p, eps)) with an independent implementation; plus the same ratio for the audit's two-block
# second moment sum_u |F(u)|^2 / (K 2^m).  (ii) (C) at rational theta: |S_n(theta)| / 2^(n-1) for denominators with powers of 2.
import math
import numpy as np
from fractions import Fraction
from oa_engine import g_auto, g_block, S_exact, IRR

def gv(A):
    return np.array([g_auto(int(x)) for x in np.atleast_1d(A)], dtype=np.int64)

def Cvec(m):
    A = np.arange(1 << m, dtype=np.int64)
    r = np.zeros_like(A)
    x = A.copy()
    for _ in range(m):
        r = (r << 1) | (x & 1); x >>= 1
    return (1 << m) - 1 - r          # reverse complement

def E(x):
    return np.exp(2j * np.pi * x)

thetas = [('3/37', 3 / 37), ('5/101', 5 / 101), ('11/257', 11 / 257), ('phi', float(IRR['phi'])), ('1/7 (major)', 1 / 7)]
ratios = []
print('(i) three-block Q / diagonal:')
for (a, L, l) in [(9, 5, 9), (10, 5, 10), (9, 7, 9)]:
    Av = np.arange(1 << L, dtype=np.int64); Hv = gv(Av); Cv = Cvec(L)
    epsv = (Av < (1 << (L - 1))).astype(np.int64); pv = Av & 1
    Cl = Cvec(l)
    for name, t in thetas:
        row = []
        for p in (0, 1):
            Au = (1 << (a - 1)) + p + 2 * np.arange(1 << (a - 2), dtype=np.int64)
            for eps in (0, 1):
                Aw = (1 - eps) * (1 << (l - 1)) + np.arange(1 << (l - 1), dtype=np.int64)
                Cw = Cl[Aw]
                Hhat = Hv + p * epsv + pv * eps
                F = (E(t * np.outer(Au, Cv).astype(np.float64)) * E(t * Hhat.astype(np.float64))[None, :]) @ E(t * np.outer(Av, Cw).astype(np.float64))
                r = float((np.abs(F) ** 2).sum()) / ((1 << L) * Au.size * Aw.size)
                row.append(r)
                if 'major' not in name:
                    ratios.append(r)
        print(f'  split {(a, L, l)} theta {name:>12s}: Q/diag over (p,eps) = ' + ' '.join(f'{r:.4f}' for r in row))
print(f'  P3 cases (48): min {min(ratios):.4f} max {max(ratios):.4f}; outside [0.85, 1.2]: {sum(not 0.85 <= r <= 1.2 for r in ratios)}')

print('(ii) two-block second moment sum_u |F(u)|^2 / (K 2^m), n = 22:')
n = 22
for a in (8, 11, 14):
    m = n - a; Ay = np.arange(1 << m, dtype=np.int64); Hy = g_block_vals = None
    Hy = np.array([g_auto(int(x)) for x in Ay], dtype=np.int64)
    Cy = Cvec(m); epsy = (Ay < (1 << (m - 1))).astype(np.int64)
    out = []
    for name, t in thetas:
        vals = []
        for p in (0, 1):
            Au = (1 << (a - 1)) + p + 2 * np.arange(1 << (a - 2), dtype=np.int64)
            b = E(t * (Hy - 1 + p * epsy).astype(np.float64))
            Fu = np.array([(E(t * (int(x) * Cy).astype(np.float64)) * b).sum() for x in Au])
            vals.append(float((np.abs(Fu) ** 2).sum()) / (Au.size * (1 << m)))
        out.append(f'{name}: {vals[0]:.3f}/{vals[1]:.3f}')
    print(f'  a = {a:2d}, m = {m:2d}: ' + '  '.join(out))

print('(iii) (C) at rational theta: |S_n(theta)| / 2^(n-1)')
rats = [Fraction(1, 4), Fraction(1, 8), Fraction(3, 8), Fraction(1, 16), Fraction(5, 32), Fraction(1, 6), Fraction(1, 12), Fraction(5, 24),
        Fraction(1, 3), Fraction(1, 7), Fraction(1, 9), Fraction(1, 5)]
hdr = ' '.join(f'{str(r):>7s}' for r in rats)
print('   n  ' + hdr)
for n in range(10, 27, 2):
    vals = g_block(n)
    print(f'  {n:2d}  ' + ' '.join(f'{abs(S_exact(vals, r)) / 2 ** (n - 1):7.4f}' for r in rats))
