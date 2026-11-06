# Compare the bilinear bound with exact |S_n(theta)| / 2^(n-1), S_n(theta) = sum over v in [2^(n-1), 2^n) of e(theta g(v)).
import numpy as np, math
from fractions import Fraction
def gtable(N):
    g = np.zeros(N, dtype=np.int64); g[0] = 1
    for v in range(1, N):
        g[v] = g[v >> 1] if v & 1 else g[v >> 1] + (v >> 1) + ((v >> 1) & 1)
    return g
def bound(n, q):
    best = float('inf')
    for a in range(2, n):
        for L in range(1, n - a):
            l = n - a - L
            if l < 1: continue
            R = 2 * (2**L / q + 1) * math.sqrt((3 + 2**(3 - a) * q) * (3 + 2**(2 - l) * q))
            best = min(best, 2**(-L / 2) * math.sqrt(1 + R))
    return best
NMAX = 24
G = gtable(1 << NMAX)
def S(n, theta):
    vals = G[1 << (n - 1): 1 << n].astype(np.float64)
    return abs(np.exp(2j * np.pi * theta * vals).sum()) / (1 << (n - 1))
viol = 0; rows = []
tests = [Fraction(1, 7), Fraction(1, 3), Fraction(2, 5), Fraction(3, 11), Fraction(5, 31), Fraction(7, 63), Fraction(11, 93),
         Fraction(1, 127), Fraction(17, 255), Fraction(100, 1023), Fraction(1, 4095), Fraction(1234, 65537), Fraction(3, 17)]
for n in (12, 16, 20, 24):
    for th in tests:
        two = 2 * th; q = two.denominator
        if q == 1: continue
        s = S(n, float(th)); b = bound(n, q)
        rows.append((n, str(th), q, s, b)); viol += (s > b)
    for th in ((math.sqrt(5) - 1) / 2, math.sqrt(2) - 1, math.pi - 3, math.e - 2):
        # best approximation of 2*theta with q <= 2^(n/2): continued fractions
        fr = Fraction(2 * th).limit_denominator(1 << (n // 2))
        q = fr.denominator; s = S(n, th); b = bound(n, q)
        rows.append((n, f'{th:.6f}', q, s, b)); viol += (s > b)
for r in rows: print(f"n={r[0]:2d} theta={r[1]:>12} q={r[2]:>6} exact={r[3]:.4f} bound={r[4]:.4f}")
print('violations:', viol)
