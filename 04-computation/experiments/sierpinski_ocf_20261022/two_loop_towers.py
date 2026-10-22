# Two-loop prediction for the ternary towers S_k = C_3^[k] = Cay(Z/3^k, unit part = 1 mod 3), k = 2..6 (n <= 729),
# against the exact counts in ternary_fast_k6.out.  The towers are circulants, so dumbbells vanish; the theta term
# survives because the outlier frequencies are additively closed (multiples of 3^j).  mac-mini-2026-10-22-sierp.
import re, math, sys
import numpy as np
from mpmath import mp, mpf, factorial
sys.path.insert(0, '.')
from two_loop_predictor import predict
mp.dps = 30
H = {}
for line in open('ternary_fast_k6.out'):
    m = re.match(r'k=(\d+) n=(\d+): H\(S_k\) = (\d+)', line)
    if m: H[int(m.group(1))] = (int(m.group(2)), int(m.group(3)))
def conn(n):
    S = []
    for z in range(1, n):
        t = z
        while t % 3 == 0: t //= 3
        if t % 3 == 1: S.append(z)
    return S
print(' k    n   exact y      y_pred     diff     n*diff   [n F8, n Theta, E0 C(L,2)]')
for k in sorted(H):
    n, h = H[k]
    if n < 9: continue
    rho, yp, f8, th, el = predict(n, conn(n))
    kap = mpf(h) * mpf(2) ** (n - 1) / factorial(n)
    y = float(n * (1 - kap / mpf(rho)))
    print(f"{k:2d} {n:4d}  {y:.6f}  {yp:.6f}  {y - yp:+.6f}  {n * (y - yp):+.4f}   [{f8:+.5f}, {th:+.5f}, {el:.5f}]", flush=True)
