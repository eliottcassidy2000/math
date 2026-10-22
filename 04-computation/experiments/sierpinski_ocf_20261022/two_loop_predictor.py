# Two-loop prediction of the self-blow-up residual y = n(1 - kappa/rho) for circulant (vertex-transitive) regular
# tournaments, from the linear-forest cluster expansion:
#   kappa = sum_L N_L/(n)_L,   log Phi = log rho (cycles, excess 0) + X_1 (excess 1) + O(n^-2),
#   y_pred = -n X_1 - E_0[C(L,2)],
#   X_1 = sum_v [2 G_vv H_vv - (3/2) G_vv^2] + sum_(u,v) G_uv (H_uv^2 - O_uv^2 - (2/3) G_uv^2),   (dumbbells = 0)
#   G = -Shat^2 (I + Shat^2)^-1, H = Shat^2 (I - Shat^2)^-1, O = Shat (I - Shat^2)^-1,  Shat = S/n,
#   E_0[C(L,2)] = (D2 + D1^2 - D1)/2, D1 = sum_k 2x^2/(1-x^4), D2 = sum_k 4x^2(1+x^4)/(1-x^4)^2, x_k = sigma_k/n.
# Kernels from two_loop_kernels.py (K8 = 8(e_a+e_b) - 12 e_a e_b; Kth = 4(e_a+e_b+e_c) - 8 e_a e_b e_c for even
# lengths, -4 e_c for (odd, odd, even c); e_a = (-1)^(a/2)).  Exact data: circulant classes n <= 25
# (tournament_pseudoprimes_20261021, lattice_bridges_20261021), R_n to 33, P_q to 31.  mac-mini-2026-10-22-sierp.
import re, glob, math, sys
import numpy as np
from math import factorial
from collections import defaultdict

def predict(n, S):
    k = np.arange(n)
    x = np.array([2 * sum(math.sin(2 * math.pi * kk * s / n) for s in S) / n for kk in k])
    lam = 1j * x
    g = x ** 2 / (1 - x ** 2); h = -x ** 2 / (1 + x ** 2); o = 1j * x / (1 + x ** 2)
    # first rows c_f(w) = (1/n) sum_k f_k w^(-k w); real parts (matrices are real)
    W = np.exp(-2j * math.pi * np.outer(np.arange(n), k) / n)
    gr = (W @ g).real / n; hr = (W @ h).real / n; orow = (W @ o).real / n
    F8 = n * (2 * gr[0] * hr[0] - 1.5 * gr[0] ** 2)
    TH = n * np.sum(gr * (hr ** 2 - orow ** 2 - (2 / 3) * gr ** 2))
    X1 = F8 + TH
    D1 = np.sum(2 * x ** 2 / (1 - x ** 4)); D2 = np.sum(4 * x ** 2 * (1 + x ** 4) / (1 - x ** 4) ** 2)
    EL = (D2 + D1 ** 2 - D1) / 2
    rho = float(np.prod((1 + x ** 2) / (1 - x ** 2)) ** 0.5)
    return rho, -n * X1 - EL, n * F8, n * TH, EL

if __name__ == '__main__':
    rows = defaultdict(dict)
    for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
        for line in open(f):
            m = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
            if m: rows[int(m.group(1))][tuple(sorted(int(x) for x in m.group(2).split(',')))] = int(m.group(4))
    for line in open('../lattice_bridges_20261021/circ25_H.txt'):
        m = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
        if m: rows[25][tuple(int(x) for x in m.group(1).split(','))] = int(m.group(2))
    print('All circulant classes: exact y = n(1 - kappa/rho) against two-loop y_pred')
    for n in sorted(rows):
        if n < 9: continue
        Y, P = [], []
        for S, H in rows[n].items():
            rho, yp, f8, th, el = predict(n, S)
            kap = H * 2 ** (n - 1) / factorial(n)
            Y.append(n * (1 - kap / rho)); P.append(yp)
        Y, P = np.array(Y), np.array(P)
        rk = lambda v: np.argsort(np.argsort(v))
        print(f" n={n:2d} ({len(Y):3d} classes): corr(y, y_pred) = {np.corrcoef(Y, P)[0, 1]:+.4f}, spearman = {np.corrcoef(rk(Y), rk(P))[0, 1]:+.4f}, "
              f"mean(y - y_pred) = {np.mean(Y - P):+.4f}, rms(y - y_pred - mean) = {np.std(Y - P):.4f}, y range [{Y.min():.3f}, {Y.max():.3f}]", flush=True)
    Hrot = {25: 2418453569285650675, 27: 426990256001863046739, 29: 87133267870433699742519,
            31: 20351875726012821643802839, 33: 5394851076280218488551446699}
    print('\nRotational R_n: n, exact y, y_pred, (y - y_pred), n*(y - y_pred)')
    for n, H in sorted(Hrot.items()):
        rho, yp, f8, th, el = predict(n, range(1, (n - 1) // 2 + 1))
        y = n * (1 - H * 2 ** (n - 1) / factorial(n) / rho)
        print(f"  n={n}: y={y:.5f} y_pred={yp:.5f} diff={y - yp:+.5f} n*diff={n * (y - yp):+.3f}  [nF8={f8:+.4f} nTheta={th:+.4f} E0C(L,2)={el:.4f}]")
    Hp = {19: 1172695746915, 23: 15760206976379349, 31: 19876729650849384218015865}
    print('\nPaley P_q (prime q): q, exact y, y_pred, diff, q*diff')
    for q, H in sorted(Hp.items()):
        sq = sorted({(t * t) % q for t in range(1, q)})
        rho, yp, f8, th, el = predict(q, sq)
        y = q * (1 - H * 2 ** (q - 1) / factorial(q) / rho)
        print(f"  q={q}: y={y:.5f} y_pred={yp:.5f} diff={y - yp:+.5f} q*diff={q * (y - yp):+.3f}  [nF8={f8:+.4f} nTheta={th:+.4f} E0C(L,2)={el:.4f}]")
    print('\nLimits: R_n at n = 2001 and Paley-like flat spectrum at q = 2003 (prime 3 mod 4):')
    for n, S in ((2001, range(1, 1001)), ):
        rho, yp, f8, th, el = predict(n, S); print(f"  R_{n}: y_pred = {yp:.6f} [nF8={f8:+.5f} nTheta={th:+.5f} E0={el:.5f}]")
