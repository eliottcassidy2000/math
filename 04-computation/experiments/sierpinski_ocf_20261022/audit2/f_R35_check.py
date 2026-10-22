# Audit2 (D), bonus: the session's new H(R_35) (sierpinski_ocf_20261022/R35.out, ie_ham64, NOT recomputed here) as an
# out-of-sample point for HYP-9365 (b): y(R_35) against y_pred(R_35), the drift of n(y - y_pred), and the residue-class
# extrapolation of y(R_n) (n = 3 mod 4) with the new point.
import re, sys
from mpmath import mp, mpf, mpc, factorial as mfac, matrix, lu_solve, re as mre
sys.path.insert(0, '.')
from e_predictor import two_loop_circ
mp.dps = 40
H = {}
for f in ('d_rot_paley_H.out', 'd_rot_paley_H_big.out'):
    for line in open(f):
        m = re.match(r'R_(\d+) H=(\d+)', line.strip())
        if m: H[int(m.group(1))] = int(m.group(2))
H[33] = 5394851076280218488551446699
m = re.search(r'H=(\d+)', open('../R35.out').read()); H[35] = int(m.group(1))
def rho(n): return (mpf(n + 1) ** n + mpf(n - 1) ** n) / (2 * mre(mpc(n, 1) ** n))
def y(n): return n * (1 - mpf(H[n]) * mpf(2) ** (n - 1) / mfac(n) / rho(n))
for n in (31, 33, 35):
    t = two_loop_circ(n, range(1, (n - 1) // 2 + 1))
    print(f"R_{n}: kappa = {mp.nstr(mpf(H[n]) * 2 ** (n - 1) / mfac(n), 12)}  y = {mp.nstr(y(n), 8)}  y_pred = {t['y_pred']:.6f}  n(y - y_pred) = {n * (float(y(n)) - t['y_pred']):+.4f}")
def extrap(ns):
    k = len(ns); A = matrix(k, k); b = matrix(k, 1)
    for i, n in enumerate(ns):
        for j in range(k): A[i, j] = mpf(1) / mpf(n) ** j
        b[i] = y(n)
    return lu_solve(A, b)[0]
sub = [n for n in sorted(H) if n % 4 == 3 and n >= 11]
print('residue 3 mod 4 fits of y(R_n), last k points:', ', '.join(f"k={k}: {mp.nstr(extrap(sub[-k:]), 6)}" for k in range(2, len(sub) + 1)))
