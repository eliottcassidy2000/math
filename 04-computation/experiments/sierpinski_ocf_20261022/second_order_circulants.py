# Second-order (three-loop) prediction on every circulant class at n = 21, 25 and on Paley P_q for larger q (the DRT
# constant c_2 = lim n^2 (kappa/e - 1 + 3/(2n))).  mac-mini-2026-10-22-sierp.
import re, glob, math, sys
import numpy as np
from collections import defaultdict
sys.path.insert(0, '.')
from second_order import predict2
from mpmath import mpf, factorial, e as E
rows = defaultdict(dict)
for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
    for line in open(f):
        m = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if m: rows[int(m.group(1))][tuple(sorted(int(x) for x in m.group(2).split(',')))] = int(m.group(4))
for line in open('../lattice_bridges_20261021/circ25_H.txt'):
    m = re.match(r'S=\{([\d,]+)\} H=(\d+)', line)
    if m: rows[25][tuple(int(x) for x in m.group(1).split(','))] = int(m.group(2))
def circ(n, S):
    A = np.zeros((n, n))
    for u in range(n):
        for s in S: A[u, (u + s) % n] = 1
    return A - A.T
for n in (21, 25):
    r1, r2 = [], []
    for S, H in rows[n].items():
        rho, first, second, X2 = predict2(circ(n, S))
        ex = float(mpf(H) * mpf(2) ** (n - 1) / factorial(n)) / rho - 1
        r1.append(ex - first); r2.append(ex - first - second)
    r1, r2 = np.array(r1), np.array(r2)
    print(f"circulants n={n} ({len(r1)} classes): residual after 1st order mean {r1.mean():+.3e} sd {r1.std():.3e} max|.| {np.abs(r1).max():.3e};"
          f"  after 2nd order mean {r2.mean():+.3e} sd {r2.std():.3e} max|.| {np.abs(r2).max():.3e}", flush=True)
print('Paley: predicted z(n) = n^2 (kappa_pred/e - 1 + 3/(2n)), first + second order')
for q in (19, 23, 31, 43, 47, 59):
    sq = {(t * t) % q for t in range(1, q)}
    rho, first, second, X2 = predict2(circ(q, sorted(sq)))
    kp = rho * (1 + first + second)
    z = q * q * (kp / float(E) - 1 + 1.5 / q)
    print(f"  q={q}: rho = {rho:.10f}, kappa_pred = {kp:.10f}, z_pred = {z:.4f}, X_2 q^2 = {X2 * q * q:+.4f}", flush=True)
