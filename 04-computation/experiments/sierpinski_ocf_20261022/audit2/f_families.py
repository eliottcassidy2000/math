# Audit2 (D): the families of HYP-9365 -- R_n, P_q (doubly regular), ternary towers S_k -- and the closed-form limits.
#   exact y = n(1 - kappa/rho) from exact H (audit2 recomputations where available, else the repo values, marked),
#   two-loop y_pred (e_predictor.two_loop_circ, FFT; dense cross-check),
#   closed forms: doubly regular y_pred(n) as an exact rational function of n (sympy), and its 1/n expansion;
#   R_n: n F8 -> -2 tan1 tanh1 - (3/2) tan^2 1, E0 -> 1 + tan^2 1 + tan1 tanh1, n Theta -> 0 (rate), rho_T(R_n) expansion;
#   towers: y_pred to k = 10, Richardson.
import sys, re
import numpy as np
import sympy as sp
from mpmath import mp, mpf, mpc, cosh, cos, tan, tanh, factorial as mfac, matrix, lu_solve, re as mre, e as me
sys.path.insert(0, '.')
from e_predictor import two_loop_circ, two_loop_dense, rho_circ_exact
from common import circulant_A, skew, rot_A, paley_A
mp.dps = 40

def parse_H(path, pat):
    d = {}
    for line in open(path):
        m = re.match(pat, line.strip())
        if m: d[int(m.group(1))] = int(m.group(2))
    return d

Hrot = parse_H('d_rot_paley_H.out', r'R_(\d+) H=(\d+)')
Hrot.update(parse_H('d_rot_paley_H_big.out', r'R_(\d+) H=(\d+)'))
src = {n: 'audit2' for n in Hrot}
repo_rot = {31: 20351875726012821643802839, 33: 5394851076280218488551446699}
for n, v in repo_rot.items():
    if n not in Hrot: Hrot[n] = v; src[n] = 'repo'
    else: assert Hrot[n] == v, (n, Hrot[n], v)
Hp = parse_H('d_rot_paley_H.out', r'P_(\d+) H=(\d+)')
Hp.update(parse_H('d_rot_paley_H_big.out', r'P_(\d+) H=(\d+)'))
srcp = {q: 'audit2' for q in Hp}
for q, v in {27: 418425945688463432031, 31: 19876729650849384218015865}.items():
    if q not in Hp: Hp[q] = v; srcp[q] = 'repo'
    else: assert Hp[q] == v

def kap(H, n): return mpf(H) * mpf(2) ** (n - 1) / mfac(n)
def rho_rot(n): return (mpf(n + 1) ** n + mpf(n - 1) ** n) / (2 * mre(mpc(n, 1) ** n))
L = cosh(1) / cos(1)

def extrap(pts):
    k = len(pts); A = matrix(k, k); b = matrix(k, 1)
    for i, (n, v) in enumerate(pts):
        for j in range(k): A[i, j] = mpf(1) / mpf(n) ** j
        b[i] = v
    return lu_solve(A, b)[0]

print('R_n: n, H source, kappa, rho_T (closed form), y, y_pred, y - y_pred, n(y - y_pred), [n F8, n Theta, E0]')
ys = []
for n in sorted(Hrot):
    if n < 9: continue
    r = rho_rot(n)
    assert abs(r - rho_circ_exact(n, range(1, (n - 1) // 2 + 1))) < mpf(10) ** -30
    k = kap(Hrot[n], n); y = float(n * (1 - k / r))
    t = two_loop_circ(n, range(1, (n - 1) // 2 + 1))
    if n <= 15:
        t2 = two_loop_dense(skew(rot_A(n))); assert abs(t2['y_pred'] - t['y_pred']) < 1e-10
    ys.append((n, y))
    print(f"  n={n:2d} [{src[n]}] kappa={mp.nstr(k, 12)} rho={mp.nstr(r, 12)} y={y:.5f} y_pred={t['y_pred']:.5f} diff={y - t['y_pred']:+.5f} "
          f"n*diff={n * (y - t['y_pred']):+.3f} [{t['nF8']:+.4f}, {t['nTH']:+.4f}, {t['E0']:.4f}]")
print('  polynomial-in-1/n extrapolation of y(R_n) on the last k points:',
      ', '.join(f"k={k}: {mp.nstr(extrap([(n, mpf(v)) for n, v in ys[-k:]]), 6)}" for k in range(2, 7)))
for r in (1, 3):
    sub = [(n, v) for n, v in ys if n % 4 == r]
    print(f"  residue n = {r} mod 4, last k points:", ', '.join(f"k={k}: {mp.nstr(extrap([(n, mpf(v)) for n, v in sub[-k:]]), 6)}" for k in range(2, len(sub) + 1)))

t1, th1 = tan(1), tanh(1)
yinf = t1 * th1 + t1 ** 2 / 2 - 1
print(f"\nR_n closed forms: y_inf = tan1 tanh1 + tan^2(1)/2 - 1 = {mp.nstr(yinf, 12)}; n F8 -> {mp.nstr(-2 * t1 * th1 - 1.5 * t1 ** 2, 10)}; "
      f"E0 -> {mp.nstr(1 + t1 ** 2 + t1 * th1, 10)}")
for n in (101, 1001, 10001, 100001):
    t = two_loop_circ(n, range(1, (n - 1) // 2 + 1))
    print(f"  n={n}: y_pred={t['y_pred']:.7f} nF8={t['nF8']:+.7f} E0={t['E0']:.7f} nTheta={t['nTH']:+.3e} (n * nTheta = {n * t['nTH']:.4f})")
print('  rho_T(R_n) = (cosh1/cos1)(1 - c_n/n): c_n =', ', '.join(f"{n}: {mp.nstr(n * (1 - rho_rot(n) / L), 8)}" for n in (33, 101, 1001, 10001)))
print('  => conditional on HYP-9365, kappa(R_n) = (cosh1/cos1)(1 - (1 + y_inf)/n + O(n^-2)), 1 + y_inf =', mp.nstr(1 + yinf, 10))

# doubly regular: exact rational y_pred(n)
n = sp.symbols('n', positive=True)
nF8 = -2 * (n - 1) / (n + 1) - sp.Rational(3, 2)
nTH = (2 * n - 2) / (n + 1) ** 2 - sp.Rational(2, 3) / n + sp.Rational(2, 3) / (n * (n - 1) ** 2)
D1 = 2 * n / (n + 1)
D2 = 2 * n / (n - 1) + 2 * n * (n - 1) / (n + 1) ** 2
E0 = (D2 + D1 ** 2 - D1) / 2
ypd = sp.simplify(-nF8 - nTH - E0)
print('\nDoubly regular: y_pred(n) =', sp.factor(ypd))
print('  expansion in 1/n:', sp.series(ypd.subs(n, 1 / sp.Symbol('u')), sp.Symbol('u'), 0, 3))
for q in (19, 23, 27, 31):
    val = float(ypd.subs(n, q))
    chk = two_loop_dense(skew(paley_A(q)))['y_pred'] if q != 27 else float('nan')
    print(f"  q={q}: closed-form y_pred = {val:.10f} (dense on cyclotomic P_q: {chk:.10f})")
print('Paley P_q: q, H source, kappa, rho, y, y_pred, q(y - y_pred), kappa/(e(1 - 3/(2q))) - 1')
yp = []
for q in sorted(Hp):
    rho = (mpf(q + 1) / (q - 1)) ** ((q - 1) // 2)
    k = kap(Hp[q], q); y = float(q * (1 - k / rho)); pred = float(ypd.subs(n, q))
    yp.append((q, y))
    print(f"  q={q} [{srcp[q]}] kappa={mp.nstr(k, 12)} y={y:.5f} y_pred={pred:.5f} q*diff={q * (y - pred):+.3f} "
          f"kappa/(e(1-3/(2q)))-1 = {mp.nstr(k / (me * (1 - mpf(3) / (2 * q))) - 1, 4)}")
print('  extrapolation of y(P_q), last k points:', ', '.join(f"k={k}: {mp.nstr(extrap([(q, mpf(v)) for q, v in yp[-k:]]), 6)}" for k in range(2, len(yp) + 1)))

# towers
def tower_conn(N):
    out = []
    for z in range(1, N):
        t = z
        while t % 3 == 0: t //= 3
        if t % 3 == 1: out.append(z)
    return out
Ht = {}
for line in open('../ternary_fast_k6.out'):
    m = re.match(r'k=(\d+) n=(\d+): H\(S_k\) = (\d+)', line)
    if m: Ht[int(m.group(2))] = int(m.group(3))
print('\nTowers S_k: n, exact y (repo H), y_pred, n(y - y_pred)')
for N in sorted(Ht):
    if N < 9: continue
    t = two_loop_circ(N, tower_conn(N))
    r = rho_circ_exact(N, tower_conn(N)) if N <= 81 else mpf(t['rho'])
    y = float(N * (1 - kap(Ht[N], N) / r))
    print(f"  n={N}: y={y:.6f} y_pred={t['y_pred']:.6f} n*diff={N * (y - t['y_pred']):+.4f} [nF8 {t['nF8']:+.5f}, nTheta {t['nTH']:+.5f}, E0 {t['E0']:.5f}]")
ypt = []
for k in range(2, 11):
    N = 3 ** k; t = two_loop_circ(N, tower_conn(N)); ypt.append((N, t['y_pred']))
print('  y_pred(S_k), k = 2..10:', ', '.join(f"{v:.7f}" for _, v in ypt))
r1 = [(3 * ypt[i + 1][1] - ypt[i][1]) / 2 for i in range(len(ypt) - 1)]
print('  Richardson (1/n, ratio 3):', ', '.join(f"{v:.7f}" for v in r1[-3:]))
