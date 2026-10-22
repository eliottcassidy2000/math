# Audit2 (F): the n mod 4 wobble of kappa(R_n) and the alternating count F(R_n,-1)/n! (HYP-9365 (e)).
#  - F(R_n,-1) for n <= 21 recomputed with hcount.c (f_fneg*.out): sign pattern, decay per vertex;
#  - w_n = (kappa(n) - smooth fit of the other residue class at n)/kappa(n), the session's definition (cubic in 1/n
#    through the 4 nearest other-class points, interpolation only), plus two alternative smooth fits (quadratic through
#    3 nearest; quartic through 5 nearest), to measure how much the ratio w_n kappa_n/(F/n!) depends on the fit.
import re
from mpmath import mp, mpf, factorial as mfac
mp.dps = 40
F = {}
for f in ('f_fneg_small.out', 'f_fneg.out'):
    for line in open(f):
        m = re.match(r'n=(\d+) F\(T,-1\)=(-?\d+)', line.strip())
        if m: F[int(m.group(1))] = int(m.group(2))
H = {}
for f in ('d_rot_paley_H.out', 'd_rot_paley_H_big.out'):
    for line in open(f):
        m = re.match(r'R_(\d+) H=(\d+)', line.strip())
        if m: H[int(m.group(1))] = int(m.group(2))
H[33] = 5394851076280218488551446699          # repo value (R33.out), not recomputed here
kap = {n: mpf(H[n]) * 2 ** (n - 1) / mfac(n) for n in H}
print('F(R_n,-1): n, value, sign, (-1)^((n+1)/2), F/n!, |F/n!| ratio to n-2')
prev = None
for n in sorted(F):
    v = mpf(F[n]) / mfac(n)
    sgn = (F[n] > 0) - (F[n] < 0)
    rat = '' if prev is None or prev == 0 else f"{mp.nstr(abs(v / prev), 5)} (per vertex {mp.nstr(abs(v / prev) ** 0.5, 5)})"
    print(f"  n={n:2d} F={F[n]:>20d} sign={sgn:+d} expected={(-1) ** ((n + 1) // 2):+d}  F/n!={mp.nstr(v, 6):>12s}  {rat}")
    prev = v

def interp(pts, x):
    tot = mpf(0)
    for i, a in enumerate(pts):
        w = mpf(1)
        for j, b in enumerate(pts):
            if j != i: w *= (mpf(1) / x - mpf(1) / b) / (mpf(1) / a - mpf(1) / b)
        tot += w * kap[a]
    return tot

ns = sorted(n for n in kap if n >= 9)
print('\nw_n and the ratio w_n kappa_n / (F(R_n,-1)/n!) for three smooth fits of the other residue class:')
for n in ns:
    if n < 13: continue
    other = [m for m in ns if m % 4 != n % 4]
    row = []
    for k in (3, 4, 5):
        near = sorted(sorted(other, key=lambda m: abs(m - n))[:k])
        if max(near) < n or min(near) > n: row.append(None); continue
        w = (kap[n] - interp(near, n)) / kap[n]
        r = w * kap[n] / (mpf(F[n]) / mfac(n)) if n in F else None
        row.append((w, r))
    cells = []
    for k, c in zip((3, 4, 5), row):
        if c is None: cells.append(f"k={k}: --"); continue
        w, r = c
        cells.append(f"k={k}: w={mp.nstr(w, 4)}" + (f", ratio={mp.nstr(r, 4)}" if r is not None else ''))
    print(f"  n={n:2d} sign(-1)^((n+1)/2)={(-1) ** ((n + 1) // 2):+d}  " + ' | '.join(cells))
