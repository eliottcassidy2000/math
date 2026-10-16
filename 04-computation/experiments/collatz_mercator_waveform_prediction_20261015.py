#!/usr/bin/env python3
"""Pre-registered out-of-sample prediction of the Collatz parity sum S(N) on [2^37, 2^40] from its Mercator waveform.

S(N) = sum_(n<=N) t(n), t(n) = (-1)^tau(n).  The amplitude note's conditional pole of F(s) = sum t(n) n^(-s) at the
zero s_0 ~ 0.99821 + 240.248i of D(s) = 1 + 2^(-s) + (3/2)^s/3 predicts S(N) ~ Re(c N^(s_0)) + smaller terms.  This
script (i) locates the zeros of D with Re s > RMIN = 0.65 and 0 < Im s < IMAX = 3000 (28 zeros; the
model was chosen on holdouts: fitted on 2^22..2^34 it predicts 2^34..2^37 to 2.05% relative rms, against 19% for the
single zero s_0 and 7.5% for the 12 zeros with Re s > 0.85), (ii) fits complex coefficients c_rho of
Re(sum_rho c_rho N^rho) by least squares to S at the grid N_j = round(2^(10+j/64)) with N in [2^FITLO, 2^37) (2^37
census extension), and (iii) writes the predicted values of S at the grid points of [2^37, 2^40) (the streaming
2^40 census was running and unread when this file was committed), plus S(2^38), S(2^39), S(2^40).
Session opus-2026-10-15-S37.  Usage: python collatz_mercator_waveform_prediction_20261015.py ext37.npz out.json
"""
import sys, json
import numpy as np
import mpmath as mp

Z = np.load(sys.argv[1])
OUT = sys.argv[2]
RMIN, IMAX, FITLO = 0.65, 3000.0, 22
grid = Z['grid'].astype(np.float64)
S = Z['classsums'].sum(axis=1).astype(np.float64)
Z.close()
th = mp.log(3) / mp.log(2)


def D(s):
    return 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3


# zeros of D: 2^-s (1 + 3^(s-1)) = -1; scan Im s on a fine grid for sign changes of the argument, refine with findroot
zeros = []
tt = np.arange(1.0, IMAX, 0.25)
for t in tt:
    for sg in (0.6, 0.75, 0.9, 0.97, 0.995):
        try:
            z = mp.findroot(D, mp.mpc(sg, t))
        except Exception:
            continue
        if RMIN < z.real < 1.0 and 0 < z.imag < IMAX and all(abs(z - w) > 1e-8 for w in zeros):
            zeros.append(z)
zeros.sort(key=lambda z: -z.real)
print('zeros of D with Re s >', RMIN, ':', [(round(float(z.real), 6), round(float(z.imag), 4)) for z in zeros], flush=True)
fit = (grid >= 2 ** FITLO) & (grid < 2 ** 37)
u = np.log(grid)
cols = []
for z in zeros:
    w = np.exp(float(z.real) * u) * np.exp(1j * float(z.imag) * u)
    cols += [w.real, -w.imag]
A = np.stack(cols, axis=1)
coef, *_ = np.linalg.lstsq(A[fit], S[fit], rcond=None)
pred_all = A @ coef
res = S[fit] - pred_all[fit]
rel_in = float(np.sqrt((res ** 2).mean()) / np.sqrt((S[fit] ** 2).mean()))
cz = [complex(coef[2 * i], coef[2 * i + 1]) for i in range(len(zeros))]
print('fitted |c_rho|:', [f'{abs(c):.4e}' for c in cz], '; in-sample relative rms residual', f'{rel_in:.4f}', flush=True)
raw = sorted(set(int(round(2 ** (10 + j / 64))) for j in range(27 * 64, 30 * 64 + 1)))
gnew = np.array([g for g in raw if 2 ** 37 <= g <= 2 ** 40], dtype=np.float64)
un = np.log(gnew)
cn = []
for z in zeros:
    w = np.exp(float(z.real) * un) * np.exp(1j * float(z.imag) * un)
    cn += [w.real, -w.imag]
pn = np.stack(cn, axis=1) @ coef
key = {int(2 ** k): float(pn[np.nonzero(gnew == 2 ** k)[0][0]]) for k in (37, 38, 39, 40)}
print('predictions: S(2^37) (census value -152262364), S(2^38), S(2^39), S(2^40):', {k: f'{v:.4e}' for k, v in key.items()}, flush=True)
json.dump(dict(zeros=[[str(z.real), str(z.imag)] for z in zeros], fit_range=[2 ** FITLO, 2 ** 37],
               coefficients=[[c.real, c.imag] for c in cz], insample_relative_rms_residual=rel_in,
               prediction_grid=[int(g) for g in gnew], prediction_S=[float(x) for x in pn],
               predicted_S_at_powers=key,
               test='PRE-REGISTERED: relative rms error of the predicted S over the grid points of [2^37, 2^40) is at most 0.03 (holdout value on 2^34..2^37: 0.0205); the single-zero model is expected near 0.19'),
          open(OUT, 'w'), indent=1)
