#!/usr/bin/env python3
"""Least-squares residues of F(s) = sum t(n) n^(-s) from S(N) on a log-uniform grid with an irrational number of points
per octave (1000 sqrt 2), where the Mercator family of zeros of D does not alias onto itself (THM-4659 statement 5),
compared with the boundary residue estimates from the exact 2^40 block moments.

Model: S(N) ~ Re sum_rho c_rho N^rho over the zeros with 1/2 < Re rho < 1 and Im rho < TMAX (below the Nyquist height
of the grid), fitted on [2^FITLO, 2^LOGMAX); |r_fit| = |c_rho| |rho| / 2.  Also a holdout (fit to 2^(LOGMAX-2), test
the last two octaves).  NUMERICAL.
Session opus-2026-10-16-S38.  Usage: python collatz_irrational_grid_residue_fit_20261016.py irr.npz zeros.json FITLO
       TMAX out.json
"""
import sys, json, math
import numpy as np

NPZ, ZJ, FITLO, TMAX, OUT = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
with np.load(NPZ) as Z:
    g = np.array(Z['grid']).astype(np.float64)
    S = np.array(Z['S']).astype(np.float64)
    per = float(Z['per'])
    logmax = int(Z['logmax'])
rows = json.load(open(ZJ))['zeros_and_residues']
ZS = [complex(r['re'], r['im']) for r in rows if 0.5 < r['re'] < 1 and r['im'] < TMAX]
BR = {complex(r['re'], r['im']): r['abs_r'] for r in rows}
nyq = per / 2 * 2 * math.pi / math.log(2)


def design(N, zs):
    lnN = np.log(N)
    cols = []
    for z in zs:
        w = np.exp(z.real * lnN) * np.exp(1j * z.imag * lnN)
        cols += [w.real, -w.imag]
    return np.stack(cols, axis=1)


def fit(mask):
    A = design(g[mask], ZS)
    sc = np.abs(A).max(axis=0)
    c, *_ = np.linalg.lstsq(A / sc, S[mask], rcond=None)
    return c / sc


hf = (g >= 2 ** FITLO) & (g < 2 ** (logmax - 2))
ht = (g >= 2 ** (logmax - 2)) & (g < 2 ** logmax)
c = fit(hf)
r = S[ht] - design(g[ht], ZS) @ c
hold = float(np.sqrt((r ** 2).mean() / (S[ht] ** 2).mean()))
full = (g >= 2 ** FITLO) & (g < 2 ** logmax)
c = fit(full)
res = S[full] - design(g[full], ZS) @ c
rs = float(np.sqrt(((res / np.sqrt(g[full])) ** 2).mean()))
out = dict(per=per, nyquist_height=nyq, n_zeros=len(ZS), fit=[2 ** FITLO, 2 ** logmax], holdout_rel=hold,
           residual_over_sqrtN=rs, rows=[])
print(f'grid density {per:.4f} per octave (Nyquist height {nyq:.1f}); {len(ZS)} zeros below {TMAX:.0f}; holdout {hold:.5f}; '
      f'full-fit residual/sqrt N {rs:.3f}', flush=True)
for i, z in enumerate(ZS):
    cc = complex(c[2 * i], c[2 * i + 1])
    rf = abs(cc) * abs(z) / 2
    out['rows'].append(dict(re=z.real, im=z.imag, abs_r_fit=rf, abs_r_boundary=BR[z], phase_fit=float(np.angle(cc * z / 2))))
fam = [row for row in out['rows'] if row['re'] > 0.98]
print('Mercator family (Re > 0.98): Im, |r_fit|, |r_boundary|')
for row in fam:
    print(f"   {row['im']:9.3f}  {row['abs_r_fit']:.5f}  {row['abs_r_boundary']:.5f}", flush=True)
big = sorted(out['rows'], key=lambda x: -x['abs_r_boundary'] * 2 / math.hypot(x['re'], x['im']) * 2 ** (38 * (x['re'] - 1)))[:8]
print('largest modes:')
for row in big:
    print(f"   {row['re']:.5f}+{row['im']:.2f}i  {row['abs_r_fit']:.5f}  {row['abs_r_boundary']:.5f}", flush=True)
json.dump(out, open(OUT, 'w'), indent=1)
