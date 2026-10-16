#!/usr/bin/env python3
"""Explicit-formula models of the Collatz parity sum S(N) on the dense 2^40 census, holdout calibration, and the
pre-registered prediction of S on the dense grid of [2^40, 2^42).

Data: the dense census through 2^40 (collatz_census_dense_resumable_20261016.py, 1024 grid points per octave; S at
every grid point) and the 64-per-octave S37 census (collatz_census_extension_2p40_20261015.npz).  Zeros of
D(s) = 1 + 2^(-s) + (3/2)^s/3 with 1/2 < Re < 1 from collatz_D_zeros_and_boundary_residues_20261016.json.
Model: S(N) ~ Re sum_rho c_rho N^rho (least squares, one complex coefficient per zero).
  * M1: the S37 model (the 28 zeros with 0.65 < Re < 1, Im < 3000) refitted on the 64-per-octave grid of [2^22, 2^40).
  * M2 candidates on the dense grid: zeros with Re > RMIN and Im < T, T <= 4641 (the Nyquist height of 1024 points
    per octave: no aliasing), fitted on [2^22, 2^37), tested on [2^37, 2^40); the candidate with the smallest holdout
    error is M2, refitted on [2^22, 2^40).
Outputs: holdout table, the chosen models, and the predicted S at every dense grid point of [2^40, 2^42) (the grid of
the resumed census run) for M1 and M2.  NUMERICAL.
Session opus-2026-10-16-S38.  Usage: python collatz_explicit_formula_dense_preregistration_20261016.py dense40.npz
       census64.npz zeros.json out.json
"""
import sys, json, math
import numpy as np

DENSE, C64, ZJ, OUT = sys.argv[1:5]
with np.load(DENSE) as Z:
    g = np.array(Z['grid']).astype(np.float64)
    S = np.array(Z['S']).astype(np.float64)
    per = int(Z['per'])
with np.load(C64) as Z:
    g64 = np.array(Z['grid']).astype(np.float64)
    S64 = np.array(Z['classsums']).sum(axis=1).astype(np.float64)
ZS = [complex(r['re'], r['im']) for r in json.load(open(ZJ))['zeros_and_residues']]


def design(N, zs):
    lnN = np.log(N)
    cols = []
    for z in zs:
        w = np.exp(z.real * lnN) * np.exp(1j * z.imag * lnN)
        cols += [w.real, -w.imag]
    return np.stack(cols, axis=1)


def fit_eval(N, Y, zs, fitmask, testmask):
    A = design(N, zs)
    # scale columns for conditioning
    sc = np.abs(A[fitmask]).max(axis=0)
    coef, *_ = np.linalg.lstsq(A[fitmask] / sc, Y[fitmask], rcond=None)
    coef = coef / sc
    r = Y - A @ coef
    rel_in = float(np.sqrt((r[fitmask] ** 2).mean() / (Y[fitmask] ** 2).mean()))
    rel_te = float(np.sqrt((r[testmask] ** 2).mean() / (Y[testmask] ** 2).mean())) if testmask.any() else None
    rs = float(np.sqrt(((r[testmask] / np.sqrt(N[testmask])) ** 2).mean())) if testmask.any() else None
    return coef, rel_in, rel_te, rs


NYQ = (per / 2) * 2 * math.pi / math.log(2)
res = dict(per=per, nyquist_height=NYQ, candidates=[])
fit = (g >= 2 ** 22) & (g < 2 ** 37)
te = (g >= 2 ** 37) & (g < 2 ** 40)
m1z = [z for z in ZS if 0.65 < z.real < 1 and z.imag < 3000]
# M1 holdout on the coarse grid (as in S37) and on the dense grid
f64 = (g64 >= 2 ** 22) & (g64 < 2 ** 37)
t64 = (g64 >= 2 ** 37) & (g64 < 2 ** 40)
_, a1, b1, c1 = fit_eval(g64, S64, m1z, f64, t64)
_, a2, b2, c2 = fit_eval(g, S, m1z, fit, te)
res['M1_holdout'] = dict(n_zeros=len(m1z), coarse_grid=[a1, b1, c1], dense_grid=[a2, b2, c2])
print(f'M1 ({len(m1z)} zeros): holdout 2^37..2^40 coarse grid {b1:.4f} (in-sample {a1:.4f}); dense grid {b2:.4f}; '
      f'residual/sqrt N {c2:.2f}', flush=True)
best = None
for T in (3000, 3600, 4200, NYQ - 1):
    for rmin in (0.5, 0.55, 0.65):
        zs = [z for z in ZS if rmin < z.real < 1 and z.imag < T]
        _, rin, rte, rs = fit_eval(g, S, zs, fit, te)
        row = dict(T=T, rmin=rmin, n_zeros=len(zs), rel_in=rin, rel_holdout=rte, holdout_residual_over_sqrtN=rs)
        res['candidates'].append(row)
        print(f'T < {T:7.1f}, Re > {rmin}: {len(zs):3d} zeros; holdout {rte:.5f} (in-sample {rin:.5f}); residual/sqrt N {rs:.3f}',
              flush=True)
        if best is None or rte < best['rel_holdout']:
            best = row
res['M2_choice'] = best
print('M2 =', best, flush=True)
# refit on [2^22, 2^40) and predict the dense grid of [2^40, 2^42)
raw = sorted(set(int(round(2 ** (10 + j / per))) for j in range(30 * per, 32 * per + 1)))
gnew = np.array([x for x in raw if 2 ** 40 <= x < 2 ** 42], dtype=np.float64)
full = (g >= 2 ** 22) & (g < 2 ** 40)
m2z = [z for z in ZS if best['rmin'] < z.real < 1 and z.imag < best['T']]
c2, *_ = fit_eval(g, S, m2z, full, np.zeros_like(full))
P2 = design(gnew, m2z) @ c2
full64 = (g64 >= 2 ** 22) & (g64 < 2 ** 40)
c1, *_ = fit_eval(g64, S64, m1z, full64, np.zeros_like(full64))
P1 = design(gnew, m1z) @ c1
res['prediction'] = dict(grid=[int(x) for x in gnew], M1=[float(x) for x in P1], M2=[float(x) for x in P2],
                         M1_zeros=[[z.real, z.imag] for z in m1z], M2_zeros=[[z.real, z.imag] for z in m2z],
                         M1_coef=[float(x) for x in c1], M2_coef=[float(x) for x in c2])
for k in (40.5, 41, 41.5):
    i = int(np.argmin(np.abs(np.log2(gnew) - k)))
    print(f'prediction at N = {int(gnew[i])} (2^{math.log2(gnew[i]):.4f}): M1 {P1[i]:.6e}, M2 {P2[i]:.6e}', flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
