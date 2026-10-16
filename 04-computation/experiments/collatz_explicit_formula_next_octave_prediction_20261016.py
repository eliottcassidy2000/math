#!/usr/bin/env python3
"""Refit the explicit-formula model M2 of the S38 pre-registration on all dense data so far and predict S on the dense
grid of the next range [2^A, 2^B) (third pre-registered test; written before the corresponding census run starts).

Inputs: dense census NPZ files of collatz_census_dense_resumable_20261016.py (consecutive ranges, same grid formula),
the zeros JSON of collatz_D_zeros_and_boundary_residues_20261016.py, and the model choice (zero set: 1/2 < Re < 1,
Im < T).  Least squares as in collatz_explicit_formula_dense_preregistration_20261016.py, fitted on [2^22, 2^A);
a holdout (fit [2^22, 2^(A-2)), test [2^(A-2), 2^A)) is reported first.
Session opus-2026-10-16-S38.  Usage: python collatz_explicit_formula_next_octave_prediction_20261016.py zeros.json T
       A B out.json dense1.npz [dense2.npz ...]
"""
import sys, json, math
import numpy as np

ZJ, T, A, Bk, OUT = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
files = sys.argv[6:]
g, S = [], []
per = None
for f in files:
    with np.load(f) as Z:
        g.append(np.array(Z['grid']).astype(np.float64))
        S.append(np.array(Z['S']).astype(np.float64))
        per = float(Z['per'])
g = np.concatenate(g)
S = np.concatenate(S)
assert np.all(np.diff(g) > 0)
ZS = [complex(r['re'], r['im']) for r in json.load(open(ZJ))['zeros_and_residues'] if 0.5 < r['re'] < 1 and r['im'] < T]


def design(N, zs):
    lnN = np.log(N)
    cols = []
    for z in zs:
        w = np.exp(z.real * lnN) * np.exp(1j * z.imag * lnN)
        cols += [w.real, -w.imag]
    return np.stack(cols, axis=1)


def fit(N, Y, mask):
    Am = design(N[mask], ZS)
    sc = np.abs(Am).max(axis=0)
    c, *_ = np.linalg.lstsq(Am / sc, Y[mask], rcond=None)
    return c / sc


hfit = (g >= 2 ** 22) & (g < 2 ** (A - 2))
hte = (g >= 2 ** (A - 2)) & (g < 2 ** A)
c = fit(g, S, hfit)
r = S[hte] - design(g[hte], ZS) @ c
hold = float(np.sqrt((r ** 2).mean() / (S[hte] ** 2).mean()))
hold_rs = float(np.sqrt(((r / np.sqrt(g[hte])) ** 2).mean()))
print(f'{len(ZS)} zeros (Im < {T:.1f}); holdout fit [2^22, 2^{A - 2}) test [2^{A - 2}, 2^{A}): relative rms {hold:.5f}, '
      f'residual/sqrt N {hold_rs:.3f}', flush=True)
full = (g >= 2 ** 22) & (g < 2 ** A)
c = fit(g, S, full)
raw = sorted(set(int(round(2 ** (10 + j / per))) for j in range(int((A - 10) * per), int((Bk - 10) * per) + 1)))
gnew = np.array([x for x in raw if 2 ** A <= x < 2 ** Bk], dtype=np.float64)
P = design(gnew, ZS) @ c
json.dump(dict(T=T, A=A, B=Bk, per=per, n_zeros=len(ZS), holdout_rel=hold, holdout_residual_over_sqrtN=hold_rs,
               grid=[int(x) for x in gnew], prediction=[float(x) for x in P], coef=[float(x) for x in c],
               zeros=[[z.real, z.imag] for z in ZS]), open(OUT, 'w'), indent=1)
for k in (A + 0.25, A + 0.5, A + 0.75):
    i = int(np.argmin(np.abs(np.log2(gnew) - k)))
    print(f'prediction at 2^{math.log2(gnew[i]):.4f}: {P[i]:.6e}', flush=True)
