#!/usr/bin/env python3
"""Evaluation of the pre-registered tests P1-P6 (collatz_census_2p40_preregistration_20261015.md) on the exact census
to 2^40 (collatz_census_extension_2p40_20261015.py), plus the full octave tables 2^14..2^40.

P1  relative rms error of the 28-zero explicit-formula prediction of S over the grid points of [2^37, 2^40) <= 0.03;
    reference models refitted on [2^22, 2^37): the single zero s_0 and the 12 zeros with Re > 0.85.
P2  E_X/X^2 <= 2 for X = 2^37, 2^38, 2^39 (exact octave sums of H^2).
P3  E_X/X^2 <= 1 for the same windows.
P4  mean of the three E_X/X^2 below 0.25.
P5  mean square of the normalized channels g_(q,r) = H_(q,r)/sqrt(6N/3^q), averaged over octaves 37..39 and the
    channels of a level, below 0.8 for each level q = 3, 4, 5.
P6  octave mean of g = H/sqrt(2N) positive in at least two of the octaves 37, 38, 39.
Session opus-2026-10-15-S37.  Usage: python collatz_census_2p40_evaluation_20261015.py census40.npz census40.json
       prediction.json out.json
"""
import sys, json
import numpy as np
import mpmath as mp

with np.load(sys.argv[1]) as Z:
    grid = np.array(Z['grid']).astype(np.int64)
    C = np.array(Z['classsums']).astype(np.float64)
    Hg = np.array(Z['H']).astype(np.float64)
J = json.load(open(sys.argv[2]))
PR = json.load(open(sys.argv[3]))
OUT = sys.argv[4]
N = grid.astype(np.float64)
u = np.log2(N)
S = C.sum(axis=1)
res = {'S_2p40': J['S'], 'H_2p40': J['H'], 'reproduces_2p37': J['reproduces_2p37']}

# P1
pg = np.array(PR['prediction_grid'], dtype=np.int64)
ps = np.array(PR['prediction_S'])
sel = (pg >= 2 ** 37) & (pg < 2 ** 40)
idx = np.searchsorted(grid, pg[sel])
assert np.array_equal(grid[idx], pg[sel])
act = S[idx]
rel = float(np.sqrt(((act - ps[sel]) ** 2).mean() / (act ** 2).mean()))
res['P1'] = dict(relative_rms_error=rel, passed=rel <= 0.03, n_points=int(sel.sum()))
for k in (37, 38, 39):
    i = int(np.nonzero(grid == 2 ** k)[0][0])
    res['P1'][f'S(2^{k})'] = dict(actual=int(S[i]), predicted=PR['predicted_S_at_powers'][str(2 ** k)])
res['P1']['S(2^40)'] = dict(actual=J['S'], predicted=PR['predicted_S_at_powers'][str(2 ** 40)])


def D(s):
    return 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3


zeros = [mp.mpc(mp.mpf(a), mp.mpf(b)) for a, b in PR['zeros']]
lnN = np.log(N)
for label, zs in (('single_zero_s0', [z for z in zeros if z.real > 0.995]), ('zeros_Re_gt_0.85', [z for z in zeros if z.real > 0.85]),
                  ('zeros_Re_gt_0.65_refit_check', zeros)):
    cols = []
    for z in zs:
        w = np.exp(float(z.real) * lnN) * np.exp(1j * float(z.imag) * lnN)
        cols += [w.real, -w.imag]
    A = np.stack(cols, axis=1)
    fit = (N >= 2 ** 22) & (N < 2 ** 37)
    te = (N >= 2 ** 37) & (N < 2 ** 40)
    coef, *_ = np.linalg.lstsq(A[fit], S[fit], rcond=None)
    r = S[te] - A[te] @ coef
    res['P1'][label] = dict(n_zeros=len(zs), relative_rms_error=float(np.sqrt((r ** 2).mean() / (S[te] ** 2).mean())))
# fit on everything to 2^40 for the record (amplitude of s_0 and the implied |r| = |c_0| |s_0| / 2)
cols = []
for z in zeros:
    w = np.exp(float(z.real) * lnN) * np.exp(1j * float(z.imag) * lnN)
    cols += [w.real, -w.imag]
A = np.stack(cols, axis=1)
fit = (N >= 2 ** 22)
coef, *_ = np.linalg.lstsq(A[fit], S[fit], rcond=None)
c0 = complex(coef[0], coef[1])
res['fit_to_2p40'] = dict(abs_c0=abs(c0), implied_abs_residue=abs(c0) * abs(complex(zeros[0])) / 2,
                          relative_rms_last_octave=float(np.sqrt(((S - A @ coef)[(N >= 2 ** 39)] ** 2).mean() / (S[N >= 2 ** 39] ** 2).mean())))

# P2-P4
ex = {}
for k in range(20, 40):
    s1, s2 = (int(x) for x in J['octave_sums_H_H2'][str(k)])
    X = 2 ** k
    ex[k] = dict(e=s2 / X ** 2, mean_part=(s1 / X) ** 2 / X, mean_over_sqrtX=s1 / X / np.sqrt(X))
res['dyadic_energy'] = {str(k): v for k, v in ex.items()}
e3 = [ex[k]['e'] for k in (37, 38, 39)]
res['P2'] = dict(values=e3, passed=all(x <= 2 for x in e3))
res['P3'] = dict(values=e3, passed=all(x <= 1 for x in e3))
res['P4'] = dict(mean=float(np.mean(e3)), passed=float(np.mean(e3)) < 0.25)

# P5 and level table
lev = {0: S[:, None]}
for q in range(1, 6):
    m = 3 ** q
    lev[q] = np.stack([C[:, a::m].sum(axis=1) for a in range(m)], axis=1)
levtab = {}
for k in range(14, 40):
    s = (u >= k) & (u < k + 1)
    row = []
    for q in range(1, 6):
        m, mp_ = 3 ** q, 3 ** (q - 1)
        G = np.stack([(3 * lev[q][s][:, r] - lev[q - 1][s][:, r % mp_]) / np.sqrt(6 * N[s] / m) for r in range(m)], axis=1)
        row.append(float((G ** 2).mean()))
    levtab[k] = row
res['level_energy_by_octave'] = {str(k): v for k, v in levtab.items()}
p5 = [float(np.mean([levtab[k][q - 1] for k in (37, 38, 39)])) for q in (3, 4, 5)]
res['P5'] = dict(level3_4_5=p5, insample_14_33=[float(np.mean([levtab[k][q - 1] for k in range(14, 34)])) for q in (3, 4, 5)],
                 passed=all(x < 0.8 for x in p5))

# P6 and profile
g = Hg / np.sqrt(2 * N)
om = {k: float(g[(u >= k) & (u < k + 1)].mean()) for k in range(14, 40)}
res['profile_octave_means'] = {str(k): v for k, v in om.items()}
pos = sum(1 for k in (37, 38, 39) if om[k] > 0)
res['P6'] = dict(values=[om[37], om[38], om[39]], positive=pos, passed=pos >= 2)

# the five HYP-9282 channels at octaves 37..39 (its pre-registration covered 2^34..2^40)
chan = {}
for (q, r) in ((3, 18), (3, 9), (5, 155), (5, 78), (5, 59)):
    m, mp_ = 3 ** q, 3 ** (q - 1)
    v = (3 * lev[q][:, r] - lev[q - 1][:, r % mp_]) / np.sqrt(6 * N / m)
    chan[f'H({q},{r})'] = {str(k): float(v[(u >= k) & (u < k + 1)].mean()) for k in range(33, 40)}
res['hyp9282_channels_33_39'] = chan
json.dump(res, open(OUT, 'w'), indent=1)
for key in ('S_2p40', 'H_2p40', 'reproduces_2p37', 'P1', 'fit_to_2p40', 'P2', 'P3', 'P4', 'P5', 'P6'):
    print(key, json.dumps(res[key]))
print('dyadic e_X:', {k: round(v['e'], 4) for k, v in ex.items()})
print('level energies by octave (levels 1..5):')
for k in range(26, 40):
    print(k, [round(x, 3) for x in levtab[k]], 'g mean', round(om[k], 3))
print('HYP-9282 channels:', json.dumps(chan))
