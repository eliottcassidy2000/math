#!/usr/bin/env python3
"""Evaluation of the S38 second-stage pre-registered tests (collatz_census_2p43_preregistration_20261016.md) on the
checked dense census of [2^40, 2^43) (collatz_census_dense_resumable_checked_20261016.py, resumed from the 2^40 state).

R1  relative rms error of the stage-2 M2 prediction (collatz_explicit_formula_prediction_2p43_20261016.json; the 69 zeros
    below the Nyquist height refitted on [2^22, 2^42)) over the 1024 grid points of [2^42, 2^43): pass if <= 0.0035.
R2  E_X/X^2 <= 2 for X = 2^42 (exact octave sum of H^2).
R3  the boundary residue estimate at rho* = 0.99999921 + 3014.0202i before the future tail (certificate JSON, moments
    through 2^43) has modulus in [0.0098901, 0.0100899] (within 1% of 0.00999).
R4  the conditional Arb certificate at rho* with N = 2^43 and the budget (a, b) = (2, 0) gives a lower bound for |r| that
    is positive.
R5  replication: the checked run reproduces the first run of [2^40, 2^42) exactly (all 2048 grid values S, H and the exact
    prefix sums of H, H^2, and every moment row).
Session opus-2026-10-16-S38.  Usage: python collatz_census_2p43_evaluation_20261016.py run43.npz run43.json
       run43_moments.npz prediction.json runB.npz runB_moments.npz cert_rhostar_2p43.json out.json
"""
import sys, json
import numpy as np

NPZ, JS, MOM, PR, BNPZ, BMOM, CERT, OUT = sys.argv[1:9]
Z = np.load(NPZ)
g = np.array(Z['grid'], dtype=np.int64)
S = np.array(Z['S'], dtype=np.float64)
J = json.load(open(JS))
P = json.load(open(PR))
res = {'S_2p43': J['S'], 'H_2p43': J['H'], 'checks': J.get('checks')}
sel = (g >= 2 ** 42) & (g < 2 ** 43)
pg = np.array(P['grid'], dtype=np.int64)
assert np.array_equal(pg, g[sel]), 'prediction grid differs from the census grid'
Sa = S[sel]
pred = np.array(P['prediction'], dtype=np.float64)
rel = float(np.sqrt(((Sa - pred) ** 2).mean() / (Sa ** 2).mean()))
res['R1'] = dict(n_points=int(sel.sum()), relative_rms_error=rel, bound=0.0035, passed=rel <= 0.0035,
                 residual_over_sqrtN=float(np.sqrt((((Sa - pred) / np.sqrt(g[sel].astype(np.float64))) ** 2).mean())))
s1, s2 = (int(x) for x in J['octave_sums_H_H2']['42'])
ex = s2 / (1 << 42) ** 2
res['R2'] = dict(e_X=ex, passed=ex <= 2)
C = json.load(open(CERT))
rb = C['residue_before_future_tail']
re_s, im_s = rb.split('] + [')
re_v = float(re_s.strip('[ ').split(' +/- ')[0])
im_v = float(im_s.strip(']j ').split(' +/- ')[0])
mod = (re_v ** 2 + im_v ** 2) ** 0.5
res['R3'] = dict(residue=[re_v, im_v], modulus=mod, window=[0.0098901, 0.0100899], passed=0.0098901 <= mod <= 0.0100899)
row = next(r for r in C['rows'] if int(r['a']) == 2 and int(r['b']) == 0)
low = row['abs_residue_lower']
low_mid = float(low.strip('[]').split(' +/- ')[0])
res['R4'] = dict(abs_residue_lower=low, abs_residue_upper=row['abs_residue_upper'], nonzero_certified=row['nonzero_certified'],
                 passed=bool(row['nonzero_certified']) and low_mid > 0)
B = np.load(BNPZ)
gb = np.array(B['grid'], dtype=np.int64)
idx = np.searchsorted(g, gb)
same_grid = bool(np.all(idx < len(g)) and np.array_equal(g[idx], gb))
same_vals = same_grid and all(np.array_equal(np.array(Z[k])[idx], np.array(B[k])) for k in ('S', 'H', 'preH1_str', 'preH2_str'))
MA = np.load(MOM)
MB = np.load(BMOM)
ma, mb = np.array(MA['moments']), np.array(MB['moments'])
same_mom = (int(MA['block0']) == int(MB['block0'])) and np.array_equal(ma[:len(mb)], mb)
res['R5'] = dict(grid_points=int(len(gb)), grid_values_identical=bool(same_vals), moment_rows=int(len(mb)),
                 moments_identical=bool(same_mom), passed=bool(same_vals and same_mom))
json.dump(res, open(OUT, 'w'), indent=1)
for k in ('R1', 'R2', 'R3', 'R4', 'R5', 'checks'):
    print(k, json.dumps(res[k]))
