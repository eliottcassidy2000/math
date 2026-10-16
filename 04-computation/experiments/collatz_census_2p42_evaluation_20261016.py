#!/usr/bin/env python3
"""Evaluation of the S38 pre-registered tests on the resumed dense census of [2^40, 2^42)
(collatz_census_2p42_preregistration_20261016.md; predictions in collatz_explicit_formula_dense_preregistration_20261016.json).

Q1  relative rms error of M1 (S37's 28-zero model refitted on the 64-point grid to 2^40) over the dense grid of
    [2^40, 2^42): pass if <= 0.03.
Q2  relative rms error of M2 (dense-grid model) over the same points: pass if <= the bound stated in the
    pre-registration (read from its JSON field 'Q2_bound').
Q3  E_X/X^2 <= 2 for X = 2^40, 2^41 (exact octave sums of H^2).
Q4  E_X/X^2 <= 1 for X = 2^40, 2^41.
Q5  mean window energy of H over window starts in [2^40, 2^41) (windows [N_j, 2 N_j) on the dense grid) below 1.
Q6  rms of the M2 residual divided by sqrt(N) over [2^40, 2^42) below the bound 'Q6_bound'.
Session opus-2026-10-16-S38.  Usage: python collatz_census_2p42_evaluation_20261016.py dense42b.npz dense42b.json
       prediction.json out.json
"""
import sys, json, math
import numpy as np

NPZ, JS, PR, OUT = sys.argv[1:5]
with np.load(NPZ) as Z:
    g = np.array(Z['grid']).astype(np.int64)
    S = np.array(Z['S']).astype(np.float64)
    H = np.array(Z['H']).astype(np.float64)
    p2 = np.array([int(x) for x in Z['preH2_str']], dtype=object)
    per = int(Z['per'])
J = json.load(open(JS))
P = json.load(open(PR))
pg = np.array(P['prediction']['grid'], dtype=np.int64)
assert np.array_equal(pg, g[(g >= 2 ** 40) & (g < 2 ** 42)]), 'prediction grid differs from the census grid'
sel = (g >= 2 ** 40) & (g < 2 ** 42)
Sa = S[sel]
res = {'n_points': int(sel.sum()), 'S_2p42': J['S'], 'H_2p42': J['H'], 'checks': J.get('checks')}
for key, name in (('M1', 'Q1'), ('M2', 'Q2')):
    pred = np.array(P['prediction'][key])
    rel = float(np.sqrt(((Sa - pred) ** 2).mean() / (Sa ** 2).mean()))
    res[name] = dict(model=key, relative_rms_error=rel)
res['Q1']['passed'] = res['Q1']['relative_rms_error'] <= 0.03
res['Q2']['bound'] = P['Q2_bound']
res['Q2']['passed'] = res['Q2']['relative_rms_error'] <= P['Q2_bound']
ex = {}
for k in (40, 41):
    s1, s2 = (int(x) for x in J['octave_sums_H_H2'][str(k)])
    ex[k] = s2 / (1 << k) ** 2
res['Q3'] = dict(values=[ex[40], ex[41]], passed=all(v <= 2 for v in ex.values()))
res['Q4'] = dict(values=[ex[40], ex[41]], passed=all(v <= 1 for v in ex.values()))
# window energies on the dense grid: windows [N_j, N_(j+per)) with N_j in [2^40, 2^41)
idx = np.nonzero(sel)[0]
gg = g[sel]
HH = H[sel]
pp = p2[sel]
w = []
for i in range(len(gg)):
    if gg[i] >= 2 ** 41:
        break
    j = i + per
    if j >= len(gg):
        break
    E = (pp[j] - int(HH[j]) ** 2) - (pp[i] - int(HH[i]) ** 2)
    w.append(float(E) / float(gg[i]) ** 2)
res['Q5'] = dict(n_windows=len(w), mean=float(np.mean(w)), max=float(np.max(w)), passed=float(np.mean(w)) < 1)
pred2 = np.array(P['prediction']['M2'])
rs = float(np.sqrt((((Sa - pred2) / np.sqrt(gg.astype(np.float64))) ** 2).mean()))
res['Q6'] = dict(rms_residual_over_sqrtN=rs, bound=P['Q6_bound'], passed=rs <= P['Q6_bound'])
res['profile'] = dict(g_mean=float((HH / np.sqrt(2 * gg)).mean()), g_rms=float(np.sqrt(((HH / np.sqrt(2 * gg)) ** 2).mean())))
json.dump(res, open(OUT, 'w'), indent=1)
for k in ('Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'profile', 'checks'):
    print(k, json.dumps(res[k]))
