#!/usr/bin/env python3
"""S38 audit, items 4, 5, 7: the pre-registered explicit-formula models and the aliasing explanation.

(1) Re-implementation of the least-squares models of collatz_explicit_formula_dense_preregistration_20261016.py
    (S(N) ~ Re sum_rho c_rho N^rho, one complex coefficient per zero) on dense40.npz and the S37 64-point census:
    the twelve-candidate holdout table, M1, the choice of M2 and the predicted S on the 2048 grid points of [2^40, 2^42),
    compared with the committed JSON.
(2) Least-squares residues |r_fit| = |c_rho| |rho|/2 of M2 refitted on [2^22, 2^40) against the numbers quoted in
    section 3 of the S38 note and against the boundary residues.
(3) Synthetic test of the aliasing explanation: a signal S_syn(N) = Re sum (2 r_hat/rho) N^rho over ALL 295 zeros
    (boundary residues at 2^40) on the same grid; the same M2 fit; fitted vs input residues; the same on grids with
    1000 and 1000 sqrt 2 points per octave; holdout drop when the zeros with 4200 < Im < 4640 enter.
Usage: python a7_prereg_and_alias.py dense40.npz census64.npz zeros.json prereg.json
"""
import sys, json, math, os
import numpy as np

DENSE, C64, ZJ, PJ = sys.argv[1:5]
HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(DENSE)
g = np.array(d['grid']).astype(np.float64)
S = np.array(d['S']).astype(np.float64)
c64 = np.load(C64)
g64 = np.array(c64['grid']).astype(np.float64)
S64 = np.array(c64['classsums']).sum(axis=1).astype(np.float64)
rows = json.load(open(ZJ))['zeros_and_residues']
ZS = np.array([complex(r['re'], r['im']) for r in rows])
RB = np.array([complex(r['r_re'], r['r_im']) for r in rows])
P = json.load(open(PJ))
out = {}


def design(N, zs):
    lnN = np.log(N)
    W = np.exp(np.outer(lnN, zs))            # N^rho
    A = np.empty((len(N), 2 * len(zs)))
    A[:, 0::2] = W.real
    A[:, 1::2] = -W.imag
    return A


def fit(N, Y, zs, mask):
    A = design(N[mask], zs)
    sc = np.abs(A).max(axis=0)
    c, *_ = np.linalg.lstsq(A / sc, Y[mask], rcond=None)
    return c / sc


def relerr(N, Y, zs, c, mask):
    r = Y[mask] - design(N[mask], zs) @ c
    return float(np.sqrt((r ** 2).mean() / (Y[mask] ** 2).mean())), float(np.sqrt(((r / np.sqrt(N[mask])) ** 2).mean()))


NYQ = 512 * 2 * math.pi / math.log(2)
fitm = (g >= 2 ** 22) & (g < 2 ** 37)
tem = (g >= 2 ** 37) & (g < 2 ** 40)
m1z = ZS[(ZS.real > 0.65) & (ZS.imag < 3000)]
f64 = (g64 >= 2 ** 22) & (g64 < 2 ** 37)
t64 = (g64 >= 2 ** 37) & (g64 < 2 ** 40)
c = fit(g64, S64, m1z, f64)
m1c = relerr(g64, S64, m1z, c, t64)
c = fit(g, S, m1z, fitm)
m1d = relerr(g, S, m1z, c, tem)
out['M1_holdout_coarse_dense'] = [m1c, m1d]
print('(1) M1 holdout coarse', m1c, 'dense', m1d, '; JSON', P['M1_holdout']['coarse_grid'][1], P['M1_holdout']['dense_grid'][1])
cand = []
for T in (3000, 3600, 4200, NYQ - 1):
    for rmin in (0.5, 0.55, 0.65):
        zs = ZS[(ZS.real > rmin) & (ZS.imag < T)]
        c = fit(g, S, zs, fitm)
        h, hr = relerr(g, S, zs, c, tem)
        cand.append([T, rmin, len(zs), h, hr])
jc = [[x['T'], x['rmin'], x['n_zeros'], x['rel_holdout'], x['holdout_residual_over_sqrtN']] for x in P['candidates']]
out['candidates_mine'] = cand
out['candidates_max_abs_diff_holdout'] = float(max(abs(a[3] - b[3]) for a, b in zip(cand, jc)))
print('(1) candidate holdouts (mine vs JSON):')
for a, b in zip(cand, jc):
    print(f'   T<{a[0]:7.1f} Re>{a[1]}: n={a[2]} holdout {a[3]:.6f} (JSON {b[3]:.6f}) resid/sqrtN {a[4]:.3f} (JSON {b[4]:.3f})')
best = min(cand, key=lambda x: x[3])
out['M2_choice_mine'] = best
gnew = np.array(P['prediction']['grid'], dtype=np.float64)
full = (g >= 2 ** 22) & (g < 2 ** 40)
m2z = ZS[(ZS.real > best[1]) & (ZS.imag < best[0])]
c2 = fit(g, S, m2z, full)
P2 = design(gnew, m2z) @ c2
full64 = (g64 >= 2 ** 22) & (g64 < 2 ** 40)
c1 = fit(g64, S64, m1z, full64)
P1 = design(gnew, m1z) @ c1
pj1, pj2 = np.array(P['prediction']['M1']), np.array(P['prediction']['M2'])
out['prediction_M1_max_rel_diff'] = float(np.max(np.abs(P1 - pj1)) / np.sqrt(np.mean(pj1 ** 2)))
out['prediction_M2_max_rel_diff'] = float(np.max(np.abs(P2 - pj2)) / np.sqrt(np.mean(pj2 ** 2)))
out['prediction_M1_M2_rel_rms_difference'] = float(np.sqrt(np.mean((pj1 - pj2) ** 2) / np.mean(pj2 ** 2)))
print('(1) predictions reproduced: M1 max |diff|/rms', out['prediction_M1_max_rel_diff'], ' M2', out['prediction_M2_max_rel_diff'],
      '; rms(M1 - M2)/rms(M2) =', out['prediction_M1_M2_rel_rms_difference'], flush=True)
# (2) residues from M2 on [2^22, 2^40)
quoted = {240.248: (0.18398, 0.18427), 131.156: (0.5065, 0.5081), 720.740: (0.0701, 0.0648), 1572.556: (0.2036, 0.2025),
          22.253: (2.570, 2.579), 611.617: (0.1361, 0.1357), 3014.020: (0.00999, 0.00712), 3494.515: (0.0128, 0.0326),
          3975.0: (0.0135, 0.0592), 4455.5: (0.0124, 0.151)}
res2 = []
for im, (qb, qf) in quoted.items():
    i = int(np.argmin(np.abs(m2z.imag - im)))
    z = m2z[i]
    cc = complex(c2[2 * i], c2[2 * i + 1])
    rf = abs(cc) * abs(z) / 2
    j = int(np.argmin(np.abs(ZS - z)))
    res2.append([float(z.real), float(z.imag), float(abs(RB[j])), rf, qb, qf])
    print(f'(2) {z.real:.5f}+{z.imag:9.3f}i  boundary {abs(RB[j]):.5f} (quoted {qb})  fit {rf:.5f} (quoted {qf})')
out['residues_two_ways'] = res2


# (3) synthetic aliasing test
def synth(N, zs=ZS, rs=RB):
    W = np.exp(np.outer(np.log(N), zs))
    return (W * (2 * rs / zs)[None, :]).real.sum(axis=1)


Ssyn = synth(g)
c2s = fit(g, Ssyn, m2z, full)
syn = []
for im in quoted:
    i = int(np.argmin(np.abs(m2z.imag - im)))
    z = m2z[i]
    j = int(np.argmin(np.abs(ZS - z)))
    rf = abs(complex(c2s[2 * i], c2s[2 * i + 1])) * abs(z) / 2
    syn.append([float(z.imag), float(abs(RB[j])), rf])
out['synthetic_dense1024_fit_vs_input'] = syn
print('(3) synthetic signal (all 295 zeros, boundary residues), 1024/octave, M2 fit on [2^22,2^40): Im, input |r|, fitted |r|')
for x in syn:
    print(f'      {x[0]:9.3f}  {x[1]:.5f}  {x[2]:.5f}')
# synthetic holdout drop
hs = []
for T in (4200, NYQ - 1):
    zs = ZS[(ZS.real > 0.5) & (ZS.imag < T)]
    c = fit(g, Ssyn, zs, fitm)
    hs.append([T, len(zs)] + list(relerr(g, Ssyn, zs, c, tem)))
out['synthetic_holdout_T4200_vs_Nyquist'] = hs
print('(3) synthetic holdout (fit 2^22..2^37, test 2^37..2^40):', hs)
# synthetic out-of-sample error of M2 on [2^40, 2^42) (aliasing and detuning only, no noise)
Snew = synth(gnew)
err = float(np.sqrt(np.mean((Snew - design(gnew, m2z) @ c2s) ** 2) / np.mean(Snew ** 2)))
out['synthetic_M2_error_on_2p40_2p42'] = err
print('(3) synthetic M2 relative rms error on [2^40, 2^42) (model-internal expectation, not the Q2 evaluation):', err)
# other grids
for per in (1000.0, 1000 * math.sqrt(2)):
    raw = sorted(set(int(round(2 ** (10 + j / per))) for j in range(0, int(30 * per) + 1)))
    gg = np.array([x for x in raw if 2 ** 22 <= x < 2 ** 40], dtype=np.float64)
    nyq = per / 2 * 2 * math.pi / math.log(2)
    zz = ZS[(ZS.real > 0.5) & (ZS.imag < nyq - 1)]
    Sg = synth(gg)
    cg = fit(gg, Sg, zz, np.ones(len(gg), bool))
    rr = []
    for im in quoted:
        i = int(np.argmin(np.abs(zz.imag - im)))
        z = zz[i]
        j = int(np.argmin(np.abs(ZS - z)))
        rr.append([float(z.imag), float(abs(RB[j])), abs(complex(cg[2 * i], cg[2 * i + 1])) * abs(z) / 2])
    out[f'synthetic_per{per:.2f}'] = dict(nyquist=nyq, n_zeros=len(zz), rows=rr)
    print(f'(3) grid {per:.2f}/octave (Nyquist {nyq:.1f}, {len(zz)} zeros): Im, input |r|, fitted |r|:')
    for x in rr:
        print(f'      {x[0]:9.3f}  {x[1]:.5f}  {x[2]:.5f}')
json.dump(out, open(os.path.join(HERE, 'a7_prereg_and_alias.json'), 'w'), indent=1)
