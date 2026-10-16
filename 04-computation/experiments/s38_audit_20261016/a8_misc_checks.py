#!/usr/bin/env python3
"""S38 audit, items 4, 5, 7: small checks.

(1) The census grid of the pre-registered evaluation is unchanged by the edit of the dense script made after the
    pre-registration (c3bfa60c5c: 'per' parsed as float, j-range int((LOGMAX - 10) * PER)) for PER = 1024, LOGMAX = 42,
    and equals the prediction grid of the committed JSON.
(2) Stability of the boundary residue estimates at all 295 zeros between N = 2^36, 2^38 and 2^40 (relative change of
    r_hat), by band of Re rho.
(3) The almost-period 2 pi 665/ln 2: which shifted zeros do not land near a zero of the strip.
(4) Phase drift accumulated by the k = 1 alias pairs on the 1024-per-octave grid over [2^22, 2^40) and [2^40, 2^42).
Usage: python a8_misc_checks.py moments.npz zeros.json prereg.json torus_zeros.npy
"""
import sys, os, json, math
import numpy as np

MOM, ZJ, PJ, TZ = sys.argv[1:5]
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
# (1)
PER_i, PER_f, LOGMAX = 1024, 1024.0, 42
old = sorted(set(int(round(2 ** (10 + j / PER_i))) for j in range(0, (LOGMAX - 10) * PER_i + 1)))
new = sorted(set(int(round(2 ** (10 + j / PER_f))) for j in range(0, int((LOGMAX - 10) * PER_f) + 1)))
P = json.load(open(PJ))
pg = P['prediction']['grid']
sel_new = [x for x in new if 2 ** 40 <= x < 2 ** 42]
out['grid_old_equals_new'] = (old == new)
out['eval_grid_equals_prediction_grid'] = (sel_new == pg)
print('(1) grid unchanged by the post-registration edit:', old == new, '; equals the prediction grid:', sel_new == pg, flush=True)
# (2)
L2, L15 = math.log(2.0), math.log(1.5)
B = 1 << 17
v = np.load(MOM)['moments']
sm_all = v[:, 0].astype(np.float64)
mu1_all = (2 * v[:, 1].astype(np.float64) - (B - 1) * sm_all) / 2
cen_all = np.arange(len(v), dtype=np.float64) * B + (B - 1) / 2
lc_all = np.log(cen_all)
rows = json.load(open(ZJ))['zeros_and_residues']
Z = np.array([complex(r['re'], r['im']) for r in rows])


def rhat(s, N):
    tN = (-1) ** (N.bit_length() - 1)
    Lc = (2 * N - 1) // 3
    Jc = Lc // 2
    jl = (Jc + 1 + B - 1) // B
    js = (Lc + 1 + B - 1) // B
    sl = slice(jl, N // B)
    val = np.exp(-s * lc_all[sl]) * (sm_all[sl] - s * mu1_all[sl] / cen_all[sl])
    last = tN * np.exp(-s * math.log(N))
    lsum = val.sum() + last
    ssum = val[js - jl:].sum() + last
    C = np.exp(s * math.log(2 / 3)) * (ssum + np.exp(-s * L2) * lsum)
    Dp = -L2 * np.exp(-s * L2) + L15 * np.exp(s * L15) / 3
    return -np.exp(s * L15) * C / Dp


r36 = np.array([rhat(z, 1 << 36) for z in Z])
r38 = np.array([rhat(z, 1 << 38) for z in Z])
r40 = np.array([rhat(z, 1 << 40) for z in Z])
ch1 = np.abs(r40 - r38) / np.abs(r40)
ch2 = np.abs(r38 - r36) / np.abs(r38)
bands = {}
for lo_, hi_ in ((0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 0.98), (0.98, 1.0)):
    m = (Z.real > lo_) & (Z.real <= hi_)
    bands[f'{lo_}-{hi_}'] = dict(n=int(m.sum()), median_rel_change_38_40=float(np.median(ch1[m])),
                                 frac_change_below_10pct=float(np.mean(ch1[m] < 0.1)),
                                 median_rel_change_36_38=float(np.median(ch2[m])))
out['residue_stability_by_band'] = bands
print('(2) relative change of r_hat between N = 2^38 and 2^40 (and 2^36 -> 2^38), by band of Re rho:')
for k, x in bands.items():
    print('   ', k, x)
imax = int(np.argmax(np.abs(r40)))
out['largest_abs_rhat_history'] = [[float(abs(r36[imax])), float(abs(r38[imax])), float(abs(r40[imax]))], [float(Z[imax].real), float(Z[imax].imag)]]
print('    largest |r_hat| (at', Z[imax], '): |r_hat| at 2^36, 2^38, 2^40 =', out['largest_abs_rhat_history'][0], flush=True)
# (3)
Zt = np.load(TZ)
P665 = 2 * math.pi * 665 / L2
far = []
for z in Zt:
    w = z + 1j * P665
    if w.imag < 20000:
        dmin = np.abs(Zt - w).min()
        if dmin > 0.01:
            far.append([float(z.real), float(z.imag), float(dmin)])
out['shift_far_cases'] = far
print('(3) zeros whose shift by 2 pi 665/ln 2 is not within 0.01 of a zero in 1/2 < Re < 1:', far, flush=True)
# (4) alias phase drift
th = math.log(3) / L2
kap = L2 - math.log(1.5)
det = (2048 * th - 3246) * L2 / (2 * kap)      # leading-order detuning, cycles per octave
out['alias_detuning_cycles_per_octave_leading'] = det
out['alias_phase_drift_2p22_2p40_degrees'] = det * 18 * 360
out['alias_phase_drift_2p40_2p42_degrees'] = det * 2 * 360
print(f'(4) k = 1 alias pairs on the 1024 grid: detuning {det:.5f} cycles/octave (leading order): {det * 18 * 360:.1f} deg '
      f'over [2^22, 2^40), {det * 2 * 360:.1f} deg over [2^40, 2^42)', flush=True)
json.dump(out, open(os.path.join(HERE, 'a8_misc_checks.json'), 'w'), indent=1)
