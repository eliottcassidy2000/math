#!/usr/bin/env python3
"""Analysis of the 2^37 census extension: pre-registered tests of HYP-9282 and HYP-9283.  Session opus-2026-10-15-S37.

(1) Channels g_(q,r)(N) = H_(q,r)(N)/sqrt(6N/3^q), H_(q,r) = 3 C(r mod 3^q) - C(r mod 3^(q-1)), q = 1..5 (363 channels),
    octave means over the 64 log-uniform grid points N_j = round(2^(10+j/64)) of each octave [2^k, 2^(k+1)).
    In-sample octaves 14..32 reproduce the HYP-9282 table; out-of-sample octaves 34, 35, 36 test the pre-registered
    signs and the +-0.6 bands; all channels with in-sample |t| > 5 are followed out of sample.
(2) Octave-to-octave autocorrelation of the channel octave means (pooled by level), and t-values corrected for it.
(3) Window energies e = E/Y^2 of H over [N_j, N_(j+64)), coarse-mean and refinement parts, for window starts up to
    2^36; profile g(u) = H(2^u)/sqrt(2^(u+1)): moments, autocorrelation, two-octave regression slope.
(4) The HEURISTIC Gaussian AR(8) surrogate of HYP-9283, refitted exactly as in collatz_octave_memory_analysis_20261015.py
    on the in-sample windows (window ends <= 2^34), simulated over the out-of-sample window starts (2^33, 2^36]:
    predictive distribution of the maximum window energy, compared with the observed maximum.
Usage: python collatz_census_extension_analysis_20261015.py census_ext.npz out.json
"""
import sys, json, math
import numpy as np

Z = np.load(sys.argv[1])
OUT = sys.argv[2]
grid = Z['grid'].astype(np.int64)
C = Z['classsums'].astype(np.float64)
N = grid.astype(np.float64)
u = np.log2(N)
res = {'logmax': int(Z['logmax'])}

# ---------- (1) channels ----------
lev = {0: C.sum(axis=1)[:, None]}
for q in range(1, 6):
    m = 3 ** q
    lev[q] = np.stack([C[:, a::m].sum(axis=1) for a in range(m)], axis=1)
names, cols = [], []
for q in range(1, 6):
    m, mp = 3 ** q, 3 ** (q - 1)
    for r in range(m):
        names.append((q, r))
        cols.append((3 * lev[q][:, r] - lev[q - 1][:, r % mp]) / np.sqrt(6 * N / m))
G = np.stack(cols, axis=1)
KS = list(range(14, int(Z['logmax'])))
OM = np.stack([G[(u >= k) & (u < k + 1)].mean(axis=0) for k in KS])          # octaves x channels
nper = [int(((u >= k) & (u < k + 1)).sum()) for k in KS]
ins = [i for i, k in enumerate(KS) if k <= 32]
outk = [34, 35, 36]
outs = [KS.index(k) for k in outk if k in KS]
mi = OM[ins].mean(axis=0)
si = OM[ins].std(axis=0, ddof=1)
ti = mi / (si / math.sqrt(len(ins)))
kk = np.array([KS[i] for i in ins], dtype=np.float64)
drift = np.array([np.polyfit(kk, OM[ins, c], 1)[0] for c in range(OM.shape[1])])
res['samples_per_octave'] = dict(zip([str(k) for k in KS], nper))
idx = {nm: i for i, nm in enumerate(names)}
table = ['H(1,2)', 'H(1,1)', 'H(2,3)', 'H(2,6)', 'H(3,9)', 'H(3,18)', 'H(3,16)', 'H(4,51)', 'H(4,55)', 'H(4,22)',
         'H(5,78)', 'H(5,155)', 'H(5,59)', 'H(5,79)']
rows = []
for lab in table:
    q, r = map(int, lab[2:-1].split(','))
    c = idx[(q, r)]
    rows.append(dict(channel=lab, insample_mean=float(mi[c]), insample_sd=float(si[c]), insample_t=float(ti[c]),
                     insample_sign_octaves=int((np.sign(OM[ins, c]) == np.sign(mi[c])).sum()),
                     drift=float(drift[c]), oct33=float(OM[KS.index(33), c]) if 33 in KS else None,
                     out_octaves={str(KS[i]): float(OM[i, c]) for i in outs},
                     out_mean=float(OM[outs, c].mean()) if outs else None))
res['table'] = rows
pre = [((3, 18), -1, -2.44), ((3, 9), +1, 2.03), ((5, 155), -1, -1.42), ((5, 78), +1, 1.55), ((5, 59), +1, 0.92)]
prer = []
for (q, r), sg, centre in pre:
    c = idx[(q, r)]
    vals = [float(OM[i, c]) for i in outs]
    av = float(np.mean(vals))
    prer.append(dict(channel=f'H({q},{r})', predicted_sign=sg, out_octave_means=dict(zip([str(KS[i]) for i in outs], vals)),
                     signs_kept=int(sum(1 for v in vals if np.sign(v) == sg)), n_out=len(vals), out_average=av,
                     band=[centre - 0.6, centre + 0.6], in_band=bool(abs(av - centre) <= 0.6)))
res['preregistered'] = prer
strong = [c for c in range(len(names)) if abs(ti[c]) > 5]
very = [c for c in range(len(names)) if abs(ti[c]) > 8]
def follow(sel):
    if not outs or not sel:
        return None
    mo = OM[outs][:, sel].mean(axis=0)
    keep_each = [int((np.sign(OM[i, sel]) == np.sign(mi[sel])).sum()) for i in outs]
    slope = float((mo * mi[sel]).sum() / (mi[sel] ** 2).sum())
    return dict(n=len(sel), sign_kept_out_mean=int((np.sign(mo) == np.sign(mi[sel])).sum()),
                sign_kept_per_octave=dict(zip([str(KS[i]) for i in outs], keep_each)),
                slope_out_on_in_through_origin=slope, corr_out_in=float(np.corrcoef(mo, mi[sel])[0, 1]))
res['follow_t_above_5'] = follow(strong)
res['follow_t_above_8'] = follow(very)
allc = list(range(len(names)))
res['follow_all_channels'] = follow(allc)
# ---------- (2) octave autocorrelation of the octave means ----------
def pooled_acf(Mx, lags=(1, 2, 3)):
    X = Mx - Mx.mean(axis=0)
    v = (X ** 2).mean(axis=0)
    out = {}
    for L in lags:
        out[L] = float(np.mean((X[L:] * X[:-L]).mean(axis=0) / v))
    return out
acfs = {}
for q in range(1, 6):
    sel = [c for c, (qq, r) in enumerate(names) if qq == q]
    acfs[q] = pooled_acf(OM[:, sel])
res['octave_mean_acf_by_level'] = acfs
tcorr = np.zeros(len(names))
for c, (q, r) in enumerate(names):
    rho = acfs[q]
    infl = 1 + 2 * sum(max(rho[L], 0) * (1 - L / len(ins)) for L in (1, 2, 3))
    tcorr[c] = ti[c] / math.sqrt(infl)
res['t_corrected_counts'] = dict(above5=int((np.abs(tcorr) > 5).sum()), above8=int((np.abs(tcorr) > 8).sum()),
                                 above5_raw=len(strong), above8_raw=len(very))
# full-range (14..36) statistics for the pre-registered channels
full = []
for (q, r), sg, centre in pre:
    c = idx[(q, r)]
    v = OM[:, c]
    full.append(dict(channel=f'H({q},{r})', mean_14_36=float(v.mean()), sd=float(v.std(ddof=1)),
                     t_raw=float(v.mean() / (v.std(ddof=1) / math.sqrt(len(v)))),
                     drift_14_36=float(np.polyfit(np.array(KS, dtype=float), v, 1)[0])))
res['preregistered_full_range'] = full

# ---------- (3) window energies and profile of H ----------
H = Z['H'].astype(np.float64)
p1 = Z['preH1'].astype(np.float64)
p2 = Z['preH2'].astype(np.float64)
jj = np.arange(len(grid) - 64)
Y = N[jj]
L = N[jj + 64] - N[jj]
E = (p2[jj + 64] - H[jj + 64] ** 2) - (p2[jj] - H[jj] ** 2)
SH = (p1[jj + 64] - H[jj + 64]) - (p1[jj] - H[jj])
e = E / Y ** 2
meanpart = (SH / L) ** 2 * L / Y ** 2
uY = np.log2(Y)
def wstats(mask):
    return dict(n=int(mask.sum()), mean=float(e[mask].mean()), max=float(e[mask].max()),
                argmax_log2=float(uY[mask][np.argmax(e[mask])]), frac_mean_part=float(meanpart[mask].sum() / e[mask].sum()),
                q50=float(np.quantile(e[mask], 0.5)), q99=float(np.quantile(e[mask], 0.99)))
res['window_energy_insample'] = wstats((uY >= 14) & (uY + 1 <= 34 + 1e-9))
res['window_energy_out'] = wstats(uY > 33 + 1e-9)
dy = {}
for k in range(14, int(Z['logmax'])):
    j = int(np.nonzero(grid == (1 << k))[0][0])
    if j + 64 < len(grid) or (1 << (k + 1)) < (1 << int(Z['logmax'])):
        if j < len(jj):
            dy[str(k)] = dict(e=float(e[j]), mean_part=float(meanpart[j]), refinement=float(e[j] - meanpart[j]))
res['dyadic_windows'] = dy
g = H / np.sqrt(2 * N)
sel = u >= 14
gs = g[sel]
res['profile'] = dict(mean=float(gs.mean()), sd=float(gs.std()), kurtosis=float(((gs - gs.mean()) ** 4).mean() / gs.var() ** 2),
                      max_abs=float(np.abs(gs).max()), argmax_abs_log2=float(u[sel][np.argmax(np.abs(gs))]),
                      out_mean=float(g[u >= 34].mean()), out_sd=float(g[u >= 34].std()))
x = gs - gs.mean()
res['profile_acf'] = {f'{L / 64:.3f}': float((x[L:] * x[:-L]).mean() / x.var()) for L in (16, 32, 64, 128, 192)}
res['two_octave_slope'] = float(np.polyfit(gs[:-128], gs[128:], 1)[0])
res['two_octave_slope_insample'] = float(np.polyfit(g[sel & (u < 32)][:-128], g[sel & (u < 32)][128:], 1)[0])

# ---------- (4) AR surrogate refit and predictive check ----------
fitmask = (uY >= 14) & (jj + 64 <= int(np.nonzero(grid == (1 << 34))[0][0]))
mw = (SH / L / np.sqrt(Y))[fitmask]
ew = e[fitmask]
mp = meanpart[fitmask]
sub = 16
ms = mw[::sub]
p_ar = 8
Xa = np.column_stack([ms[p_ar - i - 1:len(ms) - i - 1] for i in range(p_ar)])
yv = ms[p_ar:]
A = np.column_stack([Xa, np.ones(len(yv))])
coef, *_ = np.linalg.lstsq(A, yv, rcond=None)
sig = (yv - A @ coef).std()
ref = ew - mp
rng = np.random.default_rng(2026)
def simulate(Koct, nsim=4000):
    steps = int(Koct * 4)
    mxs, mxds = [], []
    for _ in range(nsim):
        h = list(ms[-p_ar:])
        mx = mxd = 0.0
        for t in range(steps):
            nxt = coef[-1] + sum(coef[i] * h[-1 - i] for i in range(p_ar)) + rng.normal(0, sig)
            h.append(nxt)
            ee = nxt * nxt + ref[rng.integers(0, len(ref))]
            mx = max(mx, ee)
            if t % 4 == 3:
                mxd = max(mxd, ee)
        mxs.append(mx)
        mxds.append(mxd)
    return np.array(mxs), np.array(mxds)
mxs, mxds = simulate(3)
obs_mask = (uY > 33 + 1e-9) & (np.abs((uY * 4) - np.round(uY * 4)) < 1e-9)    # quarter-octave starts as in the surrogate
obs_max = float(e[obs_mask].max()) if obs_mask.any() else None
dmask = (np.abs(uY - np.round(uY)) < 1e-9) & (uY > 33 + 1e-9)
obs_dmax = float(e[dmask].max()) if dmask.any() else None
res['surrogate_check_window_starts_2^33_to_2^36'] = dict(
    ar_coefficients=[float(c) for c in coef], innovation_sd=float(sig),
    predicted_median_max=float(np.median(mxs)), predicted_q10_q90=[float(np.quantile(mxs, 0.1)), float(np.quantile(mxs, 0.9))],
    predicted_prob_max_above_2=float((mxs > 2).mean()), observed_max_quarter_starts=obs_max,
    predictive_upper_tail_prob_of_observed=float((mxs >= obs_max).mean()) if obs_max is not None else None,
    predicted_median_dyadic_max=float(np.median(mxds)), observed_dyadic_max=obs_dmax,
    predictive_upper_tail_prob_dyadic=float((mxds >= obs_dmax).mean()) if obs_dmax is not None else None)
json.dump(res, open(OUT, 'w'), indent=1)
for key in ['preregistered', 'follow_t_above_5', 'follow_t_above_8', 'follow_all_channels', 'octave_mean_acf_by_level',
            't_corrected_counts', 'preregistered_full_range', 'window_energy_insample', 'window_energy_out', 'profile',
            'profile_acf', 'two_octave_slope', 'two_octave_slope_insample', 'surrogate_check_window_starts_2^33_to_2^36']:
    print(key, json.dumps(res[key]))
print('table:')
for r in rows:
    print(r['channel'], f"{r['insample_mean']:+.3f} {r['insample_sd']:.3f} {r['insample_t']:+.1f} {r['insample_sign_octaves']}/19 "
          f"drift {r['drift']:+.3f} | 33: {r['oct33']:+.3f} | out", {k: round(v, 3) for k, v in r['out_octaves'].items()})
print('dyadic windows:', {k: round(v['e'], 3) for k, v in dy.items()})
