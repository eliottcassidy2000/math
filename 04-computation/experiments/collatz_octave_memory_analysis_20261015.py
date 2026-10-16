#!/usr/bin/env python3
"""Analysis of the exact octave census (collatz_octave_memory_census_20261015.py).  Session opus-2026-10-15-S37.

(A) Per dyadic octave X = 2^k: normalised energy e_X = E_X/X^2 of H = H_(1,2), its mean part m^2 and Haar refinement
    part, the level-2 increment energy delta_X = sum_(X<=N<2X) (H(4N)-H(N))^2/X^2, the cross term
    c_X = sum H(N)(H(4N)-H(N))/X^2, kappa = c/sqrt(e delta), and the exact two-octave energy transport
    e_(4X) = (e_X + 2 c_X + delta_X)/4 + O(X^(-1/2))  (Proposition B of the S37 note), plus the 4-adic descent bound
    sqrt(e_(4X)) <= (sqrt(e_X) + sqrt(delta_X))/2 + sqrt(3/X).
(B) On the grid N_j = round(2^(10+j/64)): normalised profiles g_c(u) = H_c(N)/sqrt(N w_q) (w_1 = 2, w_2 = 2/3) for the
    level-1/2 channels; moving-window octave means m(u) and energies e(u) from prefix sums (64 phases per octave);
    autocorrelation of g_top and m in u; periodogram of g_top(u) and of S(N)/N (the Mercator line at 26.5 cycles per
    octave is the positive control); cross-channel correlations.
(C) A stationary Gaussian AR surrogate for the window energy e(u), fitted on u >= 16, simulated forward to estimate
    P(e(u) > 2 for some window start in [2^34, 2^K]).  This is a HEURISTIC model, not a statement about Collatz.
Usage: python collatz_octave_memory_analysis_20261015.py census.json census.npz out.json
"""
import sys, json, math
import numpy as np

J = json.load(open(sys.argv[1]))
Z = np.load(sys.argv[2])
OUT = sys.argv[3]
res = {}

# ---------- (A) octave table ----------
oct_rows = {r['k']: r for r in J['octaves']}
ks = sorted(k for k in oct_rows if k >= 10)
tab = []
for k in ks:
    r = oct_rows[k]
    e = r['energy_over_X2'][2]
    row = dict(k=k, e=e, mean_part=r['top_mean_part'], refinement=r['top_innovation_total'],
               refinement_top3=sum(r['top_innovation_by_level'][-3:]) if len(r['top_innovation_by_level']) >= 3 else None)
    if 'delta_energy_over_X2' in r:
        dl = r['delta_energy_over_X2']
        cx = r['H_delta_cross_over_X2']
        row.update(delta=dl, cross=cx, kappa=cx / math.sqrt(e * dl) if e > 0 and dl > 0 else None)
        if k + 2 in oct_rows:
            e4 = oct_rows[k + 2]['energy_over_X2'][2]
            pred = (e + 2 * cx + dl) / 4
            row.update(e_4X=e4, transport_prediction=pred, transport_rel_err=(e4 - pred) / e4,
                       descent_bound=((math.sqrt(e) + math.sqrt(dl)) / 2 + math.sqrt(3 / 2 ** k)) ** 2)
    tab.append(row)
res['octaves'] = tab

# ---------- (B) grid profiles ----------
grid = Z['grid'].astype(np.float64)
S = Z['S'].astype(np.float64)
C3 = Z['C3'].astype(np.float64)
C9 = Z['C9'].astype(np.float64)
u = np.log2(grid)
Hc = np.zeros((len(grid), 12))
for r in range(3):
    Hc[:, r] = 3 * C3[:, r] - S
for r in range(9):
    Hc[:, 3 + r] = 3 * C9[:, r] - C3[:, r % 3]
w = np.array([2.0] * 3 + [2.0 / 3] * 9)
g = Hc / np.sqrt(grid[:, None] * w[None, :])
sel = u >= 14
gt = g[sel, 2]
ut = u[sel]
res['g_top_stats'] = dict(n=int(sel.sum()), mean=float(gt.mean()), var=float(gt.var()), max_abs=float(np.abs(gt).max()),
                          kurtosis=float(((gt - gt.mean()) ** 4).mean() / gt.var() ** 2))
step = 1 / 64


def acf(x, maxlag):
    x = x - x.mean()
    v = (x * x).mean()
    return [float((x[:len(x) - l] * x[l:]).mean() / v) for l in range(maxlag + 1)]


A = acf(gt, 4 * 64)
res['g_top_acf_by_octave_lag'] = {f'{l / 64:.4f}': A[l] for l in [1, 8, 16, 32, 48, 64, 96, 128, 160, 192, 256]}
# channel correlations at equal u
cc = np.corrcoef(g[sel].T)
res['corr_level1'] = [[float(cc[i, j]) for j in range(3)] for i in range(3)]
res['var_by_channel'] = [float(g[sel, c].var()) for c in range(12)]
# two-octave memory at the grid: g(u+2) vs g(u)/2
lag2 = 128
x0 = gt[:-lag2]
x2 = gt[lag2:]
res['two_octave_regression'] = dict(slope=float(np.polyfit(x0, x2, 1)[0]), corr=float(np.corrcoef(x0, x2)[0, 1]))
# periodogram (uniform in u)
def periodogram(x, dt):
    x = (x - x.mean()) * np.hanning(len(x))
    F = np.fft.rfft(x)
    f = np.fft.rfftfreq(len(x), dt)
    P = np.abs(F) ** 2
    return f, P
fq, Pg = periodogram(gt, step)
sn = (S / grid)[sel]
fq2, Ps = periodogram(sn, step)
def top_peaks(f, P, kk=6, fmin=0.3):
    m = f >= fmin
    idx = np.argsort(P[m])[::-1][:kk]
    return [(float(f[m][i]), float(P[m][i] / P[m].sum())) for i in idx]
res['periodogram_g_top_peaks'] = top_peaks(fq, Pg)
res['periodogram_S_over_N_peaks'] = top_peaks(fq2, Ps)
res['periodogram_g_top_power_fraction_below_1_cycle_per_octave'] = float(Pg[(fq > 0) & (fq < 1)].sum() / Pg[fq > 0].sum())
# moving-window octave means and energies from prefix sums: window [N_j, N_(j+64))
ps = Z['prefix_sum'].astype(np.float64)
pq = Z['prefix_sq'].astype(np.float64)
jj = np.arange(len(grid) - 64)
Y = grid[jj]
L = grid[jj + 64] - grid[jj]
mu = (ps[jj + 64, 2] - ps[jj, 2]) / L
ew = (pq[jj + 64, 2] - pq[jj, 2]) / Y ** 2
mw = mu / np.sqrt(Y)
uu = np.log2(Y)
keep = uu >= 14
mw = mw[keep]; ew = ew[keep]; uu = uu[keep]; L = L[keep]; Y = Y[keep]
meanpart = (L / Y) * mw ** 2
res['window_energy'] = dict(n=int(len(ew)), mean=float(ew.mean()), max=float(ew.max()), argmax_log2=float(uu[np.argmax(ew)]),
                            frac_mean_part=float(meanpart.sum() / ew.sum()), quantiles={q: float(np.quantile(ew, q)) for q in (0.5, 0.9, 0.99)})
res['window_mean_acf'] = {f'{l / 64:.4f}': v for l, v in zip([1, 16, 32, 64, 128, 192], [acf(mw, 192)[i] for i in [1, 16, 32, 64, 128, 192]])}

# ---------- (C) Gaussian AR surrogate for m(u) at quarter-octave resolution ----------
sub = 16
ms = mw[::sub]          # quarter-octave samples
p_ar = 8
X = np.column_stack([ms[p_ar - i - 1:len(ms) - i - 1] for i in range(p_ar)])
yv = ms[p_ar:]
coef, *_ = np.linalg.lstsq(np.column_stack([X, np.ones(len(yv))]), yv, rcond=None)
resid = yv - np.column_stack([X, np.ones(len(yv))]) @ coef
sig = resid.std()
ref = ew - meanpart          # refinement part of window energy
rng = np.random.default_rng(2026)
def simulate(Koct, nsim=4000):
    steps = int(Koct * 4)
    exceed = 0
    exceed_dyadic = 0
    maxima = []
    for _ in range(nsim):
        h = list(ms[-p_ar:])
        mx = 0.0
        mxd = 0.0
        for t in range(steps):
            nxt = coef[-1] + sum(coef[i] * h[-1 - i] for i in range(p_ar)) + rng.normal(0, sig)
            h.append(nxt)
            e = nxt * nxt + ref[rng.integers(0, len(ref))]
            mx = max(mx, e)
            if t % 4 == 3:          # windows starting at integer octaves (dyadic X)
                mxd = max(mxd, e)
        maxima.append(mx)
        exceed += mx > 2
        exceed_dyadic += mxd > 2
    return exceed / nsim, exceed_dyadic / nsim, float(np.median(maxima))
sims = {}
for Kfut in (6, 16, 30, 66, 166):
    pr, prd, med = simulate(Kfut)
    sims[f'{34 + Kfut}'] = dict(prob_some_window_energy_above_2=pr, prob_some_dyadic_window_energy_above_2=prd,
                                median_max_window_energy=med)
res['AR_surrogate'] = dict(order=p_ar, resolution_octaves=0.25, innovation_sd=float(sig),
                           stationary_var_m=float(np.var(ms)), simulation_from_2_34_to_2_K=sims,
                           note='HEURISTIC surrogate fitted to the finite census; not a model-free prediction')
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps({k: res[k] for k in res if k != 'octaves'}, indent=1))
print('octave table (k, e, mean part, refinement, delta, kappa, e_4X, prediction, descent bound):')
for r in tab:
    print(r['k'], round(r['e'], 4), round(r['mean_part'], 4), round(r['refinement'], 4),
          round(r.get('delta', float('nan')), 4), round(r['kappa'], 3) if r.get('kappa') is not None else None,
          round(r.get('e_4X', float('nan')), 4), round(r.get('transport_prediction', float('nan')), 4),
          round(r.get('descent_bound', float('nan')), 4))
