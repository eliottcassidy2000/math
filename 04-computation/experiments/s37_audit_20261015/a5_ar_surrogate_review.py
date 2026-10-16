#!/usr/bin/env python3
"""S37 audit, part 5: HYP-9283(A).  (1) Recompute the quoted profile statistics of g(u) = H(2^u)/sqrt(2^(u+1)) and the
window energies from the AUDIT census (a2_census34.npz).  (2) Re-fit the session's Gaussian AR(8) surrogate on the
quarter-octave window means, check stationarity (AR roots), and (3) run a sensitivity analysis of the exceedance
probabilities P(some dyadic window E_X > 2X^2 for X in [2^34, 2^K)) with vectorised simulations:
AR order, bias removed, refinement fixed at its mean, innovation sd +-20 %, fit on u >= 20 only, and an iid-Gaussian
(no memory) control.  All probabilities are Monte Carlo with nsim = 4000 (s.e. <= 0.8 %).
"""
import os, json, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results'))
Z = np.load(os.path.join(HERE, 'a2_census34.npz'))
AN = json.load(open(os.path.join(RES, 'collatz_octave_memory_analysis_20261015.json')))
grid = Z['grid'].astype(np.float64)
C243 = Z['C243'].astype(np.float64)
S = C243.sum(axis=1)
C2 = C243[:, 2::3].sum(axis=1)
H = 3 * C2 - S
u = np.log2(grid)
sel = u >= 14
g = (H / np.sqrt(2 * grid))[sel]
print(f'g(u), u in [14,34]: n = {len(g)}, mean {g.mean():.4f}, sd {g.std():.4f}, kurtosis '
      f'{((g - g.mean()) ** 4).mean() / g.var() ** 2:.3f}, max |g| {np.abs(g).max():.4f}')


def acf(x, lags):
    x = x - x.mean()
    v = (x * x).mean()
    return [float((x[:len(x) - l] * x[l:]).mean() / v) for l in lags]


lags = [1, 8, 16, 32, 48, 64, 128]
print('ACF at lags (octaves) 1/64, 1/8, 1/4, 1/2, 3/4, 1, 2:', [round(v, 3) for v in acf(g, lags)])
x0, x2 = g[:-128], g[128:]
print(f'two-octave regression slope {np.polyfit(x0, x2, 1)[0]:.4f}, corr {np.corrcoef(x0, x2)[0, 1]:.4f}')
# rough uncertainty of the kurtosis: effective number of independent samples ~ (range in octaves)/(decorrelation ~0.75)
neff = 20 / 0.75
print(f'kurtosis standard error for ~{neff:.0f} effective samples: ~sqrt(24/neff) = {math.sqrt(24 / neff):.2f}')
xx = (g - g.mean()) * np.hanning(len(g))
P = np.abs(np.fft.rfft(xx)) ** 2
f = np.fft.rfftfreq(len(g), 1 / 64)
print(f'periodogram power fraction below 1 cycle/octave: {P[(f > 0) & (f < 1)].sum() / P[f > 0].sum():.4f}; '
      f'frequency resolution {1 / 20:.3f} cycles/octave')
ps = Z['prefix_sum'].astype(np.float64)
pq = Z['prefix_sq']
jj = np.arange(len(grid) - 64)
Y = grid[jj]
L = grid[jj + 64] - grid[jj]
mu = (ps[jj + 64] - ps[jj]) / L
ew = (pq[jj + 64] - pq[jj]) / Y ** 2
mw = mu / np.sqrt(Y)
uu = np.log2(Y)
keep = uu >= 14
mw, ew, uu, L, Y = mw[keep], ew[keep], uu[keep], L[keep], Y[keep]
meanpart = (L / Y) * mw ** 2
print(f'window energies: n = {len(ew)}, mean {ew.mean():.4f}, 99% quantile {np.quantile(ew, 0.99):.4f}, max {ew.max():.4f} '
      f'at log2 X = {uu[np.argmax(ew)]:.3f}; mean-part fraction {meanpart.sum() / ew.sum():.4f}')
print(f'session analysis JSON: {AN["window_energy"]["mean"]:.4f}, {AN["window_energy"]["quantiles"]["0.99"]:.4f}, '
      f'{AN["window_energy"]["max"]:.4f}, {AN["window_energy"]["frac_mean_part"]:.4f}')
dy = np.isclose(uu, np.round(uu), atol=1e-9)
print(f'dyadic windows only (k = 14..33): max e = {ew[dy].max():.4f}; count {int(dy.sum())}')

ms_full = mw[::16]
ref_full = ew - meanpart


def fit_ar(ms, p, intercept=True):
    X = np.column_stack([ms[p - i - 1:len(ms) - i - 1] for i in range(p)])
    yv = ms[p:]
    A = np.column_stack([X, np.ones(len(yv))]) if intercept else X
    coef, *_ = np.linalg.lstsq(A, yv, rcond=None)
    resid = yv - A @ coef
    if not intercept:
        coef = np.append(coef, 0.0)
    return coef, resid.std()


def simulate(ms, coef, sig, ref, Koct, nsim=4000, seed=1, fixed_ref=False, start=None):
    rng = np.random.default_rng(seed)
    p = len(coef) - 1
    h = np.tile((ms[-p:] if start is None else start)[::-1], (nsim, 1))   # h[:, 0] most recent
    steps = int(Koct * 4)
    mx = np.zeros(nsim)
    mxd = np.zeros(nsim)
    mx_after = np.zeros(nsim)
    rbar = ref.mean()
    for t in range(steps):
        nxt = coef[-1] + h @ coef[:p] + rng.normal(0, sig, nsim)
        h = np.column_stack([nxt, h[:, :-1]])
        e = nxt * nxt + (rbar if fixed_ref else ref[rng.integers(0, len(ref), nsim)])
        mx = np.maximum(mx, e)
        if t >= 3:
            mx_after = np.maximum(mx_after, e)
        if t % 4 == 3:
            mxd = np.maximum(mxd, e)
    return float((mx > 2).mean()), float((mxd > 2).mean()), float((mx_after > 2).mean()), float(np.median(mx))


coef8, sig8 = fit_ar(ms_full, 8)
roots = np.roots(np.concatenate([[1.0], -coef8[:8]]))
print(f'AR(8) fit: coefficients {np.round(coef8[:8], 3).tolist()}, intercept {coef8[8]:.4f}, innovation sd {sig8:.4f} '
      f'(session {AN["AR_surrogate"]["innovation_sd"]:.4f})')
print(f'  largest AR root modulus {np.abs(roots).max():.4f} (stationary iff < 1); implied mean '
      f'{coef8[8] / (1 - coef8[:8].sum()):.4f} vs sample mean {ms_full.mean():.4f}')
# implied stationary variance via long simulation
rng = np.random.default_rng(5)
z = np.zeros(200000)
for i in range(8, len(z)):
    z[i] = coef8[8] + coef8[:8] @ z[i - 8:i][::-1] + rng.normal(0, sig8)
print(f'  implied stationary variance of m {z[1000:].var():.4f} vs sample variance {ms_full.var():.4f} '
      f'(n = {len(ms_full)} quarter-octave samples, {len(ms_full) / 4:.1f} octaves)')

print('sensitivity of P(dyadic E_X > 2X^2 for some X in [2^34, 2^K)), K = 40, 64, 200 '
      '(and all-phase windows starting at or after 2^34, K = 64):')
configs = []
for p in (2, 4, 8, 12):
    c, s = fit_ar(ms_full, p)
    configs.append((f'AR({p}) bias incl.', c, s, ref_full, False, ms_full))
c, s = fit_ar(ms_full - ms_full.mean(), 8, intercept=False)
configs.append(('AR(8) bias removed', c, s, ref_full, False, ms_full - ms_full.mean()))
configs.append(('AR(8) refinement fixed', coef8, sig8, ref_full, True, ms_full))
configs.append(('AR(8) sd x 0.8', coef8, 0.8 * sig8, ref_full, False, ms_full))
configs.append(('AR(8) sd x 1.2', coef8, 1.2 * sig8, ref_full, False, ms_full))
half = uu[::16] >= 20
c, s = fit_ar(ms_full[half], 8)
configs.append(('AR(8) fit u >= 20', c, s, ref_full, False, ms_full))
# iid control: dyadic window means iid Gaussian with the sample mean/variance of dyadic m
c = np.zeros(2)
c[-1] = ms_full.mean()
configs.append(('iid Gaussian (no memory)', np.array([0.0, ms_full.mean()]), ms_full.std(), ref_full, False, ms_full))
for name, c, s, ref, fixed, ms in configs:
    out = []
    for K in (6, 30, 166):
        pr, prd, pra, med = simulate(ms, c, s, ref, K, fixed_ref=fixed)
        out.append((prd, pra, med))
    rts = np.roots(np.concatenate([[1.0], -c[:-1]])) if len(c) > 2 else np.array([0.0])
    print(f'  {name:26s} root {np.abs(rts).max():.3f}: K=40 {100 * out[0][0]:5.1f}%  K=64 {100 * out[1][0]:5.1f}%  '
          f'K=200 {100 * out[2][0]:5.1f}%  | all-phase K=64 {100 * out[1][1]:5.1f}% | median max (K=64) {out[1][2]:.2f}')
print('session JSON (AR(8)): K=40 5.7%, K=64 30.7%, K=200 87.7%; all-phase K=64 46.6% (includes three windows '
      'starting at 2^33.25..2^33.75, inside the census)')
