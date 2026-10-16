#!/usr/bin/env python3
"""S37 audit 2, script b4: independent recomputation of the HYP-9282 out-of-sample test from the 2^37 census NPZ.

Channels g_(q,r)(N) = H_(q,r)(N)/sqrt(6N/3^q), H_(q,r) = 3 C(r mod 3^q) - C(r mod 3^(q-1)), q = 1..5, from the class
sums modulo 243 at the grid points N_j = round(2^(10+j/64)) (the grid itself was verified exactly by b3b).  Octave
means = plain averages over the grid points with 2^k <= N_j < 2^(k+1).  In-sample octaves 14..32 (19), out of sample
34..36 (octave 33 reported).  Also: naive t-values, sign persistence, slopes through the origin of out-of-sample on
in-sample means (all / |t| > 5 / |t| > 8), pooled octave-mean autocorrelation by level (14..36), AR(1)-corrected
t-values (effective n = n (1 - rho)/(1 + rho)), and a check of the pre-registration text in the HYP-9282 history.
Usage: python b4_hyp9282_recompute.py census_ext.npz
"""
import sys, math
import numpy as np

Z = np.load(sys.argv[1])
grid = Z['grid'].astype(np.int64)
C = Z['classsums'].astype(np.int64)
u = np.log2(grid.astype(np.float64))
N = grid.astype(np.float64)
names, vals = [], []
for q in range(1, 6):
    m, mp_ = 3 ** q, 3 ** (q - 1)
    # C(r mod 3^q) = sum of the mod-243 classes congruent to r mod 3^q
    Cq = np.stack([C[:, r::m].sum(axis=1) for r in range(m)], axis=1)
    Cp = np.stack([C[:, r::mp_].sum(axis=1) for r in range(mp_)], axis=1) if q > 1 else C.sum(axis=1)[:, None]
    for r in range(m):
        Hqr = 3 * Cq[:, r] - Cp[:, r % mp_]
        names.append((q, r))
        vals.append(Hqr / np.sqrt(6 * N / m))
G = np.stack(vals, axis=1)
lev = np.array([q for q, r in names])
KS = list(range(14, 37))
cnt = [int(((u >= k) & (u < k + 1)).sum()) for k in KS]
OM = np.stack([G[(u >= k) & (u < k + 1)].mean(axis=0) for k in KS])
print(f'grid points per octave 14..36: min {min(cnt)}, max {max(cnt)}')
ins = [KS.index(k) for k in range(14, 33)]
oo = [KS.index(k) for k in (34, 35, 36)]
mi = OM[ins].mean(axis=0)
si = OM[ins].std(axis=0, ddof=1)
ti = mi / (si / math.sqrt(len(ins)))
idx = {nm: i for i, nm in enumerate(names)}
print(f'naive |t| > 5: {(abs(ti) > 5).sum()}, > 8: {(abs(ti) > 8).sum()}, max |t| {abs(ti).max():.2f}; '
      f'one sign in all 19 in-sample octaves: {int(((OM[ins] > 0).all(0) | (OM[ins] < 0).all(0)).sum())}')
order = np.argsort(-abs(ti))
print('top 6 channels by naive |t|: ' + ', '.join(f'H{names[c]} t={ti[c]:+.1f}' for c in order[:6]))
pre = [((3, 18), -1, -2.44), ((3, 9), +1, 2.03), ((5, 155), -1, -1.42), ((5, 78), +1, 1.55), ((5, 59), +1, 0.92)]
print('pre-registered channels: in-sample mean | octave 33 | 34 35 36 | out mean | signs kept | in band +-0.6')
for (q, r), sg, centre in pre:
    c = idx[(q, r)]
    o = OM[oo, c]
    print(f'  H({q},{r}): {mi[c]:+.3f} (t {ti[c]:+.1f}) | {OM[KS.index(33), c]:+.3f} | '
          f'{o[0]:+.3f} {o[1]:+.3f} {o[2]:+.3f} | {o.mean():+.3f} | {int((np.sign(o) == sg).sum())}/3 | '
          f'{abs(o.mean() - centre) <= 0.6}')


def follow(sel):
    mo = OM[oo][:, sel].mean(axis=0)
    return (len(sel), int((np.sign(mo) == np.sign(mi[sel])).sum()), float((mo * mi[sel]).sum() / (mi[sel] ** 2).sum()),
            float(np.corrcoef(mo, mi[sel])[0, 1]))
for lab, sel in (('all', np.arange(len(names))), ('|t|>5', np.nonzero(abs(ti) > 5)[0]), ('|t|>8', np.nonzero(abs(ti) > 8)[0])):
    n, kept, slope, corr = follow(sel)
    print(f'  follow {lab}: n = {n}, sign of out-of-sample mean kept {kept}, slope {slope:.3f}, corr {corr:.3f}')
X = OM - OM.mean(axis=0)
v = (X ** 2).mean(axis=0)
for L in (1, 2):
    acf = [float(np.mean((X[L:, lev == q] * X[:-L, lev == q]).mean(axis=0) / v[lev == q])) for q in range(1, 6)]
    print(f'  pooled octave-mean ACF lag {L} by level 1..5 (octaves 14..36): ' + ', '.join(f'{a:.3f}' for a in acf))
# alternative estimator: per-channel lag-1 correlation with in-sample octaves only (unbiased-ish), median by level
r1 = []
for c in range(len(names)):
    x = OM[ins, c] - OM[ins, c].mean()
    r1.append((x[1:] * x[:-1]).sum() / (x ** 2).sum())
r1 = np.array(r1)
print('  median per-channel lag-1 correlation of in-sample octave means by level: ' +
      ', '.join(f'{np.median(r1[lev == q]):.3f}' for q in range(1, 6)))
for rho in (0.63, 0.65):
    neff = len(ins) * (1 - rho) / (1 + rho)
    tc = mi / (si / math.sqrt(neff))
    print(f'  AR(1)-corrected (rho = {rho}, n_eff = {neff:.2f}): |t| > 5: {(abs(tc) > 5).sum()}, > 3.5: {(abs(tc) > 3.5).sum()}')
# H = H_(1,2) octave means 33..36 and in-sample mean
c = idx[(1, 2)]
print(f'  H_(1,2) octave means 33..36: ' + ', '.join(f'{OM[KS.index(k), c]:+.3f}' for k in (33, 34, 35, 36)) +
      f'; in-sample mean {mi[c]:+.3f}')
