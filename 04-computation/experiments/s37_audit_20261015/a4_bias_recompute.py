#!/usr/bin/env python3
"""S37 audit, part 4: HYP-9282 residue-class biases recomputed from the audit census (a2_census34.npz: class sums
mod 243 at the grid round(2^(10+j/64)) through 2^34, computed independently of the session's code).

1. All 363 channels g_(q,r) = H_(q,r)/sqrt(6N/3^q), octave means over octaves 14..32 (the session's range) and
   14..33 (one more octave); compare with collatz_residue_bias_deep_levels_20261015.json.
2. Statistics: |t| counts; lag-1 autocorrelation of octave means; AR(1)-corrected effective sample size; the
   t_18 null tail; dependence (sibling sums vanish identically); relative drifts and drift t-values.
3. D_j splits; the averaged refinement consistency check of HYP-9282.
4. Mercator line (26.5 cycles/octave, real part 0.99821) in S versus in the wavelet channels of levels 1..3:
   least-squares amplitude of N^sigma0 (a cos + b sin)(gamma ln N) in each channel, relative to S.
"""
import json, math, os
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.normpath(os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results'))
Z = np.load(os.path.join(HERE, 'a2_census34.npz'))
SB = json.load(open(os.path.join(RES, 'collatz_residue_bias_deep_levels_20261015.json')))
AN = json.load(open(os.path.join(RES, 'collatz_octave_memory_analysis_20261015.json')))
grid = Z['grid'].astype(np.float64)
C243 = Z['C243'].astype(np.float64)
u = np.log2(grid)


def classes(q):
    m = 3 ** q
    return np.stack([C243[:, a::m].sum(axis=1) for a in range(m)], axis=1)


Cq = {q: classes(q) for q in range(0, 6)}


def channel(q, r):
    return 3 * Cq[q][:, r] - Cq[q - 1][:, r % 3 ** (q - 1)]


def octave_means(x, ks):
    return np.array([x[(u >= k) & (u < k + 1)].mean() for k in ks])


def tail_t(x, nu):
    # two-sided P(|T_nu| > x)
    return float(mp.betainc(nu / 2, 0.5, 0, nu / (nu + x * x), regularized=True))


def stats(om, ks):
    n = len(om)
    m = om.mean()
    sd = om.std(ddof=1)
    t = m / (sd / math.sqrt(n))
    c = om - m
    r1 = float((c[:-1] * c[1:]).sum() / (c * c).sum())
    coef = np.polyfit(ks, om, 1)
    resid = om - np.polyval(coef, ks)
    se_slope = math.sqrt((resid ** 2).sum() / (n - 2) / ((np.array(ks) - np.mean(ks)) ** 2).sum())
    return m, sd, t, int((om > 0).sum()), float(coef[0]), r1, float(coef[0] / se_slope)


for (ks, label) in [(list(range(14, 33)), 'octaves 14..32 (session range, grid <= 2^33)'),
                    (list(range(14, 34)), 'octaves 14..33 (audit extension to 2^34)')]:
    rows = []
    for q in range(1, 6):
        wq = 6.0 / 3 ** q
        for r in range(3 ** q):
            g = channel(q, r) / np.sqrt(grid * wq)
            om = octave_means(g, ks)
            rows.append((q, r) + stats(om, ks) + (om,))
    n = len(ks)
    ts = np.array([x[4] for x in rows])
    print(f'== {label}: {len(rows)} channels; |t| > 5: {int((np.abs(ts) > 5).sum())}, |t| > 8: {int((np.abs(ts) > 8).sum())}')
    if n == 19:
        sess = {(c['q'], c['r']): c for c in SB['channels']}
        md = max(max(abs(x[2] - sess[(x[0], x[1])]['mean']), abs(x[4] - sess[(x[0], x[1])]['t']),
                     abs(x[6] - sess[(x[0], x[1])]['drift'])) for x in rows)
        mp_ = max(abs(x[5] - sess[(x[0], x[1])]['positive_octaves']) for x in rows)
        print(f'   vs session JSON (mean, t, drift): max abs difference {md:.2e}; positive-octave counts max diff {mp_}')
    r1s = np.array([x[7] for x in rows])
    print(f'   lag-1 autocorrelation of octave means across channels: median {np.median(r1s):+.3f}, mean {r1s.mean():+.3f}, '
          f'quartiles {np.quantile(r1s, 0.25):+.3f} / {np.quantile(r1s, 0.75):+.3f}')
    # AR(1)-corrected t: n_eff = n (1 - rho)/(1 + rho) with rho = max(pooled median, 0)
    rho = max(float(np.median(r1s)), 0.0)
    neff = n * (1 - rho) / (1 + rho)
    tc = ts * math.sqrt(neff / n)
    print(f'   pooled rho = {rho:.3f} -> n_eff = {neff:.1f}; corrected |t| > 5: {int((np.abs(tc) > 5).sum())}, '
          f'> 3.5: {int((np.abs(tc) > 3.5).sum())}')
    p5 = tail_t(5.0, n - 1)
    print(f'   iid null: P(|T_{n - 1}| > 5) = {p5:.2e} -> expected {len(rows) * p5:.3f} of {len(rows)}; '
          f'with n_eff: P(|T_{neff - 1:.1f}| > 5/sqrt(n/n_eff) = {5 * math.sqrt(neff / n):.2f}) = '
          f'{tail_t(5 * math.sqrt(neff / n), neff - 1):.2e} -> expected {len(rows) * tail_t(5 * math.sqrt(neff / n), neff - 1):.2f}')
    # sibling dependence: in each sibling triple the three channels sum to zero identically
    nfree = sum(2 * 3 ** (q - 1) for q in range(1, 6))
    print(f'   dependence: sibling triples sum to 0 identically; free channels {nfree} of {len(rows)}')
    top = sorted(rows, key=lambda x: -abs(x[4]))
    print('   channel, mean, sd, t, pos, drift, rel drift %/octave, drift t, lag-1 acf of octave means')
    want = [(1, 2), (1, 1), (2, 3), (2, 6), (3, 9), (3, 18), (3, 16), (4, 51), (4, 55), (4, 22), (5, 78), (5, 155), (5, 59), (5, 79)]
    for x in rows:
        if (x[0], x[1]) in want:
            print(f'   H_({x[0]},{x[1]:3d}) {x[2]:+.3f} {x[3]:.3f} {x[4]:+6.1f} {x[5]:2d} {x[6]:+.4f} {100 * x[6] / abs(x[2]):+6.1f} '
                  f'{x[8]:+5.1f} {x[7]:+.2f}')
    strong = [x for x in rows if abs(x[4]) > 8]
    big = [(f'H_({x[0]},{x[1]})', round(100 * x[6] / abs(x[2]), 1)) for x in strong if abs(x[6] / x[2]) > 0.025]
    toward0 = [x for x in strong if x[6] * x[2] < 0]
    print(f'   among the {len(strong)} channels with |t| > 8: relative drift > 2.5%/octave in {len(big)}: {big}')
    print(f'   among them drifting toward zero: {len(toward0)}; drift |t| > 3: '
          f'{[f"H_({x[0]},{x[1]})" for x in strong if abs(x[8]) > 3]}')
    tab = [x for x in rows if (x[0], x[1]) in want]
    big_t = [(f'H_({x[0]},{x[1]})', round(100 * x[6] / abs(x[2]), 1)) for x in tab if abs(x[6] / x[2]) > 0.025]
    print(f'   HYP-9282 table rows with relative drift > 2.5%/octave: {big_t}')
    # opposite signs of H_(3,18), H_(3,9) in every octave
    o918 = [x[9] for x in rows if (x[0], x[1]) in [(3, 9), (3, 18)]]
    print(f'   H_(3,9) > 0 and H_(3,18) < 0 in every octave: {bool(np.all(o918[0] > 0) and np.all(o918[1] < 0))}')
    # D_j splits
    for j in range(4):
        m1 = 3 ** (j + 1)
        D = Cq[j + 1][:, (3 ** j) % m1] - Cq[j + 1][:, (2 * 3 ** j) % m1]
        st = stats(octave_means(D / np.sqrt(grid), ks), ks)
        print(f'   D_{j}: mean {st[0]:+.4f}, t {st[2]:+.1f}, positive {st[3]}/{n}, drift {st[4]:+.4f}, lag-1 acf {st[5]:+.2f}')

# averaged refinement check (octave means of g over 14..33 of the audit census)
ks = list(range(14, 34))
gm = {}
for q in (1, 2):
    for r in range(3 ** q):
        gm[(q, r)] = octave_means(channel(q, r) / np.sqrt(grid * 6.0 / 3 ** q), ks).mean()
print('averaged refinement m_(1,r) vs -(m_(1,e(r)) + m_(2,o(r)))/sqrt(2), octaves 14..33:')
for r, e, o in [(2, 1, 8), (1, 2, 2), (0, 0, 5)]:
    print(f'   r = {r}: {gm[(1, r)]:+.3f} vs {-(gm[(1, e)] + gm[(2, o)]) / math.sqrt(2):+.3f}')

# correlation of adjacent octave means implied by the point ACF of g_top (from the session analysis JSON)
acf = AN['g_top_acf_by_octave_lag']
print(f"point ACF of g_top (session): {acf}")
print(f"window-mean ACF at lag 1 octave (session analysis JSON, adjacent one-octave windows): {AN['window_mean_acf']['1.0000']:+.3f}")
gt = channel(1, 2) / np.sqrt(2 * grid)
sel = u >= 14
x = gt[sel] - gt[sel].mean()
lags = np.arange(0, 129)
rho_pt = np.array([(x[:len(x) - l] * x[l:]).mean() / (x * x).mean() for l in lags])
# Corr(mean over [0,1), mean over [1,2)) for a stationary process with ACF rho(tau), tau in 1/64 steps
num = sum(rho_pt[abs(l)] * (64 - abs(l - 64)) for l in range(1, 128))
den = sum(rho_pt[abs(l)] * (64 - abs(l)) for l in range(-63, 64))
print(f'implied correlation of adjacent octave means of g_top from its own point ACF: {num / den:+.3f}')

# Mercator line in channels: fit N^sigma0 (a cos + b sin)(gamma ln N) + c N^sigma0 + d sqrt(N) on u in [20, 34]
s0 = complex('0.9982092207652568136445639331675041376460+240.2478301406912443253007903875287137674j')
sel = u >= 20
Nn = grid[sel]
ph = s0.imag * np.log(Nn)
Ns = Nn ** s0.real
Bm = np.column_stack([Ns * np.cos(ph), Ns * np.sin(ph), Ns, np.sqrt(Nn)])


def merc_amp(y):
    coef, *_ = np.linalg.lstsq(Bm, y, rcond=None)
    return math.hypot(coef[0], coef[1])


aS = merc_amp(Cq[0][sel, 0])
print(f'Mercator amplitude in S(N)/N^sigma0 (fit on u in [20,34]): {aS:.4e} (this is 2|r/s0| ~ 2*0.177/240.25 = {2 * 0.177 / 240.25:.4e})')
for q in (1, 2, 3):
    amps = []
    for r in range(3 ** q):
        amps.append(merc_amp(channel(q, r)[sel]))
    amps = np.array(amps)
    # a Haar-distributed share of S's line in one class has amplitude aS/3^q; channels cancel it exactly
    print(f'   level {q}: max channel Mercator amplitude {amps.max():.3e} = {amps.max() / aS:.2e} x (amplitude in S); '
          f'a class share would be {1 / 3 ** q:.3f} x')
