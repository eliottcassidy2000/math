#!/usr/bin/env python3
"""S37 audit, part 4b: sign persistence of HYP-9282 channels against an autocorrelated zero-mean null, and the
first octave beyond the session's mod-243 data ([2^33, 2^34), available in the audit census) for the five
pre-registered channels.  Null model: per channel, zero-mean stationary Gaussian AR(1) octave means with lag-1
correlation rho in {0.3 (H's adjacent-octave value), 0.5, 0.63 (median over channels)}; P(19 octaves of one sign).
"""
import os, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(HERE, 'a2_census34.npz'))
grid = Z['grid'].astype(np.float64)
C243 = Z['C243'].astype(np.float64)
u = np.log2(grid)
Cq = {q: np.stack([C243[:, a::3 ** q].sum(axis=1) for a in range(3 ** q)], axis=1) for q in range(6)}


def g(q, r):
    return (3 * Cq[q][:, r] - Cq[q - 1][:, r % 3 ** (q - 1)]) / np.sqrt(grid * 6.0 / 3 ** q)


ks = list(range(14, 33))
allsame = 0
for q in range(1, 6):
    for r in range(3 ** q):
        x = g(q, r)
        om = np.array([x[(u >= k) & (u < k + 1)].mean() for k in ks])
        allsame += int(np.all(om > 0) or np.all(om < 0))
print(f'channels (of 363) whose 19 octave means 14..32 all have one sign: {allsame}')
rng = np.random.default_rng(9282)
for rho in (0.0, 0.3, 0.5, 0.63, 0.75):
    nsim = 400000
    e = rng.standard_normal((nsim, 19))
    x = np.empty_like(e)
    x[:, 0] = e[:, 0]
    for i in range(1, 19):
        x[:, i] = rho * x[:, i - 1] + math.sqrt(1 - rho * rho) * e[:, i]
    p = float(np.mean(np.all(x > 0, axis=1) | np.all(x < 0, axis=1)))
    print(f'  AR(1) zero-mean null, rho = {rho:.2f}: P(19 same-sign octaves) = {p:.2e} -> expected {363 * p:.2f} of 363 '
          f'(242 free channels: {242 * p:.2f})')
pre = [(3, 18, '-'), (3, 9, '+'), (5, 155, '-'), (5, 78, '+'), (5, 59, '+')]
print('pre-registered channels: octave means at octaves 31, 32 (session data) and 33 (first octave beyond the '
      'session mod-243 census; inside its level-1/2 census range):')
for q, r, sgn in pre:
    x = g(q, r)
    vals = [x[(u >= k) & (u < k + 1)].mean() for k in (31, 32, 33)]
    ok = (vals[2] > 0) == (sgn == '+')
    print(f'  H_({q},{r}): {vals[0]:+.3f} {vals[1]:+.3f} | {vals[2]:+.3f}  predicted sign {sgn}: '
          f'{"kept" if ok else "FLIPPED"}')
