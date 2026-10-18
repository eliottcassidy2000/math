#!/usr/bin/env python3
"""Decay law of THM-4667's ideal-ladder cross terms E_d - E_(d-1), d = 1..12 (values from S41 audit a4,
04-computation/experiments/s41_audit_20261017/a4_ideal_ladder.out, quoted to the printed 6 digits).
Prints consecutive ratios, local power exponents alpha_d = -log(c_(d+1)/c_d)/log((d+1)/d), least-squares fits of
log|c_d| against d (geometric) and against log d (power) on d = 6..12, and the positivity budget: since E_d >= 0 and
every cross term is negative, sum_(d > 12) |c_d| <= E_12 = 1.479358, so the local exponent cannot settle at <= 1."""
import math
import numpy as np
E = [4.857142857, 3.388402513, 2.575381952, 2.137429991, 1.915078985, 1.778418992, 1.693329901, 1.634125798,
     1.590908284, 1.556295497, 1.527640161, 1.501753864, 1.479358488]
c = [E[d] - E[d - 1] for d in range(1, 13)]
print('cross terms d = 1..12:', ['%.6f' % x for x in c])
print('ratios c_(d+1)/c_d:', ['%.4f' % (c[i + 1] / c[i]) for i in range(11)])
alpha = [-math.log(c[i + 1] / c[i]) / math.log((i + 2) / (i + 1)) for i in range(11)]
print('local power exponents alpha_d (d -> d+1):', ['%.2f' % a for a in alpha])
d = np.arange(6, 13); y = np.log(-np.array(c[5:12]))
g = np.polyfit(d, y, 1); p = np.polyfit(np.log(d), y, 1)
rg = np.sqrt(np.mean((np.polyval(g, d) - y) ** 2)); rp = np.sqrt(np.mean((np.polyval(p, np.log(d)) - y) ** 2))
print('d = 6..12 fits: geometric ratio %.4f (rms %.4f) ; power exponent %.3f (rms %.4f)' % (math.exp(g[0]), rg, -p[0], rp))
print('positivity budget: E_12 = %.6f bounds the whole tail sum_(d>12) |c_d|; a c/d tail with c = %.3f would exceed it'
      % (E[12], -c[11] * 12))
