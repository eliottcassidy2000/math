#!/usr/bin/env python3
"""Decay law of THM-4667's ideal-ladder cross terms E_d - E_(d-1), d = 1..12 (values from S41 audit a4,
04-computation/experiments/s41_audit_20261017/a4_ideal_ladder.out, quoted to the printed 6 digits).
Prints consecutive ratios, local power exponents alpha_d = -log(c_(d+1)/c_d)/log((d+1)/d), least-squares fits of
log|c_d| against d (geometric) and against log d (power) on d = 6..12, and the positivity budget: since E_d >= 0,
sum_(d > 12) c_d >= -E_12 = -1.479358, so the cross terms cannot be eventually <= -C d^(-alpha) with C > 0 and
alpha <= 1 (their sign at all depths is not known, THM-4667). Values for d <= 6 are exact, d = 7..12 float64 with a
rigorous tail bound (audit G1, G1-6)."""
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
print('positivity budget: since E_d >= 0, sum_(d>12) c_d >= -E_12 = -%.6f; so no eventual bound c_d <= -C/d^alpha with alpha <= 1 (e.g. C = %.3f, alpha = 1)'
      % (E[12], -c[11] * 12))
d8 = np.arange(8, 13); y8 = np.log(-np.array(c[7:12]))
g8 = np.polyfit(d8, y8, 1); p8 = np.polyfit(np.log(d8), y8, 1)
print('d = 8..12 only: geometric ratio %.4f (rms %.4f) ; power exponent %.3f (rms %.4f)  [close: the power-law reading is HEURISTIC]'
      % (math.exp(g8[0]), np.sqrt(np.mean((np.polyval(g8, d8) - y8) ** 2)), -p8[0], np.sqrt(np.mean((np.polyval(p8, np.log(d8)) - y8) ** 2))))
