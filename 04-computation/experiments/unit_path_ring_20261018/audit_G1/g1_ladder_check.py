#!/usr/bin/env python3
"""audit G1 (independent): ladder numbers of the results note sections 5-6, recomputed from the S41 a4 printout
(04-computation/experiments/s41_audit_20261017/a4_ideal_ladder.out) and the clock-graded output
(04-computation/experiments/unit_path_ring_20261018/clock_graded_ladder_d8.out). Parses the files; no session code."""
import re, math, os
import numpy as np
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
a4 = open(os.path.join(HERE, '..', '..', 's41_audit_20261017', 'a4_ideal_ladder.out')).read()
E = {}
for m in re.finditer(r'\(b\) d = +(\d+): E ([0-9.]+)', a4): E[int(m.group(1))] = float(m.group(2))
ex = {int(m.group(1)): Fr(m.group(3)) for m in re.finditer(r'\(a\) d = (\d+):.*?\n\s+exact: gamma = (\S+), E = (\S+)', a4)}
print('E_d (a4, part b):', E)
print('exact E_d available for d =', sorted(ex), ':', {d: float(v) for d, v in ex.items()})
c = {d: E[d] - E[d - 1] for d in range(1, 13)}
print('cross terms:', {d: round(v, 6) for d, v in c.items()})
print('ratios c_(d+1)/c_d:', [round(c[d + 1] / c[d], 4) for d in range(1, 12)])
print('mean ratio for d <= 4 (c2/c1, c3/c2, c4/c3):', round(np.mean([c[2] / c[1], c[3] / c[2], c[4] / c[3]]), 4))
alpha = {d: -math.log(c[d + 1] / c[d]) / math.log((d + 1) / d) for d in range(1, 12)}
print('local exponents alpha_d (d -> d+1):', {d: round(a, 3) for d, a in alpha.items()})
print('  quoted in the note for d = 6..11: 2.4, 2.4, 1.9, 1.8, 1.1, 1.7 ; recomputed:', [round(alpha[d], 1) for d in range(6, 12)])
ds = np.arange(6, 13); y = np.log([-c[d] for d in ds])
g = np.polyfit(ds, y, 1); p = np.polyfit(np.log(ds), y, 1)
rg = math.sqrt(np.mean((np.polyval(g, ds) - y) ** 2)); rp = math.sqrt(np.mean((np.polyval(p, np.log(ds)) - y) ** 2))
print('fits on d = 6..12: geometric ratio %.4f rms %.4f ; power exponent %.3f rms %.4f' % (math.exp(g[0]), rg, -p[0], rp))
# leave-one-out robustness of the power exponent
loo = []
for k in ds:
    m = ds != k; loo.append(-np.polyfit(np.log(ds[m]), y[m], 1)[0])
print('power exponent, leave-one-out over d = 6..12: %.3f .. %.3f' % (min(loo), max(loo)))
# fits on d = 8..12 only
ds2 = np.arange(8, 13); y2 = np.log([-c[d] for d in ds2])
g2 = np.polyfit(ds2, y2, 1); p2 = np.polyfit(np.log(ds2), y2, 1)
rg2 = math.sqrt(np.mean((np.polyval(g2, ds2) - y2) ** 2)); rp2 = math.sqrt(np.mean((np.polyval(p2, np.log(ds2)) - y2) ** 2))
print('fits on d = 8..12: geometric ratio %.4f rms %.4f ; power exponent %.3f rms %.4f' % (math.exp(g2[0]), rg2, -p2[0], rp2))
print('E_12 =', E[12], '; sum of the observed cross terms d = 1..12 =', round(sum(c.values()), 6), '= E_12 - E_0 =', round(E[12] - E[0], 6))
# clock-graded rows
cg = open(os.path.join(HERE, '..', 'clock_graded_ladder_d8.out')).read()
rows = {}
for m in re.finditer(r'd = (\d+): E_d = ([0-9.]+) ; cross term E_d - E_\(d-1\) = ([-0-9.]+) ; lag 0: ([-+0-9.]+) lag 1: ([-+0-9.]+) lag 2: ([-+0-9.]+) lag 3: ([-+0-9.]+)', cg):
    rows[int(m.group(1))] = tuple(float(m.group(k)) for k in range(2, 8))
for d, r in sorted(rows.items()):
    print('clock-graded d = %d: E %.12f cross %.6f lags %+.6f %+.6f %+.6f %+.6f ; a4 E_d = %s ; diff %.1e%s'
          % (d, r[0], r[1], r[2], r[3], r[4], r[5], E.get(d), abs(r[0] - E[d]) if d in E else float('nan'),
             (' ; exact diff %.1e' % abs(r[0] - float(ex[d]))) if d in ex else ''))
print('rows present in clock_graded_ladder_d8.out: d =', sorted(rows), '(the d = 8 row is', 'present)' if 8 in rows else 'NOT present at audit time)')
for d in (6, 7):
    r = rows[d]
    print('d = %d: sum of the four printed lags = %+.6f vs cross %+.6f (remaining lags >= 4: %+.6f)' % (d, sum(r[2:6]), r[1], r[1] - sum(r[2:6])))
print('P3 of the pre-registration is implied by P1 and the known cross_8 = -0.043218: [-0.043218 + 0.098, -0.043218 + 0.106] =',
      [round(-0.043218 + 0.098, 6), round(-0.043218 + 0.106, 6)])
# full-precision comparison against the a4 json (positive-gauge float values, exact rationals for d <= 3)
import json
js = json.load(open(os.path.join(HERE, '..', '..', 's41_audit_20261017', 'a4_ideal_ladder.json')))
Ej = {d: row['E'] for d, row in enumerate(js['positive_gauge'])}
for d, r in sorted(rows.items()):
    print('d = %d: |clock-graded E_d - a4 json E_d| = %.1e' % (d, abs(r[0] - Ej[d])))
