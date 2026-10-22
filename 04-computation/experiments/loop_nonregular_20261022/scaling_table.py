"""Scaling of the second-order errors with n (HYP-9369 tables).  Session opus-2026-10-22-S54.

From budget{n}.txt (exact forest numbers + Session30 predictor) and s30_pred20.txt:
  y-units:   n^4 RMSE(y - y2), n^4 scatter(y - y2), n^2 mean(y - y2);
  log units: n^3 mean(A - W1 - W2), n^3 mean(B - B_f), corr of the two errors, n^3 mean(joint), n^5 sd(joint),
             joint = log(kappa/rho) - W1 - W2 - B_f.
Usage: python scaling_table.py"""
import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    rows = [l.split() for l in open(path) if l.strip()]
    head = rows[0]
    return {k: np.array([float(r[i]) for r in rows[1:]]) for i, k in enumerate(head) if k not in ('code', 'H')}


print('  n   count   n^4 RMSE   n^4 scatter   n^2 offset | n^3 (A-W)   n^3 (B-Bf)   corr     n^3 joint   n^5 sd(joint)')
for n, f in ((8, 'budget8.txt'), (10, 'budget10.txt'), (12, 'budget12.txt'), (14, 'budget14.txt'), (16, 'budget16.txt'),
             (20, 's30_pred20.txt')):
    d = load(os.path.join(HERE, f))
    r = d['y'] - d['y2']
    line = f'{n:3d} {len(r):7d} {n ** 4 * np.sqrt(np.mean(r ** 2)):10.1f} {n ** 4 * r.std():12.1f} {n ** 2 * r.mean():12.3f}'
    if 'A' in d:
        dA = d['A'] - d['W1'] - d['W2']; dB = d['B'] - d['Bf']; J = dA + dB
        line += (f' | {n ** 3 * dA.mean():9.1f} {n ** 3 * dB.mean():11.1f} {np.corrcoef(dA, dB)[0, 1]:9.4f} '
                 f'{n ** 3 * J.mean():10.2f} {n ** 5 * J.std():12.1f}')
    print(line)
