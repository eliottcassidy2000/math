"""Output tables for the almost-regular loop predictors.  Session opus-2026-10-22-S54.

Reads this lane's predictions pred{n}.txt (predict_nonregular.py: y_mm = regular two-loop formula applied blindly,
y_mm_t = + path trees, y1/y2 = this lane's measured-order first/second order) and the Session30 decorated-core predictor
s30_pred{n}.txt (crosscheck_s30.py: y1/y2 with the weight grading w = e + leaves), merged by code.
For each n prints corr, Spearman, scatter (sd of y - pred), mean offset, MAE, RMSE, and (for the Session30 pair) how
many tournaments the second order improves.  y = n (1 - kappa/rho).  Usage: python summarize_tables.py n [n ...]"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def read(path):
    rows = [l.split() for l in open(path) if l.strip()]
    head = rows[0]
    return {r[0]: dict(zip(head, r)) for r in rows[1:]}


def rank(v):
    return np.argsort(np.argsort(v))


for n in map(int, sys.argv[1:]):
    mine = read(os.path.join(HERE, f'pred{n}.txt')) if os.path.exists(os.path.join(HERE, f'pred{n}.txt')) else {}
    s30 = read(os.path.join(HERE, f's30_pred{n}.txt')) if os.path.exists(os.path.join(HERE, f's30_pred{n}.txt')) else {}
    codes = [c for c in s30 if c in mine] if mine and s30 else list(s30 or mine)
    src = s30 if s30 else mine
    y = np.array([float(src[c]['y']) for c in codes])
    cols = []
    if mine:
        cols += [('regular two-loop formula (before)', np.array([float(mine[c]['y_mm']) for c in codes])),
                 ('+ path trees (mac-mini t)', np.array([float(mine[c]['y_mm_t']) for c in codes])),
                 ('measured-order 2nd (this lane)', np.array([float(mine[c]['y2']) for c in codes]))]
    if s30:
        cols += [('weight-1 predictor (Session30)', np.array([float(s30[c]['y1']) for c in codes])),
                 ('weight-2 predictor (Session30)', np.array([float(s30[c]['y2']) for c in codes]))]
    print(f'n={n}: {len(codes)} tournaments; y mean {y.mean():.4f}, sd {y.std():.4f}, range [{y.min():.4f}, {y.max():.4f}]')
    print(f'  {"predictor":36s} {"corr":>7s} {"Spear":>7s} {"scatter":>8s} {"offset":>8s} {"MAE":>7s} {"RMSE":>7s}')
    for label, p in cols:
        r = y - p
        print(f'  {label:36s} {np.corrcoef(y, p)[0, 1]:7.4f} {np.corrcoef(rank(y), rank(p))[0, 1]:7.4f} {r.std():8.4f} '
              f'{r.mean():+8.4f} {np.abs(r).mean():7.4f} {np.sqrt((r ** 2).mean()):7.4f}')
    if s30:
        e1 = np.abs(y - cols[-2][1]); e2 = np.abs(y - cols[-1][1])
        print(f'  Session30 second order improves |y - pred| on {int((e2 < e1).sum())} of {len(codes)}; '
              f'n * RMSE: first {n * np.sqrt(np.mean(e1 ** 2)):.4f}, second {n * np.sqrt(np.mean(e2 ** 2)):.4f}; '
              f'n^2 * RMSE second {n * n * np.sqrt(np.mean(e2 ** 2)):.3f}; n^2 * mean(y - y2) {n * n * np.mean(y - cols[-1][1]):+.3f}')
