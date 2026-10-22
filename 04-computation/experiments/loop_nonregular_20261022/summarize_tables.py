"""Output tables for the almost-regular loop predictors.  Session opus-2026-10-22-S54.

Reads this lane's predictions pred{n}.txt or base{n}.txt (predict_nonregular.py or baselines.py: y_mm = regular two-loop formula applied blindly,
y_mm_t = + path trees, y2 = this lane's measured-order second order) and the Session30 decorated-core predictor
(weight grading w = e + leaves) from budget{n}.txt (budget_cw.py, which also carries the exact z = 1 budget) or
s30_pred{n}.txt (crosscheck_s30.py), merged by code.  y = n (1 - kappa/rho).
Table 1: corr, Spearman, scatter (sd of y - pred), mean offset, MAE, RMSE; how often the second order improves.
Table 2 (budget files only): exact log(kappa/rho) = A + B against the weight-graded truncation W1 + W2 of the cluster
part and the second-order inverse-factorial (finite-size) factor Bf, in units of 1/n.
Usage: python summarize_tables.py n [n ...]"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def read(path):
    rows = [l.split() for l in open(path) if l.strip()]
    if len(rows) < 2 or rows[0][0] != 'code':
        return {}
    head = rows[0][:8] if rows[0][3] == 'y_mm' else rows[0]          # pred files: the X(...) headers contain spaces
    return {r[0]: dict(zip(head, r)) for r in rows[1:]}


def rank(v):
    return np.argsort(np.argsort(v))


for n in map(int, sys.argv[1:]):
    f = lambda name: os.path.join(HERE, name)
    mine = read(f(f'pred{n}.txt')) if os.path.exists(f(f'pred{n}.txt')) else {}
    if not mine and os.path.exists(f(f'base{n}.txt')):
        mine = read(f(f'base{n}.txt'))                                # baselines.py: no measured-order y2
    budget = os.path.exists(f(f'budget{n}.txt'))
    s30 = read(f(f'budget{n}.txt')) if budget else (read(f(f's30_pred{n}.txt')) if os.path.exists(f(f's30_pred{n}.txt')) else {})
    codes = [c for c in s30 if c in mine] if (mine and s30) else list(s30 or mine)
    src = s30 if s30 else mine
    y = np.array([float(src[c]['y']) for c in codes])
    cols = []
    if mine:
        cols += [('regular two-loop formula (before)', np.array([float(mine[c]['y_mm']) for c in codes])),
                 ('+ path trees (mac-mini t)', np.array([float(mine[c]['y_mm_t']) for c in codes])),
                 ] + ([('measured-order 2nd (this lane)', np.array([float(mine[c]['y2']) for c in codes]))]
                      if 'y2' in mine[codes[0]] else [])
    if s30:
        cols += [('weight-1 predictor (Session30)', np.array([float(s30[c]['y1']) for c in codes])),
                 ('weight-2 predictor (Session30)', np.array([float(s30[c]['y2']) for c in codes]))]
    print(f'n={n}: {len(codes)} tournaments; exact y mean {y.mean():.4f}, sd {y.std():.4f}, range [{y.min():.4f}, {y.max():.4f}]')
    print(f'  {"predictor":36s} {"corr":>7s} {"Spear":>7s} {"scatter":>8s} {"offset":>8s} {"MAE":>7s} {"RMSE":>7s}')
    for label, p in cols:
        r = y - p
        print(f'  {label:36s} {np.corrcoef(y, p)[0, 1]:7.4f} {np.corrcoef(rank(y), rank(p))[0, 1]:7.4f} {r.std():8.4f} '
              f'{r.mean():+8.4f} {np.abs(r).mean():7.4f} {np.sqrt((r ** 2).mean()):7.4f}')
    if s30:
        e1 = np.abs(y - cols[-2][1]); e2 = np.abs(y - cols[-1][1])
        print(f'  Session30 second order improves |y - pred| on {int((e2 < e1).sum())} of {len(codes)}; '
              f'n * RMSE: first {n * np.sqrt(np.mean(e1 ** 2)):.4f}, second {n * np.sqrt(np.mean(e2 ** 2)):.4f}; '
              f'n^2 mean(y - y2) {n * n * np.mean(y - cols[-1][1]):+.3f}, n^2 scatter(y - y2) {n * n * np.std(y - cols[-1][1]):.3f}')
    if budget:
        g = lambda k: np.array([float(s30[c][k]) for c in codes])
        A, B, W1, W2, Bf = g('A'), g('B'), g('W1'), g('W2'), g('Bf')
        dA1, dA2, dB, J = n * (A - W1), n * (A - W1 - W2), n * (B - Bf), n * (A + B - W1 - W2 - Bf)
        print(f'  budget (units 1/n): n(A - W1) {dA1.mean():+.4f}; n(A - W1 - W2) {dA2.mean():+.4f} (sd {dA2.std():.4f}); '
              f'n(B - Bf) {dB.mean():+.4f} (sd {dB.std():.4f}); joint n(logk/rho - W1 - W2 - Bf) {J.mean():+.4f} (sd {J.std():.4f}); '
              f'n^2 x joint {n * J.mean():+.3f}; corr(dA2, dB) {np.corrcoef(dA2, dB)[0, 1]:+.4f}')
