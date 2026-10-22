"""The exact z = 1 budget for the 13 doubly regular 23-tournament classes (HYP-9369 section 2 on a regular family).
Session opus-2026-10-22-S54.  Inputs (mac-mini-2026-10-22-loop4): ../sierpinski_ocf_20261022/drt23_classes.json (codes)
and drt23_forest_exact.json (exact N_L).  Predictor: Session30's predict.py (unchanged).
log(kappa/rho) = A + B exactly; compares A with W1 + W2 and B with the second-order inverse-factorial factor B_f.
Usage: python budget_drt23.py"""
import sys, os, json, math
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'sierpinski_ocf_20261022')
ENGINE = os.path.normpath(os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results', 'session30_nonregular_loops_20261022',
                                       'session30-decorated-engine'))
sys.path.insert(0, ENGINE)
from predict import predict          # Session30 predictor (unchanged)

n = 23
cls = {c['label']: c for c in json.load(open(os.path.join(SRC, 'drt23_classes.json')))}
iu = np.triu_indices(n, 1)
rows = []
for d in json.load(open(os.path.join(SRC, 'drt23_forest_exact.json'))):
    code = cls[d['label']]['code']
    A_ = np.zeros((n, n), dtype=int); b = np.array([1 if ch == '1' else 0 for ch in code])
    A_[iu[0], iu[1]] = b; A_[iu[1], iu[0]] = 1 - b
    S = A_ - A_.T
    H = int(d['H']); N = [int(x) for x in d['N']]
    kap = Fraction(H * 2 ** (n - 1), math.factorial(n))
    assert sum(Fraction(N[L], math.prod(range(n - L + 1, n + 1))) for L in range(len(N))) == kap
    Phi1 = sum(Fraction(N[L], n ** L) for L in range(len(N)))
    p = predict(S)
    x2 = np.linalg.eigvalsh(1j * S / n) ** 2
    d1 = sum(float(np.sum(x2 ** m)) * (2 * m) / m for m in range(1, 151, 2))
    A = math.log(Phi1) - math.log(p['rho']); B = math.log(kap / Phi1)
    Bf = p['P1'] / n + ((p['D2_A1'] + (2 * d1 - 1) * p['D_A1']) / 2) / n + (p['P2'] - p['P1'] ** 2 / 2) / n ** 2
    y = n * (1 - float(kap) / p['rho'])
    rows.append((d['label'], y, p['y1'], p['y2'], A - p['A1'] - p['A2'], B - Bf, A + B - p['A1'] - p['A2'] - Bf))
    print(f"{d['label']:6s}: y {y:.6f}, y2 {p['y2']:.6f}; n^3 (A - W1 - W2) {n ** 3 * rows[-1][4]:+.2f}, "
          f"n^3 (B - B_f) {n ** 3 * rows[-1][5]:+.2f}, n^3 joint {n ** 3 * rows[-1][6]:+.4f}", flush=True)
R = np.array([r[1:] for r in rows])
print(f'n = 23, {len(rows)} doubly regular classes: mean n^3 (A - W1 - W2) {n ** 3 * R[:, 3].mean():+.2f}, '
      f'mean n^3 (B - B_f) {n ** 3 * R[:, 4].mean():+.2f}, mean n^3 joint {n ** 3 * R[:, 5].mean():+.4f}; '
      f'mean y - y2 {np.mean(R[:, 0] - R[:, 2]):+.2e}')
