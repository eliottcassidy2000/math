# Audit (HYP-9371 section 8): Delta X_4(1) (the five-loop term at z = 1 used in the strict n^-4 column) recomputed in
# the engine's MATRIX mode (actual 23 x 23 propagators, no doubly regular reduction) for D23_0, D23_1, D23_2 against
# P_23, compared with drt23_five_loop.json (DRT mode): D23_0 +3.544140e-06, D23_1 - D23_2 = -7.314e-09.
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.dirname(HERE)
sys.path.insert(0, D)
import loop_engine5 as LE5
from drt_switch import code_to_S
cls = {c['label']: c for c in json.load(open(os.path.join(D, 'drt23_classes.json')))}
eng = json.load(open(os.path.join(D, 'drt23_five_loop.json')))
X = {}
for lab in ('P_23', 'D23_0', 'D23_1', 'D23_2'):
    t0 = time.time()
    S = code_to_S(cls[lab]['code'])
    X[lab] = LE5.X_general(S / 23, 4, np.array([1.0]))[0].real
    print(f"{lab}: X_4(1) matrix mode {X[lab]:+.15e} ({time.time() - t0:.0f}s)", flush=True)
for lab in ('D23_0', 'D23_1', 'D23_2'):
    k = eng['labels'].index(lab)
    print(f"{lab}: Delta X_4(1) matrix mode {X[lab] - X['P_23']:+.6e}, DRT mode {eng['X4'][k]['dvalue']:+.6e}")
k1, k2 = eng['labels'].index('D23_1'), eng['labels'].index('D23_2')
print(f"D23_1 - D23_2: matrix mode {X['D23_1'] - X['D23_2']:+.4e}, DRT mode {eng['X4'][k1]['dvalue'] - eng['X4'][k2]['dvalue']:+.4e}")
