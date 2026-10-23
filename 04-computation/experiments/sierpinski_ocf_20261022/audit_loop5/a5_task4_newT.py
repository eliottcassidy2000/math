# Audit (HYP-9371 section 8), task 4: the engine's X_4 z-series (loop_engine5.series_general, all 1076 excess-4 cores)
# on tournaments NOT used in excess4_validate.py, against the exact blow-up scaling audit2/b_mscaling.py.
#   * a regular 9-tournament, gentourng class CLASS9 (checked non-isomorphic to R_9 and to class 2 of
#     excess4_validate.py by canonical labelling), to z^LMAX9 (z^12 brings in the cubic 8-vertex cores with injective
#     maps available, which n = 7 cannot test);
#   * optionally R_9 again with the aliasing-free contour (the session attributed R_9's 2.6e-9 to aliasing but
#     re-ran only R_7).
# Contour K = 64, r = 0.35.  Usage: python3 a5_task4_newT.py CLASS9 LMAX9 [R9]
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.dirname(HERE)
sys.path.insert(0, D); sys.path.insert(0, os.path.join(D, 'audit2'))
import loop_engine5 as LE5
from b_mscaling import mscaling
from common import rot_A, skew, regular_classes, code_to_A
import subprocess

def canon(A):
    n = A.shape[0]
    bits = ''.join(str(int(A[i, j])) for i in range(n) for j in range(n))
    # digraph6-like string for labelg
    b = [int(A[i, j]) for i in range(n) for j in range(n)]
    while len(b) % 6: b.append(0)
    s = '&' + chr(63 + n) + ''.join(chr(63 + int(''.join(map(str, b[k:k + 6])), 2)) for k in range(0, len(b), 6))
    return subprocess.run(['labelg', '-q'], input=s + '\n', capture_output=True, text=True).stdout.strip()

cl9 = int(sys.argv[1]); Lmax = int(sys.argv[2]); doR9 = len(sys.argv) > 3
c9 = regular_classes(9)
A = code_to_A(c9[cl9])
cR = canon(rot_A(9)); c2 = canon(code_to_A(c9[2])); cN = canon(A)
print(f"gentourng regular 9-tournaments: {len(c9)} classes; class {cl9}: distinct from R_9 {cN != cR}, from class 2 {cN != c2}; "
      f"score check {sorted(set(A.sum(1)))}", flush=True)
cases = [(f'N_9 class {cl9}', skew(A), Lmax)]
if doR9: cases.append(('R_9', skew(rot_A(9)), 10))
for name, S, Lm in cases:
    n = S.shape[0]; Sh = S / n
    t0 = time.time(); X = mscaling(S, Lm, verbose=False); tm = time.time() - t0
    t0 = time.time(); co = LE5.series_general(Sh, 4, Lm, K=64, r=0.35); te = time.time() - t0
    print(f"{name} (n = {n}): exact blow-up scaling {tm:.0f}s, engine (K = 64, r = 0.35) {te:.0f}s", flush=True)
    worst_rel = 0.0; worst_abs = 0.0
    for L in range(2, Lm + 1, 2):
        ex = X[L].get(4, 0); d = co[L] - float(ex)
        if ex != 0: worst_rel = max(worst_rel, abs(d) / abs(float(ex)))
        else: worst_abs = max(worst_abs, abs(d))
        print(f"   L={L:2d}: engine {co[L]:+.15e}  exact {float(ex):+.15e} ({ex})  diff {d:+.1e}", flush=True)
    print(f"   max relative error over nonzero orders {worst_rel:.1e}; max |engine| at vanishing orders {worst_abs:.1e}", flush=True)
