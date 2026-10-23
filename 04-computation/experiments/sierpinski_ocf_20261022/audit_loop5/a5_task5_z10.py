# Audit (HYP-9371 section 8), task 5: the z^10 identity through an independent code path, in exact arithmetic.
# Part A (validation): for regular 7- and 9-tournaments, the first-principles z^10 coefficients X_e^(10), e = 1..5
#   (a5_z10_lib: own kernel enumerator a5_kernel.c + exact integer hom numbers, own core-list checks in task 1),
#   against the exact blow-up scaling audit2/b_mscaling.py (no core theory), as exact rationals.
# Part B (n = 23): for all 14 doubly regular 23-tournaments of the session (12 annealed classes, P_23, doubled Paley),
#   exact Delta X_e^(10) (relative to P_23) for e = 0..5 and the exact Delta c_10 from the forest numbers; checks
#     Delta X_0^(10) = Delta X_1^(10) = Delta X_2^(10) = 0 (universality),
#     Delta c_10 = Delta X_3^(10) + Delta X_4^(10) + Delta X_5^(10) EXACTLY (as rationals),
#   and compares with the engine values in drt23_five_loop.json (floating point).
# Usage: python3 a5_task5_z10.py KERNEL_BINARY [A|B|AB]
import os, sys, json, time
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(D, "audit2")); sys.path.insert(0, D)
import a5_z10_lib as Z
from b_mscaling import mscaling
from common import rot_A, paley_A, skew, regular_classes, code_to_A

binary = sys.argv[1]; parts = sys.argv[2] if len(sys.argv) > 2 else 'AB'
L = 10
t0 = time.time()
tab = Z.kernel_table(binary, L, (1, 2, 3, 4, 5), cache_fn=os.path.join(HERE, 'a5_task5_K_L10.json'))
print(f"kernel table at L = {L}: " + ', '.join(f"excess {e}: {len(r)} cores, {sum(len(x[4]) for x in r)} length vectors, "
      f"{sum(1 for x in r for _, K in x[4] if K)} nonzero K" for e, r in tab.items()) + f" ({time.time() - t0:.0f}s)", flush=True)

def X10(S):
    H = Z.Homs(S, L)
    return {e: Z.X_coeff(tab[e], H, L) for e in tab}

if 'A' in parts:
    cases = [('R_7', skew(rot_A(7))), ('P_7', skew(paley_A(7)))]
    c7 = regular_classes(7); c9 = regular_classes(9)
    cases += [(f'n7 class {k}', skew(code_to_A(c))) for k, c in enumerate(c7)]
    cases += [('R_9', skew(rot_A(9)))] + [(f'n9 class {k}', skew(code_to_A(c9[k]))) for k in (2, 5, 11)]
    for name, S in cases:
        t1 = time.time()
        ex = mscaling(S, L, verbose=False)[L]
        mine = X10(S)
        res = ' '.join(f"e={e}:{'=' if mine[e] == ex.get(e, 0) else 'DIFF'}" for e in sorted(mine))
        print(f"   {name} (n = {S.shape[0]}): own X_e^(10) vs exact blow-up scaling: {res}; X_4^(10) = {float(mine[4]):+.13e}, "
              f"X_5^(10) = {float(mine[5]):+.13e} ({time.time() - t1:.0f}s)", flush=True)

if 'B' in parts:
    from drt_switch import code_to_S
    cls = json.load(open(os.path.join(D, 'drt23_classes.json'))); dbl = json.load(open(os.path.join(D, 'drt23_doubled.json')))
    labels = [c['label'] for c in cls] + ['D23_dbl']
    codes = [c['code'] for c in cls] + [dbl['code']]
    TT5 = [c['TT5'] for c in cls] + [dbl['TT5']]
    Ns = {r['label']: r['N'] for r in json.load(open(os.path.join(D, 'drt23_forest_exact.json')))}; Ns['D23_dbl'] = dbl['N']
    eng = json.load(open(os.path.join(D, 'drt23_five_loop.json')))
    assert eng['labels'] == labels
    n = 23
    X = {}; c10 = {}
    for lab, code in zip(labels, codes):
        t1 = time.time()
        S = code_to_S(code)
        assert np.array_equal(S @ S.T, n * np.eye(n, dtype=int) - np.ones((n, n), dtype=int))
        X[lab] = X10(S)
        # X_0^(10) exactly: -tr(Shat^10)/5 (THM-4757 (3): only odd m = L/2 contribute, -tr(Shat^(2m))/m)
        P10 = np.linalg.matrix_power(S.astype(object), 10)
        X[lab][0] = Fraction(-int(np.trace(P10)), 5 * n ** 10)
        c10[lab] = Z.logphi_coeffs(Ns[lab], n, 10)[10]
        print(f"   {lab}: X_e^(10) computed ({time.time() - t1:.0f}s)", flush=True)
    P = 'P_23'
    print("\nExact z^10 differences relative to P_23 (own first-principles code):")
    print(f"{'class':8s} {'TT5':>5s} {'dX0':>4s} {'dX1':>4s} {'dX2':>4s} {'Delta c_10 (exact)':>20s} {'dX3':>14s} {'dX4':>15s} {'dX5':>15s} "
          f"{'identity':>9s} {'max |own - engine| (dX3, dX4, dX5)':>36s}")
    allok = True
    for k, lab in enumerate(labels):
        if lab == P: continue
        d = {e: X[lab][e] - X[P][e] for e in range(6)}
        dc = c10[lab] - c10[P]
        ident = dc == d[3] + d[4] + d[5]
        univ = all(d[e] == 0 for e in (0, 1, 2))
        allok &= ident and univ
        ed = [abs(float(d[e]) - eng[f'X{e}'][k]['dcoeffs'][10]) for e in (3, 4, 5)]
        print(f"{lab:8s} {TT5[k]:5d} {'0' if d[0] == 0 else 'NZ':>4s} {'0' if d[1] == 0 else 'NZ':>4s} {'0' if d[2] == 0 else 'NZ':>4s} "
              f"{float(dc):+20.12e} {float(d[3]):+14.7e} {float(d[4]):+15.9e} {float(d[5]):+15.9e} {'EXACT' if ident else 'FAIL':>9s} "
              f"{max(ed):.1e} ({ed[0]:.1e}, {ed[1]:.1e}, {ed[2]:.1e})")
    print(f"all 13 classes: universality of excess 0-2 and the exact identity Delta c_10 = dX3 + dX4 + dX5: {'HOLD' if allok else 'FAIL'}")
    i1, i2 = labels.index('D23_1'), labels.index('D23_2')
    print(f"D23_2 - D23_1: Delta c_10 {float(c10['D23_2'] - c10['D23_1']):+.6e}; dX3 {float(X['D23_2'][3] - X['D23_1'][3]):+.1e}, "
          f"dX4 {float(X['D23_2'][4] - X['D23_1'][4]):+.6e}, dX5 {float(X['D23_2'][5] - X['D23_1'][5]):+.1e}")
    # store exact values
    json.dump({lab: {str(e): [X[lab][e].numerator, X[lab][e].denominator] for e in range(6)} | {'c10': [c10[lab].numerator, c10[lab].denominator]}
               for lab in labels}, open(os.path.join(HERE, 'a5_task5_z10_exact.json'), 'w'))
