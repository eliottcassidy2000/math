# Anatomy by forest length (exact) of (i) the D23_1 / D23_2 split (equal 5-census) and (ii) Paley's offset from the
# 13-class TT5 line, in units of Delta log H:  Delta kappa/kappa = sum_L g(L) Delta phi_L / kappa, g(L) = n^L/(n)_L.
# Output: a5_task6_split_anatomy.out.  Run from audit_loop5/.
import json
from fractions import Fraction
from math import factorial
import numpy as np
n = 23
cls = json.load(open('../drt23_classes.json')); dbl = json.load(open('../drt23_doubled.json'))
N = {r['label']: [int(x) for x in r['N']] for r in json.load(open('../drt23_forest_exact.json'))}; N['D23_dbl'] = [int(x) for x in dbl['N']]
H = {c['label']: int(c['H']) for c in cls}; H['D23_dbl'] = int(dbl['H'])
TT5 = {c['label']: c['TT5'] for c in cls}; TT5['D23_dbl'] = dbl['TT5']
kap = {lab: Fraction(H[lab] * 2 ** (n - 1), factorial(n)) for lab in H}
g = lambda L: Fraction(n ** L * factorial(n - L), factorial(n))
def sl(lab, L, ref='P_23'):
    return float(g(L) * Fraction(N[lab][L] - N[ref][L], n ** L) / kap[ref])
Ls = list(range(8, 23, 2))
tot = sum(sl('D23_1', L, 'D23_2') for L in range(0, n, 2))
print("D23_1 - D23_2 (relative to kappa(D23_2)): total %+.4e;" % tot, ', '.join(f"L={L}: {sl('D23_1', L, 'D23_2') / tot:+.3f}" for L in Ls[:5]),
      f", L>=18: {sum(sl('D23_1', L, 'D23_2') for L in Ls[5:]) / tot:+.3f}")
labs = [l for l in H if l != 'P_23']
tt = np.array([TT5[l] for l in labs], float); A = np.column_stack([np.ones(len(tt)), tt])
offs = {}
for L in Ls:
    y = np.array([sl(l, L) for l in labs]); (a, b), *_ = np.linalg.lstsq(A, y, rcond=None); offs[L] = -(a + b * 1518)
T = sum(offs.values())
print(f"Paley offset (13-class line) in Delta kappa/kappa: total {T:+.4e}; " + ', '.join(f"L={L}: {offs[L] / T:+.3f}" for L in Ls[:4]) +
      f", L>=16: {sum(offs[L] for L in Ls[4:]) / T:+.3f}")
