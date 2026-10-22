# Five-loop predictions for the doubly regular 23-tournaments against exact H (reads drt23_five_loop.json).
# Strict 1/n expansion of Delta log H = log H(D) - log H(P_23) (all lower orders are universal, HYP-9371 section 2):
#   n^-3 : Delta X_3(1);
#   n^-4 : Delta X_3(1) + Delta X_4(1) + dressing, dressing = sum_j [z^j] Delta X_3 (j D_1 + C(j,2))/n, D_1 = 2n/(n+1)
#          (the first-order finite-size response of the four-loop term; this is the complete n^-4 difference).
# Laws over the 13 non-Paley classes (12 annealed + doubled Paley): slope per transitive 5-set, Paley offset (Paley's
# value 0 minus the line at TT5 = 1518), max residual.  Per coefficient: exact Delta c_L against delta_L + Delta X_4^(L)
# (+ Delta X_5^(10) at L = 10, complete there).  Output: drt23_five_loop_analysis.out.  mac-mini-2026-10-22-loop4.
import json
from fractions import Fraction
import numpy as np
from mpmath import mp, mpf, log
mp.dps = 40
n = 23
d = json.load(open('drt23_five_loop.json'))
labels, TT5 = d['labels'], d['TT5']
iP = labels.index('P_23')
cls = json.load(open('drt23_classes.json')); dbl = json.load(open('drt23_doubled.json'))
H = {c['label']: int(c['H']) for c in cls}; H['D23_dbl'] = int(dbl['H'])
D1 = 2 * n / (n + 1)
rows = []
for k, lab in enumerate(labels):
    if k == iP: continue
    x3, x4, x5 = d['X3'][k], d['X4'][k], d['X5'][k]
    dl = np.array(x3['dcoeffs'])
    dress = float(np.sum(dl * np.array([(j * D1 + j * (j - 1) / 2) / n for j in range(len(dl))])))
    exact = float(log(mpf(H[lab]) / mpf(H['P_23'])))
    rows.append(dict(label=lab, TT5=TT5[k], exact=exact, s3=x3['dvalue'], x4=x4['dvalue'], dress=dress,
                     s4=x3['dvalue'] + x4['dvalue'] + dress, x4byV=x4['byV'],
                     dc3=x3['dcoeffs'], dc4=x4['dcoeffs'], dc5=x5['dcoeffs']))
print("Per class (relative to P_23).  s3 = Delta X_3(1);  s4 = s3 + Delta X_4(1) + dressing (complete n^-4)")
print(f"{'class':8s} {'TT5':>5s} {'exact':>12s} {'s3':>12s} {'Delta X_4(1)':>13s} {'dressing':>12s} {'s4':>12s} {'s4/exact':>9s}")
for r in rows:
    print(f"{r['label']:8s} {r['TT5']:5d} {r['exact']:+.5e} {r['s3']:+.5e} {r['x4']:+.6e} {r['dress']:+.5e} {r['s4']:+.5e} {r['s4'] / r['exact']:9.3f}")
tt = np.array([r['TT5'] for r in rows], float)
A = np.column_stack([np.ones(len(tt)), tt])
print("\nLaws over the 13 non-Paley classes: slope per TT5, Paley offset (0 minus the line at TT5 = 1518), max residual")
fits = {}
for name in ('exact', 's3', 'x4', 'dress', 's4'):
    y = np.array([r[name] for r in rows]); (a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    res = y - A @ [a, b]; fits[name] = (b, -(a + b * 1518), res)
    print(f"   {name:6s}: slope {b:+.4e}, Paley offset {-(a + b * 1518):+.4e}, max residual {np.max(np.abs(res)):.2e}")
r_ex, r_s4 = fits['exact'][2], fits['s4'][2]
print(f"   non-TT5 residuals: corr(exact, s4) = {np.corrcoef(r_ex, r_s4)[0, 1]:+.4f}; "
      f"max |exact - s4| after removing both lines {np.max(np.abs(r_ex - r_s4)):.2e}")
i1, i2 = [i for i, r in enumerate(rows) if r['label'] in ('D23_1', 'D23_2')]
print(f"   D23_1 vs D23_2 (equal 5-census): exact difference {rows[i1]['exact'] - rows[i2]['exact']:+.3e}, "
      f"Delta X_4(1) difference {rows[i1]['x4'] - rows[i2]['x4']:+.3e}")
print("\nDelta X_4(1) by core vertex count (class D23_0): " + ', '.join(f"V={V}: {v:+.3e}" for V, v in sorted(rows[0]['x4byV'].items(), key=lambda t: int(t[0]))))

# per coefficient: exact Delta c_L from forest numbers
ex = {r['label']: [int(x) for x in r['N']] for r in json.load(open('drt23_forest_exact.json'))}
ex['D23_dbl'] = [int(x) for x in dbl['N']]
def logc(N, LM=20):
    phi = [Fraction(N[L], n ** L) for L in range(LM + 1)]
    c = [Fraction(0)] * (LM + 1)
    for L in range(1, LM + 1):
        s = L * phi[L]
        for j in range(1, L): s -= j * c[j] * phi[L - j]
        c[L] = s / L
    return c
cP = logc(ex['P_23'])
print("\nPer coefficient: exact Delta c_L against loops through five (excess 5 complete at L = 10 only), class D23_0 and the TT5 slopes")
for L in (8, 10, 12, 14, 16):
    exs = np.array([float(logc(ex[r['label']])[L] - cP[L]) for r in rows])
    pred = np.array([r['dc3'][L] + r['dc4'][L] + (r['dc5'][L] if L == 10 else 0) for r in rows])
    (a1, b1), *_ = np.linalg.lstsq(A, exs, rcond=None); (a2, b2), *_ = np.linalg.lstsq(A, pred, rcond=None)
    print(f"   L = {L:2d}: D23_0 exact {exs[0]:+.6e}, loops 3-{5 if L == 10 else 4} {pred[0]:+.6e}; max |exact - loops| over classes "
          f"{np.max(np.abs(exs - pred)):.2e}; slope exact {b1:+.3e} vs loops {b2:+.3e}; Paley offset exact {-(a1 + b1 * 1518):+.3e} vs loops {-(a2 + b2 * 1518):+.3e}")
json.dump({'rows': rows}, open('drt23_five_loop_analysis.json', 'w'), indent=1)

# Forest-length slices fixed by loops <= 5 alone: Delta phi_8 = delta_8, Delta phi_10 = Delta c_10 + phi_2 Delta c_8
# with Delta c_10 = delta_10 + Delta X_4^(10) + Delta X_5^(10) (all from the engine) and phi_2 = c_2 = (n-1)/n
# (universal); weights g(L) = n^L/(n)_L, normalised by the universal doubly regular curve
# kappa = e(1 - 3/(2n) + 25/(8 n^2)) (HYP-9367), no exact input.  Compared with the exact slices and the exact total.
from math import e as E_, factorial as fact
g = lambda L: n ** L / (fact(n) / fact(n - L))
kap_u = E_ * (1 - 3 / (2 * n) + 25 / (8 * n * n))
phi2 = (n - 1) / n
print(f"\nSlices L = 8, 10 from loops <= 5 alone (universal kappa = {kap_u:.6f}) against the exact H difference")
print(f"{'class':8s} {'loops: L=8':>12s} {'loops: L=10':>12s} {'sum':>12s} {'exact total':>12s} {'share':>7s}")
for r in rows:
    d8 = r['dc3'][8]; d10 = r['dc3'][10] + r['dc4'][10] + r['dc5'][10]
    s8 = g(8) * d8 / kap_u; s10 = g(10) * (d10 + phi2 * d8) / kap_u
    r['slice8'], r['slice10'] = s8, s10
    print(f"{r['label']:8s} {s8:+.5e} {s10:+.5e} {s8 + s10:+.5e} {r['exact']:+.5e} {(s8 + s10) / r['exact']:7.3f}")
y = np.array([r['slice8'] + r['slice10'] for r in rows]); (a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
print(f"   loop slices L <= 10: TT5 slope {b:+.4e} (exact {fits['exact'][0]:+.4e}), Paley offset {-(a + b * 1518):+.4e} (exact {fits['exact'][1]:+.4e})")
json.dump({'rows': rows}, open('drt23_five_loop_analysis.json', 'w'), indent=1)
