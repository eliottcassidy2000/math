# Audit (HYP-9371 section 8), task 6: every number of section 8 recomputed from the stored engine values
# (drt23_five_loop.json), the exact H values and forest numbers (drt23_classes.json, drt23_doubled.json,
# drt23_forest_exact.json), with own code (exact rationals where possible).  Also: fits with the two converse
# duplicates removed (11 distinct classes up to converse), exact kappa instead of the universal curve, and the
# per-class ratios quoted in the z^10 bullet.
import os, sys, json
from fractions import Fraction
from math import factorial, comb
import numpy as np
from mpmath import mp, mpf, log, e as E_
mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from a5_z10_lib import logphi_coeffs

n = 23
d = json.load(open(os.path.join(D, 'drt23_five_loop.json')))
labels, TT5 = d['labels'], d['TT5']
iP = labels.index('P_23')
cls = json.load(open(os.path.join(D, 'drt23_classes.json'))); dbl = json.load(open(os.path.join(D, 'drt23_doubled.json')))
H = {c['label']: int(c['H']) for c in cls}; H['D23_dbl'] = int(dbl['H'])
Ns = {r['label']: [int(x) for x in r['N']] for r in json.load(open(os.path.join(D, 'drt23_forest_exact.json')))}
Ns['D23_dbl'] = [int(x) for x in dbl['N']]
conv = {c['label']: c.get('converse_class') for c in cls}
print("converse_class field:", {k: v for k, v in conv.items() if v is not None})

# exact kappa and the forest identity
for lab in ('P_23', 'D23_0'):
    kap = Fraction(H[lab] * 2 ** (n - 1), factorial(n))
    kap2 = sum(Fraction(Ns[lab][L], factorial(n) // factorial(n - L)) for L in range(n))
    print(f"{lab}: kappa = H 2^22/23! = {float(kap):.9f}; forest identity sum N_L/(n)_L equal: {kap == kap2}")
kapP = float(Fraction(H['P_23'] * 2 ** (n - 1), factorial(n)))
kap_u = float(E_ * (1 - mpf(3) / (2 * n) + mpf(25) / (8 * n * n)))
print(f"universal kappa e(1 - 3/(2n) + 25/(8n^2)) = {kap_u:.9f}; exact kappa(P_23) = {kapP:.9f}; ratio - 1 = {kap_u / kapP - 1:+.2e}")

D1 = 2 * n / (n + 1)
rows = []
cP = logphi_coeffs(Ns['P_23'], n, 16)
for k, lab in enumerate(labels):
    if k == iP: continue
    x3, x4, x5 = d['X3'][k], d['X4'][k], d['X5'][k]
    dl = np.array(x3['dcoeffs'])
    dress = float(sum(dl[j] * (j * D1 + comb(j, 2)) / n for j in range(len(dl))))
    exact = float(log(mpf(H[lab]) / mpf(H['P_23'])))
    c = logphi_coeffs(Ns[lab], n, 16)
    dc = {L: c[L] - cP[L] for L in (8, 10, 12, 14, 16)}
    g = lambda L: Fraction(n ** L * factorial(n - L), factorial(n))
    phi2 = Fraction(Ns['P_23'][2], n ** 2)
    d8 = x3['dcoeffs'][8]; d10 = x3['dcoeffs'][10] + x4['dcoeffs'][10] + x5['dcoeffs'][10]
    sl8 = float(g(8)) * d8; sl10 = float(g(10)) * (d10 + float(phi2) * d8)
    # exact slices from the forest numbers: g(L) Delta phi_L / kappa(P)
    dphi = {L: Fraction(Ns[lab][L] - Ns['P_23'][L], n ** L) for L in range(0, n, 2)}
    ex_sl = {L: float(g(L) * dphi[L]) / kapP for L in (8, 10)}
    rows.append(dict(label=lab, TT5=TT5[k], exact=exact, s3=x3['dvalue'], x4=x4['dvalue'], dress=dress,
                     s4=x3['dvalue'] + x4['dvalue'] + dress, dc=dc, d3_10=x3['dcoeffs'][10], d4_10=x4['dcoeffs'][10],
                     d5_10=x5['dcoeffs'][10], d3_8=d8, slices_u=(sl8 + sl10) / kap_u, slices_x=(sl8 + sl10) / kapP,
                     ex_sl=ex_sl, phi2=phi2, x4byV=x4['byV']))
print(f"phi_2 = N_2/n^2 = {rows[0]['phi2']} (= (n-1)/n: {rows[0]['phi2'] == Fraction(n - 1, n)})")

print("\nTable of section 8 (exact, s3, s4, loop slices with universal / exact kappa; shares):")
for r in rows:
    print(f"{r['label']:8s} TT5 {r['TT5']} exact {r['exact']:+.5e}  s3 {r['s3']:+.3e} ({r['s3'] / r['exact']:.4f})  "
          f"s4 {r['s4']:+.4e} ({r['s4'] / r['exact']:.4f})  slices(univ) {r['slices_u']:+.4e} ({r['slices_u'] / r['exact']:.4f})  "
          f"slices(exact kappa) {r['slices_x'] / r['exact']:.4f}  exact L=8,10 slices {r['ex_sl'][8] + r['ex_sl'][10]:+.4e} "
          f"({(r['ex_sl'][8] + r['ex_sl'][10]) / r['exact']:.4f})")
print(f"ranges: s3/exact {min(r['s3'] / r['exact'] for r in rows):.4f}..{max(r['s3'] / r['exact'] for r in rows):.4f}; "
      f"s4/exact {min(r['s4'] / r['exact'] for r in rows):.4f}..{max(r['s4'] / r['exact'] for r in rows):.4f}; "
      f"slices/exact {min(r['slices_u'] / r['exact'] for r in rows):.4f}..{max(r['slices_u'] / r['exact'] for r in rows):.4f}")
print(f"loop slices vs exact slices (universal kappa): max relative difference "
      f"{max(abs(r['slices_u'] / (r['ex_sl'][8] + r['ex_sl'][10]) - 1) for r in rows):.2e}; with exact kappa "
      f"{max(abs(r['slices_x'] / (r['ex_sl'][8] + r['ex_sl'][10]) - 1) for r in rows):.2e}")

def fit(y, keep):
    tt = np.array([rows[i]['TT5'] for i in keep], float); yy = np.array([y[i] for i in keep])
    A = np.column_stack([np.ones(len(tt)), tt]); (a, b), *_ = np.linalg.lstsq(A, yy, rcond=None)
    res = yy - A @ [a, b]
    return b, -(a + b * 1518), res
allk = list(range(len(rows)))
dedup = [i for i, r in enumerate(rows) if r['label'] not in ('D23_10', 'D23_11')]
print("\nLaws (slope per TT5, Paley offset, max residual), 13 classes | 11 classes (converse duplicates D23_10, D23_11 removed):")
fits = {}
for name in ('exact', 's3', 'x4', 'dress', 's4', 'slices_u'):
    y = [r[name] for r in rows]
    f13 = fit(y, allk); f11 = fit(y, dedup); fits[name] = (f13, f11)
    print(f"   {name:8s}: slope {f13[0]:+.5e}, offset {f13[1]:+.5e}, max res {np.max(np.abs(f13[2])):.3e} | "
          f"slope {f11[0]:+.5e}, offset {f11[1]:+.5e}, max res {np.max(np.abs(f11[2])):.3e}")
for name in ('s4', 'slices_u'):
    for j, tag in ((0, '13'), (1, '11')):
        print(f"   {name} / exact ({tag} classes): slope {fits[name][j][0] / fits['exact'][j][0]:.4f}, offset {fits[name][j][1] / fits['exact'][j][1]:.4f}")
for j, tag, keep in ((0, '13', allk), (1, '11', dedup)):
    rex, rs4 = fits['exact'][j][2], fits['s4'][j][2]
    cc = np.corrcoef(rex, rs4)[0, 1]
    beta = float(np.dot(rex, rs4) / np.dot(rex, rex))
    print(f"   non-TT5 residuals ({tag} classes): corr(exact, s4) = {cc:+.6f}; regression s4-res on exact-res {beta:.4f}; "
          f"rms ratio {np.sqrt(np.mean(rs4 ** 2) / np.mean(rex ** 2)):.4f}; max-abs ratio {np.max(np.abs(rs4)) / np.max(np.abs(rex)):.4f}")
i1 = [i for i, r in enumerate(rows) if r['label'] == 'D23_1'][0]; i2 = [i for i, r in enumerate(rows) if r['label'] == 'D23_2'][0]
print(f"   D23_1 - D23_2: exact {rows[i1]['exact'] - rows[i2]['exact']:+.4e}; s4 {rows[i1]['s4'] - rows[i2]['s4']:+.4e}; "
      f"Delta c_10 split {float(rows[i2]['dc'][10] - rows[i1]['dc'][10]):+.4e}")

print("\nz^10 bullet: ratios per class")
for r in rows:
    print(f"   {r['label']:8s}: dX4^(10)/delta_10 = {r['d4_10'] / r['d3_10']:+.3f}; dX5^(10)/dX4^(10) = {r['d5_10'] / r['d4_10']:+.4f}")
# Paley offset per coefficient (exact), fits over 13 classes
print("\nPaley offset of exact Delta c_L (13-class line): " + ', '.join(
    f"L={L}: {fit([float(r['dc'][L]) for r in rows], allk)[1]:+.4e}" for L in (8, 10, 12, 14, 16)))
print(f"D23_0: Delta c_10 = {float(rows[0]['dc'][10]):+.10e}; delta_10 + dX4 + dX5 = {rows[0]['d3_10'] + rows[0]['d4_10'] + rows[0]['d5_10']:+.10e}")
print("Delta X_4(1) by V (D23_0): " + ', '.join(f"V={V}: {v:+.3e}" for V, v in sorted(rows[0]['x4byV'].items(), key=lambda t: int(t[0]))))

# power counting
q = np.array([19, 23, 31, 43.]); y = np.array([6340.7731, 6055.0577, 5739.4560, 5512.5521])
c1 = np.polyfit(1 / q[-2:], y[-2:], 1); c2 = np.polyfit(1 / q, y, 2); c1all = np.polyfit(1 / q, y, 1)
c2last3 = np.polyfit(1 / q[1:], y[1:], 2)
print(f"\nq^4 X_4(P_q) extrapolations: linear last two {c1[-1]:.1f}; quadratic all four {c2[-1]:.1f}; linear all four {c1all[-1]:.1f}; "
      f"quadratic last three {c2last3[-1]:.1f}")
x = d['X4'][iP]['value']; print(f"P_23: X_4(1) = {x:+.7e}, n^4 X_4 = {x * n ** 4:.4f}")
