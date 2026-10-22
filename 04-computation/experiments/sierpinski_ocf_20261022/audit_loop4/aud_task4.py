# AUDIT task 4/5: recompute the numbers of HYP-9371 Sections 3, 4 (and 6) from the AUDIT's own exact forest numbers
# (aud_task1_N.json for n = 23; own run of aud_fwdpoly for the n = 19 pair), plus the engine series stored in
# drt23_four_loop.json (four-loop delta_L; the engine itself is checked in tasks 2/3).
import json, sys, subprocess
from fractions import Fraction
from math import factorial, comb
import numpy as np
import mpmath as mp
mp.mp.dps = 60
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE + '/audit_loop4')
from aud_task1 import F_mod, lagrange_coeffs, crt2, forest_N, P1, P2, code_adj, tt

own = json.load(open(f'{HERE}/audit_loop4/aud_task1_N.json'))
cls = json.load(open(f'{HERE}/drt23_classes.json'))
fl = json.load(open(f'{HERE}/drt23_four_loop.json'))
rows = {r['class']: r for r in fl['rows']}
n = 23
N = {int(k): [int(x) for x in v['N']] for k, v in own.items()}
lab = {int(k): v['label'] for k, v in own.items()}
TT5 = {int(k): v['TT5'] for k, v in own.items()}
H = {c['class']: int(c['H']) for c in cls}
iP = [k for k in lab if lab[k] == 'P_23'][0]
nonP = [k for k in sorted(lab) if k != iP]

def logser(Nl, n, LM):
    phi = [Fraction(Nl[L], n ** L) if L < n else Fraction(0) for L in range(LM + 1)]
    c = [Fraction(0)] * (LM + 1)
    for L in range(1, LM + 1):
        s = L * phi[L]
        for k in range(1, L): s -= k * c[k] * phi[L - k]
        c[L] = s / L
    return phi, c
def gL(n, L): return Fraction(n ** L, factorial(n) // factorial(n - L))
def fit(x, y):
    A = np.column_stack([np.ones(len(x)), np.array(x, float)]); (a, b), *_ = np.linalg.lstsq(A, np.array(y, float), rcond=None)
    return a, b, float(np.max(np.abs(np.array(y, float) - A @ [a, b])))

LM = 30
ser = {k: logser(N[k], n, LM) for k in N}
phiP, cP = ser[iP]
kapP = sum(phiP[L] * gL(n, L) for L in range(n))
assert kapP == Fraction(H[iP] * 2 ** (n - 1), factorial(n))
print("== Section 3: exact log H differences, strict / dressed (engine values from drt23_four_loop.json) ==")
D1 = 2 * n / (n + 1)
ex = {k: float(mp.log(mp.mpf(H[k]) / H[iP])) for k in nonP}
for k in nonP:
    r = rows[k]
    co = r['coeffs']
    dressed24 = sum((1 + (j * D1 + j * (j - 1) / 2) / n) * co[j] for j in range(len(co)))
    print(f"   {lab[k]:7s} TT5 {TT5[k]}: exact {ex[k]:+.5e} (json {r['exact']:+.5e}); strict {r['strict']:+.4e}; dressed {r['dressed']:+.4e} "
          f"(recomputed from 24 stored coefficients {dressed24:+.4e}); exact/strict {ex[k] / r['strict']:.3f}; strict/exact {r['strict'] / ex[k]:.3f}; dressed/exact {r['dressed'] / ex[k]:.3f}")
for name, ys in (('exact', [ex[k] for k in nonP]), ('strict', [rows[k]['strict'] for k in nonP]), ('dressed', [rows[k]['dressed'] for k in nonP])):
    a, b, mr = fit([TT5[k] for k in nonP], ys)
    print(f"   TT5 law {name:8s}: slope {b:.5e} (x H(P_23) = {b * H[iP]:.5e}), line at 1518 {a + b * 1518:+.4e}, max residual {mr:.2e}")
ae, be, _ = fit([TT5[k] for k in nonP], [ex[k] for k in nonP])
_, bs, _ = fit([TT5[k] for k in nonP], [rows[k]['strict'] for k in nonP])
_, bd, _ = fit([TT5[k] for k in nonP], [rows[k]['dressed'] for k in nonP])
print(f"   strict/exact slope {bs / be:.4f}, dressed/exact slope {bd / be:.4f}, dressed/strict slope {bd / bs:.4f}; "
      f"factor 1 + (8 D_1 + 28)/n = {1 + (8 * D1 + 28) / n:.4f}")
relH = {k: float(Fraction(H[k], H[iP]) - 1) for k in nonP}
print(f"   D23_1 - D23_2 relative H difference: {(H[1] - H[2]) / H[iP]:+.3e}; exact H monotone in TT5: "
      f"{all(H[a] <= H[b] for a in nonP for b in nonP if TT5[a] < TT5[b])}")

print("\n== Section 4: exact Delta c_L (own N) vs four-loop delta_L (engine) ==")
for k in nonP:
    dc = [ser[k][1][L] - cP[L] for L in range(LM + 1)]
    d = rows[k]['coeffs']
    rel8 = abs(float(dc[8]) - d[8]) / abs(float(dc[8]))
    print(f"   {lab[k]:7s}: Delta c_8 = {float(dc[8]):+.12e} = Delta N_8/23^8 exactly: {dc[8] == Fraction(N[k][8] - N[iP][8], n ** 8)}; engine delta_8 {d[8]:+.12e} (rel diff {rel8:.1e}); "
          f"(Delta c_10 - delta_10)/delta_8 = {(float(dc[10]) - d[10]) / d[8]:+.4f}; delta_10/delta_8 = {d[10] / d[8]:+.4f}; Delta c_6 = {float(dc[6]):.1e}")
k = 0
dc0 = [ser[k][1][L] - cP[L] for L in range(LM + 1)]
print("   D23_0 Delta c_L: " + ', '.join(f"L{L} {float(dc0[L]):+.4e}" for L in range(8, 17, 2)) + "; ratios " + ', '.join(f"{float(dc0[L + 2] / dc0[L]):+.3f}" for L in range(8, 15, 2)))
print("   D23_0 delta_L : " + ', '.join(f"L{L} {rows[0]['coeffs'][L]:+.4e}" for L in range(8, 17, 2)) + f"; Delta c_10 - delta_10 = {float(dc0[10]) - rows[0]['coeffs'][10]:+.4e}")
for kk in (iP, 0):
    phi = ser[kk][0]
    coef = [mp.mpf(phi[2 * m].numerator) / phi[2 * m].denominator for m in range(12)]
    roots = mp.polyroots(coef[::-1], maxsteps=400, extraprec=400)
    r0 = min(roots, key=lambda r: abs(r))
    print(f"   {lab[kk]}: zero of Phi(u) (u = z^2) nearest 0: {mp.nstr(r0, 8)}; 1/u0 = {mp.nstr(1 / r0, 6)}; radius in z = {mp.nstr(mp.sqrt(abs(r0)), 6)}")
w = [sum(phiP[M] * gL(n, M + j) for M in range(n - j)) / kapP for j in range(n)]
print(f"   exact response weights: w_8 = {float(w[8]):.4g}, w_10 = {float(w[10]):.4g}, w_22 = {float(w[22]):.4g}")
four_ex = sum(float(w[j]) * rows[0]['coeffs'][j] for j in range(n))
print(f"   four-loop Delta X_3 through exact P_23 weights (D23_0): {four_ex:+.4e} vs exact {ex[0]:+.4e} (ratio {four_ex / ex[0]:.1f})")

def perL(k, n, NP, NK, kap):
    return [float((Fraction(NK[L] - NP[L], n ** L)) * gL(n, L) / kap) for L in range(n)]
print("\n   shares of Delta H/H by forest length (own N):")
shares8 = []
for k in nonP:
    pl = perL(k, n, N[iP], N[k], kapP); tot = sum(pl)
    assert abs(tot - relH[k]) < 1e-15
    shares8.append(pl[8] / tot)
    if k in (0, 1, 4):
        print(f"   {lab[k]:7s}: L8 {pl[8] / tot:.1%}, L10 {pl[10] / tot:.1%}, L12 {pl[12] / tot:.1%}, L14 {pl[14] / tot:.2%}, L>=16 {sum(pl[16:]) / tot:.2%}")
print(f"   L = 8 share over the 12 classes: {min(shares8):.1%} .. {max(shares8):.1%}")
print("   Paley offset by forest length (fit a + b TT5 of g(L) Delta phi_L/kappa over 12 non-Paley classes; offset = 0 - line(1518)):")
offs = {}
for L in range(8, n, 2):
    a, b, mr = fit([TT5[k] for k in nonP], [perL(k, n, N[iP], N[k], kapP)[L] for k in nonP])
    offs[L] = -(a + b * 1518)
tot = sum(offs.values())
print("   " + ', '.join(f"L{L}: {offs[L]:+.3e} ({offs[L] / tot:.1%})" for L in offs) + f"; total {tot:+.4e}")

print("\n== n = 19: own forest numbers for P_19 and D_19.1 ==")
codes19 = {}
for line in open(f'{HERE}/drt_results.out'):
    t = line.split()
    if len(t) > 6 and t[0] == '19': codes19[t[1]] = (int(t[2]), t[6])
n19 = 19; R = {}
for name, (Hh, code) in codes19.items():
    v1 = F_mod(code, P1, range(n19)); v2 = F_mod(code, P2, range(n19))
    f = [crt2(a, b) for a, b in zip(lagrange_coeffs([v1[x] for x in range(n19)], P1), lagrange_coeffs([v2[x] for x in range(n19)], P2))]
    assert f[0] == Hh and sum(f) == factorial(n19) and f == f[::-1]
    Nn = forest_N(f, n19); nn, out = code_adj(code)
    R[name] = dict(H=Hh, N=Nn, TT5=tt(5, (1 << n19) - 1, out))
    kap = sum(Fraction(Nn[L], factorial(n19) // factorial(n19 - L)) for L in range(n19))
    print(f"   {name}: H {Hh}, own TT5 {R[name]['TT5']}, N_8 = {Nn[8]}, N_10 = {Nn[10]}, forest identity exact: {kap == Fraction(Hh * 2 ** 18, factorial(n19))}")
P19, D19 = R['P_19'], R['D_19.1']
dT = D19['TT5'] - P19['TT5']
print(f"   Delta TT5 = {dT}; Delta N_2,4,6 = {[D19['N'][L] - P19['N'][L] for L in (2, 4, 6)]}; Delta N_8 = {D19['N'][8] - P19['N'][8]} = 7168 x {Fraction(D19['N'][8] - P19['N'][8], 7168)}")
kap19 = sum(Fraction(P19['N'][L], factorial(n19) // factorial(n19 - L)) for L in range(n19))
pl = perL(None, n19, P19['N'], D19['N'], kap19); tot = sum(pl)
print(f"   exact H(D_19)/H(P_19) - 1 = {float(Fraction(D19['H'], P19['H']) - 1):+.4e} (sum of slices {tot:+.4e}); shares: L8 {pl[8] / tot:.1%}, L10 {pl[10] / tot:.1%}, "
      f"L12 {pl[12] / tot:.2%}, L14 {pl[14] / tot:.2%}, L>=16 {sum(pl[16:]) / tot:.2%}")
