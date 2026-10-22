# Anatomy of the Hamiltonian-path differences among doubly regular 23-tournaments, from the exact forest numbers
# (drt23_forest_exact.json) and the four-loop engine (drt23_four_loop.json).
#   c_L(T) = [z^L] log Phi_T(z) exactly (Phi_T = sum_L N_L (z/n)^L), Delta c_L = c_L(D) - c_L(P_23).
#   Theory: Delta c_L = sum_(e >= 3) Delta X_e^(L); at L = 8 only the four-loop term can differ (excess >= 4 at
#   L = 8 has <= 4 blocks, excess <= 2 reduces to <= 4-vertex cores), so Delta c_8 = delta_8 must hold exactly.
#   Response: kappa = sum_(L <= n-1) phi_L g(L), g(L) = n^L/(n)_L, so to first order in Delta log Phi
#       Delta log kappa = sum_j w_j Delta c_j,   w_j = sum_(M <= n-1-j) phi_M(P_23) g(M + j) / kappa(P_23),
#   exact weights from the exact P_23 measure (no loop approximation).  Output: drt23_anatomy.out.
# mac-mini-2026-10-22-loop4.
import json
from fractions import Fraction
from math import factorial
import numpy as np
from mpmath import mp, mpf, log
mp.dps = 50

n = 23; LM = 30
F = {r['class']: r for r in json.load(open('drt23_forest_exact.json'))}
fl = json.load(open('drt23_four_loop.json'))
rows = {r['class']: r for r in fl['rows']}
iP = [k for k, r in rows.items() if r['label'] == 'P_23'][0]

def series(N):
    phi = [Fraction(int(N[L]), n ** L) if L < n else Fraction(0) for L in range(LM + 1)]
    assert phi[0] == 1
    c = [Fraction(0)] * (LM + 1)
    for L in range(1, LM + 1):
        s = L * phi[L]
        for k in range(1, L): s -= k * c[k] * phi[L - k]
        c[L] = s / L
    return phi, c

g = [Fraction(n ** L, factorial(n) // factorial(n - L)) for L in range(n)]
ser = {k: series(F[k]['N']) for k in F}
phiP, cP = ser[iP]
kapP = sum(phiP[L] * g[L] for L in range(n))
assert kapP == Fraction(int(F[iP]['H']) * 2 ** (n - 1), factorial(n))
w = [sum(phiP[M] * g[M + j] for M in range(n - j)) / kapP for j in range(n)]

print("P_23 exact measure: phi_L = N_L/23^L for even L <= 22:")
print('   ' + ' '.join(f"{float(phiP[L]):.4e}" for L in range(0, n, 2)))
print("exact log Phi coefficients c_L(P_23):", ' '.join(f"c{L}={float(cP[L]):+.4e}" for L in range(2, 21, 2)))
print("four-loop X_3^(L)(P_23) for comparison:", ' '.join(f"{x:+.4e}" for x in fl['X3_P23_coeffs'][6:17:2]), "(L = 6..16)")
print("exact response weights w_j (P_23 measure):", ' '.join(f"w{j}={float(w[j]):.4g}" for j in range(0, n, 2)))
D1 = 2 * n / (n + 1)
print(f"   strict next-order weights 1 + (j D_1 + j(j-1)/2)/n, D_1 = 2n/(n+1): " + ' '.join(f"w{j}={1 + (j * D1 + j * (j - 1) / 2) / n:.4g}" for j in range(0, 17, 2)))

print("\nPer class: Delta c_L exact vs four-loop delta_L (L = 8, 10, 12, 14), and the response decomposition")
tab = []
for k in sorted(F):
    if k == iP: continue
    phi, c = ser[k]
    dc = [c[L] - cP[L] for L in range(LM + 1)]
    delta = rows[k]['coeffs']
    exact = float(log(mpf(int(F[k]['H'])) / mpf(int(F[iP]['H']))))
    lin = float(sum(w[j] * dc[j] for j in range(n)))
    # per forest length: Delta H/H = sum_L g(L) Delta phi_L / kappa_P; four-loop: Delta phi_L = [z^L] Phi_P Delta X_3
    dphi = [float((phi[L] - phiP[L]) * g[L] / kapP) for L in range(n)]
    dphi4 = [float(g[L] / kapP) * sum(delta[j] * float(phiP[L - j]) for j in range(L + 1)) for L in range(n)]
    parts = {j: float(w[j] * dc[j]) for j in range(0, n, 2)}
    four_exactw = sum(float(w[j]) * delta[j] for j in range(min(n, len(delta))))
    four_strict = sum(delta[j] for j in range(len(delta)))
    four_next = sum((1 + (j * D1 + j * (j - 1) / 2) / n) * delta[j] for j in range(len(delta)))
    rec = {'class': k, 'label': rows[k]['label'], 'TT5': rows[k]['TT5'], 'exact': float(exact), 'linear_response': float(lin),
           'dc': [float(x) for x in dc[:LM + 1]], 'delta': delta[:LM + 1], 'parts': parts,
           'four_strict': four_strict, 'four_next': four_next, 'four_exactw': four_exactw,
           'per_L_exact': dphi, 'per_L_four': dphi4, 'relH': float(Fraction(int(F[k]['H']), int(F[iP]['H'])) - 1)}
    tab.append(rec)
    print(f"{rows[k]['label']:7s} TT5 {rows[k]['TT5']}: exact Delta log H {float(exact):+.6e}; linear response sum_j w_j Delta c_j {float(lin):+.6e} (diff {float(lin - exact):+.1e})")
    print("        Delta c_L exact : " + ' '.join(f"L{L}:{float(dc[L]):+.5e}" for L in range(6, 17, 2)))
    print("        four-loop delta : " + ' '.join(f"L{L}:{delta[L]:+.5e}" for L in range(6, 17, 2)))
    print("        w_j Delta c_j   : " + ' '.join(f"j{j}:{parts[j]:+.3e}" for j in range(8, n, 2)))
    print(f"        four-loop: strict {four_strict:+.4e}, next-order weights {four_next:+.4e}, exact weights {four_exactw:+.4e}")
    print("        per L, g(L) Delta phi_L/kappa exact : " + ' '.join(f"L{L}:{dphi[L]:+.3e}" for L in range(8, n, 2)) + f"  (sum {sum(dphi):+.4e} = H_D/H_P - 1 = {rec['relH']:+.4e})")
    print("        per L, four-loop [z^L] Phi_P dX_3   : " + ' '.join(f"L{L}:{dphi4[L]:+.3e}" for L in range(8, n, 2)) + f"  (sum {sum(dphi4):+.4e})")

tt = np.array([r['TT5'] for r in tab], float)
print("\nLinear laws in TT5 over the 12 non-Paley classes (slope per transitive 5-set; max residual; value at Paley's TT5 = 1518):")
for r in tab:
    r['four_perL_sum'] = sum(r['per_L_four']); r['exact_upto12'] = sum(r['per_L_exact'][:13]); r['four_upto12'] = sum(r['per_L_four'][:13])
    r['exact_14up'] = sum(r['per_L_exact'][14:])
for name in ('exact', 'linear_response', 'four_strict', 'four_next', 'four_exactw', 'four_perL_sum', 'exact_upto12', 'four_upto12', 'exact_14up'):
    y = np.array([r[name] for r in tab])
    A = np.column_stack([np.ones(len(tt)), tt]); (a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    print(f"   {name:16s}: slope {b:+.4e}, max residual {np.max(np.abs(y - A @ [a, b])):.2e}, line at TT5 = 1518: {a + b * 1518:+.3e}")
for L in (8, 10, 12, 14, 16):
    y = np.array([r['dc'][L] for r in tab]); A = np.column_stack([np.ones(len(tt)), tt]); (a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    yd = np.array([r['delta'][L] for r in tab]); (ad, bd), *_ = np.linalg.lstsq(A, yd, rcond=None)
    print(f"   Delta c_{L:<2d} exact: slope {b:+.4e}, max residual {np.max(np.abs(y - A @ [a, b])):.2e}, line at 1518 {a + b * 1518:+.2e};   four-loop delta_{L}: slope {bd:+.4e}")
json.dump({'w': [float(x) for x in w], 'rows': tab}, open('drt23_anatomy.json', 'w'), indent=1)

# Where do Paley's offset and the non-TT5 scatter live?  Per forest length L, fit g(L) Delta phi_L / kappa over the
# 12 non-Paley classes by a + b TT5 and report the Paley offset (its value, 0, minus the line at TT5 = 1518).
print("\nPer forest length: TT5 law over non-Paley classes and Paley's offset from it (relative H units)")
A = np.column_stack([np.ones(len(tt)), tt])
for L in range(8, n, 2):
    y = np.array([r['per_L_exact'][L] for r in tab]); (a, b), *_ = np.linalg.lstsq(A, y, rcond=None)
    print(f"   L = {L:2d}: slope {b:+.4e}, max residual {np.max(np.abs(y - A @ [a, b])):.2e}, Paley offset {-(a + b * 1518):+.3e}")
# The alternation Delta c_(L+2)/Delta c_L ~ -0.4: zeros of the forest polynomial Phi(u), u = z^2 (degree 11 in u).
for k in (iP, 0):
    phi = ser[k][0]
    coef = [float(phi[2 * m]) for m in range(12)]
    roots = np.roots(coef[::-1])
    r0 = roots[np.argmin(np.abs(roots))]
    print(f"{rows[k]['label']}: zero of Phi(u) nearest 0: u = {r0:.6f} (|u| = {abs(r0):.4f}, -1/u = {(-1 / r0).real:+.4f}); "
          f"loop terms X_e(u) are singular only at u = +-n = +-{n}")
