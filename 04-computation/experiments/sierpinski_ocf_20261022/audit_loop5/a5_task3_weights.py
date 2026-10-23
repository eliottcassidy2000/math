# Audit (HYP-9371 section 8), task 3 (numerical part): the doubly regular Ising weights of loop_engine5.
#  (a) a = u^2/(1-u^2), b = u^2/(1+u^2), g = z/(1-u^2), h = -2u^2 z/((1+u^2)(1-u^2)), u^2 = -z^2/n, against
#      loop_engine4.drt_coeffs (alpha, gamma of the five basis functions), via B0 = (a, g), B1 = (a, -g),
#      B2 = (-a, h), B3 = (b, 0), B4 = (-a, -h);
#  (b) the same five pairs from the matrices themselves: for a real DRT (P_23 and D23_0), the 5 basis matrices
#      B_k(z Shat) equal alpha_k P + gamma_k Shat (least squares + residual);
#  (c) ops_drt edge tensor against the 5-basis map: E_drt[o_t, o_h, :] = sum_b E_basis5[o_t, o_h, b] * (alpha_b, gamma_b),
#      F+/- likewise;
#  (d) symbolic check (sympy) of the m1/m2 P and Shat coefficients from the scalar resolvents at lambda = +-i n^-1/2.
import os, sys, json
import numpy as np
import sympy as sp
D = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, D)
import loop_engine4 as LE4
import loop_engine5 as LE5
from drt_switch import code_to_S

zs = [0.3, 1.0, 1.7 + 0.4j, 2.2j, 0.9 * np.exp(0.7j)]
worst = 0.0
for n in (7, 11, 19, 23, 43):
    for z in zs:
        a, b, g, h = LE5.drt_funcs(n, np.array([z]))
        mine = np.array([[a[0], g[0]], [a[0], -g[0]], [-a[0], h[0]], [b[0], 0], [-a[0], -h[0]]])
        ref = LE4.drt_coeffs(n, complex(z))
        worst = max(worst, np.max(np.abs(mine - ref)))
print(f"(a) Ising weights vs loop_engine4.drt_coeffs, n in 7..43, 5 complex z: max |diff| = {worst:.2e}")

cls = {c['label']: c for c in json.load(open(os.path.join(D, 'drt23_classes.json')))}
for lab in ('P_23', 'D23_0'):
    S = code_to_S(cls[lab]['code']); n = S.shape[0]
    assert np.array_equal(S @ S.T, n * np.eye(n, dtype=int) - np.ones((n, n), dtype=int))
    Sh = S / n; P = np.eye(n) - np.ones((n, n)) / n
    wr = 0.0
    for z in zs:
        Bm = LE4.basis(Sh.astype(complex), complex(z))
        a, b, g, h = LE5.drt_funcs(n, np.array([z]))
        pairs = [(a[0], g[0]), (a[0], -g[0]), (-a[0], h[0]), (b[0], 0), (-a[0], -h[0])]
        for k in range(5):
            wr = max(wr, np.max(np.abs(Bm[k] - (pairs[k][0] * P + pairs[k][1] * Sh))) / max(1e-300, np.max(np.abs(Bm[k]))))
    print(f"(b) {lab}: basis matrices B_k(z Shat) = alpha_k P + gamma_k Shat, max relative residual {wr:.2e}")

ops5 = LE5.ops_basis5()
wc = 0.0
for n in (7, 23):
    for z in zs:
        od = LE5.ops_drt(n, np.array([z]))
        a, b, g, h = LE5.drt_funcs(n, np.array([z]))
        C = np.array([[a[0], g[0]], [a[0], -g[0]], [-a[0], h[0]], [b[0], 0], [-a[0], -h[0]]])   # (5, 2)
        wc = max(wc, np.max(np.abs(od['E'][0] - np.einsum('ijb,bl->ijl', ops5['E'], C))))
        wc = max(wc, np.max(np.abs(od['F+'][0] - ops5['F+'] @ C)), np.max(np.abs(od['F-'][0] - ops5['F-'] @ C)))
print(f"(c) ops_drt (E, F+, F-) equal to ops_basis5 composed with the 5-basis map: max |diff| = {wc:.2e}")

# (d) symbolic
z, n = sp.symbols('z n', positive=True)
w = 1 / sp.sqrt(n)
def parts(f):
    fp, fm = f(sp.I * w * z), f(-sp.I * w * z)
    return sp.simplify((fp + fm) / 2), sp.simplify((fp - fm) / (2 * sp.I * w))
u2 = -z ** 2 / n
A_ = u2 / (1 - u2); B_ = u2 / (1 + u2); G_ = z / (1 - u2); H_ = -2 * u2 * z / ((1 + u2) * (1 - u2))
R = {1: lambda x: x / (1 - x), -1: lambda x: -x / (1 + x)}
Cm = lambda x: (1 - x ** 2) / (1 + x ** 2)
ok = True
for ot in (1, -1):
    al, ga = parts(R[ot])                                            # m1 with traversal sign s = o_t
    ok &= sp.simplify(al - A_) == 0 and sp.simplify(ga - G_ * ot) == 0
    for oh in (1, -1):
        st, sh = ot, -oh
        al, ga = parts(lambda x: -R[st](x) * R[sh](x) * Cm(x))      # m2: Q_(s_t, s_h), s_t = o_t, s_h = -o_h
        ok &= sp.simplify(al - ((B_ - A_) + (B_ + A_) * ot * oh) / 2) == 0
        ok &= sp.simplify(ga - (H_ / 2) * (ot - oh)) == 0
print(f"(d) symbolic: m1 (P a, Shat g o_t) and m2 (P ((b-a)+(b+a) o_t o_h)/2, Shat (h/2)(o_t-o_h)) from the resolvents: {'all equal' if ok else 'MISMATCH'}")
