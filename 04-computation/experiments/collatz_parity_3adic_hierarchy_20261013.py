"""
collatz_parity_3adic_hierarchy_20261013.py -- the deterministic 3-adic refinement of the twisted transfer law.
Classes r mod 3^k.  Even n = r: n/2 = r * inv(2) (permutation P).  Odd n = r: n mod 2*3^k is determined, so
(3n+1)/2 mod 3^k = c(r) is determined (function matrix Q).  For the character eps_e, eps_o on even / odd steps,
the exact bin-mean identities read m_r(u) = (eps_e/2) m_{P r}(u-1) + (eps_o/2) m_{Q r}(u+theta-1), and a mode
v e^{su} needs 1 in spec M(s), M(s) = (eps_e a/2) P + (eps_o b/2) Q, a = e^{-s}, b = e^{s(theta-1)}.
The constant vector is invariant (P 1 = Q 1 = 1), giving the scalar law (eps_e a + eps_o b)/2 = 1 at every level.
Part 1: sympy factorisation of det(I - M) for k = 1, 2 (all three parities).
Part 2: numeric scan for k = 2, 3, 4: smallest singular value of I - M(s) on a grid s = w + 2 pi i f,
        f in [0, 30], w in [-1.2, 0]; local minima below a threshold are near-roots; compare with scalar roots.
"""
import sys, numpy as np, sympy as sp
THETA = np.log2(3.0)

def PQ(k):
    m = 3 ** k
    inv2 = pow(2, -1, m)
    P = np.zeros((m, m)); Q = np.zeros((m, m))
    for r in range(m):
        P[r, (r * inv2) % m] = 1
        n = r if r % 2 == 1 else r + m          # odd representative of the class r mod 3^k
        Q[r, ((3 * n + 1) // 2) % m] = 1
    return P, Q

chars = {'tau': (-1, -1), 'omega': (1, -1), 'C': (-1, 1)}
print('# Part 1: symbolic characteristic polynomials det(I - M(a,b)), factored')
a, b = sp.symbols('a b')
for k in (1, 2):
    P, Q = PQ(k)
    for name, (ee, eo) in chars.items():
        M = sp.Matrix(P) * (ee * a / 2) + sp.Matrix(Q) * (eo * b / 2)
        det = sp.factor(sp.simplify((sp.eye(3 ** k) - M).det()))
        print(f' k={k} {name:5s}: {det}')

print('\n# Part 2: numeric near-roots of det(I - M(s)) for k = 2, 3, 4 (smallest singular value < 0.02)')
def scalar_roots(ee, eo):
    # roots of (ee e^{-s} + eo e^{s(theta-1)})/2 = 1 near the imaginary axis, by Newton on a grid
    out = []
    for f in np.arange(0.5, 30.5, 0.5):
        s = 2j * np.pi * f - 0.05
        for _ in range(60):
            F = (ee * np.exp(-s) + eo * np.exp(s * (THETA - 1))) / 2 - 1
            dF = (-ee * np.exp(-s) + eo * (THETA - 1) * np.exp(s * (THETA - 1))) / 2
            s = s - F / dF
        if -1.2 < s.real <= 0.001 and abs(F) < 1e-10:
            out.append(s)
    out = sorted(set(round(x.imag / (2 * np.pi), 3) + 1j * round(x.real, 4) for x in out), key=lambda z: z.real)
    return [(z.real, z.imag) for z in out]   # (f, Re s)
for k in (2, 3, 4):
    P, Q = PQ(k); m = 3 ** k
    I = np.eye(m)
    for name, (ee, eo) in chars.items():
        fs = np.arange(0.0, 30.0, 0.01); ws = np.arange(-1.2, 0.001, 0.01)
        best = {}
        for w in ws:
            for f in fs:
                s = w + 2j * np.pi * f
                Ms = (ee * np.exp(-s) / 2) * P + (eo * np.exp(s * (THETA - 1)) / 2) * Q
                sv = np.linalg.svd(I - Ms, compute_uv=False)[-1]
                key = round(f, 1)
                if sv < best.get(key, (1.0, 0, 0))[0]:
                    best[key] = (sv, f, w)
        near = sorted([(f, w, sv) for (sv, f, w) in best.values() if sv < 0.02], key=lambda t: t[0])
        sc = scalar_roots(ee, eo)
        print(f' k={k} {name:5s}: near-roots (f, Re s, sv): ' + ' '.join(f'({f:.2f},{w:+.2f},{sv:.3f})' for f, w, sv in near[:12]))
        print(f'            scalar roots (f, Re s): ' + ' '.join(f'({f:.2f},{w:+.2f})' for f, w in sc[:12]), flush=True)
print('done')
