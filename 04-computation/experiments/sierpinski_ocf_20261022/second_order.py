# Second-order (three-loop) prediction of kappa/rho for regular tournaments (THM-4757 / HYP-9365 framework):
#   kappa = sum_L phi_L g_n(L),  g_n(L) = n^L/(n)_L = exp(a(L)/n + b(L)/n^2 + ...),  a = C(L,2), b = L(L-1)(2L-1)/12,
#   log Phi(z) = X_0(z) + X_1(z) + X_2(z) + ...,
#   kappa/rho = 1 + [X_1 + E_0 a/n] + [X_2 + X_1^2/2 + X_1 E_0 a/n + delta/n + (E_0 a^2/2 + E_0 b)/n^2] + O(n^-3),
# where E_0 is the one-loop measure (cumulants D_k = (z d/dz)^k X_0 at z = 1) and delta is the first-order shift of
# E[a] from X_1 (cumulants of X_1).  X_1 by the closed formula, X_2 by the loop engine (15 three-loop cores).
# mac-mini-2026-10-22-sierp.
import sys, math
import numpy as np
sys.path.insert(0, '.')
from loop_engine import cores, decompositions_fast, automorphisms, propagator2, hom_c
import loop_engine_fast as LEF

def X1_of(Sh, z):
    n = Sh.shape[0]; I = np.eye(n); X = z * Sh; X2 = X @ X
    G = -X2 @ np.linalg.inv(I + X2); H = X2 @ np.linalg.inv(I - X2); O = X @ np.linalg.inv(I - X2)
    gd, hd = np.diag(G), np.diag(H); d = hd - gd
    return np.sum(2 * gd * hd - 1.5 * gd ** 2) + np.sum(G * (H * H - O * O - (2 / 3) * G * G)) - d @ G @ d

def X0_of(x, z):
    return 0.5 * np.sum(np.log((1 + (z * x) ** 2) / (1 - (z * x) ** 2)))

def cumulants(f, kmax=4, r=0.05, K=32):
    """D_k = (z d/dz)^k f at z = 1 = d^k/dtheta^k f(e^theta) at 0, by a Cauchy sum in theta"""
    th = r * np.exp(2j * np.pi * np.arange(K) / K)
    vals = np.array([f(np.exp(t)) for t in th])
    return [ (math.factorial(k) * np.mean(vals * th ** (-k))).real for k in range(kmax + 1)]

_DECS = None
def X2_of(Sh, z=1.0):
    global _DECS
    cs = cores(2)
    if _DECS is None: _DECS = [decompositions_fast(V, e) for V, e in cs]
    tot = 0
    for (V, edges), dec in zip(cs, _DECS):
        G = automorphisms(V, edges); pc = {}; sub = 0
        for props, w in dec.items():
            if w == 0: continue
            mats = []
            for dsc in props:
                if dsc not in pc: pc[dsc] = propagator2(dsc, Sh.astype(complex), z)
                mats.append(pc[dsc])
            sub += w * hom_c(V, edges, mats)
        tot += sub / G
    return tot.real

def predict2(S):
    n = S.shape[0]; Sh = S / n
    x = np.linalg.eigvals(Sh).imag
    D0 = cumulants(lambda z: X0_of(x, z))
    D1 = cumulants(lambda z: X1_of(Sh.astype(complex), z))
    k1, k2, k3, k4 = D0[1], D0[2], D0[3], D0[4]
    EL = k1; EL2 = k2 + k1 ** 2; EL3 = k3 + 3 * k2 * k1 + k1 ** 3
    EL4 = k4 + 4 * k3 * k1 + 3 * k2 ** 2 + 6 * k2 * k1 ** 2 + k1 ** 4
    Ea = (EL2 - EL) / 2
    Ea2 = (EL4 - 2 * EL3 + EL2) / 4
    Eb = (2 * EL3 - 3 * EL2 + EL) / 12
    delta = (D1[2] + 2 * k1 * D1[1] - D1[1]) / 2
    X1 = D1[0]; X2 = float(np.real(LEF.XB(Sh, 2)))
    first = X1 + Ea / n
    second = X2 + X1 ** 2 / 2 + X1 * Ea / n + delta / n + (Ea2 / 2 + Eb) / n ** 2
    rho = float(np.sqrt(np.prod((1 + x ** 2) / (1 - x ** 2))))
    return rho, first, second, X2

if __name__ == '__main__':
    from mpmath import mpf, factorial
    from cluster_check import tournament
    tests = [('paley', 7, 189), ('paley', 11, 95095), ('paley', 19, 1172695746915), ('paley', 23, 15760206976379349),
             ('paley', 31, 19876729650849384218015865),
             ('rot', 11, 93027), ('rot', 15, 198464295), ('rot', 19, 1184212824763), ('rot', 23, 16011537490557279),
             ('rot', 27, 426990256001863046739), ('rot', 31, 20351875726012821643802839), ('rot', 35, 1610822294345031078388661430315)]
    for kind, n, H in tests:
        A = tournament(kind, n); S = A - A.T
        rho, first, second, X2 = predict2(S)
        kap = float(mpf(H) * mpf(2) ** (n - 1) / factorial(n))
        exact = kap / rho - 1
        print(f"{kind} n={n}: exact kappa/rho - 1 = {exact:+.6e}; 1st order {first:+.6e} (resid {exact - first:+.3e}, n^2 resid {n * n * (exact - first):+.3f}); "
              f"+2nd order {first + second:+.6e} (resid {exact - first - second:+.3e}, n^3 resid {n ** 3 * (exact - first - second):+.3f}); X_2 = {X2:+.4e}", flush=True)
