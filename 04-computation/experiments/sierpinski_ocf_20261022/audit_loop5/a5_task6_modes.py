# Audit (HYP-9371 section 8), checks of engine modes used in section 8 that the session's outputs do not document:
#  (a) "In matrix mode, X_1 ... X_3 equal loop_engine4 to 1e-16 on a non-vertex-transitive regular 11-tournament"
#      (no output file in the repo): loop_engine5.X_general vs loop_engine4.X_e at three complex z, on the random
#      regular 11-tournament of loop_engine.py's self-test (seed 5, 4000 cyclic-triangle reversals of R_11);
#  (b) the DRT two-label mode used for q^4 X_4(P_q) and drt23_five_loop.py, checked in ABSOLUTE value (the z^10
#      identity only tests class differences): X_4 of P_7 (a DRT) from ops_drt + {P, Shat} homomorphism tensors,
#      z-series by Cauchy sums (K = 64, r = 0.35), against the exact blow-up scaling to z^14; and X_4(1) of P_7, P_11
#      in drt mode against matrix (general) mode at z = 1.
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.dirname(HERE)
sys.path.insert(0, D); sys.path.insert(0, os.path.join(D, 'audit2'))
import loop_engine4 as LE4
import loop_engine5 as LE5
from b_mscaling import mscaling
from common import rot_A, paley_A, skew

# (a)
rng = np.random.default_rng(5); n = 11; d = 5
A = np.array([[1 if 0 < (j - i) % n <= d else 0 for j in range(n)] for i in range(n)])
for _ in range(4000):
    u, v, w = rng.choice(n, 3, replace=False)
    if A[u, v] and A[v, w] and A[w, u]: A[u, v] = A[v, w] = A[w, u] = 0; A[v, u] = A[w, v] = A[u, w] = 1
S = (A - A.T); Sh = S / n
S2d = np.diag(np.linalg.matrix_power(S, 2)); S4d = np.diag(np.linalg.matrix_power(S, 4))
print(f"(a) random regular 11-tournament: scores {sorted(set(A.sum(1)))}, diag S^4 values {sorted(set(S4d.tolist()))} "
      f"(non-constant => not vertex-transitive)", flush=True)
zs = np.array([0.4, 0.7 * np.exp(0.9j), 1.0 + 0.3j])
for e in (1, 2, 3):
    t0 = time.time()
    v5 = LE5.X_general(Sh, e, zs)
    v4 = np.array([LE4.X_e(Sh, e, complex(z)) for z in zs])
    print(f"   X_{e}: |matrix mode - loop_engine4| = {np.max(np.abs(v5 - v4)):.1e} (|X| up to {np.max(np.abs(v4)):.2e}) ({time.time() - t0:.0f}s)", flush=True)

# (b) DRT mode in absolute value
def drt_X(S, excess, zs):
    n = S.shape[0]
    ops = LE5.ops_drt(n, zs)
    B = np.array([np.eye(n) - np.ones((n, n)) / n, S / n])
    tot = np.zeros(len(zs), dtype=complex)
    for V, edges, aut in LE5.load_cores(excess):
        W2 = LE5.SignCore(V, edges, aut).contract('drt', ops).reshape(len(zs), -1)
        tot += W2 @ LE4.Network(V, edges).contract(B).reshape(-1) / aut
    return tot

S7 = skew(paley_A(7)); assert LE4.check_drt(S7)
K, r = 64, 0.35
zc = r * np.exp(2j * np.pi * np.arange(K) / K)
t0 = time.time()
vals = drt_X(S7, 4, np.concatenate([zc, [1.0]]))
co = [(np.mean(vals[:K] * np.exp(-2j * np.pi * np.arange(K) * L / K)) / r ** L).real for L in range(15)]
X = mscaling(S7, 14, verbose=False)
print(f"(b) P_7, X_4 in DRT two-label mode (all 1076 cores, {time.time() - t0:.0f}s) against exact blow-up scaling:", flush=True)
for L in range(2, 15, 2):
    ex = float(X[L].get(4, 0))
    print(f"   L={L:2d}: drt mode {co[L]:+.14e}  exact {ex:+.14e}  diff {co[L] - ex:+.1e}", flush=True)
t0 = time.time(); g7 = LE5.X_general(S7 / 7, 4, np.array([1.0]))[0]
print(f"   X_4(1) of P_7: drt mode {vals[K].real:+.15e}, matrix mode {g7.real:+.15e}, diff {vals[K].real - g7.real:+.1e} ({time.time() - t0:.0f}s)", flush=True)
S11 = skew(paley_A(11)); assert LE4.check_drt(S11)
t0 = time.time(); d11 = drt_X(S11, 4, np.array([1.0]))[0]; t1 = time.time()
g11 = LE5.X_general(S11 / 11, 4, np.array([1.0]))[0]
print(f"   X_4(1) of P_11: drt mode {d11.real:+.15e}, matrix mode {g11.real:+.15e}, diff {d11.real - g11.real:+.1e} "
      f"(q^4 X_4 = {11 ** 4 * d11.real:.4f}; {t1 - t0:.0f}s + {time.time() - t1:.0f}s)", flush=True)
