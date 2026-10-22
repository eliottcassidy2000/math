# Validation of the four-loop (excess-3) engine loop_engine4.py: z-series coefficients of X_3 (Cauchy sums of the
# tensor engine) against two independent exact computations that use no core/kernel theory:
#   (M) blow-up scaling, audit2/b_mscaling.py: log Psi(Shat (x) J_m/m) = sum_e m^(-e) X_e(Shat), exact rationals;
#   (B) set-partition brute force, audit2/b_partitions.py (cluster_check.py's enumeration with pruning), exact
#       rationals, at L = 6, 8 (the task's "lowest orders").
# Tournaments: all three regular 7-tournaments (R_7, P_7, the non-circulant one) to z^12 (M); R_9 and two
# non-vertex-transitive regular 9-tournaments to z^10 (M) and z^8 (B); P_11 and a non-vertex-transitive regular
# 11-tournament to z^8 (M).  Also: X_1, X_2 of the new engine against the closed form and the old engine.
# Output: excess3_validate.out.  mac-mini-2026-10-22-loop4.
import sys, time
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'audit2')
import loop_engine4 as LE4
import loop_engine_fast as LEF
from second_order import X1_of
from b_mscaling import mscaling
from b_partitions import cluster_sums
from common import rot_A, paley_A, skew, regular_classes, code_to_A

def vertex_transitive_hint(S):
    """False if diag(S^4) or diag(S^6) is non-constant (then T is certainly not vertex-transitive)"""
    S4 = np.linalg.matrix_power(S.astype(np.int64), 4); S6 = S4 @ S @ S
    return len(set(np.diag(S4))) == 1 and len(set(np.diag(S6))) == 1

cases = [('R_7', skew(rot_A(7)), 12, False), ('P_7', skew(paley_A(7)), 12, False)]
for c in regular_classes(7):
    S = skew(code_to_A(c))
    if not vertex_transitive_hint(S): cases.append(('N_7 (non-circulant)', S, 12, False))
cases.append(('R_9', skew(rot_A(9)), 10, True))
nvt = [c for c in regular_classes(9) if not vertex_transitive_hint(skew(code_to_A(c)))]
cases.append(('N_9a (gentourng class %d)' % regular_classes(9).index(nvt[0]), skew(code_to_A(nvt[0])), 10, True))
cases.append(('N_9b (gentourng class %d)' % regular_classes(9).index(nvt[-1]), skew(code_to_A(nvt[-1])), 10, True))
cases.append(('P_11', skew(paley_A(11)), 8, False))
rng = np.random.default_rng(5); n = 11
A = np.array([[1 if 0 < (j - i) % n <= 5 else 0 for j in range(n)] for i in range(n)])
for _ in range(4000):
    u, v, w = rng.choice(n, 3, replace=False)
    if A[u, v] and A[v, w] and A[w, u]: A[u, v] = A[v, w] = A[w, u] = 0; A[v, u] = A[w, v] = A[u, w] = 1
cases.append(('N_11 (3-cycle-reversal chain, seed 5)', skew(A), 8, False))

worst = 0.0
for name, S, Lmax, brute in cases:
    n = S.shape[0]; Sh = S / n
    assert all(S.sum(axis=1) == 0)
    t0 = time.time(); X = mscaling(S, Lmax, verbose=False); tm = time.time() - t0
    t0 = time.time(); co = LE4.series(Sh, 3, Lmax, K=32, r=0.5); te = time.time() - t0
    bf = None
    if brute:
        t0 = time.time(); bf = cluster_sums(S, 8, regular=True, emax=3, odd=False, Lmin=6, verbose=False); tb = time.time() - t0
    x1 = LE4.X_e(Sh, 1, 0.8); x2 = LE4.X_e(Sh, 2, 0.8)
    print(f"{name} (n = {n}, vertex-transitive hint {vertex_transitive_hint(S)}): engine {te:.0f}s, m-scaling {tm:.0f}s"
          + (f", partitions {tb:.0f}s" if brute else ''), flush=True)
    print(f"   X_1(0.8) engine {x1:+.13e} closed form {X1_of(Sh, 0.8):+.13e};  X_2(0.8) engine {x2:+.13e} old engine {LEF.XB(Sh, 2, 0.8):+.13e}")
    for L in range(2, Lmax + 1, 2):
        ex = float(X[L].get(3, 0))
        d = co[L] - ex; worst = max(worst, abs(d))
        line = f"   L={L:2d}: [z^L] X_3 engine {co[L]:+.13e}   blow-up scaling (exact) {ex:+.13e}   diff {d:+.1e}"
        if bf is not None and L in bf:
            b = float(bf[L].get(3, 0)); worst = max(worst, abs(co[L] - b))
            line += f"   set partitions (exact) {b:+.13e} (== m-scaling: {bf[L].get(3, 0) == X[L].get(3, 0)})"
        print(line, flush=True)
print(f"max |engine - exact| over all cases and orders: {worst:.1e}")
