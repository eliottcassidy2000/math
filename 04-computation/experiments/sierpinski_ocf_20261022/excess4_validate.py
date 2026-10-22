# Validation of the five-loop (excess-4) term computed by the sign-network engine loop_engine5.py (1076 cores,
# cores_e4.json): z-series coefficients of X_4 (Cauchy sums) against the exact rational blow-up scaling
# log Psi(Shat (x) J_m/m) = sum_e m^(-e) X_e(Shat) (audit2/b_mscaling.py), which uses no core theory.
# Also: the same engine's X_3 against the exact values (a check of the engine at a level validated before), and the
# engine's weight tensors at excess <= 3 equal decomp.c's (checked separately, see sign_network_basis5_check.out).
# Tournaments: R_7, P_7, the non-circulant regular 7-tournament (z^6..z^14); R_9 and a non-vertex-transitive regular
# 9-tournament (z^6..z^10).  Output: excess4_validate.out.  mac-mini-2026-10-22-loop4.
import sys, time
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'audit2')
import loop_engine5 as LE5
from b_mscaling import mscaling
from common import rot_A, paley_A, skew, regular_classes, code_to_A

def vt_hint(S):
    S4 = np.linalg.matrix_power(S.astype(np.int64), 4); S6 = S4 @ S @ S
    return len(set(np.diag(S4))) == 1 and len(set(np.diag(S6))) == 1

cases = [('R_7', skew(rot_A(7)), 14), ('P_7', skew(paley_A(7)), 14)]
for c in regular_classes(7):
    S = skew(code_to_A(c))
    if not vt_hint(S): cases.append(('N_7 (non-circulant)', S, 14))
cases.append(('R_9', skew(rot_A(9)), 10))
nvt = [c for c in regular_classes(9) if not vt_hint(skew(code_to_A(c)))]
cases.append((f'N_9 (gentourng class {regular_classes(9).index(nvt[0])})', skew(code_to_A(nvt[0])), 10))

worst = {3: 0.0, 4: 0.0}
for name, S, Lmax in cases:
    n = S.shape[0]; Sh = S / n
    t0 = time.time(); X = mscaling(S, Lmax, verbose=False); tm = time.time() - t0
    line = f"{name} (n = {n}): m-scaling {tm:.0f}s"
    co = {}
    for e in (3, 4):
        t0 = time.time(); co[e] = LE5.series_general(Sh, e, Lmax, K=32, r=0.5); line += f", engine excess {e} {time.time() - t0:.0f}s"
    print(line, flush=True)
    for L in range(2, Lmax + 1, 2):
        parts = []
        for e in (3, 4):
            ex = float(X[L].get(e, 0)); d = co[e][L] - ex
            worst[e] = max(worst[e], abs(d))
            parts.append(f"X_{e}: engine {co[e][L]:+.13e} exact {ex:+.13e} diff {d:+.1e}")
        print(f"   L={L:2d}: " + ' | '.join(parts), flush=True)
print(f"max |engine - exact|: excess 3 {worst[3]:.1e}, excess 4 {worst[4]:.1e}")
