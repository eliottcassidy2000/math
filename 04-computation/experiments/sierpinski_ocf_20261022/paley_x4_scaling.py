# Five-loop power counting on Paley tournaments: q^4 X_4(P_q) at z = 1 (sign-network engine loop_engine5.py, DRT
# two-label mode, all 1076 excess-4 cores), next to q X_1, q^2 X_2, q^3 X_3 from paley_x3_scaling.out.
# Output: paley_x4_scaling.out.  mac-mini-2026-10-22-loop4.
import sys, time
import numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, 'audit2')
import loop_engine4 as LE4
import loop_engine5 as LE5
from common import paley_A, skew
res = []
for q in (19, 23, 31, 43):          # q = 47 exceeds the contraction planner's 2^27-entry memory cap
    t = time.time()
    S = skew(paley_A(q)); assert LE4.check_drt(S)
    ops = LE5.ops_drt(q, np.array([1.0]))
    B = np.array([np.eye(q) - np.ones((q, q)) / q, S / q])
    x4 = 0.0
    for V, edges, aut in LE5.load_cores(4):
        W2 = LE5.SignCore(V, edges, aut).contract('drt', ops).reshape(-1)
        x4 += (W2 @ LE4.Network(V, edges).contract(B).reshape(-1)).real / aut
    res.append((q, q ** 4 * x4))
    print(f"P_{q}: q^4 X_4 = {q ** 4 * x4:+.4f}  ({time.time() - t:.0f}s)", flush=True)
qs = np.array([r[0] for r in res], float); ys = np.array([r[1] for r in res])
for deg in (1, 2):
    c = np.polyfit(1 / qs[-(deg + 1):] if deg == 1 else 1 / qs, ys[-(deg + 1):] if deg == 1 else ys, deg)
    print(f"fit of degree {deg} in 1/q ({'last 2 points' if deg == 1 else 'all points'}): limit {c[-1]:+.1f}")
print("loop coefficients of P_q at z = 1 (limits): q X_1 -> -3.5, q^2 X_2 -> +29, q^3 X_3 -> about -341, q^4 X_4 -> see fits;"
      " alternating in sign, ratios about 8.3, 11.8, 14-15 (asymptotic series)")
