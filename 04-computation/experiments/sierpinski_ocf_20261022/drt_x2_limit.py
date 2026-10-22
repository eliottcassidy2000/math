# Limit of q^2 X_2 (three-loop cluster sum) for Paley P_q, q = 3 mod 4 prime; with the exact limits of the other
# second-order terms for doubly regular tournaments (X_1^2/2 -> 6.125, X_1 E_0 a/n -> -10.5, delta/n -> -49,
# E_0 a^2/2 -> 21, E_0 b -> 29/6, from rho/e and y_pred: 5/3) this gives c_2 = lim q^2 X_2 - 25.875.
import sys, time
import numpy as np
sys.path.insert(0, '.')
import loop_engine_fast as LEF
for q in [int(a) for a in sys.argv[1:]]:
    sq = {(t * t) % q for t in range(1, q)}
    A = np.array([[1.0 if (j - i) % q in sq else 0.0 for j in range(q)] for i in range(q)])
    Sh = (A - A.T) / q
    t0 = time.time(); x2 = float(np.real(LEF.XB(Sh, 2)))
    print(f"q={q}: q^2 X_2 = {q * q * x2:.10f}  ({time.time() - t0:.1f}s)", flush=True)
