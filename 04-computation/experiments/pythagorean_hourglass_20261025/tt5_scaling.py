# Size of the arithmetic part of TT5(P_p): high-precision least-squares degree-5 polynomial over all primes, then
# the scaling exponent of the residual (the residual contains the oscillating arithmetic term up to its projection).
import re, math
import mpmath as mp
mp.mp.dps = 60
rows = [(int(m.group(1)), int(m.group(2))) for m in (re.match(r'p=(\d+) TT5=(\d+)', l) for l in open('paley_tt5_data.out')) if m]
P = [mp.mpf(p) for p, t in rows]; T = [mp.mpf(t) for p, t in rows]
# scaled variable u = p / 1000 to condition the normal equations
U = [p / 1000 for p in P]
A = mp.matrix([[u ** k for k in range(6)] for u in U]); b = mp.matrix(T)
x = mp.lu_solve(A.T * A, A.T * b)
res = [T[i] - sum(x[k] * U[i] ** k for k in range(6)) for i in range(len(rows))]
import numpy as np
lp = np.array([math.log(float(p)) for p in P]); lr = np.array([math.log(abs(float(r)) + 1e-30) for r in res])
sel = lp > math.log(300)
slope = np.polyfit(lp[sel], lr[sel], 1)[0]
print(f"polynomial main term (in u = p/1000): {[mp.nstr(c, 10) for c in x]}")
print(f"leading coefficient per p^5: {mp.nstr(x[5] / 1000 ** 5, 10)} (1/1024 = {1 / 1024:.10f})")
print(f"residual |TT5 - poly| grows like p^{slope:.3f} (fit on p > 300); typical sizes:",
      [(int(P[i]), f"{float(res[i]):.3g}") for i in range(0, len(rows), 15)])
print("residual / p^3.5 at large p:", [f"{float(res[i]) / float(P[i]) ** 3.5:+.4f}" for i in range(len(rows) - 12, len(rows))])
