# AUDIT Section 1 power counting: q^3 X_3(P_q) values of paley_x3_scaling.out, (i) the polynomial extrapolations in
# 1/q are re-fitted; (ii) test of the stated reason for the 0.03 fit residual ("TT5(P_q) carries arithmetic
# fluctuations"): add the regressor r(q) = 7168 TT5(P_q) / q^5 (the z^8 TT5-dependence of q^3 X_3, from
# N_8 = U_8 + 7168 TT5) and see whether the residual collapses.  TT5(P_q) computed here (own recursion).
import re, sys
import numpy as np
sys.path.insert(0, '.')
from aud_task1 import tt
vals = {}
for line in open('../paley_x3_scaling.out'):
    m = re.match(r'P_(\d+): .* q\^3 X_3 = ([+-][\d.]+)', line)
    if m: vals[int(m.group(1))] = float(m.group(2))
qs = sorted(vals); y = np.array([vals[q] for q in qs]); x = 1 / np.array(qs, float)
def TT5_paley(q):
    sq = {(t * t) % q for t in range(1, q)}
    out = [sum(1 << j for j in range(q) if (j - i) % q in sq) for i in range(q)]
    return tt(5, (1 << q) - 1, out)
T5 = np.array([TT5_paley(q) for q in qs], float)
print("q:      ", qs); print("q^3 X_3:", list(y)); print("TT5(P_q):", [int(t) for t in T5])
r = 7168 * T5 / np.array(qs, float) ** 5
for deg in (3, 4, 5):
    A = np.column_stack([x ** k for k in range(deg + 1)])
    c, *_ = np.linalg.lstsq(A, y, rcond=None); res0 = np.max(np.abs(y - A @ c))
    A2 = np.column_stack([A, r]); c2, *_ = np.linalg.lstsq(A2, y, rcond=None); res2 = np.max(np.abs(y - A2 @ c2))
    A3 = np.column_stack([A, r, r * x]); c3, *_ = np.linalg.lstsq(A3, y, rcond=None); res3 = np.max(np.abs(y - A3 @ c3))
    print(f"degree {deg}: limit {c[0]:+.3f}, max resid {res0:.1e} | + r(q): limit {c2[0]:+.3f}, coeff of r {c2[-1]:+.3f}, max resid {res2:.1e} "
          f"| + r, r/q: limit {c3[0]:+.3f}, max resid {res3:.1e}")
# fluctuation of TT5 itself about a smooth polynomial in q (degree 5)
A = np.column_stack([np.array(qs, float) ** k for k in range(6)])
c, *_ = np.linalg.lstsq(A, T5, rcond=None)
print("TT5(P_q) minus its degree-5 polynomial fit in q:", np.round(T5 - A @ c, 1))
# implied limit of q^3 X_3 = (limit of the TT5-corrected fit) + beta * lim r(q), lim r = 7168 lim TT5/q^5;
# quasi-randomness of Paley tournaments gives lim TT5/q^5 = 5!/(5! 2^10) = 1/1024, i.e. lim r = 7.
q = np.array(qs, float)
for deg in (2, 3, 4):
    A = np.column_stack([x ** k for k in range(deg + 1)])
    c, *_ = np.linalg.lstsq(A, T5 / q ** 5, rcond=None)
    print(f"TT5/q^5 polynomial fit in 1/q, degree {deg}: limit {c[0]:.4e} (1/1024 = {1 / 1024:.4e}), max resid {np.max(np.abs(T5 / q ** 5 - A @ c)):.1e}")
for deg in (3, 4, 5):
    A = np.column_stack([x ** k for k in range(deg + 1)])
    A2 = np.column_stack([A, r]); c2, *_ = np.linalg.lstsq(A2, y, rcond=None)
    print(f"degree {deg} + r(q): implied lim q^3 X_3 = {c2[0]:+.3f} + {c2[-1]:.3f} x 7 = {c2[0] + 7 * c2[-1]:+.3f}")
