"""Audit: non-Heegner CM points with conductor divisible by N (D(tau) = -3, D(N tau) = -3 N^2 with h <= 2):
predicted by hand m_5 = -108 +- 48 sqrt5 (D(5tau) = -75) and m_7 = -112 +- 24 sqrt21 (D(7tau) = -147); irrational."""
import sys, os
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a3_cm_classification import m_value, best_rep
mp.mp.dps = 50
rho = mp.mpc(-0.5, mp.sqrt(3) / 2)
for N, pred in ((5, (-108, 48, 5)), (7, (-112, 24, 21))):
    vals = [m_value(best_rep(rho, N), N)] + [m_value(best_rep(-1 / (rho + k), N), N) for k in range(N)]
    a, b, d = pred
    targets = [a + b * mp.sqrt(d), a - b * mp.sqrt(d)]
    print(f'N={N}: m_N at tau=rho and its {N} coset images:', [mp.nstr(v.real, 12) for v in vals])
    print(f'      predicted {a} +- {b} sqrt{d} =', [mp.nstr(t, 12) for t in targets],
          ' both present:', all(min(abs(v - t) for v in vals) < 1e-30 for t in targets))
