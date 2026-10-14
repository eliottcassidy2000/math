"""
hofstadter_q_moebius_iteration_20261013.py -- the generation-to-generation iteration of the lag-1 correlation of the Q
jitter is the Moebius map F(r) = (r - 1)/(3 - r) (from the Yule-Walker equations of THM-4630 with c = r_prev/2:
r_new = (2c - 1)/(3 - 2c)), a degree-one rational function on the Riemann sphere with fixed points 1 -+ sqrt 2, the
attracting one 1 - sqrt 2 with multiplier F'(1 - sqrt 2) = 2/(2 + sqrt 2)^2 = 3 - 2 sqrt 2 = (sqrt 2 - 1)^2 = 0.1716, matrix
[[1, -1], [-1, 3]] of determinant 2 and eigenvalues 2 +- sqrt 2 (Codex, golden_nested_trees_20261013c.md section 6, derived
there as a two-step Pell mechanism).  Test: rho_k = lag-1 autocorrelation of d over generation k, for k = 4..20, against
the iteration rho_(k+1) = F(rho_k) started from the measured rho_4, and the geometric approach |rho_k - (1 - sqrt 2)|.
"""
import numpy as np
from numba import njit

@njit(cache=True)
def hofq(N):
    q = np.zeros(N + 1, dtype=np.int32); q[1] = 1; q[2] = 1
    for n in range(3, N + 1):
        i = n - q[n - 1]; j = n - q[n - 2]
        if i < 1 or j < 1:
            q[0] = n; break
        q[n] = q[i] + q[j]
    return q

q = hofq(1 << 22)
F = lambda r: (r - 1) / (3 - r)
rs = 1 - np.sqrt(2)
print("F fixed points:", np.roots([1, -2, -1]), " F'(1 - sqrt2) =", 2 / (2 + np.sqrt(2)) ** 2, " = 3 - 2 sqrt2 =", 3 - 2 * np.sqrt(2))
print(" k | measured rho_k   F(rho_(k-1))   |rho_k - (1-sqrt2)|   ratio of successive distances   (geometric rate predicted 0.1716)")
prev = None; dprev = None
for k in range(4, 21):
    lo, hi = 1 << k, 1 << (k + 1)
    d = q[lo:hi] / np.arange(lo, hi, dtype=np.float64) - 0.5
    d = d - d.mean()
    r = float((d[:-1] * d[1:]).mean() / (d * d).mean())
    dist = abs(r - rs)
    line = f" {k:2d} | {r:+.5f}        {('%+.5f' % F(prev)) if prev is not None else '   -   '}        {dist:.5f}"
    if dprev:
        line += f"               {dist / dprev:.3f}"
    print(line)
    prev = r; dprev = dist
print("done")
