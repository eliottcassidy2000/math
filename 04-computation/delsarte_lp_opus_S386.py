# opus-2026-07-17-S386 -- HYP-7720: THE DELSARTE LP FOR LRC(14).
#
# THE SETUP.  Let P >= 0 be a trig polynomial, P = |sum_k a_k e(kt)|^2.  Then
#     int P * 1_{U}  =  int P  -  int P * 1_{union D_v}
#                    >= int P  -  sum_v int P * 1_{D_v}.
# Since 1_{D_v}(t) = h(v t) with h the arc indicator, its Fourier coefficients
# live on multiples of v, giving
#     int P * 1_{D_v} = sum_n Phat(n v) hhat(n).
# With Phat(m) = sum_k a_k a_{k+m} (real autocorrelation) this is a QUADRATIC
# FORM a^T M a, where M is the symmetric TOEPLITZ matrix
#     M[k,l] = c_{l-k},   c_d = sum over {v in V : v | d} of hhat(d/v),
# and c_0 = sum_v hhat(0) = 13 * 2*lam = 13/7.
# So: if lambda_min(M) < 1 = Phat(0), the uncovered set has positive measure and
# LRC(14) holds for V.  A free optimisation over P confronting the 13/7 term --
# exactly the shape of a Delsarte bound, and never tried in this project.
from math import sin, pi, gcd
import numpy as np
from fractions import Fraction as F
LAM = 1.0/14
def hhat(n):
    if n == 0: return 2*LAM
    return sin(2*pi*n*LAM)/(pi*n)

def toeplitz_M(V, K):
    c = np.zeros(K+1)
    for d in range(K+1):
        s = 0.0
        for v in V:
            if d % v == 0: s += hhat(d // v)
        c[d] = s
    M = np.empty((K+1, K+1))
    for k in range(K+1):
        for l in range(K+1):
            M[k, l] = c[abs(l-k)]
    return M, c

def covering_multiplicity_min(V, grid=200000):
    """min over t of #{v : ||v t|| < lam} -- the true covering multiplicity."""
    ts = np.arange(grid) / grid
    m = np.zeros(grid, dtype=int)
    for v in V:
        r = (v*ts) % 1.0
        m += (np.minimum(r, 1-r) < LAM).astype(int)
    return m.min()

print("(1) lambda_min(M) vs K, and the true covering multiplicity")
print("    criterion: lambda_min(M) < 1  =>  uncovered set nonempty")
fams = [("{1,...,13} (tight)", list(range(1,14))),
        ("{1..11,13,24}",      [1,2,3,4,5,6,7,8,9,10,11,13,24]),
        ("random spread",      [3,7,11,19,23,31,37,41,47,53,59,61,67]),
        ("AP d=8",             [1+8*i for i in range(13)])]
for nm, V in fams:
    row = []
    for K in (30, 80, 160):
        M, c = toeplitz_M(V, K)
        row.append(np.linalg.eigvalsh(M).min())
    mult = covering_multiplicity_min(V)
    print(f"    {nm:22s} lam_min: K=30 {row[0]:+.4f}  K=80 {row[1]:+.4f}  K=160 {row[2]:+.4f}"
          f"   | true min multiplicity = {mult}")

print()
print("(2) WHY -- the Toeplitz SYMBOL.  c_d = sum_{v|d} hhat(d/v), so")
print("      f(theta) = sum_d c_d e(d theta) = sum_v sum_n hhat(n) e(n v theta)")
print("               = sum_v h(v theta) = THE COVERING MULTIPLICITY FUNCTION.")
print("    As K -> infinity, lambda_min(Toeplitz) -> min_theta f(theta), i.e. the")
print("    minimum number of arcs covering a point.  So the LP criterion")
print("    'lambda_min < 1' is EXACTLY 'some point is covered 0 times'.")
print("    Check: does lambda_min approach the integer multiplicity?")
for nm, V in fams:
    M, c = toeplitz_M(V, 160)
    lm = np.linalg.eigvalsh(M).min()
    print(f"    {nm:22s} lambda_min(K=160) = {lm:+.5f}   true multiplicity = {covering_multiplicity_min(V)}")
