"""Numerical check that the core sum (5) of THM-4845 converges away from 4^k (meromorphic continuation)."""
from audit_gf import partition_core_kernel
from fractions import Fraction
from math import comb
import mpmath as mp
mp.mp.dps = 40
Kp = partition_core_kernel(90)
def F(y, mmax):
    tot = mp.mpc(0)
    for (m, t), K in Kp.items():
        if m > mmax: continue
        k = (m - t) // 2
        kappa = mp.mpf(K) / mp.mpf(2) ** comb(m, 2)
        tot += kappa * y ** m / (1 - y / mp.mpf(4) ** k) ** (m + 1)
    return tot
for y in (mp.mpc(20, 0), mp.mpc(-50, 7), mp.mpc(1000, 1), mp.mpc(3, 2)):
    print('y =', y, ' m<=40:', mp.nstr(F(y, 40), 20), ' m<=60:', mp.nstr(F(y, 60), 20), ' m<=90:', mp.nstr(F(y, 90), 20))
