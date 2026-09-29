"""Universality test of the 'growing remnant': the negative family to n=NMAX on m in [1,MMAX]; normalised window energy
E_n = 3^n sum_{m=1}^{MMAX} |v_n(m)|^2, the argmax path, and mu_hat_n(1) at the exits; for each predicted ridge arrival
(n_arr) the energy trajectory over the 60 levels before it: do all arriving packets grow?"""
import sys, math, numpy as np
sys.path.insert(0, '.')
from jseries import phases_level
NMAX = int(sys.argv[1]); MMAX = int(sys.argv[2]); A = 60
w = 2.0**(-np.arange(1, A+1))
prev_lo = -NMAX*A - MMAX - A; prev = np.ones(-prev_lo+1, dtype=np.complex128)
E = np.zeros(NMAX+1); PK = np.zeros(NMAX+1); PM = np.zeros(NMAX+1, dtype=int); MU1 = np.zeros(NMAX+1)
for n in range(1, NMAX+1):
    lo = -(NMAX-n)*A - MMAX
    ph = phases_level(n, lo-A, 0)
    P = ph * prev[(lo-A)-prev_lo : (0-prev_lo)+1]
    cur = np.zeros(-lo+1, dtype=np.complex128)
    for a in range(1, A+1):
        cur += w[a-1]*P[A-a : A-a+(-lo+1)]
    prev, prev_lo = cur, lo
    vv = 3.0**n * np.abs(cur[(0-lo)-MMAX:(0-lo)])**2      # m = MMAX..1
    E[n] = vv.sum(); PK[n] = vv.max(); PM[n] = MMAX - int(vv.argmax()); MU1[n] = 3.0**(n/2)*abs(cur[0-lo])
print(f"NMAX={NMAX} MMAX={MMAX}: normalised energy E_n of the window m in [1,{MMAX}] (mean {E[20:].mean():.1f}, median {np.median(E[20:]):.1f}), and 3^(n/2)|mu_hat_n(1)|")
print("  n: E_n  peak m  3^n|v|^2 peak   3^(n/2)|mu1|")
for n in range(20, NMAX+1, 5):
    print(f"  {n:4d}: {E[n]:8.1f}  {PM[n]:3d}  {PK[n]:7.1f}   {MU1[n]:6.2f}")
print("\nlocal maxima of E_n above 3x the median:", [(n, round(E[n])) for n in range(21, NMAX) if E[n] > 3*np.median(E[20:]) and E[n] >= E[n-1] and E[n] >= E[n+1]])
print("local maxima of 3^(n/2)|mu1| >= 2:", [(n, round(MU1[n],1)) for n in range(2, NMAX) if MU1[n] >= 2 and MU1[n] >= MU1[n-1] and MU1[n] >= MU1[n+1]])
np.save('arrivals_E.npy', np.vstack([E, PK, PM, MU1]))
