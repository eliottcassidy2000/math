"""Truncation control for mu_hat_n(1): the a priori error h 2^-A exceeds 3^-(n/2) beyond n ~ 1.26 A, so the values are
certified beyond that only by agreement between different truncations. Compare A = 60, 100 (N = 600) with the A = 120 run."""
import sys, math, numpy as np
sys.path.insert(0, '.')
from jseries import phases_level
def mu1_run(N, A):
    w = 2.0**(-np.arange(1, A+1)); prev_lo = -N*A - A; prev = np.ones(-prev_lo+1, dtype=np.complex128)
    mu1 = np.zeros(N+1, dtype=np.complex128); mu1[0] = 1
    for n in range(1, N+1):
        lo = -(N-n)*A; ph = phases_level(n, lo-A, 0); P = ph*prev[(lo-A)-prev_lo:(0-prev_lo)+1]
        cur = np.zeros(-lo+1, dtype=np.complex128)
        for a in range(1, A+1): cur += w[a-1]*P[A-a:A-a+(-lo+1)]
        prev, prev_lo = cur, lo; mu1[n] = cur[-lo]
    return mu1
m60 = mu1_run(600, 60); m100 = mu1_run(600, 100); m120 = np.load('mu1_N1200.npy')[:601]
for n in (100, 131, 200, 300, 400, 500, 600):
    print(f"n={n}: |mu1| A=60 {abs(m60[n]):.6e}  A=100 {abs(m100[n]):.6e}  A=120 {abs(m120[n]):.6e};  rel diff 60 vs 120 {abs(m60[n]-m120[n])/abs(m120[n]):.1e}, 100 vs 120 {abs(m100[n]-m120[n])/abs(m120[n]):.1e};  a priori bound n 2^-60 / |mu1| = {n*2.0**-60/abs(m120[n]):.1e}")
print("max relative difference A=60 vs A=120 over n<=600:", max(abs(m60[n]-m120[n])/abs(m120[n]) for n in range(1,601)))
print("max relative difference A=100 vs A=120 over n<=600:", max(abs(m100[n]-m120[n])/abs(m120[n]) for n in range(1,601)))
