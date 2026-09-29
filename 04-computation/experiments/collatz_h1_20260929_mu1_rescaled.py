"""Rescaled negative-window recursion: v_n(k) = 3^(n/2) f_n(k) (no underflow), v_n = sqrt3 * sum_a 2^-a omega_n(k-a) v_{n-1}(k-a).
Reports 3^(n/2)|mu_hat_n(1)| = |v_n(0)| to N, the sup of |v_n(0)| (3^(-1/2)/rho)^n for rho = 0.585, 0.58, and the local maxima."""
import sys, math, numpy as np
sys.path.insert(0, '.')
from jseries import phases_level
N = int(sys.argv[1]); A = int(sys.argv[2])
w = math.sqrt(3.0) * 2.0**(-np.arange(1, A+1)); prev_lo = -N*A - A; prev = np.ones(-prev_lo+1, dtype=np.complex128)
v0 = np.zeros(N+1); v0[0] = 1
for n in range(1, N+1):
    lo = -(N-n)*A; ph = phases_level(n, lo-A, 0); P = ph*prev[(lo-A)-prev_lo:(0-prev_lo)+1]
    cur = np.zeros(-lo+1, dtype=np.complex128)
    for a in range(1, A+1): cur += w[a-1]*P[A-a:A-a+(-lo+1)]
    prev, prev_lo = cur, lo; v0[n] = abs(cur[-lo])
np.save(f'v0_N{N}.npy', v0)
print(f"N={N} A={A}: 3^(n/2)|mu_hat_n(1)|: median {np.median(v0[1:]):.4f}, max {v0[1:].max():.3f} at n={v0[1:].argmax()}")
for rho in (0.585, 0.58, 0.5774):
    r = v0[1:] * (3**-0.5/rho)**np.arange(1, N+1)
    print(f"  sup_n |mu_hat_n(1)|/{rho}^n = {r.max():.3f} at n={r.argmax()+1}; over n>=200: {r[199:].max():.3f} at n={r[199:].argmax()+200}; over n>=1200: {r[1199:].max():.3e} at n={r[1199:].argmax()+1200}")
print("  local maxima of 3^(n/2)|mu_hat_n(1)| >= 1.0:", [(n, round(v0[n],2)) for n in range(2, N) if v0[n] >= 1.0 and v0[n] >= v0[n-1] and v0[n] >= v0[n+1]])
print("  window medians of 3^(n/2)|mu_hat_n(1)| by 250-blocks:", [f"{np.median(v0[a:a+250]):.3f}" for a in range(1, N, 250)])
for (a_, b_) in ((200, N), (N//2, N), (1200, N)):
    if b_ > a_ + 50:
        ns = np.arange(a_, b_+1); sl = np.polyfit(ns, np.log(v0[a_:b_+1] * 3.0**(-ns/2) if False else np.log(v0[a_:b_+1]) - ns*math.log(3)/2), 1)[0]
        print(f"  least-squares rate of |mu_hat_n(1)| over n={a_}..{b_}: {math.exp(sl):.4f}")
