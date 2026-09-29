"""mu_hat_n(1) = E e(Y_n/3^n) for n <= N (the H1 sequence), by the negative-window recursion; saves mu1.npy;
reports sup_n rho^-n |mu_hat_n(1)| for rho in (0.585, 0.58, 0.5774), the local maxima of 3^(n/2)|mu_hat_n(1)| above 2,
and the predicted ridge arrivals n ~ k + (Q - r_k)/1.3 for the deepest seeds (from ridge_inventory)."""
import sys, math, numpy as np
sys.path.insert(0, '.')
from jseries import phases_level
N = int(sys.argv[1]); A = int(sys.argv[2])
LOG23 = math.log2(3.0)
w = 2.0**(-np.arange(1, A+1))
prev_lo = -N*A - A; prev = np.ones(-prev_lo + 1, dtype=np.complex128)
mu1 = np.zeros(N+1, dtype=np.complex128); mu1[0] = 1
for n in range(1, N+1):
    lo = -(N-n)*A
    ph = phases_level(n, lo - A, 0)
    P = ph * prev[(lo-A)-prev_lo : (0-prev_lo)+1]
    cur = np.zeros(-lo+1, dtype=np.complex128)
    for a in range(1, A+1):
        cur += w[a-1]*P[A-a : A-a+(-lo+1)]
    prev, prev_lo = cur, lo
    mu1[n] = cur[-lo]
np.save(f'mu1_N{N}.npy', mu1)
v = np.abs(mu1[1:]) * 3.0**(np.arange(1, N+1)/2)
print(f"N={N} A={A}: 3^(n/2)|mu_hat_n(1)|: median {np.median(v):.3f}, mean {v.mean():.3f}, max {v.max():.3f} at n={v.argmax()+1}")
for rho in (0.585, 0.58, 0.5774, 0.575):
    r = np.abs(mu1[1:]) / rho**np.arange(1, N+1)
    print(f"  sup_n |mu_hat_n(1)| / {rho}^n = {r.max():.3f} at n={r.argmax()+1};  over n>=200: {r[199:].max():.3f} at n={r[199:].argmax()+200}")
print("  local maxima of 3^(n/2)|mu_hat_n(1)| >= 2.5 (n, value):")
loc = [(n, v[n-1]) for n in range(2, N) if v[n-1] >= 2.5 and v[n-1] >= v[n-2] and v[n-1] >= v[n]]
print("   ", " ".join(f"({n},{x:.1f})" for n, x in loc))
print("  predicted arrivals of the deepest seeds (u,sign,Q,depth k): n_arr ~ k + (Q - (k log2 3 - 6 - log2 u))/1.3:")
seeds = [(55,'+',423,15),(41,'-',271,12),(7,'+',335,8),(13,'-',154,7),(1,'-',162,5),(1,'-',486,6),(1,'+',243,6),(1,'+',729,7),(1,'-',1458,7),(1,'+',2187,8),(11,'-',689,8),(35,'-',69,8),(5,'-',463,7),(1,'+',567,5),(1,'-',324,5),(1,'+',405,5)]
for (u,sg,Q,k) in sorted(seeds, key=lambda t: t[2]):
    r_k = k*LOG23 - 6 - math.log2(u)
    n_arr = k + (Q - r_k)/1.3
    if n_arr <= N + 40:
        lo_, hi_ = max(1, int(n_arr) - 25), min(N, int(n_arr) + 25)
        seg = v[lo_-1:hi_]
        print(f"   u={u:3d}{sg} Q={Q:4d} k={k:2d}: n_arr ~ {n_arr:6.0f};  max of 3^(n/2)|mu_hat_n(1)| in [{lo_},{hi_}] = {seg.max():.2f} at n={lo_+seg.argmax()}")
# windows of the least-squares rate
for (a_, b_) in ((20, N), (N//2, N), (3*N//4, N)):
    ns = np.arange(a_, b_+1); sl = np.polyfit(ns, np.log(np.abs(mu1[a_:b_+1])), 1)[0]
    print(f"  least-squares rate over n={a_}..{b_}: {math.exp(sl):.4f}")
