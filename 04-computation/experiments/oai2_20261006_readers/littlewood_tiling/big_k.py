import numpy as np, time, sys
from tower_merit import tower, autocorr_int, supnorm
from switched_tower import rs_seq
def stats(X, over):
    N = X.shape[1]; bs = max(1, (1 << 22) // N)
    s2 = []
    for s in range(0, N, bs):
        C = autocorr_int(X[s:s+bs]); s2.append((C[:, 1:] ** 2).sum(axis=1))
    s2 = np.concatenate(s2); F = N * N / (2.0 * s2)
    mx = []
    b2 = max(1, bs // over)
    for s in range(0, N, b2):
        ub, m, mn = supnorm(X[s:s+b2], over); mx.append(ub)
    mx = np.concatenate(mx) / np.sqrt(N)
    return F, mx, int(s2.sum())
for k in [int(a) for a in sys.argv[1:]]:
    t = time.time()
    H = tower(k); N = H.shape[0]
    F, mx, tot = stats(H, 8)
    print(f"tower    k={k} N={N}: rows Fmax {F.max():.4f} Fmean {F.mean():.4f} Fmin {F.min():.6f}; best-row sup/sqrtN (ub) {mx.min():.4f}; sum_rows sumC^2 = {tot}  [{time.time()-t:.0f}s]", flush=True)
    a = rs_seq(k)
    H *= a[:, None]; H *= a[None, :]
    F, mx, tot = stats(H, 8)
    print(f"switched k={k} N={N}: rows Fmin {F.min():.4f} Fmean {F.mean():.4f} Fmax {F.max():.4f}; worst-row sup/sqrtN (ub) {mx.max():.4f}  [{time.time()-t:.0f}s]", flush=True)
    del H
