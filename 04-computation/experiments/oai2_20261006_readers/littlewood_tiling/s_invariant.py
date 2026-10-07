"""Check the invariant-interval argument: s(x) = S(x)/||x||_4^4 with S(x) = sum_{u=1}^{m-1} C_u C_{m-u}.
Pure Walsh products: s in [0, 1/4] (PROVED by f([-1/4,1/4]) = [0,1/4], f(t) = (1+4t)/(6+8t)).
Tower: rows r_i and flipped rows r_i - 2 e_i at each level: report min/max s, and min over rows of R_{k}/R_{k-2}
along the tower recursion (R = ||.||_4^4 / N^2)."""
import numpy as np
from tower_merit import tower, autocorr_int
# Walsh, all patterns to k=20
R = np.array([1.0]); S = np.array([0.0])  # normalized: R, s
for k in range(1, 21):
    t_plus, t_minus = S, -S
    R = np.concatenate([R * (1.5 + 2 * t_plus), R * (1.5 + 2 * t_minus)])
    S = np.concatenate([(1 + 4 * t_plus) / (6 + 8 * t_plus), (1 + 4 * t_minus) / (6 + 8 * t_minus)])
    if k in (5, 10, 15, 20):
        print(f"Walsh k={k}: s in [{S.min():.6f}, {S.max():.6f}]   min R = {R.min():.4f}  (1.5^(k/2) = {1.5**(k/2):.4f})")
def stats(X):
    C = autocorr_int(X); m = X.shape[1]
    L4 = m * m + 2 * (C[:, 1:] ** 2).sum(axis=1)
    Sx = np.array([int(np.dot(C[i, 1:m], C[i, m-1:0:-1])) for i in range(X.shape[0])]) if m > 1 else np.zeros(X.shape[0], dtype=np.int64)
    return L4 / m**2, Sx / L4
prevR = {}
for k in range(1, 13):
    H = tower(k).astype(np.int8); N = H.shape[0]
    Hf = H.copy(); Hf[np.arange(N), np.arange(N)] = -1   # r_i - 2 e_i  (= -c_i)
    R1, s1 = stats(H); R2, s2 = stats(Hf)
    line = f"tower k={k:2d}: rows s in [{s1.min():+.4f},{s1.max():+.4f}], flipped s in [{s2.min():+.4f},{s2.max():+.4f}], min R rows {R1.min():.4f}"
    if k - 2 in prevR:
        # row i at level k comes from row (i mod 2^(k-2)) [or its flip] at level k-2: compare R(row) to min of the two ancestors
        Rp, Rpf = prevR[k - 2]
        anc = np.arange(N) % (N // 4)
        ratio = R1 / np.maximum(Rp[anc], Rpf[anc])
        line += f", min_i R_k(i)/max(R_(k-2) ancestors) = {ratio.min():.4f}"
    prevR[k] = (R1, R2)
    print(line)
