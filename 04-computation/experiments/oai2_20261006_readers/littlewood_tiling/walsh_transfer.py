"""Exact L4 norms of all Walsh/Morse products prod_m (1 + eps_m z^{2^m}) via the 2x2 transfer
v=(||y||_4^4, S(y)) -> A_eps v, A_eps = [[6, 8 eps],[1, 4 eps]], v_0 = (1, 0).  (PROVED identity.)
Check against direct FFT for small k; report max merit factor over all 2^k sign patterns (= over Sylvester rows)."""
import numpy as np
from tower_merit import autocorr_int
def direct(eps):
    P = np.array([1], dtype=np.int64)
    for e in eps:
        P = np.concatenate([P, e * P])
    C = autocorr_int(P[None, :].astype(np.int8))[0]
    return int(len(P) ** 2 + 2 * (C[1:] ** 2).sum()), int(sum(C[u] * C[len(P) - u] for u in range(1, len(P))))
# check identity
import itertools
for k in range(1, 9):
    for eps in itertools.product([1, -1], repeat=k):
        v = np.array([1, 0], dtype=object)
        for e in eps:
            v = np.array([6 * v[0] + 8 * e * v[1], v[0] + 4 * e * v[1]], dtype=object)
        L4, S = direct(eps)
        assert v[0] == L4 and v[1] == S, (eps, v, L4, S)
print("transfer identity verified for all sign patterns, k <= 8")
# all sign patterns up to k=22 (float64 fine for ratios)
R = np.array([1.0]); Sg = np.array([0.0])
for k in range(1, 23):
    Rp = np.concatenate([(6 * R + 8 * Sg) / 4, (6 * R - 8 * Sg) / 4])
    Sp = np.concatenate([(R + 4 * Sg) / 4, (R - 4 * Sg) / 4])
    R, Sg = Rp, Sp
    Fmax = 1.0 / (R.min() - 1.0)
    tm = None
    if k % 2 == 0 or k > 0:
        pass
    print(f"k={k:2d} N=2^{k}: max over Walsh products F = {Fmax:.6f}   ratio to (4/(1+sqrt17))^k: {Fmax / (4/(1+17**0.5))**k:.4f}")
