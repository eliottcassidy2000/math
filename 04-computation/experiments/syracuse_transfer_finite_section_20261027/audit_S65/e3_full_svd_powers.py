"""Auditor's closed form for the full singular spectrum of P_K^n (from the exact Gram identity
(P^n)^T P^n = sum_m w_m(n) E_m, w_m(n) = binom(K-m-1, n-1) 2^(m-K) for m >= 2):
  1 (once) and sqrt(P(Bin(K-m, 1/2) >= n)) with multiplicity 2^(m-2), m = 2..K.   Float SVD check, all n, K <= 11."""
import math
import numpy as np
from e1_float import P_K
worst = 0.0
for K in range(3, 12):
    P = P_K(K)
    Pn = np.eye(len(P))
    for n in range(1, K):
        Pn = Pn @ P
        s = np.sort(np.linalg.svd(Pn, compute_uv=False))[::-1]
        pred = [1.0]
        for m in range(2, K + 1):
            t = sum(math.comb(K - m, j) for j in range(n, K - m + 1)) / 2 ** (K - m)
            pred += [math.sqrt(t)] * (2 ** (m - 2))
        pred = np.sort(np.array(pred))[::-1]
        worst = max(worst, float(np.max(np.abs(s - pred))))
print(f'full singular spectrum of P_K^n matches the closed form for K = 3..11, n = 1..K-1: max deviation {worst:.1e}')
