"""Audit A: best Walsh (Sylvester) row merit factor via exact transfer recursion over all sign patterns.
Walsh rows of length 2^k = prod_j (1 + eps_j z^{2^{j-1}}); (n4, S) -> (6 n4 + 8 eps S, n4 + 4 eps S), start (1, 0)."""
import numpy as np, math
lam = 4 / (1 + math.sqrt(17))
n4 = np.array([1], dtype=object); S = np.array([0], dtype=object)
for k in range(1, 23):
    n4n = np.concatenate([6 * n4 + 8 * S, 6 * n4 - 8 * S])
    Sn = np.concatenate([n4 + 4 * S, n4 - 4 * S])
    n4, S = n4n, Sn
    # dedupe states to keep it small
    st = np.unique(np.stack([n4.astype(float), S.astype(float)]), axis=1)
    pairs = sorted(set(zip(n4.tolist(), S.tolist())))
    n4 = np.array([p[0] for p in pairs], dtype=object); S = np.array([p[1] for p in pairs], dtype=object)
    N = 2 ** k
    best_n4 = min(n4)
    F = N * N / (best_n4 - N * N) if best_n4 > N * N else math.inf
    if k >= 4:
        print(f"k={k:2d}: states {len(n4):6d}, best F = {F:.6f}, F/lam^k = {F / lam ** k:.5f}", flush=True)
