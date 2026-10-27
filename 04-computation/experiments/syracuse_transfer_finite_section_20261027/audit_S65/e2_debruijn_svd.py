"""Side check: singular values of the T-map Ulam (de Bruijn) matrix on Z/2^K (exact structure + float SVD):
P P^T is block-diagonal with 2x2 averaging blocks on the pairs {a, a'} with T(a) = T(a') mod 2^(K-1), so the singular
values are 1 (multiplicity 2^(K-1)) and 0 (multiplicity 2^(K-1)): sigma_2 = 1 exactly, singular gap 0 for every K."""
import numpy as np
for K in range(2, 12):
    n = 1 << K
    P = np.zeros((n, n))
    for a in range(n):
        t = a // 2 if a % 2 == 0 else (3 * a + 1) // 2
        b = t % (n // 2)
        P[a, b] += 0.5
        P[a, b + n // 2] += 0.5
    s = np.sort(np.linalg.svd(P, compute_uv=False))[::-1]
    n1 = int(np.sum(np.abs(s - 1) < 1e-9)); n0 = int(np.sum(s < 1e-9))
    print(f'K = {K:2d}: singular values: {n1} equal to 1, {n0} equal to 0 (of {n}); sigma_2 = {s[1]:.12f}')
