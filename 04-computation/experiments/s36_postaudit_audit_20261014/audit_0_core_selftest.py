"""Self-test of the referee core: brute force over permutations vs explicit DFS (H, C, Nd, L identical), Landau score
test vs Warshall closure for strongness, triangle test vs distinct-score test for transitivity."""
import time
import numpy as np
from pa_core import (random_tournament, brute_force_perms, enumerate_paths_dfs, strong_landau, strong_warshall,
                     transitive_by_triangles, count_paths_dfs)

rng = np.random.default_rng(20261014)
t0 = time.time()
bad = 0; tested = 0
for n in range(3, 10):
    for trial in range(40 if n <= 8 else 10):
        A = random_tournament(n, rng)
        H1, C1, N1, L1 = brute_force_perms(n, A)
        C2 = np.zeros((n, n), dtype=np.int64); N2 = np.zeros((n, n + 1), dtype=np.int64); L2 = np.zeros((n, n), dtype=np.int64)
        H2 = enumerate_paths_dfs(n, A, C2, N2, L2)
        H3 = count_paths_dfs(n, A)
        tested += 1
        if not (H1 == H2 == H3 and (C1 == C2).all() and (N1 == N2).all() and (L1 == L2).all() and H1 % 2 == 1):
            bad += 1
print(f'brute force vs DFS: {tested} random tournaments n = 3..9, mismatches (or even H): {bad}  [{time.time() - t0:.1f} s]')

bad_s = 0; bad_t = 0; ns = 0
for trial in range(20000):
    n = int(rng.integers(3, 13))
    A = random_tournament(n, rng)
    if strong_landau(n, A) != strong_warshall(n, A):
        bad_s += 1
    scores = sorted(int(A[i].sum()) for i in range(n))
    if transitive_by_triangles(n, A) != (scores == list(range(n))):
        bad_t += 1
    ns += 1
# transitive tournaments and near-transitive ones (where strongness is delicate)
for n in range(3, 13):
    for lo in range(n):
        for hi in range(lo + 1, n):
            A = np.zeros((n, n), dtype=np.int8)
            for a in range(n):
                for b in range(a):
                    A[a, b] = 1
            A[hi, lo] = 0; A[lo, hi] = 1
            ns += 1
            if strong_landau(n, A) != strong_warshall(n, A):
                bad_s += 1
            if strong_landau(n, A) != (lo == 0 and hi == n - 1):
                bad_s += 1
print(f'strongness Landau vs Warshall (and transitive+1 strong iff extreme arc): {ns} tests, mismatches {bad_s}; '
      f'transitivity triangle test vs distinct scores mismatches {bad_t}')
