# Audit2 (C): two random (non-vertex-transitive) regular 15-tournaments (3-cycle-reversal chain from R_15): brute-force
# excess-0/1 cluster sums (set partitions) against the closed formulas, exact rationals, through z^8.
import sys, random
import numpy as np
sys.path.insert(0, '.')
from b_partitions import cluster_sums
from common import rot_A, skew, X0_series, X1_series, int_powers
random.seed(15)
n = 15
A = rot_A(n).copy()
for k in range(2):
    for step in range(3000):
        a, b, c = random.sample(range(n), 3)
        if A[a, b] and A[b, c] and A[c, a]:
            A[a, b] = A[b, c] = A[c, a] = 0; A[b, a] = A[c, b] = A[a, c] = 1
    S = skew(A)
    d6 = sorted(set(int(int_powers(S, 6)[6][v, v]) for v in range(n)))
    cs = cluster_sums(S, 8, regular=True, emax=1, odd=False, Lmin=2, verbose=False)
    x0 = X0_series(S, 8); x1 = X1_series(S, 8)
    res = [(L, cs[L].get(0, 0) == x0[L], cs[L].get(1, 0) == x1[L]) for L in (2, 4, 6, 8)]
    print(f"random regular 15-tournament #{k} (diag S^6 takes {len(d6)} values): " +
          ', '.join(f"z^{L}: X0 {'ok' if a else 'FAIL'}, X1 {'ok' if b else 'FAIL'}" for L, a, b in res), flush=True)
