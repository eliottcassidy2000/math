"""finite_blocks_trend.py - misses of the L-bounded block semigroups by dyadic block (follow-up of addendum E)."""
import sys
import numpy as np
from chain_semigroup import g_list
from finite_blocks import blocks, closure

K = int(sys.argv[1]) if len(sys.argv) > 1 else 17
X = 1 << K
LS = [int(t) for t in sys.argv[2:]] or [3, 4, 5]
for L in LS:
    G = g_list(1 << (L + 1))
    bl = blocks(L, G)
    arr = [np.array([b[i] for b in bl], dtype=np.int64) for i in range(7)]
    hit = np.zeros(X + 1, dtype=np.uint8)
    table = np.full(1 << 27, -1, dtype=np.int64)
    stK = np.zeros(1 << 26, dtype=np.int64)
    states = closure(X, arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], hit, table, stK)
    odd = np.arange(1, X + 1, 2)
    miss = odd[hit[odd] == 0]
    row = []
    for k in range(8, K):
        a, b = 1 << k, 1 << (k + 1)
        row.append("%d:%d" % (k, int(np.sum((miss >= a) & (miss < b)))))
    print("L = %d (%d blocks), below 2^%d: %d live states, %d misses; misses in [2^k, 2^(k+1)), k = 8..%d: %s; largest miss %d"
          % (L, len(bl), K, states, len(miss), K - 1, " ".join(row), int(miss[-1])))
    sys.stdout.flush()
