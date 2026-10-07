import numpy as np
from tower_merit import tower
from tower_aggregate import fwht
for k in range(1, 12):
    H = tower(k); N = H.shape[0]
    W = fwht(H)
    nz = (W != 0).sum(axis=1)
    vals = set()
    for i in range(N):
        vals |= set(np.unique(np.abs(W[i][W[i] != 0])).tolist())
    # for 4-term rows check xor of the support indices
    xorok = True; supports = []
    for i in range(N):
        sup = np.nonzero(W[i])[0]
        if len(sup) == 4:
            x = 0
            for s in sup: x ^= int(s)
            if x != 0: xorok = False
    cnt = {int(c): int((nz == c).sum()) for c in np.unique(nz)}
    print(k, N, "support sizes:", cnt, "abs values / N:", sorted(v / N for v in vals), "4-supports xor to 0:", xorok)
