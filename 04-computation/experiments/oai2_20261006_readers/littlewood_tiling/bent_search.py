import numpy as np, itertools
from tower_merit import tower
rng = np.random.default_rng(1)
for m in [4, 6, 8]:
    N = 1 << m
    H = tower(m).astype(np.int64)
    bits = ((np.arange(N)[:, None] >> np.arange(m)[None, :]) & 1)
    pairs = [(a, b) for a in range(m) for b in range(a + 1, m)]
    Pm = np.stack([bits[:, a] * bits[:, b] for a, b in pairs], axis=1)  # N x P
    best = -1; bestq = None
    P = len(pairs)
    if P <= 15:
        Q = np.array(list(itertools.product([0, 1], repeat=P)), dtype=np.int64)
    else:
        Q = rng.integers(0, 2, size=(200000, P))
    for s in range(0, Q.shape[0], 4096):
        q = Q[s:s+4096]
        D = (-1) ** ((Pm @ q.T) % 2)          # N x batch
        L = rng.integers(0, 2, size=(m, D.shape[1]))
        D2 = D * (-1) ** ((bits @ L) % 2)      # random affine twist too
        for DD in (D, D2):
            Y = np.abs(H @ DD).min(axis=0)
            i = int(Y.argmax())
            if Y[i] > best:
                best = int(Y[i]); bestq = q[i]
    print(f"m={m} N={N}: best over quadratic(+affine) sign vectors of min_i |(H d)_i| / sqrtN = {best/np.sqrt(N):.3f}")
