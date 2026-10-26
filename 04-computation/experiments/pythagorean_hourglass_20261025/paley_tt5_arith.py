# Arithmetic in the 5-census of Paley tournaments: TT5(P_p) for primes p = 3 (mod 4) up to Pmax, fast O(m^3) count
# (TT5 = p * TT4(T[QR]); TT4(X) = sum_w t3(N+_X(w)); t3(Y) = sum_u C(d_Y(u), 2)), then the residual against the best
# quartic in p and its correlation with candidate CM data (p = x^2 + 2y^2, x^2 + 7y^2, x^2 + xy + 2y^2 ...).
import sys, math
import numpy as np
from sympy import primerange
Pmax = int(sys.argv[1])
rows = []
for p in primerange(19, Pmax):
    if p % 4 != 3: continue
    sq = sorted({(x * x) % p for x in range(1, p)})
    sqs = set(sq)
    Q = np.array(sq)
    # adjacency restricted to QR: a -> b iff b - a in QR
    D = (Q[None, :] - Q[:, None]) % p
    A = np.isin(D, sq).astype(np.int64)                # m x m
    m = len(sq)
    # t3 of N+(w) for every w: within Y = N+(w), d_Y(u) = #{v in Y : u -> v} = (A[u] * A[w]).sum over v
    # TT4(QR) = sum_w sum_{u in N+(w)} C(d_{N+(w)}(u), 2)
    M = A @ A.T                                         # M[u, w] = #{v : u -> v, w -> v}
    tot = 0
    for w in range(m):
        Y = np.nonzero(A[w])[0]
        d = M[Y, w]                                     # d_Y(u) for u in Y: v ranges over N+(u) cap N+(w)
        tot += int((d * (d - 1) // 2).sum())
    rows.append((p, p * tot))
    print(f"p={p} TT5={p * tot}", flush=True)
