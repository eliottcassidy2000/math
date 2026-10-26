# Seven-vertex level of the Paley census: TT7(P_p) = p * TT6(T[QR]), TT6(X) = sum_w TT5(N+(w)), with the chain counter
# TT5(Y) = sum_u sum_z C(d, 2) (source u, second z, d = common out-degree), so TT7 counts chains w -> u -> z -> {pair}
# inside the QR-tournament.  Usage: paley_tt7.py Pmin Pmax  (primes p = 3 mod 4).  mac-mini-2026-10-25-pyth.
import sys
import numpy as np
from sympy import primerange

def tt7(p):
    sq = sorted({(x * x) % p for x in range(1, p)})
    Q = np.array(sq)
    A = np.isin((Q[None, :] - Q[:, None]) % p, sq).astype(np.int64)    # tournament on QR (vertex 0 removed)
    m = len(Q); tot = 0
    for w in range(m):
        Y = np.nonzero(A[w])[0]
        AY = A[np.ix_(Y, Y)]
        for u in range(len(Y)):
            Z = np.nonzero(AY[u])[0]
            AZ = AY[np.ix_(Z, Z)]
            MZ = AZ @ AZ.T
            for z in range(len(Z)):
                W = np.nonzero(AZ[z])[0]
                d = MZ[W, z]
                tot += int((d * (d - 1) // 2).sum())
    return p * tot

def brute_tt7(p):
    """independent check for small p: count transitive 7-sets by extending transitive chains in a fixed order"""
    sq = {(x * x) % p for x in range(1, p)}
    beats = [[((y - x) % p) in sq for y in range(p)] for x in range(p)]
    count = 0
    def ext(chain, cands):
        nonlocal count
        if len(chain) == 7: count += 1; return
        for v in cands:
            ext(chain + [v], [c for c in cands if beats[v][c]])
    ext([], list(range(p)))
    return count

if __name__ == '__main__':
    if sys.argv[1] == 'check':
        for p in (7, 11, 19, 23, 31, 43):
            a, b = tt7(p), brute_tt7(p)
            print(f"p={p}: chain TT7 = {a}, brute force = {b}, {'OK' if a == b else 'MISMATCH'}", flush=True)
        sys.exit()
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    for p in primerange(lo, hi):
        if p % 4 == 3:
            print(f"p={p} TT7={tt7(p)}", flush=True)
