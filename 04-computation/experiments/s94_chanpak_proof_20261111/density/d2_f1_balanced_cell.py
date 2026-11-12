# S94 density lane. The balanced cell of F1: N(h1,h2) = eta(h1) eta(h2) + 2^sigma (2h1+1)(2h2+1), h1 <= h2 both of
# n bits (top bit 1).  Exact integer statistics: distinct values, collision energy, same-row collisions
# (same h1, different h2), share of coprime states, and the second-moment constant
#   K = E * (number of odd integers in [min N, max N]) / M^2     (E = sum_N r(N)^2, M = number of pairs).
import sys, numpy as np
from math import gcd
def build_q(n):
    N = 1 << n
    q = np.zeros(N, np.int64)
    for h in range(1, N):
        q[h] = q[h >> 1] + ((h & 1) * (((h >> 1) + 1) >> 1))
    return q
sigma = int(sys.argv[1]); ns = [int(a) for a in sys.argv[2:]]
for n in ns:
    q = build_q(n)
    lo, hi = 1 << (n - 1), 1 << n
    h = np.arange(lo, hi, dtype=np.int64)
    x = 2 * h + 1 - 2 * q[lo:hi]; y = 2 * h + 1
    cop = np.array([gcd(int(a), int(b)) == 1 for a, b in zip(x, y)])
    # all ordered pairs matrix
    Nmat = np.outer(x, x) + (1 << sigma) * np.outer(y, y)
    # same-row collisions: for each h1 the number of ordered pairs h2 != h2' with equal N
    rowcoll = 0; rowcoll_cop = 0
    for i in range(hi - lo):
        u, c = np.unique(Nmat[i], return_counts=True)
        k = int((c * (c - 1)).sum())
        rowcoll += k
        if cop[i]: rowcoll_cop += k
    iu = np.triu_indices(hi - lo)
    vals = Nmat[iu]
    M = vals.size
    u, c = np.unique(vals, return_counts=True)
    E = int((c.astype(np.int64) ** 2).sum())
    span_odd = (int(vals.max()) - int(vals.min())) // 2 + 1
    K = E * span_odd / M ** 2
    print("sigma=%d n=%d pairs M=%d distinct=%d (%.4f of M; %.4f of the odd numbers of the span) E=%d E/M=%.4f K=%.4f max r=%d"
          % (sigma, n, M, u.size, u.size / M, u.size / span_odd, E, E / M, K, int(c.max())))
    print("    coprime states %d of %d; same-row ordered collisions: all rows %d, coprime rows %d"
          % (int(cop.sum()), hi - lo, rowcoll, rowcoll_cop))
    sys.stdout.flush()
