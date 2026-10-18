"""Audit out-of-sample test of HYP-9331 on AESZ 36 (non-hypergeometric):
D = theta^4 - 2^4 z (2theta+1)^2 (3theta^2+3theta+1) + 2^9 z^2 (2theta+1)^2 (2theta+3)^2,
A_n = C(2n,n)^2 sum_k C(n,k) C(2k,k) C(2n-2k,n-k) (checked against D below);
BEKK arXiv:2204.06506 Table 1: rational rank-2 attractor z0 = -1/2^6 -> m = -64; conifolds m = 64, 128.
Prediction of HYP-9331: (B1') survivors = {-64, 64, 128}."""
import sys, os
from math import comb
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a4_detector_scan import tables, PRIMES
from a7_slopes_valuations import analyse_point
M = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
E = [sum(comb(n, k) * comb(2*k, k) * comb(2*n - 2*k, n - k) for k in range(n + 1)) for n in range(420)]
A = [comb(2*n, n) ** 2 * E[n] for n in range(420)]
bad = [n for n in range(2, 420) if n**4 * A[n] != 16*(2*n-1)**2*(3*n*n-3*n+1)*A[n-1] - 512*(2*n-3)**2*(2*n-1)**2*A[n-2]]
print('A =', A[:5], '; satisfies the AESZ 36 operator for 2 <= n < 420:', not bad, bad[:3], '; n=1:', A[1] == 16*1*1*A[0])
print('Lucas congruence fails at', [p for p in PRIMES if any((A[p+j] - A[1]*A[j]) % p for j in range(p))])
ms = np.concatenate([np.arange(-M, 0), np.arange(1, M + 1)]).astype(np.int64)
aB = np.ones(ms.size, bool); aO = np.ones(ms.size, bool)
for p in PRIMES:
    pB, pO = tables(A, p)
    idx = ms % (p * p)
    aB &= pB[idx]; aO &= pO[idx]
print(f"(B1') survivors 0<|m|<={M}:", ms[aB].tolist(), "  (B1'_ord):", ms[aO].tolist())
for m in (-64, 64, 128, 5):
    rows = analyse_point(A, m, PRIMES)
    print(f'm={m}: ordinary v_p {dict(sorted(Counter(r["v"] for r in rows if r["ordinary"]).items()))}; '
          f'non-ordinary v_p {dict(sorted(Counter(r["v"] for r in rows if not r["ordinary"]).items()))}')
