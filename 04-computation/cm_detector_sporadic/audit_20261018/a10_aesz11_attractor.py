"""Audit out-of-sample test of HYP-9331's conjecture on AESZ 11 = 4F3(1/4,1/3,2/3,3/4), u_k = (4k)!(3k)!/((2k)! k!^5)
(holomorphic period in the AESZ variable z), which has a rational rank-2 attractor point z0 = -1/(2^4 3^3) = -1/432
(Boenisch-Elmi-Kashani-Poor-Klemm arXiv:2204.06506, Table 1) and conifold z = 1/1728.
Prediction of HYP-9331: (B1') survivors = {-432, 1728} (attractor + singular point); v_p = 2 or 3 at -432.
Scan 0 < |m| <= M with the residue tables of a4_detector_scan.py; valuations via a7_slopes_valuations.analyse_point."""
import sys, os
from math import factorial
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a4_detector_scan import tables, PRIMES
from a7_slopes_valuations import analyse_point
M = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
u = [factorial(4*k) * factorial(3*k) // (factorial(2*k) * factorial(k) ** 5) for k in range(420)]
print('u =', u[:5])
lucas_bad = [p for p in PRIMES if any((u[p+j] - u[1]*u[j]) % p for j in range(p))]
print('Lucas congruence fails at', lucas_bad)
ms = np.concatenate([np.arange(-M, 0), np.arange(1, M + 1)]).astype(np.int64)
aB = np.ones(ms.size, bool); aO = np.ones(ms.size, bool)
for p in PRIMES:
    pB, pO = tables(u, p)
    idx = ms % (p * p)
    aB &= pB[idx]; aO &= pO[idx]
print(f"(B1') survivors 0<|m|<={M}:", ms[aB].tolist(), "  (B1'_ord):", ms[aO].tolist())
for m in (-432, 1728, 5, 7):
    rows = analyse_point(u, m, PRIMES)
    print(f'm={m}: ordinary v_p {dict(sorted(Counter(r["v"] for r in rows if r["ordinary"]).items()))}; '
          f'non-ordinary v_p {dict(sorted(Counter(r["v"] for r in rows if not r["ordinary"]).items()))}')
