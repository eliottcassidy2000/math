#!/usr/bin/env python3
"""S43 referee: THM-4669 at chosen Frobenius points (p, r), exact big integers (parallel helper for audit_thm4669 part 2).
Usage: python frob_cases.py p1:r1 p2:r2 ...
"""
import sys, time
from audit_thm4669 import s7, apery, vp_int

for arg in sys.argv[1:]:
    p, r = map(int, arg.split(':'))
    t = time.time()
    n = (p ** r - 1) // 2
    d = s7(n) - (-1) ** n * apery(n)
    print('[2] p = %d, r = %d, n = %d: mod p^2 %s, mod p^3 %s, v_p(diff) = %d  (%.1fs)'
          % (p, r, n, d % (p * p) == 0, d % p ** 3 == 0, vp_int(d, p), time.time() - t), flush=True)
