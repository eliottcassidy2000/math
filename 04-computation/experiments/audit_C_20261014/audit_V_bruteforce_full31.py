#!/usr/bin/env python3
"""Audit C: fully unnormalised injective 6-tuple enumeration at q = 25, 27, 31 (complements audit_V_bruteforce.py)."""
import sys, time
sys.path.insert(0, '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_C_20261014')
from audit_V_bruteforce import brute, prime_tables, field_tables
from audit_common import GF
t0 = time.time()
for p, r in ((5, 2), (3, 3), (31, 1)):
    q = p ** r
    if r == 1:
        sub, chi = prime_tables(q)
    else:
        sub, chi = field_tables(GF(p, r))
    V = brute(q, sub, chi, 'full')
    print(f"q = {q}: V_q by full unnormalised enumeration = {V}  [{time.time() - t0:.1f}s]", flush=True)
