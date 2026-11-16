#!/usr/bin/env python3
"""nu_factor.py q imax: nu(2^i q) = number of primes of order exactly 2^i q, by factoring the cyclotomic value Phi_d(2)
with sympy (ECM/pollard inside factorint), one line per level, flushing as it goes.  Records each factorisation."""
import sys, time
from sympy import cyclotomic_poly, factorint, n_order, isprime
q = int(sys.argv[1]); imax = int(sys.argv[2])
tot = 0
for i in range(imax + 1):
    d = 2**i * q; t = time.time()
    val = int(cyclotomic_poly(d, 2))
    f = factorint(val)
    prims = [p for p in f if n_order(2, p) == d]
    tot += len(prims)
    print(f"q={q} i={i} d={d} Phi_d(2)={val if val < 10**30 else str(val)[:12]+'...'} ({len(str(val))} digits) factors={sorted(f.items())} nu={len(prims)} cumulative c({i},{q})={tot} {'>= q: covers at a=' + str(i) if tot >= q else ''} ({time.time()-t:.0f}s)", flush=True)
