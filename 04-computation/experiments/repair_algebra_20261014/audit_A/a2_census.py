#!/usr/bin/env python3
"""Audit A (independent), THM-4640 (7) / HYP-9270 (a): census N(E) of primitive values at -1.

Different algorithm from the audited lexsort census.  S_e = set of values K(t) over all increasing t with exponents
<= e (t = () included; K(()) = -1, K(t.e) = 3K(t) + 2^e).  The values with a representation of last exponent exactly e
form T_e = 3 S_(e-1) + 2^e.  A value is primitive at level E (two representations with distinct last exponents, all
exponents <= E) iff it lies in T_e and in S_(e-1) for some e <= E.  So
    P_E = P_(E-1) u (T_E n S_(E-1)),   S_E = S_(E-1) u T_E,   N(E) = |P_E|.
All arrays are sorted unique int64 (max |K| < 3^27 fits).  Usage: python3 a2_census.py [EMAX=25]
"""
import sys, time, math
import numpy as np

EMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 25
t0 = time.time()
S = np.array([-1], dtype=np.int64)
P = np.array([], dtype=np.int64)
rows = []
for e in range(1, EMAX + 1):
    T = 3 * S + (1 << e)                       # sorted, unique
    idx = np.searchsorted(S, T)
    idx[idx == len(S)] = len(S) - 1
    hit = T[S[idx] == T]                        # T_e n S_(e-1)
    del idx
    newvals = np.setdiff1d(hit, P, assume_unique=True)
    P = np.union1d(P, hit)
    rows.append((e, len(newvals), len(P), len(S) + len(T)))
    if e < EMAX:
        S = np.union1d(S, T)
    del T
    N = len(P)
    print(f"E={e:2d} new={len(newvals):7d} N(E)={N:8d}  N sqrt(E)/2^E={N*math.sqrt(e)/2**e:.6f}  "
          f"|S_E|={len(S) if e < EMAX else '-'}  [{time.time()-t0:.1f}s]", flush=True)
print("local exponent beta(E) with N(E) ~ C 2^E E^-beta, from consecutive cumulative N:")
for (e1, _, N1, _), (e2, _, N2, _) in zip(rows, rows[1:]):
    if e2 >= 14:
        beta = -math.log(N2 / (2 * N1)) / math.log(e2 / e1)
        print(f"   E={e2:2d}: N(E)/N(E-1) = {N2/N1:.5f}, beta_local = {beta:.3f}")
# smallest primitive values
print("smallest primitive values:", sorted(P.tolist(), key=abs)[:12])
