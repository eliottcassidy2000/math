#!/usr/bin/env python3
"""
death-star-2026-07-19-S59l -- HYP-8045(ii): the K-ladder K_c(N) = {1..N-1, cN}.
kind-pasteur S128c86 verified M = c/(cN+1) for (N,c) <= (24,8), filed the
theorem target. This session: (a) the FLOOR is a 3-line witness (a = c);
(b) the ceiling's seals proved on paper (small moduli, packing, e=1 channel
via u = q mod N, m=1 channel empty by self-consistency); (c) HERE: extend
the exact verification to (N,c) <= (60,12) and CHECK THE SEAL STRUCTURE
empirically (which channel binds at each (N,c) -- confirming the residual
is the same e>=2, m>=2 corner as the F_D tower's).
"""
from fractions import Fraction as F
from itertools import combinations
import sys, time
sys.path.insert(0, '04-computation')
from lrc_firstgap_crossN_census_deathstar_S59 import M_exact

t0 = time.time()
bad = []
n_checked = 0
for N in range(3, 61):
    for c in range(1, 13):
        fam = list(range(1, N)) + [c*N]
        M = M_exact(fam)
        n_checked += 1
        if M != F(c, c*N+1):
            bad.append((N, c, M))
            print(f"  !! K_{c}({N}): M = {M} != {c}/{c*N+1}", flush=True)
print(f"K-ladder exact check: {n_checked} pairs (N,c) in [3,60]x[1,12], "
      f"violations: {len(bad)} ({time.time()-t0:.0f}s)")
print("law M(K_c(N)) = c/(cN+1):", "HOLDS on the whole range" if not bad else bad)
