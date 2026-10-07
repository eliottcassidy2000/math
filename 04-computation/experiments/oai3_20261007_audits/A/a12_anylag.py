#!/usr/bin/env python3
"""Audit A, item 10 (THM-4556 Update 2): Haar model of Mersenne switching in the Terras clock.
X Haar on 1 + 8Z_2 (precision P bits); for odd D <= DMAX the lag-D chain from (D, 3^D - 1) on (p, q_D),
q_D = 2X 3^-D - 1, driven by q_D's Terras parities. Lag-1 non-merge share q1(T) and any-lag share q(T)."""
import random, math, sys
from a3_drift import step
P = int(sys.argv[1]) if len(sys.argv) > 1 else 3200
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
DMAX = 31
TMAX = P - 64
rng = random.Random(4556)
mod = 1 << P
inv3 = pow(3, -1, mod)
t1s = []; tas = []
for s in range(NS):
    X = (rng.getrandbits(P) & ~7) | 1
    chains = {}
    for D in range(1, DMAX + 1, 2):
        q = (2 * X * pow(inv3, D, mod) - 1) % mod
        chains[D] = [D, 3 ** D - 1, q]
    t1 = ta = None
    for t in range(TMAX):
        for D in list(chains):
            k, N, q = chains[D]
            if k == 0 and N == 0:
                if ta is None: ta = t
                if D == 1: t1 = t
                del chains[D]; continue
            k, N = step(k, N, q & 1)
            q = (q >> 1) if q % 2 == 0 else ((3 * q + 1) >> 1)
            chains[D] = [k, N, q]
        if t1 is not None and ta is not None: break
        if ta is not None and 1 not in chains: break
        if ta is not None:  # only lag 1 still matters
            for D in list(chains):
                if D != 1: del chains[D]
    t1s.append(t1); tas.append(ta)
print(f"Haar model, {NS} samples, precision {P} bits, odd D <= {DMAX}:")
for T in (100, 200, 400, 800, 1600, 3200, 6400):
    if T > TMAX: break
    q1 = sum(1 for x in t1s if x is None or x > T) / NS
    qa = sum(1 for x in tas if x is None or x > T) / NS
    print(f"  T={T:5d}: lag-1 q1={q1:.4f} (sqrtT q1={math.sqrt(T)*q1:6.2f})   any-lag q={qa:.4f} (sqrtT q={math.sqrt(T)*qa:6.2f})")
