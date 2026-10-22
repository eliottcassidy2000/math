#!/usr/bin/env python3
"""Independent pure-Python cross-check of the DP audit for ranks <= NMAX (explicit multisets, exact Fractions,
printed formulas incl. the 1/3 in c3)."""
import sys
from fractions import Fraction as Fr
NMAX = int(sys.argv[1])
def cA(k): return (Fr(k*(k+1)), Fr(-k*(k+1)*(2*k-1)), Fr(4,3)*k*(k+1)*(5*k**2-7*k+3), Fr(-2)*k*(k+1)*(14*k**3-37*k**2+39*k-15))
def cD(k): return (Fr(2*k*(k-1)), Fr(-2*k*(k-1)*(4*k-7)), Fr(8,3)*k*(k-1)*(20*k**2-74*k+69), Fr(-4)*k*(k-1)*(112*k**3-672*k**2+1378*k-963))
cE = {6: (72, -1512, 50976, -2011824), 7: (126, -4158, 223272, -14196924), 8: (240, -13680, 1281600, -143445600)}
T = [('A', k, cA(k)) for k in range(1, NMAX+1)] + [('D', k, cD(k)) for k in range(4, NMAX+1)] + [('E', k, tuple(map(Fr, cE[k]))) for k in (6,7,8) if k <= NMAX]
T.sort(key=lambda x: (x[1], x[0]))
from collections import defaultdict
sigs = [defaultdict(list) for _ in range(NMAX+1)]
def rec(i, n, s, ms):
    # enumerate multisets with types index >= i (non-increasing order of choice)
    sigs[n][s].append(tuple(ms))
    for j in range(i, len(T)):
        t, k, v = T[j]
        if n + k > NMAX: continue
        rec(j, n + k, tuple(a+b for a, b in zip(s, v)), ms + [f"{t}{k}"])
rec(0, 0, (Fr(0),)*4, [])
tot = 0
for n in range(1, NMAX+1):
    cnt = sum(len(v) for v in sigs[n].values()); tot += cnt
    coll = [v for v in sigs[n].values() if len(v) > 1]
    if coll: print("rank", n, "collisions:", coll)
print(f"ranks 1..{NMAX}: {tot} multisets enumerated explicitly; collisions found at ranks:",
      [n for n in range(1, NMAX+1) if any(len(v) > 1 for v in sigs[n].values())])
