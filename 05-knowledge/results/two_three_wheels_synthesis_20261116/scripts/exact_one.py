#!/usr/bin/env python3
"""exact_one.py N limit [--brier]: decide a single period with a large node budget (same search as wheels3 / brier_sheets)."""
import sys, time
from sympy import factorint, n_order, divisors
N = int(sys.argv[1]); limit = int(float(sys.argv[2])); brier = '--brier' in sys.argv
mods = []
for d in divisors(N):
    if d == 1: continue
    for p in factorint(2**d - 1):
        if n_order(2, p) == d: mods.append(d)
mods = sorted(mods)
if brier: mods = [d for d in mods if d % 2 == 0]
full = (1 << N) - 1
def cls(d, r):
    m = 0
    for x in range(r % d, N, d): m |= 1 << x
    return m
used = [False] * len(mods); nodes = [0]; t0 = time.time()
def rec(cm, cp):
    nodes[0] += 1
    if cm == full and (not brier or cp == full): return True
    if nodes[0] > limit: raise TimeoutError
    if cm != full: x = (~cm & full); x = (x & -x).bit_length() - 1; sheet = 0
    else: x = (~cp & full); x = (x & -x).bit_length() - 1; sheet = 1
    tried = set()
    for i, d in enumerate(mods):
        if used[i] or d in tried: continue
        tried.add(d); used[i] = True
        r = x if sheet == 0 else x - d // 2
        if rec(cm | cls(d, r), (cp | cls(d, r + d // 2)) if brier else 0): return True
        used[i] = False
    return False
try:
    res = rec(0, 0); print(f"N={N} brier={brier} moduli={mods} result={'YES' if res else 'no'} nodes={nodes[0]} time={time.time()-t0:.0f}s")
except TimeoutError:
    print(f"N={N} brier={brier} undecided after {limit} nodes ({time.time()-t0:.0f}s)")
