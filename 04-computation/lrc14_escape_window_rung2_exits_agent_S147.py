#!/usr/bin/env python3
"""
D3: the SECOND RUNG of the window ladder.
Re-collect the two-far [14,48]-residue escapers up to B=20000 (all 78
complements), integer-verify a sample against the direct witness scan (numpy
cross-check), then compute every survivor's FIRST EXIT q* >= 49 and its duty/
rung admissibility.  Output: the rung-2 exit histogram = the [49, ?] window.
"""
import time, random
from math import gcd
from itertools import combinations
from functools import reduce
import numpy as np

T0 = time.time()
random.seed(2)

def tq(q): return -(-3*q // 41)
def dist(v, j, q):
    r = (v*j) % q
    return min(r, q-r)
def has_witness(W, q):
    t = tq(q)
    for j in range(1, q//2 + 1):
        if all(dist(v, j, q) >= t for v in W):
            return j
    return None
def first_exit(W, qlo, qhi):
    for q in range(qlo, qhi+1):
        j = has_witness(W, q)
        if j is not None:
            return q, j
    return None, None
def rung_ok(s):
    return 41*(s//14 + 1) < 3*s

def escape_bool(S, q):
    t = tq(q)
    A = [j for j in range(1, q) if all(dist(s, j, q) >= t for s in S)]
    if not A:
        return np.ones((q, q), dtype=bool)
    u = np.arange(q)
    esc = np.ones((q, q), dtype=bool)
    for j in A:
        r = (u*j) % q
        cu = np.minimum(r, q-r) >= t
        esc &= ~np.outer(cu, cu)
    return esc

B = 20000
QS = list(range(14, 49))
yidx = {q: (np.arange(B+1) % q) for q in QS}
survivors = []
for ab in combinations(range(1, 14), 2):
    S = [v for v in range(1, 14) if v not in ab]
    tables = {q: escape_bool(S, q) for q in QS}
    order = sorted(QS, key=lambda q: tables[q].mean())
    for x in range(14, B):
        bad = False
        for q in order[:6]:
            if not tables[q][x % q].any():
                bad = True; break
        if bad: continue
        mask = np.ones(B+1, dtype=bool)
        mask[:x+1] = False
        for q in order:
            mask &= tables[q][x % q][yidx[q]]
            if not mask.any(): break
        else:
            for y in np.nonzero(mask)[0]:
                survivors.append((ab, x, int(y)))
print(f"collected {len(survivors)} residue escapers of [14,48], "
      f"two-far, heights <= {B}  [{time.time()-T0:.0f}s]")

# integer cross-check on a random sample of 60
sample = random.sample(survivors, min(60, len(survivors)))
bad = 0
for (ab, x, y) in sample:
    S = [v for v in range(1, 14) if v not in ab]
    W = sorted(S + [x, y])
    if any(has_witness(W, q) is not None for q in range(14, 49)):
        bad += 1
print(f"integer re-verification of 60 random escapers: {60-bad}/60 confirmed "
      f"(numpy tables sound: {'YES' if bad==0 else 'NO'})")

# first exits and admissibility
hist = {}
noexit = []
minh, maxh = None, None
duty_ok_n = rung_ok_n = 0
lowest = []
for (ab, x, y) in survivors:
    S = [v for v in range(1, 14) if v not in ab]
    W = sorted(S + [x, y])
    q, j = first_exit(W, 49, 400)
    key = q if q is not None else 'NONE<=400'
    hist[key] = hist.get(key, 0) + 1
    if q is None:
        noexit.append((ab, x, y))
    duty = all(any(v % d == 0 for v in W) for d in range(2, 14)) and \
           not any(v % 14 == 0 for v in W)
    rung = any(rung_ok(a+b) for a, b in combinations(W, 2))
    duty_ok_n += duty; rung_ok_n += rung
    h = y
    if minh is None or h < minh: minh, argmin = h, (ab, x, y)
    lowest.append((y, x, ab))
lowest.sort()
print(f"\nheights: min max(x,y) = {lowest[0][0]} at removed={lowest[0][2]}, "
      f"pair=({lowest[0][1]},{lowest[0][0]})")
print(f"the 10 lowest escapers: {[(ab,x,y) for (y,x,ab) in lowest[:10]]}")
print(f"duty-admissible: {duty_ok_n}/{len(survivors)}; "
      f"rung-admissible: {rung_ok_n}/{len(survivors)}")
print(f"\nRUNG-2 FIRST-EXIT HISTOGRAM (q* >= 49):")
for k in sorted(hist, key=lambda z: (isinstance(z, str), z)):
    print(f"  {k}: {hist[k]}")
if noexit:
    print(f"deep escapers with NO exit <= 400 (need higher scan): {noexit[:10]}")
print(f"\nTOTAL {time.time()-T0:.0f}s")
