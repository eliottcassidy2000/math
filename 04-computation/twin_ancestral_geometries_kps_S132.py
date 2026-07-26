#!/usr/bin/env python3
r"""
twin_ancestral_geometries_kps_S132.py
(kind-pasteur-2026-07-26-S132; companion to HYP-9026)

Two canonical ancestral geometries on the twin ranks K = A002822
(typing per MISTAKE-268), both rooted in the parent-fibre law of
THM-2422/2443, inspired by the Farey/Stern-Brocot mediant reading:

  STAIRCASE (min-partner, THM-2422 canonical): parent(k) = k - a*(k),
  a*(k) = min partner.  Measured: depth grows like a constant
  fraction of the rank INDEX (max 118,503 at the census top; the
  mean partner E[a*] ~ 150 at 10^7 scale), and the a* histogram is
  the FIRST-HIT transform of the pair-correlation harmonics: the
  leading partners are 5, 7, 30, 70, 40 -- 70 = 2*5*7 outranking
  every smaller non-smooth value because it is boosted at both
  p = 5 and p = 7.

  MEDIANT (balanced): parent_bal(k) = k - a_bal, a_bal = max
  {a <= k/2 : a, k-a in K}.  Measured: exists for EVERY census rank
  k >= 3 (no failures), found within 26.6 candidates on average,
  and gives a genuine Stern-Brocot-like tree: max depth 28 vs
  log2(max k) = 24; per-window mean depth = log2 k + excess with
  excess ~ 3.4 falling monotonically (1.38 -> 1.16 as ratio).

The contrast is the point: the staircase reads the harmonic
(singular-series) layer, the mediant tree reads the entropic
(density) layer; each balanced split is a windowed Seymour
certificate (u + v = w with v <= k/2 < w).

Checks (hard failures):
  * no orphan k >= 3 in either tree;
  * staircase a*-histogram top-2 = {5, 7} and 70 in top-4
    (harmonic ranking);
  * balanced max depth <= log2(max k) + 8;
  * balanced ratio mean_depth/log2(k) decreasing across windows.

Universe: centers <= 1e8. Reproduction:
python 04-computation/twin_ancestral_geometries_kps_S132.py
"""
import numpy as np
from collections import Counter
import bisect
import math

LIMIT = 100_000_000

def fail(msg):
    raise SystemExit("CHECK FAILED: " + msg)

sieve = np.ones(LIMIT + 3, dtype=bool)
sieve[:2] = False
for p in range(2, int((LIMIT + 2) ** 0.5) + 1):
    if sieve[p]:
        sieve[p * p:: p] = False
mid = (np.where(sieve[:-2] & sieve[2:])[0] + 1).astype(np.int64)
K = (mid[mid >= 6] // 6)
Kl = K.tolist()
Kset = set(Kl)
print(f"universe: |K| = {len(Kl)}, max k = {Kl[-1]}")

# ---------- staircase ----------
astar = {}
for k in Kl:
    if k < 3:
        continue
    for a in Kl:
        if 2 * a >= k:
            break
        if (k - a) in Kset:
            astar[k] = a
            break
orph = [k for k in Kl if k >= 3 and k not in astar]
if orph:
    fail(f"staircase orphans {orph[:5]}")
parent_s = {k: k - a for k, a in astar.items()}

depth_s = {1: 0, 2: 0}
def dget(par, dep, k):
    stack = []
    while k not in dep:
        stack.append(k)
        k = par[k]
    d = dep[k]
    for kk in reversed(stack):
        d += 1
        dep[kk] = d
    return d

maxd_s = 0; arg_s = None
for k in Kl:
    if k < 3:
        continue
    d = dget(parent_s, depth_s, k)
    if d > maxd_s:
        maxd_s, arg_s = d, k
print(f"staircase: max depth {maxd_s} at k={arg_s}")
ah = Counter(astar.values())
top = [a for a, _ in ah.most_common(4)]
print(f"staircase a* top-4: {top}")
if set(top[:2]) != {5, 7} or 70 not in top:
    fail(f"harmonic first-hit ranking broken: {top}")
print("staircase window means (depth, index-fraction):")
idx = {k: i for i, k in enumerate(Kl)}
for j in range(8, 25, 4):
    lo, hi = 2 ** j, 2 ** (j + 2)
    sel = [k for k in Kl if lo <= k < hi and k in depth_s]
    if len(sel) < 50:
        continue
    md = np.mean([depth_s[k] for k in sel])
    mi = np.mean([idx[k] for k in sel])
    print(f"  [2^{j},2^{j+2}): depth {md:9.1f}   depth/index {md/mi:.3f}")

# ---------- mediant ----------
parent_b = {}
scan_tot = 0
for k in Kl:
    if k < 3:
        continue
    i = bisect.bisect_right(Kl, k // 2) - 1
    found = None
    while i >= 0:
        a = Kl[i]
        scan_tot += 1
        if (k - a) in Kset:
            found = a
            break
        i -= 1
    if found is None:
        fail(f"no balanced parent for {k}")
    parent_b[k] = k - found
print(f"\nmediant: mean scans {scan_tot/len(parent_b):.1f}")

depth_b = {1: 0, 2: 0}
maxd_b = 0; arg_b = None
for k in Kl:
    if k < 3:
        continue
    d = dget(parent_b, depth_b, k)
    if d > maxd_b:
        maxd_b, arg_b = d, k
lim = math.log2(Kl[-1]) + 8
print(f"mediant: max depth {maxd_b} at k={arg_b} (bound log2+8 = {lim:.1f})")
if maxd_b > lim:
    fail("mediant depth exceeds log2+8")
print("mediant window means (depth, ratio to log2 k):")
prev = None
mono = True
for j in range(8, 25, 2):
    lo, hi = 2 ** j, 2 ** (j + 2)
    sel = [k for k in Kl if lo <= k < hi and k in depth_b]
    if len(sel) < 50:
        continue
    md = float(np.mean([depth_b[k] for k in sel]))
    r = md / (j + 1)
    print(f"  [2^{j},2^{j+2}): depth {md:6.2f}   ratio {r:.3f}")
    if prev is not None and r > prev + 1e-9:
        mono = False
    prev = r
if not mono:
    fail("mediant ratio not monotone decreasing")

print("\nALL CHECKS PASSED")
