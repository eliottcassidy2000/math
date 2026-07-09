"""
mac-mini-2026-07-09-S61 -- close the DISSOCIATED branch (L<=k-6) of THM-527-A via ARC-COUNT.

The finite-Vmax glue needs only EXISTENCE of a good period (any j), not small j. The arc-count
pigeonhole gives it: #{good grid j} >= rho*.Vmax - #arcs(Good_E) > 0  <=>  #arcs < rho*.Vmax.
Since spread <= Vmax, this holds if  c := #arcs/spread < rho*.  kps-S90: c is governed by
longest-AP L -- small for low L (dissociated), large for near-AP. Since dissociated L<=k-6 <= 7 <= 8
for all k<=13, the arc-count SHOULD close it.

This computes c = #arcs/spread and rho* bucketed by longest-AP L, and checks c < rho* for L <= k-6.
"""
import numpy as np
from math import gcd
from functools import reduce
import random
random.seed(6112)

def arc_count_and_rho(E, GRID):
    x = (np.arange(GRID)+0.5)/GRID
    Ea = np.array(sorted(E), float)
    ph = np.mod(np.outer(x, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0]+1-ph[:, -1])[:, None]], axis=1)
    good = (g.max(axis=1) > 1/7).astype(int)
    tr = int(np.sum((good - np.roll(good, 1)) == 1))
    return tr, good.mean()
def prim(E):
    E = sorted(E); return len(E) >= 2 and reduce(gcd, [E[i+1]-E[i] for i in range(len(E)-1)]) == 1
def longest_ap(E):
    S = set(E); best = 2; E = sorted(E)
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            d = E[j]-E[i]; L = 2; nx = E[j]+d
            while nx in S: L += 1; nx += d
            bk = E[i]-d
            while bk in S: L += 1; bk -= d
            best = max(best, L)
    return best

print("Arc-count for the dissociated branch: c=#arcs/spread and rho*, bucketed by longest-AP L\n")
for k in (11, 13):
    diss = k-6
    buckets = {}   # L -> (max c, min rho*, count)
    for _ in range(2500):
        s = random.choice([150, 300, 600])
        mid = sorted(random.sample(range(1, s), k-2)); E = [0]+mid+[s]
        if len(set(E)) != k or not prim(E): continue
        GRID = 80*s
        L = longest_ap(E)
        na, rho = arc_count_and_rho(E, GRID)
        c = na/s
        if L not in buckets: buckets[L] = [0, 1.0, 0]
        buckets[L][0] = max(buckets[L][0], c)
        buckets[L][1] = min(buckets[L][1], rho)
        buckets[L][2] += 1
    print(f"k={k} (dissociated = L <= k-6 = {diss}):")
    print(f"   {'L':>3} {'max c=#arcs/spread':>18} {'min rho*':>9} {'c<rho*?':>8} {'n':>5}")
    for L in sorted(buckets):
        mc, mr, n = buckets[L]
        tag = "DISS" if L <= diss else "near-AP"
        ok = "YES" if mc < mr else "NO"
        print(f"   {L:>3} {mc:>18.4f} {mr:>9.4f} {ok:>8} {n:>5}  [{tag}]")
    # the dissociated verdict
    dmax_c = max((buckets[L][0] for L in buckets if L <= diss), default=0)
    dmin_rho = min((buckets[L][1] for L in buckets if L <= diss), default=1)
    print(f"   => DISSOCIATED (L<={diss}): max c = {dmax_c:.4f}, min rho* = {dmin_rho:.4f}  "
          f"{'=> c < rho*, ARC-COUNT CLOSES the dissociated branch' if dmax_c < dmin_rho else '=> FAILS'}\n")
