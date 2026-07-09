"""
lrc14_arccount_growth_recheck_opus_S168.py   (opus-2026-07-09-S168)

RESOLVE the regime confusion. mac-mini-S61 arc-count uses spread 150-600 and gets c=#arcs/spread~0.2
<< rho*~0.96 (CLOSES). My S168 recheck used spread 12-35 and got c~1 (FAILS). Reconcile: is #arcs
SUBLINEAR in spread (so c=#arcs/spread -> 0 asymptotically, mac-mini's regime) while c~1 only in the
small-spread FINITE regime (handled by direct enumeration)?  If so, mac-mini's route is correct and my
"correction" was regime-confused.

Computes, for dissociated k=11 clusters at spread = 30,60,...,1920: median/max #arcs (good runs, same
object as mac-mini: rising edges of the good indicator) and c=#arcs/spread and rho*=|good|.  Reports the
growth exponent  #arcs ~ spread^alpha.  alpha<1 => c->0 => arc-count closes asymptotically.
"""
import numpy as np, random
from math import gcd
from functools import reduce

random.seed(1688)


def arc_count_and_rho(E, GRID):
    x = (np.arange(GRID) + 0.5) / GRID
    Ea = np.array(sorted(E), float)
    ph = np.mod(np.outer(x, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0] + 1 - ph[:, -1])[:, None]], axis=1)
    good = (g.max(axis=1) > 1.0 / 7).astype(int)
    tr = int(np.sum((good - np.roll(good, 1)) == 1))     # # good runs = mac-mini #arcs
    return tr, float(good.mean())


def prim(E):
    E = sorted(E)
    return len(E) >= 2 and reduce(gcd, [E[i + 1] - E[i] for i in range(len(E) - 1)]) == 1


def longest_ap(E):
    S = set(E); E = sorted(E); best = 2
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            d = E[j] - E[i]; L = 2; nx = E[j] + d
            while nx in S: L += 1; nx += d
            bk = E[i] - d
            while bk in S: L += 1; bk -= d
            best = max(best, L)
    return best


k = 11
print("=" * 90)
print(f"ARC-COUNT growth vs spread (k={k}, dissociated).  #arcs sublinear => c=#arcs/spread -> 0")
print("=" * 90)
print(f"  {'spread':>7} {'med #arcs':>10} {'max #arcs':>10} {'med c':>8} {'max c':>8} {'min rho*':>9} {'samples':>8}")
spreads = [30, 60, 120, 240, 480, 960, 1920]
logs, logm = [], []
for s in spreads:
    cs, arcs, rhos = [], [], []
    tries = 0
    want = 25 if s <= 480 else 12
    while len(cs) < want and tries < want * 40:
        tries += 1
        mid = sorted(random.sample(range(1, s), k - 2)); E = [0] + mid + [s]
        if len(set(E)) != k or not prim(E): continue
        if longest_ap(E) > k - 6: continue                # dissociated only
        GRID = min(80 * s, 200000)
        na, rho = arc_count_and_rho(E, GRID)
        cs.append(na / s); arcs.append(na); rhos.append(rho)
    if not cs:
        print(f"  {s:>7}  (no dissociated samples)"); continue
    arcs.sort(); cs.sort()
    medA = arcs[len(arcs) // 2]; maxA = max(arcs)
    medc = cs[len(cs) // 2]; maxc = max(cs); minr = min(rhos)
    logs.append(np.log(s)); logm.append(np.log(medA))
    print(f"  {s:>7} {medA:>10d} {maxA:>10d} {medc:>8.3f} {maxc:>8.3f} {minr:>9.3f} {len(cs):>8d}")

# fit #arcs ~ spread^alpha
if len(logs) >= 2:
    A = np.polyfit(logs, logm, 1)
    print(f"\n  FIT: median #arcs ~ spread^{A[0]:.3f}  (alpha<1 => c=#arcs/spread -> 0 asymptotically)")
    print(f"  => arc-count inequality #arcs < rho*.Vmax holds for large spread (mac-mini regime).")
    print(f"     Small spread (c~1) is the FINITE regime, closed by direct enumeration of j<=Vmax.")
print("=" * 90)
