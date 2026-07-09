"""
lrc14_arccount_apriori_target_opus_S169.py   (opus-2026-07-09-S169)

Pin down the EXACT a-priori target for bounded-arc-count on the DISSOCIATED branch (longest-AP<=k-6),
the only branch that invokes arc-count (near-AP goes to klein LEM-012).  For dissociated E across a
spread sweep, measure:
  #arcs (good runs) ; rho*=meas(good) ; D3(E) (dilation-invariant density-floor lower bound on rho*) ;
  and the three candidate closing inequalities
     (i)  #arcs < rho* * (spread+1)      [pigeonhole at the binding V=spread+1 -> good period EXISTS]
     (ii) c = #arcs/spread < D3(E)        [mac-mini-S62 clean form, D3 exact + dilation-invariant]
     (iii) min good-arc length * spread   [is each arc >= const/spread? -> #arcs <= rho*/min-len]
Also the a-priori GAP: trivial bound #arcs<=2*sum|e_i-e_j| (=O(k^2 spread)) vs the true #arcs.
Reports margins so the fleet has a crisp, de-risked target for the a-priori proof.
"""
import numpy as np, random
from math import gcd
from functools import reduce

random.seed(1690)


def arc_stats(E, GRID):
    x = (np.arange(GRID) + 0.5) / GRID
    Ea = np.array(sorted(E), float)
    ph = np.mod(np.outer(x, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0] + 1 - ph[:, -1])[:, None]], axis=1)
    good = (g.max(axis=1) > 1.0 / 7).astype(int)
    up = (good - np.roll(good, 1)) == 1
    narcs = int(up.sum())
    rho = float(good.mean())
    # min good-run length (in units of 1/GRID -> fraction), circular
    minlen = 1.0
    if narcs > 0:
        idx = np.where(up)[0]
        gd = good.copy()
        for s in idx:
            t = s; ln = 0
            while gd[t % GRID] == 1:
                ln += 1; t += 1
                if ln > GRID: break
            minlen = min(minlen, ln / GRID)
    return narcs, rho, minlen


def W_moments(E, NG):
    """m1,m2,m3 = E[W^p], W = uncovered measure (continuous x, fine grid). dilation-invariant."""
    x = (np.arange(NG) + 0.5) / NG
    Ea = np.array(sorted(E), float)
    ph = np.mod(np.outer(x, Ea), 1.0); ph.sort(axis=1)
    g = np.concatenate([np.diff(ph, axis=1), (ph[:, 0] + 1 - ph[:, -1])[:, None]], axis=1)
    excess = np.clip(g - 1.0 / 7, 0, None)          # (gap - 1/7)_+
    W = excess.sum(axis=1)
    return W.mean(), (W ** 2).mean(), (W ** 3).mean()


def D3(E, NG=None):
    if NG is None:
        NG = max(20000, 300 * max(E))
    m1, m2, m3 = W_moments(E, NG)
    M = 6.0 / 7
    denom = m2 - m3 / M
    if denom <= 1e-15:
        return m1 / M
    return m1 / M + (m1 - m2 / M) ** 2 / denom


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


def prim(E):
    E = sorted(E)
    return reduce(gcd, [E[i + 1] - E[i] for i in range(len(E) - 1)]) == 1


def sum_abs_diff(E):
    return sum(abs(a - b) for i, a in enumerate(E) for b in E[i + 1:])


print("=" * 100)
print("A-PRIORI TARGET for bounded-arc-count on the DISSOCIATED branch (longest-AP<=k-6)")
print("=" * 100)
for k in (11, 13):
    diss = k - 6
    print(f"\n  k={k} (dissociated = longest-AP <= {diss}):")
    print(f"  {'spread':>7} {'#arcs':>6} {'rho*':>6} {'D3':>6} {'c=#/spr':>8} "
          f"{'(i)#<rho(spr+1)':>15} {'(ii)c<D3':>9} {'(iii)minL*spr':>13} {'trivBnd/#arcs':>13}")
    for s in [40, 80, 160, 320, 640, 1280]:
        rows = []
        tries = 0
        want = 20 if s <= 320 else 10
        while len(rows) < want and tries < want * 60:
            tries += 1
            mid = sorted(random.sample(range(1, s), k - 2)); E = [0] + mid + [s]
            if len(set(E)) != k or not prim(E) or longest_ap(E) > diss:
                continue
            GRID = min(120 * s, 300000)
            na, rho, minlen = arc_stats(E, GRID)
            d3 = D3(E)
            rows.append((na, rho, d3, minlen, sum_abs_diff(E)))
        if not rows:
            print(f"  {s:>7}  (no dissociated samples)"); continue
        # report the WORST case for each inequality (max c, min margins)
        na_arr = np.array([r[0] for r in rows]); rho_arr = np.array([r[1] for r in rows])
        d3_arr = np.array([r[2] for r in rows]); minl_arr = np.array([r[3] for r in rows])
        triv_arr = np.array([2 * r[4] for r in rows])
        c_arr = na_arr / s
        # worst (max c) representative
        w = int(np.argmax(c_arr))
        i_ok = np.all(na_arr < rho_arr * (s + 1))
        ii_ok = np.all(c_arr < d3_arr)
        minL_spr = float(minl_arr.min() * s)             # smallest arc * spread (worst)
        triv_ratio = float((triv_arr / np.maximum(na_arr, 1)).mean())
        print(f"  {s:>7} {na_arr[w]:>6d} {rho_arr[w]:>6.3f} {d3_arr[w]:>6.3f} {c_arr[w]:>8.4f} "
              f"{('YES' if i_ok else 'NO'):>15} {('YES' if ii_ok else 'NO'):>9} "
              f"{minL_spr:>13.3f} {triv_ratio:>13.0f}x")
print()
print("  READING: (i) pigeonhole #arcs<rho*(spread+1) => good period EXISTS -- the target.")
print("           (ii) mac-mini c<D3 (cleaner, exact+dilation-invariant). (iii) min-arc*spread:")
print("           if bounded below by const>0 then #arcs<=rho*/minL=O(spread) with GOOD constant.")
print("           trivBnd/#arcs = how loose the trivial O(k^2 spread) a-priori bound is (the GAP).")
print("=" * 100)
