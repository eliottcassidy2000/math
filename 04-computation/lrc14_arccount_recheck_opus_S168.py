"""
lrc14_arccount_recheck_opus_S168.py   (opus-2026-07-09-S168, recheck)

CAREFUL recheck of the arc-count inequality, correcting S168 v1 (which over-counted arcs via grid
noise and used the wrong rho*).  The good-period existence pigeonhole (mac-mini route c):
  a good grid point {j/V} exists in the GOOD set G={x: maxgap{frac(e_i x)}>1/7} once  #g < mu*V,
  #g = number of maximal intervals of G, mu = |G| (the ACTUAL good measure = density floor).
In the hard region V in (spread, 7spread/6], the binding case V->spread gives  #g/spread < mu.

This computes, for clusters bucketed by longest-AP L: mu (good measure, fine grid), #g (good-interval
count, NOISE-FILTERED by a min-length), and the ratio R(L) = max #g/(mu*spread).  R<1 => arc-count
closes.  Reports honestly whether R<1 holds and for which L (correcting the S167 rho*=D3_inf claim).
"""
import sys, random


def good_stats(E, ngrid=None, minlen_cells=3):
    """mu = |G|, #g = # good intervals (merged across gaps < minlen_cells cells)."""
    spread = max(E); k = len(E)
    if ngrid is None:
        ngrid = max(6000, 200 * spread)
    thr = 1.0 / 7
    good = bytearray(ngrid)
    gm_sum = 0
    for t in range(ngrid):
        x = (t + 0.5) / ngrid
        ph = sorted((e * x) % 1.0 for e in E)
        mg = 1.0 - ph[-1] + ph[0]
        for i in range(k - 1):
            g = ph[i + 1] - ph[i]
            if g > mg: mg = g
        if mg > thr:
            good[t] = 1; gm_sum += 1
    mu = gm_sum / ngrid
    # count good intervals, merging good runs separated by < minlen_cells bad cells (noise filter),
    # and ignoring good runs shorter than minlen_cells (spurious threshold-touching)
    # circular:
    runs = []
    t = 0
    while t < ngrid:
        if good[t]:
            s = t
            while t < ngrid and good[t]:
                t += 1
            runs.append((s, t))  # [s, t)
        else:
            t += 1
    # circular merge of first/last
    if runs and runs[0][0] == 0 and runs[-1][1] == ngrid and len(runs) > 1:
        runs[0] = (runs[-1][0] - ngrid, runs[0][1]); runs.pop()
    # filter tiny runs, merge across tiny bad gaps
    runs = [(s, e) for (s, e) in runs if e - s >= minlen_cells]
    # merge runs separated by < minlen_cells
    merged = []
    for r in sorted(runs):
        if merged and r[0] - merged[-1][1] < minlen_cells:
            merged[-1] = (merged[-1][0], r[1])
        else:
            merged.append(list(r))
    ng = len(merged)
    return mu, ng


def longest_ap(E):
    S = set(E); E = sorted(E); best = 1
    for i, a in enumerate(E):
        for b in E[i + 1:]:
            d = b - a
            if a - d in S: continue
            L = 2; x = b + d
            while x in S: L += 1; x += d
            best = max(best, L)
    return best


def main():
    print("=" * 92)
    print("ARC-COUNT recheck: R(L)=max #g/(mu*spread) < 1 ?  (#g good intervals, mu good measure)")
    print("  good period exists iff #g < mu*V; hard region V~spread => need #g/spread < mu i.e. R<1")
    print("=" * 92)
    r = random.Random(1688)
    for k in (11, 13):
        buck = {}   # L -> (max R, sample mu, sample #g, count)
        cands = []
        for _ in range(400):
            spread = r.randint(k, 32)
            E = sorted(set([0] + r.sample(range(1, spread), k - 2) + [spread]))
            if len(E) == k: cands.append(E)
        for d in range(1, 4):
            base = [i * d for i in range(k)]
            for jd in range(1, k):
                for dl in (-1, 1):
                    E = sorted(set(x + (dl if idx == jd else 0) for idx, x in enumerate(base)))
                    if len(E) == k and min(E) == 0: cands.append(E)
        for E in cands:
            L = longest_ap(E); spread = max(E)
            mu, ng = good_stats(E)
            if mu < 1e-6: continue
            R = ng / (mu * spread)
            b = buck.setdefault(L, [0.0, mu, ng, spread, 0])
            if R > b[0]: b[0] = R; b[1] = mu; b[2] = ng; b[3] = spread
            b[4] += 1
        print(f"\n  k={k}:  L | max R=#g/(mu*spread) | (mu, #g, spread at max) | R<1?")
        for L in sorted(buck):
            b = buck[L]
            print(f"      {L:2d} |   {b[0]:.3f}              | (mu={b[1]:.3f}, #g={b[2]}, spr={b[3]}) "
                  f"| {'YES' if b[0] < 1 else 'NO'}  (n={b[4]})")
    print()
    print("  READING: if R(L)<1 across L, arc-count closes (good period exists); the honest a-priori")
    print("  needs #g<=(bound) and mu>=(floor). Corrects S167's rho*=D3_inf conflation if R>=1 anywhere.")
    print("=" * 92)


if __name__ == "__main__":
    main()
