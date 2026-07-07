"""
lrc_exact_rhostar_hardcores_opus_S137.py

COUPLING THE EXACT ENGINE TO THE INTERSECTED LEDGERS (owner worklist item 1) and aiming it
at monad-S3's TWO HARD CORES (owner item: "the two hard cores as the concrete targets for
opus's avoidance-kernel machinery").

Object (THM-530): G_P = {x : ||p x|| >= 1/14 for all p in P} — an explicit finite union of
rational intervals; rho*_theta(P,E) = meas(G_P ∩ {maxgap of config(E,x) > theta}).
The order-cell engine (opus-S136) computes the superlevel of the max-of-affines per subcell;
intersecting each superlevel interval with the G_P interval list gives EXACT rho*.
Every previous rho* for general (P,E) was grid-numeric (mac-mini-S42: 80k grids).

VALIDATION TARGETS (exact):
  * meas(G_P) size-minima (THM-530): 6/7, 66/91, 55/91, 1979/4004, 2243/5880, ... ,
    m_P = 14249/252252 at P = {1,2,3,5,7,8,9,11,12,13}.
  * mac-mini-S42 probe min: rho*({1,5,9,11,13}, AP-co-offset k=8) ~ 0.3390 (80k grid).
  * THM-530 union-bound constant 1891/5880 consistency.

THE TWO HARD CORES (monad-S3 window factoring, HYP-4847): after the rational-window
clearing, the mixed-shape tight-cluster domain of the k=13 floor leg collapses to
  CORE-8:  P = {9,10,11,12,13}, k = |E| = 8   (E co-offset, 0 in E)
  CORE-9:  P = {10,11,12,13},   k = |E| = 9
in the regime spread(E) > 108 with multi-q CRT-blocking. This script: exact meas(G_P) for
both cores; exact rho* on structured CRT-blocked + spread shapes; exact adversarial descent
minimizing rho* over E (the load-bearing quantity, now exact); quasi-independence ratio
R = rho*/(meas(G_P)*mu) per family.
"""
from fractions import Fraction as F
from math import gcd
import random, time, sys

sys.path.insert(0, ".")
from lrc_exact_mu_ordercells_opus_S136 import order_cells, cell_gap_affines

THETA = F(1, 7)
THR = F(1, 14)

def GP_intervals(P, thr=THR):
    """Sorted disjoint intervals of G_P = {x in [0,1]: ||p x|| >= thr for all p in P}."""
    bad = []
    for p in P:
        w = F(thr, p)
        for m in range(0, p + 1):
            lo, hi = F(m, p) - w, F(m, p) + w
            lo = max(lo, F(0)); hi = min(hi, F(1))
            if hi > lo:
                bad.append((lo, hi))
    bad.sort()
    merged = []
    for lo, hi in bad:
        if merged and lo <= merged[-1][1]:
            if hi > merged[-1][1]:
                merged[-1] = (merged[-1][0], hi)
        else:
            merged.append((lo, hi))
    good = []
    cur = F(0)
    for lo, hi in merged:
        if lo > cur:
            good.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        good.append((cur, F(1)))
    return good

def meas(intervals):
    return sum(h - l for l, h in intervals)

def overlap(good, lo, hi):
    tot = F(0)
    for a, b in good:
        if b <= lo: continue
        if a >= hi: break
        l = max(a, lo); h = min(b, hi)
        if h > l: tot += h - l
    return tot

def rho_star_exact(P, E, theta=THETA):
    """EXACT meas(G_P ∩ {maxgap(config(E,x)) > theta}). Also returns mu exact."""
    good = GP_intervals(P)
    rho = F(0); mu = F(0)
    for a, b in order_cells(E):
        gaps = cell_gap_affines(E, a, b)
        subbp = set([a, b])
        for i in range(len(gaps)):
            ci, bi = gaps[i]
            for j in range(i + 1, len(gaps)):
                cj, bj = gaps[j]
                if ci != cj:
                    xc = (bj - bi) / (ci - cj)
                    if a < xc < b:
                        subbp.add(xc)
        for u, v in zip(sorted(subbp), sorted(subbp)[1:]):
            m2 = (u + v) / 2
            cb, bb = max(gaps, key=lambda t: t[0] * m2 + t[1])
            fu = cb * u + bb; fv = cb * v + bb
            if fu > theta and fv > theta:
                slo, shi = u, v
            elif fu <= theta and fv <= theta:
                continue
            else:
                xs = (theta - bb) / cb
                slo, shi = (u, min(v, xs)) if fu > theta else (max(u, xs), v)
            if shi > slo:
                mu += shi - slo
                rho += overlap(good, slo, shi)
    return rho, mu

def main():
    t0 = time.time()
    print("=" * 104)
    print("(0) VALIDATION")
    print("=" * 104)
    mins = {1: F(6,7), 2: F(66,91), 3: F(55,91), 4: F(1979,4004), 5: F(2243,5880)}
    # size minima: brute force P over subsets is big; spot-check the known minimizers
    known_min_P = {10: [1,2,3,5,7,8,9,11,12,13]}
    mP = meas(GP_intervals(known_min_P[10]))
    print(f"   meas(G_P) at THM-530's |P|=10 minimizer = {mP} "
          f"{'== m_P MATCH' if mP == F(14249,252252) else '*** MISMATCH ***'}")
    for psz, val in mins.items():
        # minimize over all P of that size within {1..13} (exact, feasible for small psz)
        from itertools import combinations
        best = None
        for P in combinations(range(1, 14), psz):
            m = meas(GP_intervals(P))
            if best is None or m < best[0]: best = (m, P)
        ok = "MATCH" if best[0] == val else f"*** got {best[0]} expected {val} ***"
        print(f"   min_|P|={psz} meas(G_P) = {best[0]} at {best[1]}  {ok}")
    # mac-mini-S42 probe: P={1,5,9,11,13}, E = AP co-offset k=8 = {0..7}
    r, mu = rho_star_exact([1,5,9,11,13], list(range(0, 8)))
    print(f"   rho*(P={{1,5,9,11,13}}, E={{0..7}}) = {r} = {float(r):.6f}  "
          f"[mac-mini-S42 80k-grid: 0.3390]   mu = {float(mu):.6f}")

    print()
    print("=" * 104)
    print("(1) THE TWO HARD CORES (monad-S3): exact meas(G_P), then exact rho* sweeps")
    print("=" * 104)
    CORE8 = [9, 10, 11, 12, 13]   # k = 8
    CORE9 = [10, 11, 12, 13]      # k = 9
    g8 = GP_intervals(CORE8); g9 = GP_intervals(CORE9)
    print(f"   meas(G_{{9..13}})  = {meas(g8)} = {float(meas(g8)):.6f}   [monad grid: 0.447]")
    print(f"   meas(G_{{10..13}}) = {meas(g9)} = {float(meas(g9)):.6f}   [monad grid: 0.525]")
    mPq = F(14249, 252252)

    rng = random.Random(137137)
    def sweep(P, k, label):
        print(f"\n   --- {label}: P={P}, k={k}, bar m_P = {float(mPq):.6f} ---")
        fams = {
            "AP co-offset": list(range(0, k)),
            "spread 15j": [15 * j for j in range(k)],
            "spread 14j (q=14-aligned)": [14 * j for j in range(k)],
            "CRT-block 2,3(,5)": None,
            "two-cluster 0..3,109..": list(range(0, k - 4)) + [109 + j for j in range(4)],
            "geometric-ish": [0, 1, 3, 7, 15, 31, 63, 127, 255][:k],
        }
        # CRT-blocking: cover all residues mod 2 and 3 (and 5 for k>=8) with spread > 108
        crt = [0, 15, 45, 62, 81, 100, 111, 118, 125][:k]
        fams["CRT-block 2,3(,5)"] = crt
        worst = None
        for name, E in fams.items():
            E = sorted(set(E))
            if len(E) != k: continue
            r, mu = rho_star_exact(P, E)
            Rq = float(r) / (float(meas(GP_intervals(P))) * float(mu)) if mu > 0 else float("nan")
            flag = " *** BELOW m_P ***" if r < mPq else ""
            print(f"      {name:<26} rho* = {str(r)[:24]:>24} = {float(r):.6f}  "
                  f"mu={float(mu):.4f}  R={Rq:.3f}  ({float(r)/float(mPq):.1f}x bar){flag}")
            if worst is None or r < worst[1]: worst = (name, r, E)
        # adversarial exact descent (spread regime, entries <= 260)
        cur = sorted(rng.sample(range(1, 200), k - 1)); cur = [0] + cur
        rcur, _ = rho_star_exact(P, cur)
        for it in range(60):
            cand = list(cur)
            i = rng.randrange(1, k)
            if rng.random() < 0.5:
                cand[i] = max(1, cand[i] + rng.choice([-9, -3, -1, 1, 3, 9]))
            else:
                cand[i] = rng.randrange(1, 260)
            cand = sorted(set([0] + [c for c in cand[1:] if c > 0]))
            if len(cand) != k: continue
            rc, _ = rho_star_exact(P, cand)
            if rc < rcur: cur, rcur = cand, rc
        print(f"      descent min (60 exact steps): rho* = {float(rcur):.6f} at {cur}"
              f"  ({float(rcur)/float(mPq):.1f}x bar){' *** BELOW ***' if rcur < mPq else ''}")
        if worst[1] < rcur:
            print(f"      structured worst: {worst[0]} at {float(worst[1]):.6f}")

    sweep(CORE8, 8, "CORE-8")
    sweep(CORE9, 9, "CORE-9")

    print(f"\n[total {time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
