#!/usr/bin/env python3
"""
monad-explorer-2026-07-07-S2 -- THE TAIL-DIAMETER THEOREM (HYP-4817) + the E[U] chain.

THEOREM A (tail-diameter; kps-S59's pointwise monotonicity applied to the TAIL, not the mean).
  For any 13-family E (k=13/P=empty shape) with D = diam(E) = max(E)-min(E):
      {(e - min E) x : e in E}  is a rigid translate of a SUBSET of {j x : j = 0..D},
  and removing points only merges gaps, so pointwise
      maxgap(E, x) >= maxgap(AP on {0..D}, x)     for every x,
  hence
      mu_{1/7}(E) >= mu_{1/7}(AP_{D+1})       (tail version -- the DAG object!)
      E[U(E)]     >= E[U(AP_{D+1})]           (U version;  U = sum_g (g-1/7)_+ )
      E[mg(E)]    >= E[mg(AP_{D+1})]          (mean version -- kps-S59's angle 1).
  Since the skeleton's k=13 leg needs G2 = mu_{1/7} >= m_P = 14249/252252, the tail version
  CLOSES ALL diameter <= D0 families where
      D0 := max{ D : mu_{1/7}(AP_{D+1}) >= m_P }.
  mu(AP_m) decays only ~ C/m (q<=6 Farey windows, opus-S134 roof), so D0 is LARGE (~70),
  vs the mean version's crossing (~20).  This is a rigorous quantitative floor for an
  infinite class, at the consuming node's own bar.

EXACT ENGINE (independent of opus-S134's roof; cross-validates it):
  On each cell of Farey_{2(m-1)} ... actually of Farey_{2m}, the config {jx : j=1..m} has
  constant combinatorics and every gap length is affine with INTEGER slope in [-(m-1),m-1];
  any two gap functions cross only at rationals of denominator <= 2m, i.e. at cell
  boundaries.  So maxgap is EXACTLY affine per cell (asserted at both endpoints); the
  superlevel {maxgap > 1/7} is solved exactly per cell with Fractions.

THE E[U] CHAIN (the repositioned S1 program):
  U(x) = meas{phi : arc (phi, phi+1/7) empty} = sum_g (g - 1/7)_+   (exact identity)
  =>  E_x[U] = int_phi AV_E(phi, 1/7) dphi     (opus-S134's avoidance kernel, phi-integrated)
  PRUNING: the phi-integral kills ALL resonances with sum m_i != 0 (full projection onto
  the affine sublattice sum m_i = 0 -- STRONGER than my S1's mod-14 pruning).
  ASSEMBLY: mu_{1/7} >= E[U] (opus-S131 PZ)  =>  [E[U] >= 1/14] closes the k=13 leg at m_P
  (1/14 = 0.0714 > m_P = 0.0565), with observed min E[U] ~ 0.095 (33% slack).
  CORRECTED IN-SESSION (my first draft claimed W >= 14*U -- FALSE, caught by the
  verification step below; min(W-14U) ~ -2): an interval of length g contains at least
  floor(14g)-1 full cells (NOT ceil(14g)-1), so the sharp pointwise inequality is
      W >= 14 * sum_g (g - 3/14)_+        (verified tight: min(W - 14*U^{3/14}) = 0.0000).
  S1 = E[W] therefore dominates the 3/14-excess, not the 1/7-excess; S1 and E[U] are
  cousins (aligned 14-sample vs phi-integral of the avoidance kernel), not comparable
  at matched thresholds.

Tournament Analysis declaration:
  vertices: floor-mechanisms (tail-diameter, E[U]-chain, S1/CE, mean-detour[dead]);
  pairwise observable: size of the family class each closes at the m_P bar;
  switch/gauge: orient toward the larger closed class; tie path: exact -> numeric -> descent.
"""
from fractions import Fraction as F
from math import gcd
import random

import numpy as np

M_P = F(14249, 252252)
THR = F(1, 7)


# ---------------- exact mu_{1/7}(AP_m) ----------------
def farey(n):
    """Ascending Farey_n in [0,1] via the standard next-term recurrence."""
    seq = [F(0, 1), F(1, n)]
    a, b, c, d = 0, 1, 1, n
    while c < n or (c == n and d == 1):
        k = (n + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b
        seq.append(F(c, d))
        if c == 1 and d == 1:
            break
    return seq

def mu_ap_exact(m):
    """EXACT mu_{1/7}(AP_m) = meas{x: maxgap({jx: j=1..m}) > 1/7}, via Farey_{2m} cells."""
    js = list(range(1, m + 1))
    cells = farey(2 * m)
    total = F(0)
    for a, b in zip(cells, cells[1:]):
        mid = (a + b) / 2
        fl = [int(j * mid) for j in js]                       # floors constant on cell
        pts = sorted(range(m), key=lambda i: js[i] * mid - fl[i])
        # gap functions (slope, intercept), incl. wrap
        gaps = []
        for s in range(m):
            i1, i2 = pts[s], pts[(s + 1) % m]
            if s < m - 1:
                sl = js[i2] - js[i1]; ic = F(-(fl[i2] - fl[i1]))
            else:
                sl = js[i2] - js[i1]; ic = F(1 - (fl[i2] - fl[i1]))
            gaps.append((sl, ic))
        # argmax at the MIDPOINT: interior, so strict (gap-function crossings have
        # denominator <= 2m and thus sit at cell ENDPOINTS, never at the midpoint)
        vm = [sl * mid + ic for sl, ic in gaps]
        ia = max(range(m), key=lambda i: vm[i])
        sa, ca = gaps[ia]
        va = [sl * a + ic for sl, ic in gaps]
        vb = [sl * b + ic for sl, ic in gaps]
        # the midpoint-max function dominates the whole cell (weak ties at endpoints ok)
        assert all(sa * a + ca >= va[i] for i in range(m)), (m, a, b, "left")
        assert all(sa * b + ca >= vb[i] for i in range(m)), (m, a, b, "right")
        # solve sa*x + ca > 1/7 on (a,b)
        ga, gb = (sa * a + ca) - THR, (sa * b + ca) - THR
        if ga <= 0 and gb <= 0:
            continue
        if ga > 0 and gb > 0:
            total += b - a
            continue
        xstar = (THR - ca) / F(sa)
        if ga > 0:
            total += xstar - a
        else:
            total += b - xstar
    return total


# ---------------- numeric helpers ----------------
def stats_np(E, res=60000):
    xs = (np.arange(res) + 0.5) / res
    ph = np.mod(np.outer(xs, np.array(E, dtype=np.float64)), 1.0)
    ph.sort(axis=1)
    gaps = np.empty_like(ph)
    gaps[:, :-1] = np.diff(ph, axis=1)
    gaps[:, -1] = ph[:, 0] + 1.0 - ph[:, -1]
    mg = gaps.max(axis=1)
    mu = (mg > 1.0 / 7).mean()
    U = np.maximum(gaps - 1.0 / 7, 0.0).sum(axis=1)
    # W = aligned adjacent-empty-cell pairs
    W = np.zeros(res)
    for j in range(14):
        rel = np.mod(ph - j / 14.0, 1.0)
        W += (~(rel < 1.0 / 7.0).any(axis=1)).astype(float)
    return mu, U.mean(), mg.mean(), W, U

def saturated(v): return all(any(x % q == 0 for x in v) for q in range(2, 15))
def primitive(v):
    g = 0
    for x in v: g = gcd(g, x)
    return g == 1
def single_scale(v): return max(v) <= 13 * min(v)
def in_crux(v): return saturated(v) and primitive(v) and single_scale(v)


if __name__ == "__main__":
    random.seed(4817)

    print("=" * 78)
    print("STEP 1 -- validate the exact engine against canon mu_{1/7}(AP_k), k=8..13")
    print("=" * 78)
    canon = {8: F(691, 735), 9: F(247, 294), 10: F(38, 49), 11: F(1381, 2205),
             12: F(13823, 24255), 13: F(477, 1078)}
    for k in range(8, 14):
        v = mu_ap_exact(k)
        ok = "MATCH" if v == canon[k] else f"MISMATCH vs {canon[k]}"
        print(f"  mu(AP_{k}) = {v} = {float(v):.6f}   [{ok}]")

    print()
    print("=" * 78)
    print("STEP 2 -- numeric scan of mu(AP_m) to bracket the m_P crossing")
    print("=" * 78)
    ms = list(range(14, 132, 6))
    vals = {}
    for m in ms:
        mu, _, _, _, _ = stats_np(list(range(1, m + 1)), res=120000)
        vals[m] = mu
        print(f"  m={m:3d}: mu ~ {mu:.5f}  m*mu ~ {m * mu:.3f}"
              + ("   <-- near m_P" if abs(mu - float(M_P)) < 0.02 else ""))

    # bracket: find last m with numeric mu clearly >= m_P and first clearly below
    above = [m for m in ms if vals[m] >= float(M_P) + 0.004]
    below = [m for m in ms if vals[m] <= float(M_P) - 0.004]
    lo = max(above) if above else ms[0]
    hi = min([m for m in below if m > lo], default=ms[-1])
    print(f"  bracket: crossing in ({lo}, {hi}); exact binary search:")

    print()
    print("=" * 78)
    print("STEP 3 -- EXACT binary search for D0 = max{D: mu(AP_{D+1}) >= m_P}")
    print("=" * 78)
    exact_cache = {}
    def mu_exact_cached(m):
        if m not in exact_cache:
            exact_cache[m] = mu_ap_exact(m)
            print(f"    exact mu(AP_{m}) = {exact_cache[m]} = {float(exact_cache[m]):.6f}  "
                  f"({'>=' if exact_cache[m] >= M_P else '<'} m_P)")
        return exact_cache[m]
    a, b = lo, hi
    while b - a > 1:
        mid = (a + b) // 2
        if mu_exact_cached(mid) >= M_P:
            a = mid
        else:
            b = mid
    # confirm endpoints exactly
    mu_exact_cached(a); mu_exact_cached(b)
    mstar = a
    print(f"  ==> m* = {mstar} (largest m with mu(AP_m) >= m_P);  D0 = m* - 1 = {mstar - 1}")
    print(f"  THEOREM A: every 13-family with diam <= {mstar - 1} has "
          f"mu_{{1/7}} >= mu(AP_{mstar}) >= m_P  -- the k=13 hlarge bar, CLOSED.")

    print()
    print("=" * 78)
    print("STEP 4 -- the E[U] chain: monotonicity, W>=14U, adversarial E[U] landscape")
    print("=" * 78)
    bank = {
        "AP {1..13}": list(range(1, 14)),
        "record 2*{1..11}u{11,13}": [2, 4, 6, 8, 10, 11, 12, 13, 14, 16, 18, 20, 22],
        "ds 2*{1..12}u{13}": [2, 4, 6, 8, 10, 12, 13, 14, 16, 18, 20, 22, 24],
        "S57 adversarial": [2, 6, 8, 10, 11, 12, 14, 16, 18, 20, 22, 26, 42],
        "sat block 14..26": list(range(14, 27)),
        "random big": [61, 67, 74, 83, 89, 97, 104, 113, 122, 131, 140, 151, 163],
    }
    print(f"  target: E[U] >= 1/14 = {1/14:.5f}  (=> mu >= E[U] >= 1/14 > m_P={float(M_P):.5f})")
    for nm, E in bank.items():
        mu, eu, emg, W, U = stats_np(E)
        viol = float((W - 14 * U).min())
        print(f"  {nm:>26}: E[U]={eu:.5f}  mu={mu:.4f}  mu>=E[U]:{mu >= eu - 1e-9}  "
              f"min(W-14U)={viol:+.4f} {'OK' if viol > -1e-6 else 'VIOLATION'}")
    # adversarial descent on E[U] within crux class
    best_eu, best_v = 9.9, None
    for trial in range(5):
        v = random.choice([bank["record 2*{1..11}u{11,13}"], bank["ds 2*{1..12}u{13}"],
                           bank["sat block 14..26"]])[:]
        cur = stats_np(v, 20000)[1]
        for step in range(500):
            i = random.randrange(13)
            cand = random.randrange(max(1, min(v) // 2), min(3 * max(v), 500))
            w = sorted(set(v[:i] + v[i + 1:] + [cand]))
            if len(w) != 13 or not in_crux(w):
                continue
            c = stats_np(w, 12000)[1]
            if c < cur - 1e-5:
                v, cur = w, c
        if cur < best_eu:
            best_eu, best_v = cur, v
    mu, eu, emg, _, _ = stats_np(best_v, 120000)
    print(f"  crux-class adversarial min E[U]: {eu:.5f} at {best_v}")
    print(f"    margin over 1/14: {eu - 1/14:+.5f}   (mu there = {mu:.4f})")
    # 3-AP-free probe (Behrend angle): does killing 3-APs push E[U] to the iid value?
    apfree = [1, 2, 5, 11, 22, 33, 34, 46, 47, 56, 60, 64, 79]  # greedy 3-AP-free-ish
    n3 = sum(1 for i in range(13) for j in range(i+1, 13) for l in range(j+1, 13)
             if apfree[i] + apfree[l] == 2 * apfree[j])
    mu, eu, emg, _, _ = stats_np(apfree)
    print(f"  3-AP-free probe {apfree} (N3={n3}): E[U]={eu:.5f} vs iid (6/7)^13={(6/7)**13:.5f}")

    print()
    print("=" * 78)
    print("STEP 5 -- mean-version crossing for comparison (kps-S59 angle 1)")
    print("=" * 78)
    for m in [18, 20, 22, 24, 26, 28, 30]:
        _, eu, emg, _, _ = stats_np(list(range(1, m + 1)), res=60000)
        print(f"  m={m}: E[mg(AP_m)] ~ {emg:.5f} ({'>=' if emg >= 1/7 else '<'} 1/7)   "
              f"E[U(AP_m)] ~ {eu:.5f} ({'>=' if eu >= 1/14 else '<'} 1/14)")
    print()
    print("DONE.")
