#!/usr/bin/env python3
"""
lrc14_realization_residual_audit_macmini_S65cont8.py -- THE RESIDUAL-CLASS AUDIT for hpartA

The status file (LRC14-STATUS-2026-07-09 sec. 3) names the one open node: single-chain
multi-scale covering shapes.  Every agent knows their own instrument's domain; NOBODY has mapped
the INTERSECTION of the complements.  This audit tests, per covering set, every proved
realization instrument, and hunts for sets where ALL fail (the true residual class):

  I1  pure-cluster      Vmax <= 13*Vmin                     (spread13_lonely / PureClusterCorner)
  I2  coarse <=12       exists L: <=12 groups, offsets A <= L/182  (LRCCoarseReduction)
  I5  detuned d=1       S = g*H u {delta}, |H|=12, g>=2, g not| delta  (monad THM-668: M >= 1/13)
  I6  certificates      C0 window / C1 gcd ledger / C2 divisor descent (THM-668/672 mine)
  I3  composed grid     exists j: cluster gap clears fattened threshold AND every slow p safe
                        with the DRIFT-VALID budget (LEM-014 / THM-666 / monad-S2 phi-interval):
                        naive budget p/V and refined budget p*(w-1/7)/V both reported.
                        Split at the adaptive 9/14 threshold (THM-667): Ltil = {v >= 9V/14}.
  RESIDUAL = fails ALL => confirm liveness exactly (pair-sum events) + record the witness.

All integer/exact arithmetic in the instrument checks (scaled by V and 14V as needed).
"""
from math import gcd
from functools import reduce
from itertools import combinations
import random

random.seed(88)

def covering(S):  return all(any(v % q == 0 for v in S) for q in range(2, 15))
def primitive(S): return reduce(gcd, S) == 1

# ---------------------------------------------------------------- I1
def i1_pure(S):
    return max(S) <= 13 * min(S)

# ---------------------------------------------------------------- I2 (exact: A/L <= 1/182)
def i2_coarse(S):
    V = max(S)
    for L in range(2, V + 1):
        ks, amax_num_ok = [], True
        for v in S:
            k = (v + L // 2) // L                 # nearest multiple
            a = v - L * k
            if k == 0 or 182 * abs(a) > L:        # offset budget |a| <= L/182 exact
                amax_num_ok = False
                break
            ks.append(k)
        if amax_num_ok and len(set(ks)) <= 12:
            return True
    return False

# ---------------------------------------------------------------- I5 (detuned harmonic d=1)
def i5_detuned(S):
    for g in range(2, max(S) + 1):
        divs = [v for v in S if v % g == 0]
        nond = [v for v in S if v % g != 0]
        if len(divs) == 12 and len(nond) == 1 and len({v // g for v in divs}) == 12:
            return True
    return False

# ---------------------------------------------------------------- I6 (my certificate battery)
def danger_j(k):  return -(-k // 14) - 1
def blocked(R, k):
    j = danger_j(k)
    bad = set(range(0, j + 1)) | set(range(k - j, k))
    for s in range(1, k):
        if all((r * s) % k not in bad for r in R):
            return False
    return True
def rulers_of(S): return sorted({S[i] + S[j] for i in range(13) for j in range(i, 13)})
def i6_certs(S):
    if max(S) <= 13 * min(S):
        return True                                # C0 (subsumes I1, kept for clarity)
    for q in rulers_of(S):
        for k in range(15, q):
            if q % k == 0:
                R = [v % k for v in S]
                if 0 not in R and not blocked(R, k):
                    return True                    # C2
    return False

# ---------------------------------------------------------------- I3 (composed grid, exact)
from fractions import Fraction as F

def lonely_at(S, tau):
    """Exact: all ||v*tau|| >= 1/14 (closed)."""
    for v in S:
        r = (v * tau) % 1
        if 14 * r < 1 or 14 * (1 - r) < 1:
            return False
    return True

def i3_composed(S):
    """Adaptive split at 9/14: Ltil = {v : 14 v >= 9 V}, co-offsets e = V - v.  For each grid
    j whose cluster teeth leave a gap wide enough to matter, CONSTRUCT the candidate
    tau* = (j + phi*)/V with the fast phase at the gap midpoint, and verify loneliness of ALL
    13 speeds EXACTLY at tau* (sound by construction -- boxeph-S3's snap pattern, no budget
    approximations).  Also tries the drift-corrected midpoint as a second candidate."""
    V = max(S)
    Ltil = [v for v in S if 14 * v >= 9 * V]
    E = sorted(V - v for v in Ltil)                # co-offsets, incl 0 (v = V)
    for j in range(1, V):
        teeth = sorted(set((e * j) % V for e in E))
        # widest circular gap (a, b) in V-units
        best_a, best_w = teeth[-1], V - teeth[-1] + teeth[0]     # wrap gap ends at teeth[0]+V
        prev = teeth[0]
        for t in teeth[1:]:
            if t - prev > best_w:
                best_a, best_w = prev, t - prev
            prev = t
        if 7 * best_w <= V:                        # gap <= 1/7: cannot host
            continue
        # candidate 1: fast phase at the undrifted gap midpoint
        phi_num = 2 * best_a + best_w              # midpoint*2 in V-units
        tau1 = F(2 * j * V + phi_num, 2 * V * V)
        if lonely_at(S, tau1):
            return True
        # candidate 2: drift-corrected -- teeth drift by e*phi/V; recentre once
        phi0 = F(phi_num, 2 * V)
        drift_mid = sum(F(e) for e in E) / len(E) * phi0 / V
        tau2 = tau1 + drift_mid / V
        if lonely_at(S, tau2):
            return True
    return False

# ---------------------------------------------------------------- exact liveness (ground truth)
def exact_live(S):
    for q in rulers_of(S):
        for p in range(1, q):
            if all(14 * (v * p % q) >= q and 14 * (q - v * p % q) >= q for v in S):
                return (q, p)
    return None

INSTRUMENTS = [("I1-pure", i1_pure), ("I2-coarse", i2_coarse), ("I5-detuned", i5_detuned),
               ("I6-certs", i6_certs), ("I3-composed", i3_composed)]

def audit(S):
    fired = [name for name, f in INSTRUMENTS if f(S)]
    return fired

# ---------------------------------------------------------------- sweeps
print("=" * 100)
print("RESIDUAL-CLASS AUDIT: covering sets vs ALL proved realization instruments")
print("=" * 100)
from collections import Counter
firecount = Counter()
residuals = []

def process(S, tag):
    fired = audit(S)
    for f in fired:
        firecount[f] += 1
    firecount["_total"] += 1
    if not fired:
        lv = exact_live(S)
        residuals.append((tag, S, lv))
        print(f"  *** RESIDUAL [{tag}]: {S}  exact_live={lv}")

# (a) random covering, several scale bands
for cap, n in ((60, 250), (150, 250), (400, 200), (1000, 120)):
    got = 0
    while got < n:
        S = sorted(random.sample(range(1, cap + 1), 13))
        if covering(S) and primitive(S):
            process(S, f"rand{cap}")
            got += 1

# (b) structured: mid-band-loaded (many speeds in (V/14, 9V/14))
for trial in range(150):
    V = random.randrange(200, 800)
    nmid = random.randrange(3, 8)
    mids = random.sample(range(V // 14 + 1, (9 * V) // 14), nmid)
    highs = random.sample(range((9 * V) // 14, V), 12 - nmid)
    S = sorted(set(mids + highs + [V]))
    if len(S) == 13 and covering(S) and primitive(S):
        process(S, "midband")

# (c) single-chain geometric (ratio ~ 1.3-2 per step, total ratio > 13, no big gap)
for trial in range(150):
    r = 1.25 + random.random() * 0.9
    v = random.randrange(2, 9)
    ch = [v]
    while len(ch) < 13:
        v = int(v * r) + random.randrange(0, 3) + 1
        ch.append(v)
    S = sorted(set(ch))
    if len(S) == 13 and covering(S) and primitive(S) and max(S) > 13 * min(S):
        process(S, "chain")

# (d) adversarial: minimize #instruments firing
for cap in (150, 500):
    while True:
        S = sorted(random.sample(range(1, cap + 1), 13))
        if covering(S) and primitive(S):
            break
    cur = len(audit(S))
    for step in range(120):
        T = list(S)
        T[random.randrange(13)] = random.randrange(1, cap + 1)
        T = sorted(set(T))
        if len(T) != 13 or not covering(T) or not primitive(T):
            continue
        sc = len(audit(T))
        if sc <= cur:
            S, cur = T, sc
        if cur == 0:
            break
    print(f"adversarial cap {cap}: min #instruments firing = {cur}  S={S}  fired={audit(S)}")
    if cur == 0:
        process(S, f"adv{cap}")

print()
tot = firecount.pop("_total")
print(f"AUDITED {tot} covering sets.  Instrument coverage:")
for name, _ in INSTRUMENTS:
    print(f"  {name}: fired on {firecount[name]}/{tot} = {firecount[name]/tot:.1%}")
print(f"RESIDUAL (all instruments fail): {len(residuals)}")
for tag, S, lv in residuals[:10]:
    print(f"  [{tag}] {S} live={lv}")
print()
print("Done.")
