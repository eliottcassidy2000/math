#!/usr/bin/env python3
"""
lrc_harmonic_gate_kps_S2.py -- HYP-4099 (kind-pasteur-2026-07-05-S2)

THE RATIONAL-POINT MARGIN CERTIFICATE + THE RATIO-GATE REDUCTION of
TightLooseDichotomy (klein HYP-4096 surface).

The atom (to be formalized in LRCHarmonicGate.lean):

  rational_point_margin: for integers s > 0, k, mu > 0 and any speed v with
      mu <= (v*k) mod s <= s - mu
  the point t = k/s has |v*t - m| >= mu/s for EVERY integer m.
  (Proof: |v*k - m*s| >= min(r, s-r) >= mu where r = (v*k) mod s.)

Every loose-branch witness of the dichotomy is an instance.  Two canonical
sub-grids for a base with min a, max b, s = a+b:
  k=1 (band_margin_reciprocal, kps-S1): all speeds in [a,b] qualify with mu=a.
      Margin a/s >= 2/25  <=>  2b <= 23a  (ratio <= 11.5).  THE RATIO GATE.
  k=2 (second harmonic): speeds in [a,b] AVOIDING the middle window
      (s-2a)/2 < w < (s+2a)/2 qualify with mu = 2a.  Margin 2a/s >= 2/25
      <=>  b <= 24a.  Covers the second-value witness {1..11,24} (s=25, t=2/25).

This script verifies:
  (V1) the certificate arithmetic against brute distance, random instances;
  (V2) {1..11,24}: k=2 certificate holds, margin exactly 2/25;
  (V3) the ratio gate NEVER lies: every primitive 12-set with 2b <= 23a has a
       t with min margin >= 2/25 (the reciprocal point itself);
  (V4) coverage statistics over all primitive 12-subsets of [1,N]: how many are
       closed loose by k=1, by k<=2, by k<=K harmonics of s=a+b (at margin 2/25),
       how many are tight-side (dilated-AP subsets), and what survives both
       (the hard residue of hdich after the gate).
  (V5) the surviving hard bases listed explicitly (for the next session's attack).

Owner mandate: LRC(<=13) is settled (citation); this is about the n=12 base
rigidity dichotomy inside the LRC(14) endgame.
"""

from math import gcd
from fractions import Fraction
from itertools import combinations
import random

# ---------------------------------------------------------------- the atom
def cert_margin_ok(v, k, s, mu):
    """The integer-mod certificate condition for one speed."""
    r = (v * k) % s
    return mu <= r <= s - mu

def cert_family_ok(ws, k, s, mu):
    return all(cert_margin_ok(w, k, s, mu) for w in ws)

def exact_min_dist_at(ws, t):
    """min over family of distance from w*t to nearest integer (Fraction t)."""
    best = None
    for w in ws:
        x = (w * t) % 1
        d = min(x, 1 - x)
        if best is None or d < best:
            best = d
    return best

# ---------------------------------------------------------------- V1
def verify_atom(trials=20000, seed=1):
    rng = random.Random(seed)
    bad = 0
    for _ in range(trials):
        s = rng.randint(2, 400)
        k = rng.randint(-3 * s, 3 * s)
        mu = rng.randint(1, s // 2 if s // 2 >= 1 else 1)
        v = rng.randint(-10**6, 10**6)
        if v == 0:
            continue
        if cert_margin_ok(v, k, s, mu):
            t = Fraction(k, s)
            x = (v * t) % 1
            d = min(x, 1 - x)
            if d < Fraction(mu, s):
                bad += 1
                print(f"  ATOM VIOLATION: v={v} k={k} s={s} mu={mu} dist={d}")
    print(f"[V1] rational_point_margin atom: {trials} random instances, "
          f"{bad} violations")
    return bad == 0

# ---------------------------------------------------------------- V2
def verify_second_value():
    ws = list(range(1, 12)) + [24]
    a, b = min(ws), max(ws)
    s = a + b            # 25
    mu = 2 * a           # 2
    ok = cert_family_ok(ws, 2, s, mu)
    t = Fraction(2, s)
    d = exact_min_dist_at(ws, t)
    print(f"[V2] {{1..11,24}}: k=2 certificate (s={s}, mu={mu}) holds: {ok}; "
          f"exact margin at t=2/25: {d} (= 2/25: {d == Fraction(2,25)})")
    # hole condition check: no speed strictly inside ((s-2a)/2, (s+2a)/2) = (11.5, 13.5)
    hole = [w for w in ws if 2 * w > s - 2 * a and 2 * w < s + 2 * a]
    print(f"     middle-window occupants (must be none): {hole}")
    return ok and d == Fraction(2, 25) and not hole

# ---------------------------------------------------------------- V3 + V4 + V5
def is_primitive(ws):
    g = 0
    for w in ws:
        g = gcd(g, w)
    return g == 1

def tight_side(ws):
    """Dilated-AP-subset test: all values in c*{1..12} for some c >= 2, OR
    the undilated {1..12} subset case c=1 -- the dichotomy's tight branch is
    stated with c >= 2 (c=1 undilated is subsumable by margin analysis /
    free-rider handles c>=2; klein surface takes 2 <= c).  We report both."""
    a = min(ws)
    for c in range(1, a + 1):
        if all(w % c == 0 and 1 <= w // c <= 12 for w in ws):
            return c
    return 0

def harmonics_close(ws, K, beta=Fraction(2, 25)):
    """Does some harmonic t=k/s of s=a+b certify margin >= beta for the base?"""
    a, b = min(ws), max(ws)
    s = a + b
    need = beta * s  # mu >= beta*s, mu integer: mu = ceil(beta*s)
    mu = int(need) if need == int(need) else int(need) + 1
    if mu < 1:
        mu = 1
    for k in range(1, K + 1):
        if cert_family_ok(ws, k, s, mu):
            return k
    return 0

def sweep(N=16, K=12, beta=Fraction(2, 25)):
    """All primitive 12-subsets of [1,N]."""
    total = 0
    gate1 = 0      # k=1 (ratio gate 2b <= 23a)
    gate2 = 0      # k<=2
    gateK = 0      # k<=K
    tight = 0      # dilated AP subset (c>=1)
    both = 0       # neither harmonic-closed (k<=K) nor tight
    survivors = []
    gate1_lies = 0
    for ws in combinations(range(1, N + 1), 12):
        if not is_primitive(ws):
            continue
        total += 1
        a, b = min(ws), max(ws)
        # V3: the ratio gate must never lie
        if 2 * b <= 23 * a:
            # margin at reciprocal point must be >= 2/25
            d = exact_min_dist_at(ws, Fraction(1, a + b))
            if d < beta:
                gate1_lies += 1
                print(f"  GATE LIE at {ws}: reciprocal margin {d}")
        kfound = harmonics_close(ws, K, beta)
        c = tight_side(ws)
        if kfound == 1:
            gate1 += 1
        if 1 <= kfound <= 2:
            gate2 += 1
        if kfound:
            gateK += 1
        if c:
            tight += 1
        if not kfound and not c:
            both += 1
            survivors.append(ws)
    print(f"[V3] ratio-gate lies: {gate1_lies} (must be 0)")
    print(f"[V4] primitive 12-subsets of [1,{N}]: {total}")
    print(f"     closed by k=1 harmonic (ratio gate at 2/25):  {gate1}"
          f"  ({100*gate1/total:.1f}%)")
    print(f"     closed by k<=2 harmonics:                     {gate2}"
          f"  ({100*gate2/total:.1f}%)")
    print(f"     closed by k<={K} harmonics of s=a+b:           {gateK}"
          f"  ({100*gateK/total:.1f}%)")
    print(f"     tight-side (subset of c*{{1..12}}, c>=1):       {tight}")
    print(f"     SURVIVORS (neither):                          {both}")
    return survivors

def survivors_analysis(survivors, beta=Fraction(2, 25)):
    """[V5] For each survivor: is it actually loose (M >= 2/25) via a general
    small-denominator rational point (not restricted to s=a+b)?"""
    print(f"[V5] survivor analysis (general rational points q <= 60):")
    hard = []
    for ws in survivors:
        found = None
        for q in range(2, 61):
            for k in range(1, q):
                if gcd(k, q) > 1:
                    continue
                d = exact_min_dist_at(ws, Fraction(k, q))
                if d >= beta:
                    found = (k, q, d)
                    break
            if found:
                break
        if found:
            k, q, d = found
            print(f"     {ws}: loose at t={k}/{q} (margin {d}) "
                  f"-- general grid, NOT the s=a+b harmonic")
        else:
            hard.append(ws)
            print(f"     {ws}: NO rational witness q<=60 at margin {beta} -- HARD")
    print(f"     hard survivors: {len(hard)}")
    return hard

if __name__ == "__main__":
    print("=" * 74)
    print("HYP-4099: rational-point margin certificate + ratio-gate reduction")
    print("=" * 74)
    ok1 = verify_atom()
    ok2 = verify_second_value()
    survivors = sweep(N=16, K=12)
    hard = survivors_analysis(survivors)
    print("=" * 74)
    print(f"SUMMARY: atom {'OK' if ok1 else 'FAIL'}; second-value {'OK' if ok2 else 'FAIL'}; "
          f"survivors {len(survivors)}, hard {len(hard)}")
