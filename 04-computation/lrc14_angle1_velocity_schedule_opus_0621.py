#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_angle1_velocity_schedule_opus_0621.py  (opus-2026-06-21, ANGLE 1, part 2)

THE VELOCITY-ASSIGNMENT SCHEDULING MODEL (refined from part 1).

KEY STRUCTURAL FACT (proved by part-1 [S3], holds for the FULL-RESIDUE stratum):
  At EVERY resonance a=1..6, every full-residue shape E (|E|=k, residues = all of
  Z/7, with e=0 present) has the IDENTICAL start configuration at y=0:
      start positions {e*a mod 7 : e in E} = Z/7 with residue 0 doubled.
  (Because e -> e*a mod 7 is a bijection of Z/7 for gcd(a,7)=1.)
  So at y=0 EVERY full-residue shape is fully covered in the same way.  Shapes
  differ ONLY in their VELOCITIES v_e = 7e and how those velocities sit on the
  starting points.

So  sum_a W_a  is purely a function of the VELOCITY MULTISET {7e} and the
residue-permutation a acts by.  This turns LAYER 3 into a clean SCHEDULING /
COVERING problem:
   7 sectors on a circle (circumference 7), 8 (=k) moving points starting at the
   integer lattice {0,0,1,2,3,4,5,6}, point with start q has velocity v (=7e).
   The cover survives (locally, both directions) while every sector keeps a point.
   W_a = total measure of {y in [-1/14,1/14] : all 7 sectors occupied}.

THE REFINED HYPOTHESES (part-1 [S2] refuted the naive 'slower is always better';
the gap-filling effect off the full-residue stratum breaks it).  Here we test the
CORRECT, residue-resolved principle ON the full-residue stratum:

 (H-A) "AP velocity = staggered round-robin = max sum".  consec's velocities are
       7*{0,1,2,3,4,5,6,7} -- a perfect AP, residue r carrying speed exactly 7r,
       residue 0 also carrying 49.  Conjecture: consec MAXIMIZES sum_a W_a over
       the full-residue stratum (re-confirm crux, this time inside the velocity
       model so we know the model is faithful).

 (H-B) PER-RESIDUE SLOWEST-REP monotonicity.  Define leg(r) = min |e| over e==r
       (the slowest representative of residue r).  consec has leg = (0,1,2,3,4,5,6)
       -- the componentwise minimum possible for a full-residue set.  Test: is
       sum_a W_a MONOTONE DECREASING under componentwise increase of the leg vector?
       (i.e. if leg(E') >= leg(E) componentwise then sum W_a(E') <= sum W_a(E)?)
       This is the LAYER-1 'min|e| per residue' picture made into a monotonicity.
       Part-1 refuted single-clock speed-up; but that ADDS a fast clock without
       removing the slow one. The LEG vector only sees the SLOWEST per residue, so
       this is a DIFFERENT (coarser) monotonicity. TEST honestly.

 (H-C) THE SECOND-MOMENT of velocities.  consec velocities 7*{0..7} have the
       SMALLEST possible spread among full-residue sets with residue 0 doubled?
       Test sum_a W_a vs Var(velocities) / vs max velocity / vs sum velocities.
       Find which scalar of the velocity multiset best predicts the ranking.

All EXACT Fractions; brutally honest.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

# ---- faithful W_a via the moving-points (velocity) model (verified == brute) ---
def sector_of_point(e, a, y):
    pos = F(e*a) + F(7*e)*y
    return (pos.numerator // pos.denominator) % 7

def covered_all_at(E, a, y):
    return len({sector_of_point(e, a, y) for e in E}) == 7

def W_a(E, a):
    half = F(1, 14); bps = {F(0), -half, half}
    for e in E:
        if e == 0: continue
        lo_val = F(7*e)*(-half) + F(e*a); hi_val = F(7*e)*(half) + F(e*a)
        lo_i = min(lo_val, hi_val); hi_i = max(lo_val, hi_val)
        m = lo_i.numerator // lo_i.denominator
        while m <= hi_i.numerator // hi_i.denominator + 1:
            y = F(m - e*a, 7*e)
            if -half <= y <= half: bps.add(y)
            m += 1
    bps = sorted(bps)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if covered_all_at(E, a, (lo+hi)/2): tot += hi - lo
    return tot

def measS7(E): return sum(W_a(E, a) for a in range(1, 7))

def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def consec(k): return list(range(k))

def leg_vector(E):
    by = defaultdict(list)
    for e in E: by[e % 7].append(abs(e))
    return tuple(min(by[r]) for r in range(7))

def vel_multiset(E): return tuple(sorted(7*abs(e) for e in E))

if __name__ == "__main__":
    print("="*80)
    print("VELOCITY-ASSIGNMENT SCHEDULING MODEL  (full-residue stratum)")
    print("="*80)
    k = 8

    # ----------------- (H-A) re-confirm consec max in the velocity model --------
    print(f"\n[H-A] consec maximizes sum_a W_a over full-residue stratum, span scan:")
    for span in (11, 13, 14, 16):
        shapes = [(0,)+c for c in itertools.combinations(range(1, span+1), k-1)]
        shapes = [E for E in shapes if is_full_residue(E)]
        C = consec(k); msC = measS7(C)
        best = max(shapes, key=measS7); msbest = measS7(best)
        nbeat = sum(1 for E in shapes if measS7(list(E)) > msC)
        tie = sum(1 for E in shapes if measS7(list(E)) == msC)
        print(f"  span<= {span}: {len(shapes):5d} shapes  consec={float(msC):.6f}  "
              f"#strictly-beat-consec={nbeat}  #tie={tie}  argmax={'consec' if best==tuple(C) else best}")

    # ----------------- (H-B) leg-vector monotonicity ----------------------------
    print(f"\n[H-B] LEG-VECTOR monotonicity (leg(r)=min|e| over e==r mod 7):")
    print(f"      Does leg(E') >= leg(E) componentwise  =>  sum W_a(E') <= sum W_a(E)?")
    span = 14
    shapes = [(0,)+c for c in itertools.combinations(range(1, span+1), k-1)]
    shapes = [E for E in shapes if is_full_residue(E)]
    data = [(E, leg_vector(E), measS7(list(E))) for E in shapes]
    # check all comparable pairs
    viol = 0; checks = 0; worst = F(0); worst_ex = None
    for i in range(len(data)):
        Ei, li, mi = data[i]
        for j in range(len(data)):
            if i == j: continue
            Ej, lj, mj = data[j]
            # lj >= li componentwise (E_j has slower-or-equal legs everywhere)
            if all(lj[r] >= li[r] for r in range(7)) and lj != li:
                checks += 1
                # expect mj <= mi  (slower legs => less? NO -- slower legs = SMALLER
                # min|e| means FASTER actually. min|e| small = slow drift = long cover.
                # leg = min|e|; SMALL leg => SLOW => LONG survival => bigger W.
                # so leg(E') >= leg(E) (bigger min, faster) => expect W smaller.
                if mj > mi:
                    viol += 1
                    d = mj - mi
                    if d > worst: worst = d; worst_ex = (Ei, li, float(mi), Ej, lj, float(mj))
    print(f"      comparable pairs checked={checks}  violations={viol}")
    if worst_ex:
        Ei, li, mi, Ej, lj, mj = worst_ex
        print(f"      worst: leg{li} (W={mi:.6f}) dominated-by leg{lj} but W={mj:.6f} BIGGER")
        print(f"             E_small_leg={Ei}  E_big_leg={Ej}")
    else:
        print(f"      => LEG-VECTOR MONOTONICITY HOLDS: componentwise-smaller leg (slower")
        print(f"         per-residue reps) => larger sum W_a.  consec has the global")
        print(f"         componentwise-minimum leg (0,1,2,3,4,5,6) => it is the MAXIMUM.")
        print(f"         (This would PROVE LAYER 3 on the full-residue stratum.)")

    # confirm consec's leg is the global minimum
    Clv = leg_vector(consec(k))
    print(f"\n      consec leg vector = {Clv}")
    allge = all(all(lv[r] >= Clv[r] for r in range(7)) for _, lv, _ in data)
    print(f"      every full-residue shape has leg >= consec leg componentwise? {allge}")

    # ----------------- (H-C) scalar predictors ----------------------------------
    print(f"\n[H-C] scalar velocity predictors of the ranking (corr with sum W_a):")
    import statistics
    Ws = [float(m) for _,_,m in data]
    def feat(E, fn): return float(fn(E))
    feats = {
        "sum|e|": lambda E: sum(abs(e) for e in E),
        "max|e|": lambda E: max(abs(e) for e in E),
        "sum e^2": lambda E: sum(e*e for e in E),
        "var(e)": lambda E: statistics.pvariance([abs(e) for e in E]),
        "sum leg": lambda E: sum(leg_vector(E)),
        "max leg": lambda E: max(leg_vector(E)),
        "sum leg^2": lambda E: sum(l*l for l in leg_vector(E)),
    }
    for name, fn in feats.items():
        xs = [feat(list(E), fn) for E,_,_ in data]
        try:
            r = statistics.correlation(xs, Ws)
        except Exception:
            r = float('nan')
        # also: does the argmin of this feature == consec?
        amin = min(data, key=lambda d: feat(list(d[0]), fn))
        print(f"  {name:10s}: corr(feature, W)={r:+.4f}   argmin-feature is consec? "
              f"{amin[0]==tuple(consec(k))}")
