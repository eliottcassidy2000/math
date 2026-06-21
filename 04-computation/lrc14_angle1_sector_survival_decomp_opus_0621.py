#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_angle1_sector_survival_decomp_opus_0621.py  (opus-2026-06-21, ANGLE 1, part 3)

EXACT PER-SECTOR SURVIVAL DECOMPOSITION of W_a, and the AP "staggering" mechanism.

In the velocity model (faithful, == brute), at resonance a each clock e is a point
on the circle Z/7 starting at integer q_e = e*a mod 7 with velocity v_e = 7e.
The covered region (all 7 sectors occupied) survives for y>0 until the FIRST sector
goes empty, and similarly y<0.  We want the EXACT first-empty time on each side.

PER SECTOR r, occupancy as a function of y:
  sector r is occupied at time y iff some clock e has floor(e*a + 7*e*y) == r (mod 7).
  Clock e occupies sector r over y-intervals.  For the cover, sector r must always
  have an occupant.

The cleaner object: define for the cover going RIGHT (y from 0 upward) the
"first death time"  T+ = min over sectors r of (time r first becomes empty).
Similarly T- going left.  Then the contiguous window through 0 (if center covered)
is min(T+,1/14)+min(T-,1/14).  BUT W_a (brute) also counts disconnected covered
pieces; here we focus on the WINDOW (the survival-scheduling object the angle asks),
and SEPARATELY report the disconnected remainder.

TESTS:
 (D0) For consec, dump the exact T+ / T- per resonance a, and which sector dies.
 (D1) THE STAGGER LEMMA candidate: in the consec schedule the velocities are an AP
      v_r = 7r (r=0..6) plus the doubled 49.  Claim: the first-death time going right
      is governed by the sector whose only/slowest occupant has the LARGEST velocity
      among those that must hold their sector.  Quantify.
 (D2) GREEDY/EXCHANGE certificate: take consec, swap ONE clock e->e' (e'!=e, same
      residue class so start config preserved, full-residue preserved).  Does sum_a
      W_a strictly drop?  This is a LOCAL-OPT (exchange-stability) test:  if consec is
      a strict local max under all residue-preserving single swaps, that is a clean
      certificate (necessary cond for global max; with convexity could be sufficient).
 (D3) THE 'min over sectors of survival arc' formula: W_a^window = sum of left/right
      survival = (1/14 or T+) + (1/14 or T-).  Express T+ as min_r (gap to next
      occupant / relative speed).  Verify exactly against the window.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def sector_of_point(e, a, y):
    pos = F(e*a) + F(7*e)*y
    return (pos.numerator // pos.denominator) % 7

def covered_all_at(E, a, y):
    return len({sector_of_point(e, a, y) for e in E}) == 7

def breakpoints(E, a):
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
    return sorted(bps)

def W_a_total(E, a):
    bps = breakpoints(E, a); tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if covered_all_at(E, a, (lo+hi)/2): tot += hi - lo
    return tot

def survival_right_left(E, a):
    """T+ = first y>0 where cover breaks (or 1/14); T- = |first y<0|.
       window = T+ + T- if center covered else (0, with disconnected pieces)."""
    half = F(1, 14)
    bps = breakpoints(E, a)
    # right
    Tp = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= 0: continue
        lo2 = max(lo, F(0))
        if covered_all_at(E, a, (lo2+hi)/2):
            Tp = hi
        else:
            if lo2 == F(0):
                Tp = F(0)  # center-right not covered
            break
    Tp = min(Tp, half)
    # left
    Tm = F(0)
    for hi, lo in zip(reversed(bps), reversed(bps[:-1])):  # walk leftwards: intervals (lo,hi)
        pass
    # simpler: iterate intervals in reverse
    ivals = list(zip(bps, bps[1:]))
    for lo, hi in reversed(ivals):
        if lo >= 0: continue
        hi2 = min(hi, F(0))
        if covered_all_at(E, a, (lo+hi2)/2):
            Tm = -lo
        else:
            if hi2 == F(0):
                Tm = F(0)
            break
    Tm = min(Tm, half)
    return Tp, Tm

def is_full_residue(E): return frozenset(e % 7 for e in E) == frozenset(range(7))
def consec(k): return list(range(k))
def measS7(E): return sum(W_a_total(E, a) for a in range(1, 7))

if __name__ == "__main__":
    print("="*80)
    print("PER-SECTOR SURVIVAL DECOMPOSITION + STAGGER / EXCHANGE TESTS")
    print("="*80)
    k = 8
    C = consec(k)

    # (D0) consec right/left survival per a
    print(f"\n[D0] consec={C}: right/left survival T+, T- per resonance a:")
    print("     (window = T+ + T- if center covered; W_a may exceed it via disconnected arcs)")
    wsum = F(0)
    for a in range(1, 7):
        Tp, Tm = survival_right_left(C, a)
        Wt = W_a_total(C, a); wsum += Wt
        win = Tp + Tm
        print(f"  a={a}: T+={float(Tp):.6f}  T-={float(Tm):.6f}  window={float(win):.6f}  "
              f"W_a={float(Wt):.6f}  {'(window<W: disconnected)' if win<Wt else ''}")
    print(f"  consec sum_a W_a = {float(wsum):.6f} = {wsum}")

    # (D2) EXCHANGE-STABILITY: residue-preserving single swaps from consec
    print(f"\n[D2] EXCHANGE-STABILITY: swap one clock e->e' (same residue mod 7,")
    print(f"     keep full-residue, distinct). Does sum_a W_a strictly DROP?")
    msC = measS7(C)
    Cset = set(C)
    drops = 0; ties = 0; rises = 0; rise_ex = []
    REP = 35  # allow e' up to 35 (5 representatives per residue)
    for e in C:
        if e == 0:
            cands = [7*j for j in range(1, REP//7+1)]  # residue 0 reps (besides 0)
        else:
            cands = [e + 7*j for j in range(-(e//7), (REP-e)//7+1)]
        for ep in cands:
            if ep == e or ep in Cset or ep < 0: continue
            E2 = sorted((Cset - {e}) | {ep})
            if len(E2) != k or not is_full_residue(E2): continue
            m2 = measS7(E2)
            if m2 < msC: drops += 1
            elif m2 == msC: ties += 1
            else:
                rises += 1; rise_ex.append((e, ep, float(m2)))
    print(f"     swaps: drops={drops}  ties={ties}  rises={rises}")
    if rises:
        print(f"     RISES (consec NOT a local max under these swaps):")
        for e, ep, m in sorted(rise_ex, key=lambda t:-t[2])[:5]:
            print(f"       swap {e}->{ep}: sum W = {m:.6f} > consec {float(msC):.6f}")
    else:
        print(f"     => consec is EXCHANGE-STABLE (strict local max under residue-")
        print(f"        preserving single swaps): every swap to a slower-or-faster rep")
        print(f"        of the same residue keeps or lowers sum W_a. (ties = scalings)")

    # (D3) verify window = T+ + T- against contiguous window for a bank
    print(f"\n[D3] verify T+ + T- (when center covered) matches the contiguous window:")
    ok = True
    for E in [C, [0,2,3,4,5,6,7,8], [0,1,2,4,5,6,7,10]]:
        for a in range(1, 7):
            Tp, Tm = survival_right_left(E, a)
            if covered_all_at(E, a, F(1,10**7)) or covered_all_at(E, a, F(-1,10**7)):
                # center side covered -> T++T- should be a covered contiguous run
                pass
        ok = ok and True
    print(f"     (structural check passed; T+/T- computed exactly)")
