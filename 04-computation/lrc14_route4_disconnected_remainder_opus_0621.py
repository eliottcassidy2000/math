#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_route4_disconnected_remainder_opus_0621.py  (opus-2026-06-21, ROUTE 4)

LAYER 3 via the DISCONNECTED-ARC REMAINDER.

  measS7(E) = WIN(E) + DISC(E)
    WIN  = sum_{a=1..6} (T+(a) + T-(a))   contiguous-through-center survival.
           consec UNIQUELY maximizes WIN with margin (HYP-2761, 0 ties at k=8,9,10).
    DISC = measS7 - WIN >= 0             the off-center covered arcs (disconnected pieces).

ROUTE 4 question: since consec wins WIN by a margin, LAYER 3 (consec maximizes measS7)
follows IF the disconnected remainder cannot overturn consec's WIN lead. For a rival E:

  consec total lead = measS7(consec) - measS7(E)
                    = [WIN(consec) - WIN(E)]  +  [DISC(consec) - DISC(E)]
                    =      WINlead(E)         -      DISCdeficit(E)
  where WINlead(E)    = WIN(consec) - WIN(E)   > 0   (consec wins WIN)
        DISCdeficit(E)= DISC(E) - DISC(consec)        (how much MORE disc E has)

  consec wins measS7  <=>  WINlead(E) > DISCdeficit(E)   for every rival E.
  (ties allowed: measS7 has ties, so >= and uniqueness via WIN's strictness on ties)

TESTS over the full-residue stratum (all E with 0 in E, full Z/7 residue, span<=SPAN):
 (T0) baseline consec WIN, DISC, measS7 for k=8,9,10.
 (T1) for every rival, compute WINlead and DISCdeficit; is WINlead > DISCdeficit always?
      Report the MINIMUM SLACK = min over rivals of (WINlead - DISCdeficit).
      If slack > 0 for all, ROUTE 4 closes LAYER 3 (subject to span).
 (T2) is DISC bounded by the WIN margin?  Report max DISCdeficit vs WIN margin spread.
 (T3) does consec also maximize DISC?  (would make the bound trivial)
      Report #rivals with DISC(E) > DISC(consec).
 (T4) the WORST rivals: list shapes where (WINlead - DISCdeficit) is smallest.
"""
import sys, itertools
from fractions import Fraction as F
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

def window_TpTm(E, a):
    half = F(1, 14); bps = breakpoints(E, a); ivals = list(zip(bps, bps[1:])); Tp = F(0)
    for lo, hi in ivals:
        if hi <= 0: continue
        lo2 = max(lo, F(0))
        if covered_all_at(E, a, (lo2+hi)/2): Tp = hi
        else:
            if lo2 == F(0): Tp = F(0)
            break
    Tp = min(Tp, half); Tm = F(0)
    for lo, hi in reversed(ivals):
        if lo >= 0: continue
        hi2 = min(hi, F(0))
        if covered_all_at(E, a, (lo+hi2)/2): Tm = -lo
        else:
            if hi2 == F(0): Tm = F(0)
            break
    Tm = min(Tm, half); return Tp, Tm

def is_full_residue(E): return frozenset(e % 7 for e in E) == frozenset(range(7))
def consec(k): return list(range(k))
def measS7(E): return sum(W_a_total(E, a) for a in range(1, 7))
def WIN(E): return sum(sum(window_TpTm(E, a)) for a in range(1, 7))
def DISC(E): return measS7(E) - WIN(E)

def stratum(k, span):
    """All E = {0} ∪ (k-1)-subset of {1..span}, full Z/7 residue."""
    for c in itertools.combinations(range(1, span+1), k-1):
        E = [0]+list(c)
        if is_full_residue(E):
            yield E

if __name__ == "__main__":
    SPANS = {8: 16, 9: 16, 10: 15}
    for k in [8, 9, 10]:
        span = SPANS[k]
        C = consec(k)
        mC, wC, dC = measS7(C), WIN(C), DISC(C)
        print("="*84)
        print(f"k={k}  span<={span}")
        print(f"  consec measS7 = {mC} = {float(mC):.6f}")
        print(f"  consec WIN    = {wC} = {float(wC):.6f}")
        print(f"  consec DISC   = {dC} = {float(dC):.6f}   (DISC/measS7 = {float(dC/mC):.3f})")

        rivals = [E for E in stratum(k, span) if E != C]
        n = len(rivals)
        # accumulate
        min_slack = None; min_slack_E = None
        max_discdef = F(-10); max_discdef_E = None
        min_winlead = F(10); min_winlead_E = None
        disc_gt_consec = 0
        meas_ge_consec = 0  # rivals matching or beating consec measS7
        meas_gt_consec = 0
        worst = []  # (slack, E, winlead, discdef)
        for E in rivals:
            mE, wE, dE = measS7(E), WIN(E), DISC(E)
            winlead = wC - wE         # > 0 expected
            discdef = dE - dC         # how much MORE disc E has
            slack = winlead - discdef # = mC - mE
            assert slack == mC - mE
            if min_slack is None or slack < min_slack:
                min_slack, min_slack_E = slack, E
            if discdef > max_discdef:
                max_discdef, max_discdef_E = discdef, E
            if winlead < min_winlead:
                min_winlead, min_winlead_E = winlead, E
            if dE > dC: disc_gt_consec += 1
            if mE >= mC: meas_ge_consec += 1
            if mE > mC: meas_gt_consec += 1
            worst.append((slack, tuple(E), winlead, discdef, mE, wE, dE))

        print(f"  rivals in stratum: {n}")
        print(f"  [T1] min SLACK (= min consec measS7 lead) = {min_slack} = {float(min_slack):.6f}")
        print(f"       attained at {min_slack_E}")
        print(f"       => consec measS7 lead is {'>=0 (consec at least ties)' if min_slack>=0 else 'NEGATIVE -- consec NOT max!'}")
        print(f"  [T1b] #rivals strictly beating consec measS7 = {meas_gt_consec}")
        print(f"        #rivals tying-or-beating consec measS7 = {meas_ge_consec}")
        print(f"  [T2] max DISCdeficit = {max_discdef} = {float(max_discdef):.6f}  at {max_discdef_E}")
        print(f"       min WINlead     = {min_winlead} = {float(min_winlead):.6f}  at {min_winlead_E}")
        print(f"       (ROUTE 4 closes if max DISCdeficit < min WINlead, i.e. uniform bound)")
        uniform_ok = max_discdef < min_winlead
        print(f"       uniform DISC-bound (maxDISCdef < minWINlead)? {uniform_ok}")
        print(f"  [T3] #rivals with DISC(E) > DISC(consec) = {disc_gt_consec}/{n}  "
              f"(if 0, consec also maximizes DISC -> bound trivial)")
        # [T4] worst 8 by slack
        worst.sort(key=lambda t: t[0])
        print(f"  [T4] tightest 8 rivals (smallest measS7 lead = WINlead - DISCdef):")
        for slack, E, wl, dd, mE, wE, dE in worst[:8]:
            print(f"     {str(list(E)):30s} slack={float(slack):+.6f}  "
                  f"WINlead={float(wl):+.5f} DISCdef={float(dd):+.5f}  "
                  f"measS7={float(mE):.5f}")
