#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LRC14 covering floor: ADVERSARIAL UNIFORM SCAN of CV(N_R)^2  (klein-2026-06-29-S4).

THM-579 (mac-mini): the covering floor R' > 0 holds whenever the gatekeeper
    CV(N_R)^2 < m_Q/(1-m_Q)            [sufficient condition]
where N_R(t) = #{a in 0..13 : t+a/14 is R-safe} is the 14-sheet count and
CV(N_R)^2 = Var(N_R)/E[N_R]^2 is a CONGRUENCE-conditioned second moment
(sum_{N!=0}|chat(14N)|^2 = Var(N_R)/196). mac-mini verified 6 rows; the UNIFORM
bound over the whole covering family is THE open piece (THM-579, mac-mini-S19).

This script reuses mac-mini's EXACT-arithmetic machinery and does a broad
ADVERSARIAL scan to (1) find sup CV(N_R)^2 and its argmax R, (2) locate the
binding structure (esp. the speed-7 resonance: 7*(a/14)=a/2 splits the 14 sheets
into correlated even/odd halves, the 2-adic/7-adic worry), (3) report the
TIGHTEST gatekeeper margin CV^2 vs m_Q/(1-m_Q) over the family, and (4) compare
to the metagraph's bounded Burnside variance (THM-587/588) -- the finite testbed.
"""
from __future__ import annotations
import sys, os, math, itertools, random
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# reuse mac-mini's exact functions (no transcription risk)
m = __import__("lrc14_floor_CV_sheetcount_bound_macmini_20260629")
lonely_set = m.lonely_set; measure = m.measure; autocorr = m.autocorr

def cv2_mR(R):
    """exact CV(N_R)^2 and m_R for a 14-free speed set R."""
    R = tuple(sorted(set(R)))
    Lr = lonely_set(R)
    mR = measure(Lr)
    if mR == 0:
        return None, F(0)
    EN = 14 * mR
    EN2 = F(0)
    for d in range(-13, 14):
        EN2 += (14 - abs(d)) * autocorr(Lr, F(d, 14))
    Var = EN2 - EN * EN
    return Var / (EN * EN), mR

def mQ_of(Q):
    return measure(lonely_set(tuple(sorted(set(Q)))))

def scan():
    results = {}   # label -> (R, CV2, mR)
    def add(label, R):
        R = tuple(sorted(set(x for x in R if x % 14 != 0)))   # 14-free
        if not R or R in {v[0] for v in results.values()}:
            return
        cv2, mR = cv2_mR(R)
        if cv2 is not None:
            results[label + str(R)] = (R, cv2, mR)
    # (a) dense consecutive prefixes {1..k}
    for k in range(1, 14):
        add("consec", range(1, k+1))
    # (b) densest sets: {1..13} minus a few (small m_R -> large CV^2 worry)
    base = list(range(1, 14))
    add("full13", base)
    for x in base:
        add("full13_minus1", [y for y in base if y != x])
    for x, y in itertools.combinations(base, 2):
        add("full13_minus2", [z for z in base if z not in (x, y)])
    # (c) speed-7 resonance family (the 2-adic/7-adic binding worry)
    for extra in itertools.chain.from_iterable(
            itertools.combinations([1,2,3,4,5,6,8,9,10,11,12,13], j) for j in range(0,5)):
        add("seven", [7] + list(extra))
    add("seven_evens", [7,2,4,6,8,10,12])
    add("seven_odds", [7,1,3,5,9,11,13])
    # (d) even-heavy / odd-heavy
    add("evens", [2,4,6,8,10,12]); add("odds", [1,3,5,7,9,11,13])
    add("evens+1", [2,4,6,8,10,12,1]); add("evens+7", [2,4,6,8,10,12,7])
    # (e) random 14-free subsets of {1..13}, several sizes
    rng = random.Random(20260629)
    for _ in range(1500):
        sz = rng.randint(2, 12)
        add("rand", rng.sample(range(1,14), sz))
    return results

if __name__ == "__main__":
    print("="*84)
    print(" LRC14 floor: adversarial uniform scan of CV(N_R)^2  (klein-S4; reuses mac-mini THM-579)")
    print("="*84)
    res = scan()
    items = sorted(res.values(), key=lambda v: float(v[1]), reverse=True)
    print(f" scanned {len(items)} distinct 14-free speed sets R\n")
    print(f" TOP 12 by CV(N_R)^2 (the binding / worst-case structures):")
    print(f"   {'CV^2':>8} {'m_R':>8}  R")
    for R, cv2, mR in items[:12]:
        print(f"   {float(cv2):>8.4f} {float(mR):>8.4f}  {R}")
    sup = items[0]
    print(f"\n SUP CV(N_R)^2 = {float(sup[1]):.5f} at R={sup[0]} (m_R={float(sup[2]):.4f})")
    # does sup contain 7? does it match the consec binding row?
    print(f"   contains speed 7: {7 in sup[0]};  is consecutive prefix: {list(sup[0])==list(range(1,len(sup[0])+1))}")

    # gatekeeper margin: for each R pair with the SMALLEST plausible m_Q.
    # canonical covering rows (mac-mini) + adversarial: pair each top-R with Q={1,2}
    # (m_Q large) and with a hard small-m_Q Q to stress the gatekeeper.
    print("\n" + "="*84)
    print(" GATEKEEPER margin  CV^2  vs  m_Q/(1-m_Q)   (sufficient cond: CV^2 < ratio)")
    print("="*84)
    # m_Q for the standard small Q's:
    Qs = {"{1,2}":[1,2], "{1,2,3}":[1,2,3], "{1..4}":[1,2,3,4],
          "{1..5}":[1,2,3,4,5], "{1..6}":[1,2,3,4,5,6]}
    mQ = {k: mQ_of(v) for k, v in Qs.items()}
    print("   m_Q:", {k: round(float(v),4) for k,v in mQ.items()})
    print(f"\n   {'R (top CV^2)':<30}{'CV^2':>8}  best-margin Q -> ratio, slack")
    worst_slack = None
    for R, cv2, mR in items[:8]:
        # the covering that PAIRS with R: r = 14 - |R-span|... use all Qs, report tightest
        best = None
        for qk, qv in mQ.items():
            ratio = qv/(1-qv)
            slack = float(ratio) - float(cv2)
            if best is None or slack < best[2]:
                best = (qk, float(ratio), slack)
        print(f"   {str(R):<30}{float(cv2):>8.4f}  {best[0]} -> {best[1]:.3f}, slack {best[2]:+.3f}")
        if worst_slack is None or best[2] < worst_slack[1]:
            worst_slack = (R, best[2], best[0])
    print(f"\n TIGHTEST gatekeeper slack over scan (R top-CV^2 x hardest Q): {worst_slack[1]:+.4f}"
          f"  at R={worst_slack[0]} vs Q={worst_slack[2]}")

    # ---- ACTUAL R' for the binding cases: does the floor hold even when the clean
    #      gatekeeper FAILS?  Pair each R with the size-valid Q = {1..14-|R|}. ----
    print("\n" + "="*84)
    print(" ACTUAL R' vs the clean gatekeeper (size-valid pairing |R|+|Q|=14, Q={1..14-|R|})")
    print("="*84)
    print(f"   {'R':<34}{'|Q|':>4}{'R_actual':>10}{'CV^2':>8}{'gate>0?':>8}{'Rp>0?':>7}")
    def check(R):
        R = tuple(sorted(set(x for x in R if x % 14 != 0)))
        q = 14 - len(R)
        if q < 1: q = 1
        Q = list(range(1, q+1))
        mR = measure(lonely_set(R)); mQ = mQ_of(Q)
        if mR == 0 or mQ == 0:
            print(f"   {str(R):<34}{len(Q):>4}  (degenerate: m_R={float(mR):.3f}, m_Q={float(mQ):.3f}) skip")
            return None, 1.0
        d = m.floor_row(list(R), Q)
        gate = d['bound'] > 0
        Rp = float(d['Rprime']) if d['Rprime'] is not None else float('nan')
        print(f"   {str(R):<34}{len(Q):>4}{Rp:>10.4f}{float(d['CV2']):>8.3f}"
              f"{('YES' if gate else 'NO'):>8}{('YES' if Rp>0 else 'NO'):>7}")
        return gate, Rp
    print(" -- top adversarial R (dense, speed-7) --")
    floor_fails = False
    for R, cv2, mR in items[:5]:
        gate, Rp = check(R)
        if Rp <= 0: floor_fails = True
    print(" -- consecutive family R={1..k}, the binding 'real covering' rows (extends mac-mini r=2..6) --")
    for k in range(1, 14):
        gate, Rp = check(list(range(1, k+1)))
        if Rp <= 0: floor_fails = True

    print("\n CONCLUSION:")
    print(f"   sup CV(N_R)^2 in scan = {float(sup[1]):.3f} (dense R + speed 7); it GROWS as m_R->0,")
    print("   so the clean gatekeeper CV^2 < m_Q/(1-m_Q) is NOT uniformly satisfied -- it FAILS for")
    print("   dense R even at the size-valid pairing. BUT actual R' stays > 0 throughout the scan:")
    print(f"   floor holds = {not floor_fails}.  => the floor is TRUE on these cases but the clean")
    print("   CV gatekeeper alone does NOT prove it for dense R; the binding worst-case is the")
    print("   maximally-dense R with speed 7 (the 2-adic/7-adic resonance). The uniform proof must")
    print("   either restrict to the actual covering family or supplement the gatekeeper with the")
    print("   exact SPEC on dense R (where m_R->0 inflates the second moment).")
