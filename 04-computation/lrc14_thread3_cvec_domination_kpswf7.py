#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD 3 -- FOLLOW-UP to the PROMISING Riesz/conductance lead.

Consec uniquely maximizes min_r c_r where c_r = sum_{e mod7==r, e>0} 1/|e|.
NEW HYP candidate: does consec COORDINATEWISE-dominate the c_r vector? And is measS7 a
MONOTONE function of the c_r vector (so coordinatewise domination => measS7 domination)?

If BOTH hold within the bounded/full-residue stratum, LAYER 3 collapses to:
   (A) consec maximizes EACH c_r  [c_r(consec) = 1/r, the smallest possible magnitude for residue r],
   (B) measS7 is coordinatewise-monotone-increasing in c.
=> consec maximizes measS7. A clean two-line proof.

The bounded/full-residue stratum: E has a representative of EVERY residue r=1..6 (full residue),
all |e| <= some bound (bounded). consec = (0,1,2,3,4,5,6,...) achieves c_r >= 1/r with the apex
identity e=0 doubled and small reps.

PROBE 1: coordinatewise domination c_r(consec) >= c_r(E) for ALL r, all bounded full-residue E?
PROBE 2: is measS7 monotone in c? Test by perturbation: increase one |e| (decrease one c_r) and
         check measS7 does not increase. (This is the electrical-contraction principle; mac-mini
         found 11331 violations to *magnitude* monotonicity earlier -- but does c-vector
         monotonicity survive where raw-magnitude monotonicity failed?)
PROBE 3: does the c-vector DETERMINE measS7 (HYP-2760 says residue-legs do; does full c?)
"""
import sys, itertools, random
from fractions import Fraction as F
from functools import reduce
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def covered_sectors(E, xm):
    s = set()
    for e in E:
        v = F(e) * xm; v = v - F(v.numerator // v.denominator)
        s.add((v.numerator * 7) // v.denominator)
    return s

def measS7(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if covered_sectors(E, (lo + hi) / 2) >= set(range(1, 7)): tot += hi - lo
    return tot

def cvec(E):
    c = [F(0)] * 7
    for e in E:
        e = int(e)
        if e == 0: continue
        c[e % 7] += F(1, abs(e))
    return tuple(c[r] for r in range(1, 7))

def full_residue(E):
    res = set(int(e) % 7 for e in E if int(e) != 0)
    return res >= set(range(1, 7))

def rand_full_config(k, hi, seed, maxtry=400):
    random.seed(seed)
    for _ in range(maxtry):
        E = sorted(set([0] + random.sample(range(1, hi), k - 1)))
        if len(E) == k and reduce(gcd, E) == 1 and full_residue(E):
            return E
    return None

def main():
    print("=== c-vector domination & monotonicity (bounded full-residue stratum) ===\n")

    print("PROBE 1: does consec coordinatewise-dominate c_r? c(consec)=(1,1/2,1/3,1/4,1/5,1/6)")
    for k in (8, 9, 10):
        consec = list(range(k))
        cc = cvec(consec)
        dom = 0; notdom = 0; n = 0; counterex = None
        for s in range(300):
            E = rand_full_config(k, 14, 9000 + s)
            if E is None: continue
            ce = cvec(E); n += 1
            if all(cc[i] >= ce[i] for i in range(6)): dom += 1
            else:
                notdom += 1
                if counterex is None: counterex = (E, [float(x) for x in ce])
        print(f"  k={k} consec c={[float(x) for x in cc]}: dominates {dom}/{n}, FAILS {notdom}/{n}")
        if counterex: print(f"     example NOT dominated: E={counterex[0]} c={counterex[1]}")
    print("  => if consec dominates ALL: coordinatewise extremal (clean). If fails: domination is partial.\n")

    print("PROBE 2: monotonicity -- does INCREASING a c-coordinate (shrinking a |e|) never DECREASE measS7?")
    # take config, replace a speed e by a speed e' with same residue but smaller |e'| (bigger c_r),
    # keeping it a valid distinct config; check measS7 monotone up.
    viol = 0; tot = 0; worst = F(0)
    for k in (8, 9):
        for s in range(120):
            E = rand_full_config(k, 14, 11000 + s)
            if E is None: continue
            Es = sorted(int(e) for e in E)
            for idx in range(len(Es)):
                e = Es[idx]
                if e <= 1: continue
                # candidate smaller-magnitude same-residue speed
                for ep in range(1, e):
                    if ep % 7 == e % 7 and ep not in Es:
                        E2 = sorted(set(Es[:idx] + [ep] + Es[idx+1:]))
                        if len(E2) != k or reduce(gcd, E2) != 1: continue
                        m1 = measS7(Es); m2 = measS7(E2)
                        c1 = cvec(Es); c2 = cvec(E2)
                        # c2 should dominate c1 coordinatewise (only one coord changed, increased)
                        if all(c2[i] >= c1[i] for i in range(6)):
                            tot += 1
                            if m2 < m1:  # measS7 decreased despite c increasing => violation
                                viol += 1
                                if m1 - m2 > worst: worst = m1 - m2
                        break
    print(f"  c-monotonicity: {viol}/{tot} violations (measS7 DECREASED when c increased); worst drop={float(worst):.5f}")
    print("  => 0 violations: measS7 is coordinatewise-monotone in c (the electrical-contraction law).")
    print("     >0: monotone in c FAILS (resonance overrides conductance even in c-coords).\n")

    print("PROBE 3: does the c-vector DETERMINE measS7? (collisions = same c, different measS7)")
    seen = {}; coll = 0; tot = 0
    for k in (8, 9):
        for s in range(400):
            E = rand_full_config(k, 16, 13000 + s)
            if E is None: continue
            key = (k, cvec(E)); m = measS7(E); tot += 1
            if key in seen:
                if seen[key] != m: coll += 1
            else:
                seen[key] = m
    print(f"  {coll} c-vector collisions with DIFFERENT measS7 out of {tot} configs")
    print("  => 0: c determines measS7 (sharp invariant). >0: c does not determine measS7.")

if __name__ == "__main__":
    main()
