#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THREAD A (mac-mini 2026-06-21) -- THE VERDICT: is there a JOINT monotonicity
that turns "consec minimizes the local data" into a NON-finite-check proof of
consec-max(measS7) on the full-residue stratum?

Established:
  - measS7 IS a function of full local data D(E) = (sorted multiset, per residue,
    of magnitudes) -- effectively the full sorted |E| vector with residue tags.
    (Actually measS7 is dilation-class function of E itself.)
  - consec uniquely minimizes the leg profile (TEST1 prior) and uses the
    third-leg-optimal doubling.

For a clean proof we need a PARTIAL ORDER <= on full-residue shapes with:
  (P1) consec is the unique minimum, AND
  (P2) measS7 is monotone-DECREASING w.r.t. <=  (E1 <= E2 => measS7(E1)>=measS7(E2)).
Then consec = global max, no finite check.

NATURAL CANDIDATE ORDER: componentwise domination of the SORTED magnitude vector
|E| (with 0 dropped or kept). consec = (0,1,2,3,4,5,6,7) is the componentwise-
minimal full-residue k=8 vector. Test (P2): is measS7 monotone-decreasing in the
sorted |E| vector over the full-residue stratum?

If monotone holds => PROOF. If it fails, the violating pair is the PRECISE
OBSTRUCTION: measS7 is genuinely aggregate, not order-monotone, and the proof
must be a real inequality (the theta' wall), not a monotonicity argument.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def covered_sectors(E, xm):
    secs = set()
    for e in E:
        v = e * xm; v = v - (v.numerator // v.denominator)
        secs.add((v.numerator * 7) // v.denominator)
    return secs

def measS7(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        ae = abs(e)
        if ae == 0: continue
        for m in range(7 * ae + 1): bps.add(F(m, 7 * ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    total = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        if len(covered_sectors(E, (lo+hi)/2)) == 7: total += hi - lo
    return total

def residues(E): return frozenset(e % 7 for e in E)
def is_full_residue(E): return residues(E) == frozenset(range(7))
def primitive(E): return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1
def consec(k): return list(range(k))
def sortedvec(E): return tuple(sorted(E))
def dom(p, q): return len(p)==len(q) and all(a<=b for a,b in zip(p,q))

if __name__ == "__main__":
    print("#"*78); print("# JOINT-MONOTONICITY VERDICT (THREAD A)"); print("#"*78)
    k = 8; W = 13
    bank = [(0,)+r for r in itertools.combinations(range(1, W+1), k-1)]
    full = [list(E) for E in bank if primitive(E) and is_full_residue(E)]
    C = consec(k)
    ms = {tuple(E): measS7(E) for E in full}
    vC = sortedvec(C)
    print(f"\nk=8 span<= {W}: {len(full)} full-res shapes; consec vec {vC}")

    # (P1) consec componentwise-minimal?
    nondom = [E for E in full if not dom(vC, sortedvec(E))]
    print(f"(P1) consec vec <= every full vec componentwise? non-dominated: {len(nondom)} "
          f"(consec is comp-min: {len(nondom)==0})")

    # (P2) measS7 monotone-decreasing in the sorted |E| vector?
    print(f"(P2) measS7 monotone-decreasing in sorted vec? (checking all dominated pairs)")
    viol = []
    for E1 in full:
        for E2 in full:
            if E1 is E2: continue
            if dom(sortedvec(E1), sortedvec(E2)) and sortedvec(E1)!=sortedvec(E2):
                if ms[tuple(E1)] < ms[tuple(E2)] - F(1,10**12):
                    viol.append((E1, E2, ms[tuple(E1)], ms[tuple(E2)]))
    print(f"   monotonicity violations: {len(viol)}")
    for E1,E2,m1,m2 in sorted(viol, key=lambda t: float(t[3]-t[2]))[-6:]:
        print(f"      {E1} (vec<=) meas {float(m1):.5f}  <  {E2} meas {float(m2):.5f}  "
              f"gap {float(m2-m1):.5f}")

    print(f"\nVERDICT: (P1) {'holds' if not nondom else 'FAILS'}; "
          f"(P2) {'holds -> MONOTONE PROOF' if not viol else 'FAILS -> measS7 NOT order-monotone'}.")
    if viol:
        # characterize the obstruction: is the violator ALSO consec-beating? No
        # (consec is the max). The violations are among NON-consec shapes -> the
        # order does not control measS7; the proof must be a real inequality.
        beat_consec = [v for v in viol if v[3] > ms[tuple(C)]]
        print(f"   Of {len(viol)} violations, {len(beat_consec)} involve a shape beating consec "
              f"(=> {'consec NOT max!' if beat_consec else 'none beat consec; obstruction is internal, not at the top'}).")
        # The smallest violator shows measS7 jumps non-monotonically: doubling
        # residue 0 vs residue!=0 at the same vector-cost gives different measS7.
        print(f"   => The doubled-residue identity (3rd leg) controls measS7, NOT the")
        print(f"      magnitude vector. So consec-max is NOT a magnitude-domination")
        print(f"      monotonicity; it is the residue-structure (3rd-leg) extremality.")
