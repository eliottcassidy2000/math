# -*- coding: utf-8 -*-
# klein-2026-07-11-S251: THE GLOBAL Phi-ARGMAX CENSUS (exact, exhaustive per box).
#
# CONTEXT: kps THM-701 (wide-spread recursion) + mac-mini THM-702 (explicit finite
# certificate) reduce the wide-spread wall to ONE lemma: the balanced-core base check
#     Phi(F) = p0(F) + (1/3) p1(F) <= cap_{|F|+1}   for every bounded-spread core F.
# mac-mini cont.29 REFUTED the local-compression route (Phi not locally monotone),
# so the lemma is GLOBAL.  Evidence so far: 1500 RANDOM sets (kps, sampled) + the
# consec argmax exact (mac-mini, k=8,9,10).  NOBODY has decided the global argmax.
#
# THIS SCRIPT: for each core size k, EXHAUSTIVELY enumerate ALL normalized cores
#     F = {0} u {k-1 distinct offsets in [1..W]},  gcd(F \ {0}) = 1,
# and compute Phi(F) EXACTLY (all-integer breakpoint sweep on the common denominator
# D = 7*lcm(F); sector(e*x) = floor(7*e*x) mod 7; p_j = interval-length sums in 1/D
# units).  Report: the exact global max Phi per box, its argmax, whether the argmax
# is the CONSECUTIVE core, the margin against cap_{k+1}, the p0-argmax, and the
# runner-up gap.  Convention matches kps lrc14_recursion_closure / mac-mini
# lrc14_finite_certificate exactly (sector c = [c/7,(c+1)/7); sector 0 auto by 0 in F;
# p_j = meas{x in [0,1): exactly j of ALL SEVEN sectors empty} -- with 0 in F this
# equals the inner-six convention).
#
# Exact caps (THM-532/534 via mac-mini cont.29):
#   CAP = {8: 2243/5880, 9: 1979/4004, 10: 55/91, 11: 66/91, 12: 6/7, 13: 1}.
# For k <= 6 (p0 = 0 identically: k-1 <= 5 nonzero offsets cannot hit six inner
# sectors), max Phi = (1/3) max p1 is reported as data -- the induction's small-size
# base needs cap_m for m <= 7 SPECIFIED; flagged in the session letter.

import sys
from math import gcd, comb
from fractions import Fraction as F
from itertools import combinations

CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91),
       12: F(6, 7), 13: F(1)}

def lcm(a, b):
    return a // gcd(a, b) * b

def phi_exact(E):
    """Exact (p0, p1, Phi) for offset set E (0 in E). All-integer sweep.
    Returns (p0, p1, Phi) as Fractions."""
    nz = [e for e in E if e != 0]
    L = 1
    for e in nz:
        L = lcm(L, e)
    D = 7 * L                       # common denominator: all breakpoints are m/D
    # breakpoints of sector(e*x): x = m/(7e), i.e. integer multiples of D/(7e) = L/e
    pts = set()
    for e in nz:
        step = L // e               # exact: e | L
        pts.update(range(0, D + 1, step))
    pts.add(0); pts.add(D)
    pts = sorted(pts)
    n0 = 0                          # numerator of p0 (units 1/D)
    n1 = 0                          # numerator of p1
    for t1, t2 in zip(pts, pts[1:]):
        s = t1 + t2                 # midpoint = s/(2D)
        hit = 1                     # bitmask; sector 0 auto (offset 0)
        for e in nz:
            c = (7 * e * s // (2 * D)) % 7
            hit |= 1 << c
        empty = 7 - bin(hit).count("1")
        if empty == 0:
            n0 += t2 - t1
        elif empty == 1:
            n1 += t2 - t1
    p0 = F(n0, D); p1 = F(n1, D)
    return p0, p1, p0 + p1 / 3

def census(k, W):
    """Exhaustive box: F = {0} u (k-1 offsets from [1..W]), gcd-normalized."""
    total = comb(W, k - 1)
    best_phi = None; best_set = None
    best_p0 = None; best_p0_set = None
    second_phi = None
    consec_phi = None
    n_eval = 0; n_skip = 0
    for combo in combinations(range(1, W + 1), k - 1):
        g = combo[0]
        for c in combo[1:]:
            g = gcd(g, c)
            if g == 1:
                break
        if g > 1:
            n_skip += 1
            continue
        n_eval += 1
        E = (0,) + combo
        p0, p1, phi = phi_exact(E)
        if E == tuple(range(k)):
            consec_phi = phi
        if best_phi is None or phi > best_phi:
            second_phi = best_phi
            best_phi = phi; best_set = E
        elif second_phi is None or phi > second_phi:
            if E != best_set:
                second_phi = phi
        if best_p0 is None or p0 > best_p0:
            best_p0 = p0; best_p0_set = E
    return dict(total=total, n_eval=n_eval, n_skip=n_skip,
                best_phi=best_phi, best_set=best_set, second_phi=second_phi,
                best_p0=best_p0, best_p0_set=best_p0_set, consec_phi=consec_phi)

def main():
    boxes = [(7, 24), (8, 20), (9, 17), (5, 42), (6, 30), (10, 15), (11, 14), (4, 60)]
    if len(sys.argv) > 1:  # allow single-box invocation: k W
        boxes = [(int(sys.argv[1]), int(sys.argv[2]))]
    print("THE GLOBAL Phi-ARGMAX CENSUS (exact, exhaustive per box)")
    print("Phi(F) = p0 + p1/3;  target: Phi(F) <= cap_{|F|+1} with argmax = consec")
    print("=" * 78)
    for k, W in boxes:
        r = census(k, W)
        cap = CAP.get(k + 1)
        consec = tuple(range(k))
        is_consec = (r["best_set"] == consec)
        print(f"\n|F| = k = {k}, box [1..{W}] ({r['n_eval']} normalized sets, "
              f"{r['n_skip']} gcd-skipped of C({W},{k-1}) = {r['total']}):")
        print(f"  GLOBAL max Phi = {r['best_phi']} ~ {float(r['best_phi']):.6f}")
        print(f"    argmax = {r['best_set']}  "
              f"{'== CONSEC' if is_consec else '!= consec  *** NON-CONSEC ARGMAX ***'}")
        if r["consec_phi"] is not None:
            d = r["best_phi"] - r["consec_phi"]
            print(f"    consec Phi = {r['consec_phi']} ~ {float(r['consec_phi']):.6f}"
                  f"   (max - consec = {d} ~ {float(d):.6f})")
        if r["second_phi"] is not None:
            print(f"    runner-up Phi ~ {float(r['second_phi']):.6f} "
                  f"(gap {float(r['best_phi'] - r['second_phi']):.6f})")
        if cap is not None:
            m = cap - r["best_phi"]
            print(f"  cap_{k+1} = {cap} ~ {float(cap):.6f}; "
                  f"GLOBAL margin = {m} ~ {float(m):+.6f}  "
                  f"{'OK' if m > 0 else '*** VIOLATED ***'}")
        else:
            print(f"  cap_{k+1} NOT IN LEDGER (k+1 = {k+1} < 8): max Phi reported as "
                  f"data for the small-size base spec")
        print(f"  p0-argmax: {r['best_p0_set']}  p0 = {r['best_p0']} ~ "
              f"{float(r['best_p0']):.6f}")

if __name__ == "__main__":
    main()
