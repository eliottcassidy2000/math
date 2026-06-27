#!/usr/bin/env python3
"""Three notions of 'sameness' as a fiber of invariants on the LRC lonely set (mac-mini-S62).

Owner: consider equidecomposability and equinumerosity in addition to equidistribution; find invariants
that capture the fundamental nature. Prior project frame (HYP-2187): equinumerosity = cardinal shadow
(a count); equidecomposability = a fiber retaining the predicate's invariants; H = volume, beta1 = Dehn.
HYP-2239 triune carrier = sum/product/fraction faces. Here we carry that fiber onto the LONELY SET itself.

For a speed set S (the runners), the lonely set  L(S) = {t in [0,1): ||s t|| >= 1/14 for all s}  is a
finite union of arcs.  We attach a THREE-RESOLUTION fiber:

  EQUINUMEROSITY (counting / points):   covering status mod 14 (is 0 hit), #distinct residues mod 14,
                                        and #arcs of L(S)           -- the discrete 'cardinal shadow'
  EQUIDECOMPOSABILITY (cut-and-paste):  the arc-length scissors spectrum (#distinct lengths, three-gap),
                                        and D(S) = min denominator d with some a/d in L(S)
                                        (the 'rational reassembly' denominator = scissors deficiency)
  EQUIDISTRIBUTION (limiting density):  meas(L(S)) (the limiting lonely density) + star-discrepancy proxy

CLAIM under test: the TRIPLE (covering, D(S), meas) is a finer 'fundamental-nature' invariant than any
single one, separating tight (AP/GW, M=1/14) from generic and from covering-large-apex; and D(S) (the
equidecomposability/scissors deficiency) reproduces the V*/Conjecture-7.1 constant of S61.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def norm(x):  # ||x|| distance to nearest integer, x a Fraction
    f = x - int(x)
    if f < 0: f += 1
    return min(f, 1 - f)

def lonely_arcs(S, kk1=14):
    S = sorted(set(s for s in S if s != 0))
    bp = set([F(0), F(1)])
    for s in S:
        for j in range(0, s + 1):
            for sign in (F(-1), F(1)):
                x = F(j, s) + sign * F(1, kk1 * s)
                if 0 <= x <= 1: bp.add(x)
    bp = sorted(bp); arcs = []
    for i in range(len(bp) - 1):
        a, b = bp[i], bp[i + 1]
        if b <= a: continue
        mid = (a + b) / 2
        if all(norm(s * mid) >= F(1, kk1) for s in S): arcs.append((a, b))
    return arcs

def min_witness_denom(S, kk1=14, dmax=4000):
    """D(S) = min d with some a/d satisfying ||s a/d|| >= 1/14 for all s. The scissors-reassembly deficiency."""
    for d in range(1, dmax + 1):
        for a in range(0, d):
            if all(norm(F(a * s, d)) >= F(1, kk1) for s in S):
                return d
    return None

def invariants(S, kk1=14):
    res = sorted(set(s % kk1 for s in S))
    covering = (0 in [s % kk1 for s in S])
    arcs = lonely_arcs(S, kk1)
    meas = sum(b - a for a, b in arcs)
    lengths = sorted(set(b - a for a, b in arcs))
    lmax = max((b - a for a, b in arcs), default=F(0))
    D = min_witness_denom(S, kk1)
    return {
        "covering": covering, "n_res": len(res), "n_arcs": len(arcs),
        "meas": meas, "n_lengths": len(lengths), "D": D,
        "lmax": lmax, "inv_lmax": (int(1/lmax) if lmax > 0 else None),
    }

# config families
FAM = {
    "AP {1..13} (tight)":         list(range(1, 14)),
    "GW {1..11,13,24} (tight)":   [1,2,3,4,5,6,7,8,9,10,11,13,24],
    "generic non-covering":       [1,2,3,5,8,13,21,34,55,89,100,111,122],
    "easy-cover {1..12,14}":      list(range(1,13)) + [14],
    "hard-cover {1..11,13,84}":   list(range(1,12)) + [13,84],
    "hard-cover {1..11,13,168}":  list(range(1,12)) + [13,168],
    "hard-cover {1..11,13,840}":  list(range(1,12)) + [13,840],
}
print("=" * 104)
print(" THREE-SAMENESS FIBER on the LRC lonely set  L(S)={t: ||s t||>=1/14}   (kk1=14)")
print("=" * 104)
print(f"{'config':<26} | EQUINUM            | EQUIDECOMP (scissors)              | EQUIDIST")
print(f"{'':<26} | cover #res #arcs   | #len  D(min wit.)  1/lmax(worst)   | meas(L)")
print("-" * 104)
for name, S in FAM.items():
    iv = invariants(S, kk1=14)
    Ds = str(iv["D"]) if iv["D"] is not None else ">dmax"
    il = str(iv["inv_lmax"]) if iv["inv_lmax"] is not None else "-"
    print(f"{name:<26} | {str(iv['covering']):<5} {iv['n_res']:>4} {iv['n_arcs']:>5}   "
          f"| {iv['n_lengths']:>4}  D={Ds:<8} 1/lmax={il:<8} "
          f"| {float(iv['meas']):.5f}")
print("-" * 104)
print("READINGS (corrected):")
print(" * EQUIDIST  meas(L)=0  EXACTLY characterizes the TIGHT atoms (AP, GW both -> 0 arcs, meas 0).")
print("   The density face detects the extremal locus -- this is the witness floor hitting 0.")
print(" * EQUINUM (covering) is INDEPENDENT of tightness: AP is non-covering yet tight; a dilated AP is")
print("   covering AND tight (S60). So 'covering' is a cardinal shadow, predicate-blind to tightness.")
print(" * EQUIDECOMP has TWO scissors invariants, and they are DIFFERENT: D(S)=min witness denominator")
print("   (smallest d with ANY a/d lonely -- the easiest reassembly) stays small; 1/lmax (worst arc, the")
print("   S61 V* quantity) GROWS with the apex. D(S) <= 1/lmax; the V*~200 wall is the 1/lmax face, not D.")
print(" The fiber (covering | D, 1/lmax, #len | meas) is the 'fundamental nature': no single column")
print(" separates {tight, generic, easy-cover, hard-cover}; the TRIPLE does.  meas=0 = tightness detector.")
