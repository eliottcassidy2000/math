#!/usr/bin/env python3
"""
lrc13_tight_locus_rigidity_kps_S1.py  (kind-pasteur-2026-07-05-S1, HYP-4096)

THE TIGHT-LOCUS RIGIDITY for LRC(13) (13 runners = 12 nonzero speeds).

CLAIM under test (mac-mini-S47, flagged MISTAKE-100 risk class; needed by
klein-S132 `tight_ap_free_rider` to cover the whole tight case of hcomp):

    For W a 12-element set of positive integers with gcd(W)=1,
    M(W) := max_t min_{w in W} ||w t||  =  1/13   <==>   W = {1,2,...,12}.

M(W) is computed EXACTLY: g(t)=min_w ||w t|| is piecewise linear; every local
max sits at a crossing ||w_i t||=||w_j t||, i.e. t=m/(w_i +/- w_j).  So the
global max-min is attained at some t=m/D, D a pairwise sum or difference.

is_tight(W): scan all critical t; LRC(13) guarantees M(W)>=1/13, so W is tight
iff NO critical t gives min_w||w t|| > 1/13.  Early-exit the instant one does
(the common LOOSE case) -> a large hunt is feasible.

GOAL: (a) confirm {1..12} tight; (b) HUNT for any OTHER primitive 12-set with
M=1/13 (Goddyn-Wong sporadics, structured perturbations, dilations of small
non-AP tight sets).  A counterexample REDIRECTS the free-rider route.
"""
from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce
import random

TARGET = Fraction(1, 13)

def nrm_num(a, D):
    """||a/D|| as a Fraction, a,D ints, D>0."""
    r = a % D
    return Fraction(min(r, D - r), D)

def critical_Ds(W):
    Ds = set()
    for a, b in combinations(W, 2):
        Ds.add(a + b)
        Ds.add(abs(a - b))
    Ds.discard(0)
    return Ds

def is_tight(W):
    """True iff M(W)==1/13 (given M(W)>=1/13 by LRC13). Early-exits on loose."""
    for D in critical_Ds(W):
        for m in range(1, D):
            # min_w ||w*m/D||  = min_w nrm(w*m mod D)/D
            mn = D  # track numerator of min distance *over denominator D*
            for w in W:
                r = (w * m) % D
                d = r if r <= D - r else D - r
                if d < mn:
                    mn = d
                    if mn == 0:
                        break
            # min distance = mn/D ; compare to 1/13
            if mn * 13 > D:        # mn/D > 1/13  -> loose
                return False
    return True

def M_exact(W):
    best = Fraction(0)
    best_t = None
    for D in critical_Ds(W):
        for m in range(1, D):
            mn = D
            for w in W:
                r = (w * m) % D
                d = r if r <= D - r else D - r
                if d < mn:
                    mn = d
            val = Fraction(mn, D)
            if val > best:
                best = val
                best_t = Fraction(m, D)
    return best, best_t

def is_prim(W):
    return reduce(gcd, W) == 1

def report(tag, W):
    W = tuple(sorted(set(W)))
    if len(W) != 12:
        print(f"  {tag:26s} W={W}  (not 12 distinct!) skipped")
        return None
    m, t = M_exact(W)
    if m == TARGET:
        flag = "  <<< TIGHT (=1/13)"
    elif m < TARGET:
        flag = "  [BELOW 1/13 -- IMPOSSIBLE, bug!]"
    else:
        flag = ""
    print(f"  {tag:26s} M={m} (~{float(m):.5f}) t={t}{flag}")
    return m

print("=" * 78)
print("LRC(13) TIGHT-LOCUS RIGIDITY  --  M(W)=1/13 <=> W={1..12} ?   (kps-S1)")
print("=" * 78)

print("\n[1] AP {1..12} and dilations (all should be TIGHT):")
for c in range(1, 7):
    report(f"c={c}*(1..12)", [c * j for j in range(1, 13)])

print("\n[2] Goddyn-Wong & sporadic candidates (n=13 runners, 12 speeds):")
cands = {
    "{1..11,13}":      list(range(1, 12)) + [13],
    "{1..11,24}":      list(range(1, 12)) + [24],
    "{1..10,12,13}":   list(range(1, 11)) + [12, 13],
    "{2..13}":         list(range(2, 14)),
    "{1..10,12,22}":   list(range(1, 11)) + [12, 22],
    "{1..11,23}":      list(range(1, 12)) + [23],
    "{1..9,11,12,13}": list(range(1, 10)) + [11, 12, 13],
    "{1,3,4,..,13}":   [1] + list(range(3, 14)),
    # GW template {1..n-2,n,2(n-1)} truncated / cousins for a 12-set:
    "{1..10,12,20}":   list(range(1, 11)) + [12, 20],
    "{1..10,11,21}":   list(range(1, 11)) + [11, 21],  # = {1..11,21}
}
for tag, W in cands.items():
    report(tag, W)

print("\n[3] Exhaustive: all primitive 12-subsets of {1,..,16} -- collect tight:")
tight_found, count = [], 0
for W in combinations(range(1, 17), 12):
    if not is_prim(W):
        continue
    count += 1
    if is_tight(W):
        tight_found.append(W)
print(f"  scanned {count} primitive 12-subsets of [1,16];  TIGHT: {len(tight_found)}")
for W in tight_found:
    print(f"     TIGHT: {W}   {'(= AP {1..12})' if W==tuple(range(1,13)) else '<<< NON-AP SPORADIC!'}")

print("\n[4] Randomized hunt for sporadic tight sets (entries<=60, primitive):")
random.seed(1)
sporadic, trials, scanned = set(), 300000, 0
AP = tuple(range(1, 13))
for _ in range(trials):
    hi = random.randint(13, 60)
    W = tuple(sorted(random.sample(range(1, hi + 1), 12)))
    if not is_prim(W):
        continue
    scanned += 1
    if W != AP and is_tight(W):
        sporadic.add(W)
print(f"  {scanned} primitive sets scanned;  NON-AP tight found: {len(sporadic)}")
for W in sorted(sporadic)[:25]:
    print(f"     SPORADIC TIGHT: {W}")

print("\n[5] Dilated small tight sets? test c*S for small tight S with |S|<12 padded:")
# If any non-AP 12-set is tight, its dilations are too; already covered. Instead
# check: does the AP structure survive deletion+insertion (near-AP perturbations)?
near = 0
for drop in range(1, 13):
    for ins in range(13, 40):
        W = sorted(set(range(1, 13)) - {drop} | {ins})
        if len(W) != 12 or not is_prim(W):
            continue
        if is_tight(W):
            near += 1
            print(f"     NEAR-AP TIGHT: dropped {drop}, inserted {ins}: {tuple(W)}")
print(f"  near-AP single-swap tight (besides AP itself): {near}")

print("\n[6] Second-value gap: smallest M(W)>1/13 over primitive 12-subsets [1,15]:")
above = []
for W in combinations(range(1, 16), 12):
    if not is_prim(W):
        continue
    m, _ = M_exact(W)
    if m > TARGET:
        above.append((m, W))
above.sort()
for m, W in above[:6]:
    print(f"     M={m} (~{float(m):.5f})  W={W}")
if above:
    print(f"  => spectral gap: second value ~{float(above[0][0]):.5f} vs 1/13~{float(TARGET):.5f}")

print("\nDONE.")
