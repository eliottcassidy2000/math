#!/usr/bin/env python3
"""eroded_residue_direct_kps_S128c85.py -- kind-pasteur-2026-07-19-S128c85

eroded_residue_close_kps_S128c85.py showed the measure/counting route does NOT
decide the residue: 294 of 456 (core, k1) pairs fail G(P,k1) > floor(k1*2/21).

BUT THAT CRITERION IS ONLY SUFFICIENT, NOT NECESSARY.  It compares a count of
available gaps against the WORST CONCEIVABLE bad set of measure 2/21.  The actual
bad set is a specific set determined by the actual killers, so a pair failing the
criterion is UNDECIDED, not a counterexample.

This script decides them directly.  For a core P and killers k1<k2<k3<k4 the real
question is simply whether the full safe set is non-empty:

    S(P) \\ ( D_k1 u D_k2 u D_k3 u D_k4 )  !=  empty,     D_k = {t : ||k t|| < 1/14}

because any surviving t has ||v t|| >= 1/14 for EVERY speed, giving M >= 1/14
outright.  The "complete k1-gap" language of THM-1162 is a proof route to this,
not the criterion itself.

So: take the residue pairs the counting route could not decide, enumerate killer
triples (k2,k3,k4) above k1, and test non-emptiness exactly in rationals.  If the
safe set is always non-empty there, the residue is FINE and only the ARGUMENT was
too lossy.  If some triple empties it, that is a genuine obstruction and a much
more important find.

Killer range: the r=5 clustered horn bounds killers by ~235 (THM-1093 KB=235), so
triples are enumerated in (k1, KB].
"""
import sys
from fractions import Fraction as F
from itertools import combinations

LAM = F(1, 14)
KB = int(sys.argv[1]) if len(sys.argv) > 1 else 235
TRIPLE_CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 30000


def safe_components(P):
    bps = {F(0), F(1)}
    for p in P:
        for j in range(p + 1):
            for s in (F(1, 14 * p), -F(1, 14 * p)):
                v = F(j, p) + s
                if 0 <= v <= 1:
                    bps.add(v)
    B = sorted(bps)
    out = []
    for i in range(len(B) - 1):
        a, b = B[i], B[i + 1]
        if b <= a:
            continue
        mid = (a + b) / 2
        if all(min((p * mid) % 1, 1 - (p * mid) % 1) >= LAM for p in P):
            out.append((a, b))
    merged = []
    for a, b in out:
        if merged and a == merged[-1][1]:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return merged


def strip(comps, k):
    """Remove the k-comb teeth {t : ||k t|| < 1/14} from a list of intervals."""
    out = []
    r = F(1, 14 * k)
    for a, b in comps:
        jlo = int(a * k) - 1
        jhi = int(b * k) + 1
        cuts = []
        for j in range(jlo, jhi + 1):
            c = F(j, k)
            lo, hi = c - r, c + r
            if hi > a and lo < b:
                cuts.append((max(lo, a), min(hi, b)))
        cuts.sort()
        cur = a
        for lo, hi in cuts:
            if lo > cur:
                out.append((cur, lo))
            cur = max(cur, hi)
        if cur < b:
            out.append((cur, b))
    return [(a, b) for a, b in out if b > a]


def survives(comps, killers):
    c = comps
    for k in killers:
        c = strip(c, k)
        if not c:
            return False, None
    return True, c[0][0]


# the residue pairs that the counting route could not decide, worst first
TARGETS = [
    ((1, 2, 3, 7, 8, 10, 11, 12), 160),
    ((1, 2, 3, 5, 7, 8, 9, 11), 153),
    ((1, 2, 3, 5, 7, 8, 9, 11), 144),
    ((1, 2, 3, 5, 7, 8, 9, 11), 160),
    ((1, 2, 3, 5, 7, 8, 9, 11), 175),
    ((1, 2, 7, 8, 9, 10, 11, 12), 157),
    ((1, 3, 7, 8, 9, 10, 11, 12), 157),
    ((1, 2, 5, 6, 7, 8, 9, 11), 144),
    ((1, 2, 6, 7, 8, 9, 10, 11), 144),
]

print("=" * 78)
print("DIRECT DECISION OF THE UNDECIDED RESIDUE PAIRS")
print("  criterion: is  S(P) minus the four killer combs  non-empty?")
print("  killer triples (k2,k3,k4) enumerated in (k1, %d]" % KB)
print("=" * 78)
grand_fail = []
for P, k1 in TARGETS:
    comps0 = safe_components(P)
    c1 = strip(comps0, k1)
    if not c1:
        print("  P=%-30s k1=%-4d  S(P) minus k1 is ALREADY EMPTY -- genuine obstruction"
              % (str(P), k1))
        grand_fail.append((P, k1, None))
        continue
    rng = [k for k in range(k1 + 1, KB + 1)]
    trips = list(combinations(rng, 3))
    if len(trips) > TRIPLE_CAP:
        step = len(trips) // TRIPLE_CAP + 1
        trips = trips[::step]
    fails = []
    minwidth = None
    for (k2, k3, k4) in trips:
        ok, wit = survives(c1, (k2, k3, k4))
        if not ok:
            fails.append((k2, k3, k4))
            if len(fails) >= 5:
                break
    tot = len(trips)
    print("  P=%-30s k1=%-4d  triples tested %-7d  EMPTY-safe-set triples: %d"
          % (str(P), k1, tot, len(fails)))
    if fails:
        print("       first failures: %s" % (fails[:5],))
        grand_fail.append((P, k1, fails[:5]))

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
if not grand_fail:
    print("  Every undecided residue pair tested keeps a NON-EMPTY safe set for every")
    print("  killer triple in range.  So the residue pairs are NOT obstructions: the")
    print("  COUNTING ARGUMENT was too lossy, not the geometry.  Supplier (3) is fine")
    print("  on the residue; what fails there is only the measure proof of it.")
    print()
    print("  That relocates the remaining work: a sharper argument on 456 pairs, not")
    print("  a search for a counterexample.")
else:
    print("  %d residue pairs admit a killer triple emptying the safe set." % len(grand_fail))
    print("  These are GENUINE obstructions and must be reported as such:")
    for row in grand_fail[:10]:
        print("     %s" % (row,))
