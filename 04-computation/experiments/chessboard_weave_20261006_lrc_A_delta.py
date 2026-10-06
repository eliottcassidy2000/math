#!/usr/bin/env python3
"""Task A (chessboard weave, LRC reading), 2026-10-06.  Exact arithmetic throughout.

A1  the candidate lemma behind delta() (printed proof; machine-checked in the core self-test)
A2  rider lines: delta(a,b) = floor((a+b)/2)/(a+b) for coprime a<b  (proof + exact check b<=200)
A3  tight instances delta(v) = 1/(n+1): exact census of primitive speed sets
      - via delta() itself on the "delta universes" (n=3 B=60, n=4 B=30, n=5 B=18)
      - via the faster exact test "Safe(1/(n+1)) nonempty with empty interior"
        on larger universes (n=2 B=400, n=3 B=150, n=4 B=60, n=5 B=36, n=6 B=24, n=7 B=20)
    plus the smallest non-tight values of delta (the gap above 1/(n+1)).
Reproduce:  python3 chessboard_weave_20261006_lrc_A_delta.py > chessboard_weave_20261006_lrc_A_delta.out
"""
import sys
import time
from fractions import Fraction
from math import gcd
from collections import defaultdict

sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else ".")
from chessboard_weave_20261006_lrc_core import (delta, delta_argmax, safe_components,
                                                primitive_sets, check, lcm_list)

T0 = time.time()
print("=" * 78)
print("A1. Candidate lemma (PROVED).")
print("""  f(t) = min_i ||t v_i|| is continuous, piecewise linear, every piece has slope
  +-v_i != 0.  Let t0 maximise f.  With A = active indices at t0, the right
  derivative of f is min_{i in A} g_i'(t0+) <= 0 and the left derivative is
  max_{j in A} g_j'(t0-) >= 0, so some active tent i decreases just after t0 and
  some active tent j increases just before t0.  No active tent has a zero at t0
  (f(t0) > 0).  If i = j, or if i or j has its peak at t0, then t0 is a peak
  (2k+1)/(2v) = (2k+1)/(v+v).  Otherwise near t0  g_j = t v_j - m_j and
  g_i = m_i - t v_i, and equality at t0 gives t0 = (m_i+m_j)/(v_i+v_j).
  Hence delta(v) = max f over {m/(v_i+v_j) : i <= j, 0 < m <= (v_i+v_j)/2}.
  (Difference crossings m/|v_i-v_j| need not be added: every local max is already
  a sum crossing or a peak.)  Machine cross-check against the full breakpoint
  lattice, which does include them: core self-test [2].""")

# ------------------------------------------------------------------ A2
print()
print("=" * 78)
print("A2. Two speeds (rider lines).  CLAIM: coprime 1 <= a < b  =>  delta(a,b) = floor(s/2)/s, s=a+b.")
print("""  PROOF.  gcd(a,b)=1 => the closed orbit {(ta,tb) mod 1} is the whole subgroup
  K = {(x,y) in T^2 : b x - a y in Z}  (kernel of a primitive character, connected,
  and contains the orbit circle).  ||ta||,||tb|| >= c for some t  <=>  K meets
  [c,1-c]^2  <=>  the image interval of [c,1-c]^2 under (x,y) -> bx - ay, namely
  [s c - a, b - s c] (length s(1-2c), centre (b-a)/2), contains an integer.
  s even: a,b odd, (b-a)/2 in Z, so c = 1/2 works: delta = 1/2.
  s odd : (b-a)/2 is a half-integer, need s(1-2c)/2 >= 1/2, i.e. c <= (s-1)/(2s).
  So delta(a,b) = floor(s/2)/s, depends only on the anti-diagonal s = a+b, and its
  minimum over s >= 3 is 1/3 at s = 3, i.e. ONLY the knight (1,2).  QED""")
B2 = 200
bad = 0
cnt = 0
cnt80 = 0
minval = None
minpairs = []
for b in range(2, B2 + 1):
    for a in range(1, b):
        if gcd(a, b) != 1:
            continue
        s = a + b
        d = delta((a, b))
        cnt += 1
        if b <= 80:
            cnt80 += 1
        if d != Fraction(s // 2, s):
            bad += 1
            print("  MISMATCH", (a, b), d)
        if minval is None or d < minval:
            minval, minpairs = d, [(a, b)]
        elif d == minval:
            minpairs.append((a, b))
check(bad == 0, "claim 1 mismatch")
print(f"  exact check: {cnt} coprime pairs with b <= {B2} (of which {cnt80} with b <= 80): 0 mismatches")
print(f"  minimum delta = {minval} at {minpairs}  (the knight)")
# scaling invariance spot check
for (a, b) in [(1, 2), (2, 5), (3, 7)]:
    for d_ in (2, 3, 6):
        check(delta((d_ * a, d_ * b)) == delta((a, b)), "scaling")
print("  scaling invariance delta(d v) = delta(v) spot-checked (d = 2,3,6).")
print("  argmax structure for a few riders (t in (0,1/2] attaining delta):")
for v in [(1, 2), (1, 3), (2, 3), (1, 4), (3, 4), (2, 5), (1, 6)]:
    d, pts = delta_argmax(v)
    print(f"    {v}: delta = {d}, argmax = {[str(p) for p in pts]}")

# ------------------------------------------------------------------ A2b
print()
print("  A2b. Rings of an N x N subdivided cell (N even; ring 0 = central 2x2 block, ring k =")
print("  closed squares with min(||x||,||y||) in [1/2-(k+1)/N, 1/2-k/N]).  A line enters the")
print("  CLOSED ring-0 block iff delta >= 1/2 - 1/N, so by A2 the innermost ring met by the (a,b)")
print("  rider is 0 if a+b is even and ceil(N/(2(a+b))) - 1 if a+b is odd (closed squares;")
print("  when 2(a+b) | N the line only touches the next ring at corners).")
for Nb in (8,):
    thr = Fraction(1, 2) - Fraction(1, Nb)
    miss = []
    for b in range(2, 81):
        for a in range(1, b):
            if gcd(a, b) == 1 and delta((a, b)) < thr:
                miss.append((a, b))
    print(f"  {Nb}x{Nb} board: coprime riders (b <= 80) whose line misses the closed central "
          f"2x2 block [3/8,5/8]^2: {miss}")
    check(miss == [(1, 2)], "only the knight misses the centre")
print("  => on the 8x8 cell the knight is the ONLY rider line that never enters the central")
print("     4 squares (it reaches ring 1: delta = 1/3 in [1/4, 3/8)); all others have delta >= 2/5.")

# ------------------------------------------------------------------ A3
print()
print("=" * 78)
print("A3. Tight instances: primitive (gcd 1) sets of n distinct speeds with delta = 1/(n+1).")


def tight_by_safe(v):
    """Exact: delta(v) = 1/(n+1) iff Safe(1/(n+1)) is nonempty and has empty interior
    (f is piecewise linear with nonzero slopes, so f >= c on an interval forces f > c
    somewhere).  Returns (status, safe points or None)."""
    comps, D = safe_components(v)
    if not comps:
        return "LRC-FAIL", None
    if all(lo == hi for lo, hi in comps):
        return "tight", [Fraction(lo, D) for lo, _ in comps]
    return "loose", None


def is_dilated_ap(v):
    return all(v[i] == (i + 1) * v[0] for i in range(len(v)))


# (a) delta universes: compute delta for every primitive set, list tight + smallest values
print()
print("  (a) full delta census (delta computed exactly for every primitive set):")
DELTA_U = {2: 60, 3: 60, 4: 30, 5: 18}
for n, B in DELTA_U.items():
    t1 = time.time()
    vals = defaultdict(list)
    num = 0
    for v in primitive_sets(n, B):
        vals[delta(v)].append(v)
        num += 1
    keys = sorted(vals)
    thr = Fraction(1, n + 1)
    check(keys[0] >= thr, f"LRC violated n={n}")
    tight = vals.get(thr, [])
    # cross-check with the safe-set test
    for v in tight:
        check(tight_by_safe(v)[0] == "tight", v)
    print(f"   n={n}, max speed <= {B}: {num} primitive sets ({time.time() - t1:.1f}s); "
          f"min delta = {keys[0]}; tight sets ({len(tight)}): {tight}")
    shown = 0
    print(f"     smallest delta values: ", end="")
    out = []
    for k in keys[:6]:
        L = vals[k]
        out.append(f"{k} [{len(L)} sets, e.g. {L[:3]}]")
    print("; ".join(out))

# (b) larger universes via the safe-set test
print()
print("  (b) larger universes, exact test 'Safe(1/(n+1)) nonempty with empty interior':")
SAFE_U = {2: 400, 3: 150, 4: 60, 5: 36, 6: 24, 7: 20}
all_tight = {}
for n, B in SAFE_U.items():
    t1 = time.time()
    num = 0
    tight = []
    fails = []
    for v in primitive_sets(n, B):
        num += 1
        st, pts = tight_by_safe(v)
        if st == "tight":
            tight.append((v, pts))
        elif st == "LRC-FAIL":
            fails.append(v)
    check(not fails, f"LRC fails for {fails[:5]}")
    all_tight[n] = tight
    print(f"   n={n}, max speed <= {B}: {num} primitive sets, LRC holds for all, "
          f"{len(tight)} tight ({time.time() - t1:.1f}s)")
    for v, pts in tight:
        dens = sorted({p.denominator for p in pts})
        res = sorted(x % (n + 1) for x in v)
        tag = "AP {1..n}" if is_dilated_ap(v) else "sporadic"
        print(f"      {v}  [{tag}]  safe points {[str(p) for p in pts]}  "
              f"denominators {dens}  v mod {n + 1} = {res}")

print()
print("  Summary of the denominator observation (THM-3043 R3, here re-tested on the")
print("  larger universes): every tight witness has reduced denominator n+1?")
okall = True
for n, tight in all_tight.items():
    for v, pts in tight:
        if any(p.denominator != n + 1 for p in pts):
            okall = False
            print("   EXCEPTION:", n, v, pts)
print(f"   {'yes' if okall else 'NO'} (all tight sets in the universes of (b))")
print(f"\nTotal time {time.time() - T0:.1f}s")
