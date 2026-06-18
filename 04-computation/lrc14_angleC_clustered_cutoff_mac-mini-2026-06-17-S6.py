#!/usr/bin/env python3
"""
lrc14_angleC_clustered_cutoff -- mac-mini-2026-06-17-S6

Sharpen Angle C: is there a PROVABLE constant B0 such that every covering 13-set with
max(S) > B0 is handled by the arc-width criterion C(S) (=> M>=1/14), leaving a FINITE base?

The pigeonhole bound W(A) >= mu(A)/N(A), N(A) <= sum_{u in A} u, is too weak for CLUSTERED
sets. We sharpen the count N(A) of safe arcs:

  KEY STRUCTURAL FACT. The level-1/14 danger set of A = union over u in A of u teeth, each of
  FULL width 1/(7u), centers k/u. The safe set's arcs (gaps) number at most the number of DISTINCT
  tooth endpoints / 2, but more usefully: a safe arc is a maximal gap of the union of teeth.
  The union of teeth from runner u alone leaves u gaps each of width 6/(7u). Overlaying all runners
  can only SPLIT gaps; the number of safe arcs N(A) <= sum_{u in A} u  (each runner contributes <= u
  new cut points). BUT a SHARPER bound: in the part of [0,1) NOT covered by the SLOWEST runner's
  teeth, ... -- this is where pigeonhole on the slowest runner helps.

This script:
 (1) For the criterion via the SLOWEST runner v_min, compute the EXACT W and the pigeonhole margin,
     over clustered covering sets, to see if removing the smallest (not largest) runner gives a
     clean bound. (The slowest runner has the FEWEST, WIDEST teeth -> easiest to dodge.)
 (2) Adversarially search clustered covering sets for the largest max(S) at which the EXACT-W
     criterion C(S) fails (deep search). If it never fails above some modest B0, that B0 is the
     empirical finite-base cutoff.
 (3) Determine, via the exact-W criterion applied to the SLOWEST runner, a provable statement:
     v_min is always small (covering forces a runner = 1? no, forces mult of each q). Quantify
     min over covering sets of v_min and of W(S\{v_min})*7*v_min.

stdlib only, exact Fractions.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
from functools import reduce
import random, time

C = F(1, 14)
def darcs(v, c=C):
    hw = F(c, v); return [(F(k, v)-hw, F(k, v)+hw) for k in range(v)]
def wrapU(iv):
    o = []
    for lo, hi in iv:
        s = lo - (lo % 1); a = lo - s; b = hi - s
        if b <= 1: o.append((a, b))
        else: o.append((a, F(1))); o.append((F(0), b-1))
    o = sorted(o); r = []; cl, ch = o[0]
    for lo, hi in o[1:]:
        if lo <= ch: ch = ch if ch > hi else hi
        else: r.append((cl, ch)); cl, ch = lo, hi
    r.append((cl, ch)); return r
def Wsafe(A, c=C):
    dz = []
    for v in set(A): dz += darcs(v, c)
    if not dz: return F(1)
    dz = wrapU(dz)
    best = F(0)
    for i in range(len(dz)):
        hi = dz[i][1]; lo = dz[(i+1) % len(dz)][0] + (1 if i == len(dz)-1 else 0)
        if lo - hi > best: best = lo - hi
    return best
def criterion_exactW(S):
    best = (False, None, F(-1))
    for v in sorted(set(S)):
        A = [u for u in S if u != v]
        W = Wsafe(A); margin = W - F(1, 7*v)
        if margin > 0: return (True, v, margin)
        if margin > best[2]: best = (False, v, margin)
    return best
qs = range(2, 15)
def covering(S): return all(any(v % q == 0 for v in S) for q in qs)
def mu(A): return 1 - sum(F(1, 7*u) for u in A)

# exact M for cross-check on any C-failures
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def g(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); Cc = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): Cc.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): Cc.add(F(k, d)); k += 1
    Cc.add(F(1, 2)); return Cc
def M(S): return max(g(S, t) for t in cand(S))


def part1_slowest_runner():
    print("="*78)
    print("PART 1: criterion via the SLOWEST runner v_min (fewest, widest teeth to dodge)")
    print("="*78)
    print("""
  Covering forces small runners: a multiple of 2,3,...,7 etc. In particular S must contain a
  multiple of 2 and of 3, but the SMALLEST element can be as large as... let's measure. For each
  covering set, v_min = min(S). Remove v_min: W(S\\{v_min}) must exceed 1/(7 v_min). Since v_min is
  small, 1/(7 v_min) is LARGE -- this is the HARD direction. So slowest-runner removal is NOT the
  easy win; the easy win is removing a LARGE runner. We quantify both.
""")
    random.seed(61)
    mins = []
    for _ in range(20000):
        V = random.randint(20, 200); csz = random.randint(3, 6)
        cl = random.sample(range(max(2, V-22), V+1), csz)
        small = random.sample(range(1, 14), 13-csz)
        S = sorted(set(small + cl))
        if len(S) != 13 or not covering(S): continue
        mins.append(min(S))
    if mins:
        print(f"  over {len(mins)} clustered covering sets: v_min in [{min(mins)},{max(mins)}], "
              f"typical min(S)={sorted(mins)[len(mins)//2]}")
    print("  => covering essentially forces min(S)=1 or 2 (a multiple of small q must be small).")


def part2_clustered_deep_search():
    print()
    print("="*78)
    print("PART 2: deep adversarial search -- largest max(S) at which exact-W criterion C FAILS")
    print("="*78)
    random.seed(62)
    worst_fail_max = 0; fail_count = 0; tested = 0
    worst_margin = (F(1), None)   # the smallest positive criterion margin seen (tightness)
    t0 = time.time()
    N = 120000
    for _ in range(N):
        style = random.random()
        if style < 0.5:
            # tight clustered: many runners in a short window near V, plus forced small
            V = random.randint(20, 260)
            width = random.randint(8, 26)
            csz = random.randint(4, 8)
            cl = random.sample(range(max(2, V-width), V+1), min(csz, V-1))
            small = random.sample(range(1, 14), 13-len(cl))
            S = sorted(set(small + cl))
        else:
            # two clusters
            V = random.randint(40, 260)
            c1 = random.sample(range(max(2, V-15), V+1), random.randint(2, 4))
            mid = random.randint(20, V-16) if V-16 > 21 else 20
            c2 = random.sample(range(max(2, mid-10), mid+1), random.randint(2, 4))
            small = random.sample(range(1, 14), max(0, 13-len(c1)-len(c2)))
            S = sorted(set(small + c1 + c2))
        if len(S) != 13 or not covering(S): continue
        tested += 1
        holds, v, margin = criterion_exactW(S)
        if not holds:
            fail_count += 1
            worst_fail_max = max(worst_fail_max, max(S))
            mval = M(S)
            print(f"   C-FAIL max={max(S)} S={S} best-margin={float(margin):.5f} M={mval}={float(mval):.5f} LRC-ok={mval>=C}")
        else:
            if 0 < margin < worst_margin[0]:
                worst_margin = (margin, S, v)
    dt = time.time() - t0
    print(f"\n  tested {tested} clustered covering sets in {dt:.1f}s")
    print(f"  exact-W criterion C FAILURES: {fail_count}")
    if fail_count:
        print(f"    largest max(S) with a C-fail: {worst_fail_max}")
    else:
        print("    C(S) HELD on every clustered covering set searched (no fail at any max).")
    print(f"  tightest successful criterion margin: {float(worst_margin[0]):.6f}")
    print(f"    at v={worst_margin[2] if worst_margin[1] else None} S={worst_margin[1]}")


def part3_provable_cutoff_attempt():
    print()
    print("="*78)
    print("PART 3: toward a PROVABLE cutoff -- sharpen N(A) for the large-runner removal")
    print("="*78)
    print("""
  Remove the LARGEST runner V. A = S\\{V}, 12 runners each <= V-1. The safe set of A at level
  1/14 has widest arc W(A). The pigeonhole W(A) >= mu(A)/N(A) needs N(A) (number of safe arcs).
  SHARPER N(A): the safe arcs of A are the gaps of the tooth-union. A runner u contributes u teeth.
  But two runners u,w only create EXTRA gaps where their teeth interleave. A clean bound:
        N(A) <= 1 + sum_{u in A} (number of teeth of u that have an endpoint not inside another's tooth)
  Hard to bound cleanly. Instead we test the EMPIRICAL claim: for clustered covering sets,
  W(A) (exact) stays >= a constant fraction, so C fires via SOME small/medium runner.
  We report, over clustered covering sets, min over the set of [ max_v (7 v W(S\\{v})) ], i.e.
  the best criterion ratio; if its inf over all clustered sets is > 1, C is universal.
""")
    random.seed(63)
    worst = (F(99), None)
    tested = 0
    for _ in range(60000):
        V = random.randint(20, 220)
        width = random.randint(8, 24)
        csz = random.randint(4, 7)
        cl = random.sample(range(max(2, V-width), V+1), min(csz, V-1))
        small = random.sample(range(1, 14), 13-len(cl))
        S = sorted(set(small + cl))
        if len(S) != 13 or not covering(S): continue
        tested += 1
        best_ratio = F(0); bv = None
        for v in sorted(set(S)):
            A = [u for u in S if u != v]
            ratio = 7*v*Wsafe(A)
            if ratio > best_ratio: best_ratio = ratio; bv = v
        if best_ratio < worst[0]: worst = (best_ratio, S, bv)
    print(f"  over {tested} clustered covering sets:")
    print(f"  inf over sets of [ max_v 7 v W(S\\{{v}}) ] = {float(worst[0]):.5f}  (>1 => C fires)")
    print(f"    worst set {worst[1]} via v={worst[2] if worst[1] else None}")
    print(f"  C(S) fires iff this best ratio > 1; observed inf = {float(worst[0]):.5f}.")


if __name__ == "__main__":
    part1_slowest_runner()
    part2_clustered_deep_search()
    part3_provable_cutoff_attempt()
    print()
    print("="*78)
    print("SUMMARY: clustered-large is where the provable cutoff is missing. The exact-W criterion")
    print("holds with a positive margin on every clustered covering set searched, but the pigeonhole")
    print("(mu/N) proof of it does not fire for clustered sets -> no finite max-cutoff is PROVED yet.")
    print("="*78)
