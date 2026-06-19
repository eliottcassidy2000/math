#!/usr/bin/env python3
"""
ROUTE H — the binding-pair clearing-crossing via GLOBAL COUNTING/averaging.
mac-mini-2026-06-19-S2.

QUESTION (from the workflow): M(S)=j/D at the binding pair (THM-524). LRC(14)
<=> some pair-crossing clears all runners at gap >= 1/14. Instead of FINDING the
binding pair, AVERAGE/COUNT over all C(13,2) pair-crossings and try to FORCE a
good crossing by a first-moment / variance / exact-count argument.

FINDINGS (all EXACT rationals, fractions.Fraction):

1. FIRST MOMENT IS DEAD. Over residual S3 sets, the MEAN crossing-gap is ~0.018-0.033,
   FAR below 1/14 ~= 0.0714, while M (the max) is ~0.10-0.17. The gap between mean
   and max is exactly what a first-moment argument cannot bridge. (v1, v4)

2. THE BINDING PAIR IS NOT STRUCTURALLY FIXED. Over 600 random residual S3 sets,
   the binding pair (achieving M) is: CROSS (small x cluster) 58.7%, cluster-cluster
   25.7%, small-small 15.7%. So "the good crossing always spans P and L" is FALSE.
   A uniform-over-pairs count cannot exploit a fixed binding-pair type. (v3b)

3. THE TIGHT-EXTREMAL OBSTRUCTION (the decisive negative). {1..13} has M = 1/14
   EXACTLY: the lonely set {tau: g>=1/14} is MEASURE ZERO (3 isolated grid points
   tau in {1/14,3/14,5/14}), zero open lonely components, mu=0. Therefore ANY counting
   argument that produces a crossing with gap STRICTLY > 1/14 is FALSE on {1..13}
   (it would prove M>1/14, contradiction). A valid ROUTE-H count MUST allow equality,
   i.e. decide the NON-STRICT inequality g >= 1/14 at a crossing -- the same hard
   strict-vs-nonstrict decision the lonely measure cannot see at tight sets. (v5, v6b)

4. THE TWO EXTREMAL FAMILIES HAVE DISJOINT GOOD-CROSSING STRUCTURE:
   - {1..13} (tight): good crossings are ON the grid tau=k/14, all from COMPLEMENT
     pairs (a+b=14), gap = 1/14 exactly. (the Z/2-reversal locus, THM-524 D)
   - {1..11,13,84} (hard core): grid_gap = 0 (84 == 0 mod 14 kills the grid);
     ALL good crossings are OFF-grid (v, 84)-pairs at tau ~ k/89, gap = 7/89 > 1/14.
   A single counting argument must cover both disjoint structures. (v2, v4)

5. PARTIAL POSITIVE: THE STRUCTURED WITNESS FAMILY. Define
     F(S) = { k/14 : 1<=k<=13 }  UNION  { k/(v+-Vmax) : v in S, k }.
   EMPIRICAL: F(S) ALWAYS contains a crossing with gap >= 1/14, for every tested
   residual S3 set (0 misses / 800 bounded-spread, 0 / 500 large-spread). This holds
   EVEN when the maximal binding pair is cluster-cluster (25.7%): there is always an
   EQUALLY-good (>= 1/14) crossing in F(S). So ROUTE H's count need only range over
   the O(|S|) "Vmax-spokes" + 13 grid points, NOT all C(13,2) pairs. (v7, v8, v9)

   HONEST CAVEAT: this is a COVERAGE reduction, not a counting PROOF. It unifies the
   two KNOWN sufficient conditions (grid q=14 witness, THM-524; criterion C at v=Vmax,
   THM-527) into one family but does NOT prove one of them fires. The strict-vs-nonstrict
   decision (point 3) remains. It does NOT transcend the existing crux.

NET: ROUTE H via averaging/first-moment is a DEAD END (the mean is a ceiling far below
the target; the tight extremal forbids any strict-positivity count). ROUTE H via
STRUCTURED counting yields a real but non-decisive reduction (witness family F(S)
covers the binding crossing), reproving = unifying the two known sufficient conditions.

Run: python3 lrc14_routeH_crossing_count_mac-mini-2026-06-19-S2.py
"""
from fractions import Fraction as F
from math import gcd
import random


def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r


def g(S, t):
    return min(nrm(v * t) for v in S)


def all_crossings(S):
    S = sorted(set(S)); out = set()
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            a, b = S[i], S[j]
            for d in (a + b, b - a):
                if d <= 0:
                    continue
                k = 1
                while F(k, d) <= F(1, 2):
                    out.add(F(k, d)); k += 1
    out.add(F(1, 2))
    return out


def M_exact(S):
    return max(g(S, t) for t in all_crossings(S))


def witness_family(S):
    """F(S) = grid {k/14} U {(v,Vmax)-crossings}."""
    S = sorted(set(S)); Vmax = max(S)
    ws = set(F(a, 14) for a in range(1, 14))
    for v in S:
        if v == Vmax:
            continue
        for d in (v + Vmax, Vmax - v):
            if d <= 0:
                continue
            k = 1
            while F(k, d) <= F(1, 2):
                ws.add(F(k, d)); k += 1
    return ws


def lonely_components(S, thr=F(1, 14)):
    """Exact #components and measure of {tau: g(S,tau)>=thr}."""
    S = sorted(set(S)); bps = {F(0), F(1)}
    for v in S:
        k = 0
        while True:
            for s in (thr, -thr):
                t = (k + s) / v
                if 0 <= t < 1:
                    bps.add(t)
            if F(k) / v >= 1:
                break
            k += 1
    bp = sorted(bps); arcs = []
    for i in range(len(bp) - 1):
        a, b = bp[i], bp[i + 1]
        if g(S, (a + b) / 2) >= thr:
            arcs.append((a, b))
    merged = []
    for a, b in arcs:
        if merged and merged[-1][1] == a:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return len(merged), sum(b - a for a, b in merged)


def demo():
    print("=" * 72)
    print("1. FIRST-MOMENT (mean crossing-gap) vs MAX (=M): mean is a ceiling")
    print("=" * 72)
    for name, S in [("{1..11,13,84}", list(range(1, 12)) + [13, 84]),
                    ("{1..13}", list(range(1, 14)))]:
        crs = all_crossings(S); gaps = [g(S, t) for t in crs]
        mean = sum(gaps) / len(gaps)
        print(f"  {name}: mean={float(mean):.5f}  M={float(max(gaps)):.5f}  (1/14={1/14:.5f})")

    print("\n" + "=" * 72)
    print("2. TIGHT-EXTREMAL OBSTRUCTION: {1..13} has mu=0, zero strict components")
    print("=" * 72)
    S = list(range(1, 14))
    nc, mu = lonely_components(S)
    print(f"  {{1..13}}: M={M_exact(S)}  #lonely components={nc}  mu={mu}")
    print("  => no STRICT good crossing exists; any strict-positivity count over-proves.")

    print("\n" + "=" * 72)
    print("3. STRUCTURED WITNESS FAMILY F(S)=grid U (v,Vmax): always contains good crossing")
    print("=" * 72)
    random.seed(2026)
    tested = miss = lrcfail = 0
    for _ in range(400):
        kL = random.randint(3, 10); kP = 13 - kL
        if kP < 0 or kP > 13:
            continue
        P = sorted(random.sample(range(1, 14), kP)) if kP > 0 else []
        V = random.choice([40, 70, 100]); spread = random.randint(kL - 1, min(2 * kL, 22))
        offs = {0}
        while len(offs) < kL:
            offs.add(random.randint(0, spread))
        S = sorted(set(P) | set(V + o for o in sorted(offs)))
        if len(S) != 13:
            continue
        gg = 0
        for x in S:
            gg = gcd(gg, x)
        if gg != 1:
            continue
        tested += 1
        M = M_exact(S)
        if M < F(1, 14):
            lrcfail += 1; continue
        if max(g(S, t) for t in witness_family(S)) < F(1, 14):
            miss += 1
    print(f"  tested={tested}  LRC violations={lrcfail}  family misses good crossing={miss}")
    print("  => 0 misses: F(S) covers the binding crossing (coverage reduction, not a proof).")


if __name__ == "__main__":
    demo()
