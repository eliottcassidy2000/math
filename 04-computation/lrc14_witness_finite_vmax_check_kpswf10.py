#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_witness_finite_vmax_check_kpswf10.py  (kind-pasteur 2026-06-21, THREAD B honesty)

HONEST finite-Vmax validation of the WITNESS route (maxgap>1/7) vs CRITERION route
(maxgap>2/7).  THM-527 proves CRITERION C at v=Vmax gives M>=1/14 via a finite-w0
argument.  THREAD B's claim: the WEAKER global-witness threshold (maxgap>1/7) ALSO
certifies M>=1/14, more directly.

  GLOBAL WITNESS mechanism (THM-526 LEMMA-A / direct): a point tau* is safe for ALL of
  S=P u L iff (i) tau* in G_P AND (ii) ||V0 tau*|| in (1/14,13/14) leaves the cluster teeth
  a fast-phase slot.  Concretely: at slow-time near a good x, set theta=frac(V0 tau) in the
  largest cluster-phase gap (which has width maxgap>1/7); the runner u=V0-e_i is safe iff
  theta avoids each cluster tooth -- the fast window has width w_theta=maxgap-1/7>0.

  CRITERION (via-max) needs the window to also clear the 1/14-collar on BOTH sides of the
  fast cell => width must exceed 1/7 => maxgap>2/7.

CONSTRUCT actual integer covering-type sets S = {V0 - e_i} u P-as-runners and DIRECTLY
compute the EXACT M(S)=max_tau min_v ||v tau||, then check:
  - does maxgap>2/7 region (criterion) deliver M>=1/14 ?  (THM-527 -- expect yes)
  - does the WITNESS construction (pick tau in G_P at a maxgap>1/7 slow-cell, fast theta in
    the middle of the cluster gap) deliver a point with min_v||v tau||>=1/14 ?
We verify the witness POINT is genuinely safe (min over ALL runners >= 1/14) at FINITE Vmax.

This is the load-bearing check: the witness route is only useful if maxgap>1/7 at finite
Vmax yields a REAL safe tau*, not just a slow-limit artifact.
"""
from fractions import Fraction as Fr

HALF = Fr(1, 14)
T1 = Fr(1, 7)
T2 = Fr(2, 7)


def norm(y):
    r = y % 1
    return min(r, 1 - r)


def M_exact(S):
    """EXACT M(S) = max_tau min_{v in S} ||v tau||, over the finite breakpoint set.
       The max of the min is attained at a breakpoint tau=j/(2 v) style; we scan the
       grid of all ||v tau|| extrema = {j/v} and midpoints, then take exact max of the
       (piecewise-linear) min.  For exact value we evaluate min at every candidate and at
       cell midpoints; M is the max over candidates of the local min, but since min of
       sawtooths is piecewise-linear with breakpoints at j/v and crossings, the maximum is
       attained at a crossing of two sawtooth 'tent' edges.  We compute via fine candidate
       set: all (2j+1)/(2v) (tent peaks) plus pairwise crossings."""
    S = [int(v) for v in S]
    cands = set()
    for v in S:
        for j in range(0, v):
            cands.add(Fr(2 * j + 1, 2 * v))     # peak of v's tent (||v tau||=1/2)
            cands.add(Fr(j, v))                 # zero of v's tent
    # pairwise crossings of tents: ||a tau||=||b tau|| -> a tau +- b tau = integer
    for a in S:
        for b in S:
            if a >= b:
                continue
            for sgn in (1, -1):
                denom = a + sgn * b
                if denom == 0:
                    continue
                # a tau = k - sgn*b tau + ... handle a tau - sgn b tau = integer m? Use:
                # ||a t||=||b t|| at t where a t = m1 +- (b t) ... crossings at t=m/(a+-b)
                for m in range(0, abs(denom) + 1):
                    t = Fr(m, abs(denom))
                    if Fr(0) <= t <= Fr(1):
                        cands.add(t)
    best = Fr(0)
    for t in cands:
        if t < 0 or t > 1:
            continue
        mn = min(norm(v * t) for v in S)
        if mn > best:
            best = mn
    return best


def witness_point(Pset, E, V0):
    """Construct the WITNESS tau* per the slow-fast recipe and return (tau*, min_v||v tau*||).
       Recipe: find a slow-cell x in G_P with maxgap of {frac(e_i x)} > 1/7; put the fast
       phase theta in the CENTER of that largest cluster gap; tau* = (n + theta)/V0 for the
       integer n that lands tau* in the G_P slow-cell.  We search the natural grid for x."""
    Elist = [int(e) for e in E]
    # candidate slow-cells: scan x on a fine grid; pick the one in G_P with largest maxgap
    # grid from sector + danger boundaries
    bp = {Fr(0), Fr(1)}
    for e in Elist + list(Pset):
        e = int(e)
        if e == 0:
            continue
        for t in range(0, 7 * e + 1):
            bp.add(Fr(t, 7 * e))
        for k in range(0, e + 1):
            for s in (HALF, -HALF):
                bp.add((Fr(k) + s) / e)
    xs = sorted(b for b in bp if Fr(0) <= b <= Fr(1))
    best = None
    for a, b in zip(xs, xs[1:]):
        if b <= a:
            continue
        x = (a + b) / 2
        # in G_P?
        if any(norm(int(p) * x) < HALF for p in Pset):
            continue
        ph = sorted((e * x) % 1 for e in Elist)
        # maxgap and its location
        if len(ph) <= 1:
            continue
        gaps = []
        for i in range(len(ph)):
            lo = ph[i]
            hi = ph[(i + 1) % len(ph)] + (1 if i + 1 == len(ph) else 0)
            gaps.append((hi - lo, lo, hi))
        g, lo, hi = max(gaps)
        if g > T1:
            if best is None or g > best[0]:
                center = ((lo + hi) / 2) % 1
                best = (g, x, center)
    if best is None:
        return None
    g, x, theta = best
    # tau* with frac(V0 tau*) = theta and tau* near x: tau* = (round(V0 x - theta) + theta)/V0
    n = round(float(V0 * x - theta))
    tau = (Fr(n) + theta) / V0
    return tau, x, g, theta


def main():
    print("=" * 100)
    print("FINITE-Vmax: does the WITNESS threshold (maxgap>1/7) yield a REAL safe tau* (min||v t||>=1/14)?")
    print("=" * 100)
    print("Build S = {V0 - e_i : e_i in E} u P, for moderate V0, exact M(S), and a witness point.")
    print()

    # admissible-ish configs: P u cluster L={V0-e_i}; pick P, E (consec), V0 large.
    tests = [
        ("P={1,2,3}",      [1, 2, 3],       [0, 1, 2, 3, 4],      [200, 400, 803]),
        ("P={1,2,3}",      [1, 2, 3],       list(range(7)),       [300, 700]),
        ("P={1,2,3,12}",   [1, 2, 3, 12],   list(range(9)),       [500, 1001]),
        ("P={1..5}",       [1, 2, 3, 4, 5], list(range(8)),       [400, 900]),
        ("P={1,2,3} perf", [1, 2, 3],       [0, 2, 3, 4, 5, 6, 8],[400, 803]),
    ]
    hdr = f"{'config':<18}{'V0':>6}{'k':>4}{'M(S)':>10}{'M>=1/14?':>10}{'witness min||v t*||':>22}{'wit>=1/14?':>12}"
    print(hdr)
    print("-" * len(hdr))
    all_M = True
    all_wit = True
    for name, Pset, E, V0s in tests:
        for V0 in V0s:
            L = [V0 - e for e in E]
            S = sorted(set(list(Pset) + L))
            # ensure primitive-ish and the cluster is all > 13 (large): require V0-max(E) > 13
            if min(L) <= 13:
                continue
            M = M_exact(S)
            Mok = M >= HALF
            all_M &= Mok
            wr = witness_point(Pset, E, V0)
            if wr is None:
                wstr = "no maxgap>1/7 cell"
                witval = None
                witok = False
            else:
                tau, x, g, theta = wr
                witval = min(norm(v * tau) for v in S)
                witok = witval >= HALF
            all_wit &= witok
            wvs = f"{float(witval):.5f}" if witval is not None else "n/a"
            print(f"{name:<18}{V0:>6}{len(E):>4}{float(M):>10.5f}{str(Mok):>10}{wvs:>22}{str(witok):>12}")
    print("-" * len(hdr))
    print(f"\nM(S)>=1/14 for all built sets:           {all_M}")
    print(f"witness point safe (>=1/14) for all:     {all_wit}")
    print("""
READING:
  * If the witness column is >=1/14 everywhere, the maxgap>1/7 (1/7-scale) event delivers a
    GENUINE finite-Vmax safe point -- the witness route is sound, not a slow-limit artifact.
  * Then BOTH thresholds certify LRC: the criterion (2/7, THM-527) AND the strictly-weaker
    witness (1/7).  The witness floor being larger (THREAD-B: 0.0565 vs 1/84) means the 1/7
    route is the EASIER sufficient condition.  The two routes are PARALLEL (same phi-phase,
    thresholds 1/7 vs 2/7), and the 1/7 one is the weaker/easier of the two.
""")


if __name__ == "__main__":
    main()
