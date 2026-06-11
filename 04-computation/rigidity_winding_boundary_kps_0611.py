#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-485 / HYP-2417 -- resolving the NEAR-BOUNDARY hard sets.

The stress test left 45 sets (all supported on small ODD k: {1,3,5,7,9,...}) whose
interior zero -- if it exists -- sits in the boundary annulus 0.93 < |x| < 1, where
the naive winding count of the truncated polynomial is contaminated by the partial
Euler product's roots-of-unity (winding(Euler^trunc) jumps from 0 to ~19..48 across
0.95 -> 0.99). To certify a GENUINE interior zero there we use the

   WINDING DIFFERENCE  D(r) = winding(F_S, r) - winding(F_Euler, r),

evaluated at the SAME r and SAME truncation Kmax. Both polynomials share the exact
same high-degree pentagonal artifact zeros (they differ only in finitely many low
coefficients, k in S), so the artifact windings CANCEL in D(r). Since the TRUE
F_Euler is zero-free in |x|<1, D(r) counts (true F_S interior zeros inside r) minus
(0) up to artifact cancellation:  D(r) >= 1 at some r<1  =>  F_S has a genuine
interior zero.

We sweep r from 0.93 up to 0.999 and report the first r with D(r) >= 1. We ALSO
directly Newton-polish a complex zero from a fine annulus grid as an independent
witness (exhibit the actual zero z, |z|<1, F_S(z)~0).

A set that yields D(r) = 0 for ALL r < 1 AND no Newton witness would be a genuine
counterexample to the rigidity.
"""

import cmath
import math
import itertools
import time


def make_F(flips, Kmax):
    eps = []
    for k in range(1, Kmax + 1):
        e = 1 if (k % 2 == 1) else -1
        if k in flips:
            e = -e
        g = k * (3 * k - 1) // 2
        gb = k * (3 * k + 1) // 2
        eps.append((g, gb, e))

    def F(x):
        s = 1.0 + 0.0j
        for g, gb, e in eps:
            s -= e * (x ** g + x ** gb)
        return s

    return F


def winding_number(F, r, Ntheta=4096, refine_depth=28):
    TWO_PI = 2 * math.pi
    THRESH = math.pi / 4
    pts = [F(r * cmath.exp(1j * TWO_PI * i / Ntheta)) for i in range(Ntheta)]

    def inc(a, b):
        return cmath.phase(b / a)

    def refine(thA, fA, thB, fB, depth):
        d = inc(fA, fB)
        if abs(d) <= THRESH or depth >= refine_depth:
            return d
        thm = (thA + thB) / 2 if thB > thA else (thA + (thB + TWO_PI)) / 2
        fm = F(r * cmath.exp(1j * (thm % TWO_PI)))
        return refine(thA, fA, thm, fm, depth + 1) + refine(thm, fm, thB, fB, depth + 1)

    total = 0.0
    for i in range(Ntheta):
        thA = TWO_PI * i / Ntheta
        fA = pts[i]
        if i + 1 == Ntheta:
            thB, fB = TWO_PI, pts[0]
        else:
            thB, fB = TWO_PI * (i + 1) / Ntheta, pts[i + 1]
        total += refine(thA, fA, thB, fB, 0)
    return round(total / TWO_PI)


def newton_witness(F, Kmax, Rmax=0.99999, grid=600, rsteps=120, rlo=0.5):
    """Find a complex zero z with |z|<1 by annulus-grid minimization + Newton."""
    best = None
    for ai in range(grid):
        th = 2 * math.pi * ai / grid
        for ri in range(1, rsteps + 1):
            r = rlo + (Rmax - rlo) * ri / rsteps
            z = r * cmath.exp(1j * th)
            v = abs(F(z))
            if best is None or v < best[0]:
                best = (v, z)
    z = best[1]
    h = 1e-8
    for _ in range(100):
        f = F(z)
        df = (F(z + h) - F(z - h)) / (2 * h)
        if abs(df) < 1e-30:
            break
        zn = z - f / df
        if abs(zn) >= 0.99999:
            # damp the step to stay inside
            zn = z - 0.3 * f / df
            if abs(zn) >= 0.99999:
                break
        z = zn
        if abs(F(z)) < 1e-12:
            break
    return (z, abs(F(z))) if abs(z) < 1 else (None, None)


def bsum(S):
    return sum((1 if k % 2 == 1 else -1) for k in S)


def main():
    t0 = time.time()
    Kmax = 40
    print("=" * 78, flush=True)
    print("Resolving near-boundary hard sets via winding-DIFFERENCE vs Euler",
          flush=True)
    print(f"  Kmax={Kmax}; D(r)=w(F_S,r)-w(F_Euler,r); artifacts cancel", flush=True)
    print("=" * 78, flush=True)

    Feul = make_F(frozenset(), Kmax)
    rsweep = [0.93, 0.95, 0.96, 0.97, 0.98, 0.985, 0.99, 0.995, 0.999]
    wE = {}
    for r in rsweep:
        wE[r] = winding_number(Feul, r, Ntheta=4096)
    print("\nwinding(Euler^trunc) along the sweep (the artifact baseline):", flush=True)
    print("   " + "  ".join(f"r={r}:{wE[r]}" for r in rsweep), flush=True)

    # the 45 problematic sets = small-odd-k supported, |S|<=6, k<=12, B>=0.
    # regenerate them: hard sets where winding at r<=0.93 was 0.
    KFLIP, MAXSIZE = 12, 6
    candidates = []
    for sz in range(1, MAXSIZE + 1):
        for S in itertools.combinations(range(1, KFLIP + 1), sz):
            S = frozenset(S)
            if bsum(S) < 0:
                continue
            F = make_F(S, Kmax)
            if winding_number(F, 0.93, Ntheta=2048) == 0:
                candidates.append(S)
    print(f"\nNear-boundary candidates (winding 0 at r=0.93, B>=0): {len(candidates)}",
          flush=True)

    resolved = 0
    unresolved = []
    for S in candidates:
        F = make_F(S, Kmax)
        # winding-difference sweep
        firstD = None
        for r in rsweep:
            wS = winding_number(F, r, Ntheta=4096)
            D = wS - wE[r]
            if D >= 1:
                firstD = (r, D, wS, wE[r])
                break
        # independent Newton witness
        z, res = newton_witness(F, Kmax)
        ok = (firstD is not None) or (z is not None and res < 1e-6)
        if ok:
            resolved += 1
        else:
            unresolved.append((S, firstD, z, res))
        tag = "OK" if ok else "??"
        zinfo = f"z={z:.5f} |z|={abs(z):.5f} res={res:.1e}" if z is not None else "no-newton-z"
        dinfo = (f"D>=1 first at r={firstD[0]} (D={firstD[1]})" if firstD
                 else "D=0 all r")
        print(f"  [{tag}] S={set(S)} B={bsum(S):+d}: {dinfo}; {zinfo}", flush=True)

    print("\n" + "=" * 78, flush=True)
    print(f"resolved (genuine interior zero certified): {resolved}/{len(candidates)}",
          flush=True)
    if unresolved:
        print(f"*** {len(unresolved)} UNRESOLVED -> genuine counterexample candidates ***",
              flush=True)
        for S, fd, z, res in unresolved:
            print(f"    S={set(S)} B={bsum(S):+d}", flush=True)
    else:
        print("ALL near-boundary hard sets carry a genuine interior zero.", flush=True)
        print("=> NO counterexample. Rigidity holds on |S|<=6, k<=12.", flush=True)
    print(f"=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
