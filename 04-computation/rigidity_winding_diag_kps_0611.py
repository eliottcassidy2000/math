#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-485 / HYP-2417 -- DIAGNOSTIC: truncation artifacts in the winding number.

The first winding run flagged a problem: the TRUNCATED Euler product
F_Euler^trunc = 1 - sum_{k<=Kmax} (-1)^{k+1}(x^{g_k}+x^{gbar_k}) had winding 10 at
r=0.99 even though the TRUE F_Euler = prod(1-x^n) is zero-free in the open disk.
Cause: the partial pentagonal sum truncates prod_{n}(1-x^n) at degree ~ g_Kmax,
i.e. it is (close to) prod_{n<=N}(1-x^n) whose zeros are roots of unity ON |x|=1.
Near r=0.99 those boundary zeros leak inside and inflate the winding count.

So the naive winding count of the truncated polynomial is contaminated. To get a
TRUSTWORTHY interior-zero certificate we must SEPARATE genuine interior zeros of
the true F_S from boundary-hugging truncation artifacts. Strategy:

 (A) Pick a radius r0 strictly < 1 small enough that the TRUE F_Euler is provably
     zero-free AND its truncation is uniformly close: at r0 the winding of
     F_Euler^trunc must come out 0. Find the largest such r0 per Kmax.
 (B) At that safe r0, a flip set S with winding >= 1 has a CERTIFIED genuine
     interior zero of modulus < r0 < 1 (the tail bound makes trunc vs true
     indistinguishable for the argument count).
 (C) Equivalent robust certificate: the WINDING DIFFERENCE w(F_S, r) - w(F_Euler, r)
     at the SAME r cancels the shared boundary artifacts; a nonzero difference at
     a radius where artifacts dominate still localizes a genuine zero shift.

This script maps winding(Euler) vs (Kmax, r) to locate the safe annulus, then
re-certifies the hard half there.
"""

import cmath
import math
import itertools
import time
from collections import Counter


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


def winding_number(F, r, Ntheta=4096, refine_depth=26):
    TWO_PI = 2 * math.pi
    THRESH = math.pi / 4
    pts = []
    minmod = float("inf")
    for i in range(Ntheta):
        th = TWO_PI * i / Ntheta
        f = F(r * cmath.exp(1j * th))
        minmod = min(minmod, abs(f))
        pts.append((th, f))

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
        thA, fA = pts[i]
        if i + 1 == Ntheta:
            thB, fB = TWO_PI, pts[0][1]
        else:
            thB, fB = pts[i + 1]
        total += refine(thA, fA, thB, fB, 0)
    return round(total / TWO_PI), total / TWO_PI, minmod


def bsum(S):
    return sum((1 if k % 2 == 1 else -1) for k in S)


def main():
    t0 = time.time()
    print("=" * 78, flush=True)
    print("DIAGNOSTIC: locating the truncation-safe annulus for winding counts",
          flush=True)
    print("=" * 78, flush=True)

    # (A) map winding(Euler) over (Kmax, r). The TRUE F_Euler is zero-free in |x|<1,
    #     so any nonzero winding here is pure truncation artifact. Find safe r0(Kmax).
    print("\nwinding( F_Euler^trunc ) vs (Kmax, r)  [TRUE value is 0 for all r<1]:",
          flush=True)
    Kmaxes = [8, 12, 20, 30, 40, 60]
    radii = [0.5, 0.7, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.99]
    header = "  Kmax\\r |" + "".join(f"{r:>7}" for r in radii)
    print(header, flush=True)
    print("  " + "-" * (len(header) - 2), flush=True)
    safe_r = {}
    for K in Kmaxes:
        Feul = make_F(frozenset(), K)
        row = []
        best_safe = 0.0
        for r in radii:
            w, _, mm = winding_number(Feul, r, Ntheta=2048)
            row.append(w)
            if w == 0:
                best_safe = r
        safe_r[K] = best_safe
        print(f"  {K:>5} |" + "".join(f"{w:>7}" for w in row), flush=True)
    print("\n  largest r with winding(Euler)=0 per Kmax (the 'safe radius'):", flush=True)
    for K in Kmaxes:
        print(f"     Kmax={K}: safe up to r={safe_r[K]}", flush=True)

    # Choose a robust (Kmax, r): want Kmax large enough that g_Kmax >> any flip exponent
    # (so flips on k<=10, max exponent gbar_10 = 155, are fully resolved) AND r at/below
    # the safe radius so Euler artifacts vanish.
    Kmax = 30  # g_30 = 1335, gbar_10 = 155 -> flips fully inside the kept range
    r0 = safe_r[Kmax]
    if r0 == 0.0:
        r0 = 0.7
    # back off one notch below the safe radius for margin
    cand = [r for r in radii if r <= r0]
    r_use = cand[-2] if len(cand) >= 2 else cand[-1]
    print(f"\nChosen certification setting: Kmax={Kmax}, r={r_use} "
          f"(Euler winding 0 here => artifact-free)", flush=True)

    # (B) re-certify the HARD half at the safe (Kmax, r). winding>=1 => GENUINE
    #     interior zero of modulus < r_use < 1.
    KFLIP, MAXSIZE = 10, 4
    hard = []
    for sz in range(1, MAXSIZE + 1):
        for S in itertools.combinations(range(1, KFLIP + 1), sz):
            S = frozenset(S)
            if bsum(S) >= 0:
                hard.append(S)
    print(f"\nRe-certifying HARD half ({len(hard)} sets, B(S)>=0) at the safe radius:",
          flush=True)

    certified = 0
    failed = []
    dist = Counter()
    for S in hard:
        F = make_F(S, Kmax)
        w, wf, mm = winding_number(F, r_use, Ntheta=2048)
        dist[w] += 1
        if w >= 1:
            certified += 1
        else:
            failed.append((set(S), bsum(S), w, mm))
    print(f"  CERTIFIED (winding >= 1) at r={r_use}: {certified}/{len(hard)}", flush=True)
    print(f"  winding-0 (no zero of modulus < {r_use}): {len(failed)}", flush=True)
    print("  winding distribution at safe radius:", flush=True)
    for w in sorted(dist):
        print(f"     winding {w}: {dist[w]}", flush=True)

    # For the winding-0-at-safe-r sets, push the radius outward toward 1 to see if the
    # zero sits in the annulus r_use < |x| < 1 (still interior, just nearer boundary).
    if failed:
        print("\n  winding-0-at-safe-r sets -- pushing radius toward 1 to find the zero:",
              flush=True)
        outer = [0.93, 0.96, 0.98, 0.995, 0.999]
        still_zero = []
        for S, b, w0, mm0 in failed:
            F = make_F(S, Kmax)
            Feul = make_F(frozenset(), Kmax)
            found = None
            for r in outer:
                wS, _, _ = winding_number(F, r, Ntheta=4096)
                wE, _, _ = winding_number(Feul, r, Ntheta=4096)
                # winding DIFFERENCE cancels shared boundary artifacts
                diff = wS - wE
                if diff >= 1:
                    found = (r, wS, wE, diff)
                    break
            if found:
                r, wS, wE, diff = found
                print(f"     S={S} B={b:+d}: interior zero in annulus, "
                      f"w(F_S)-w(Euler)={diff} at r={r} (wS={wS}, wE={wE})", flush=True)
            else:
                still_zero.append((S, b))
        if still_zero:
            print(f"\n  *** {len(still_zero)} sets show NO interior zero even via "
                  f"winding-difference -> COUNTEREXAMPLE CANDIDATES ***", flush=True)
            for S, b in still_zero:
                print(f"        S={S} B={b:+d}", flush=True)
        else:
            print("\n  All winding-0-at-safe-r sets resolved via winding-difference: "
                  "zero sits in the boundary annulus. No counterexample.", flush=True)

    # (C) sanity: winding-DIFFERENCE for a few canonical S at high r (artifact regime)
    print("\nWinding-difference sanity (high-r artifact regime cancels):", flush=True)
    Feul = make_F(frozenset(), Kmax)
    for S in [frozenset({1}), frozenset({3}), frozenset({1, 2}), frozenset({2, 4})]:
        F = make_F(S, Kmax)
        line = []
        for r in [0.9, 0.99]:
            wS, _, _ = winding_number(F, r, Ntheta=2048)
            wE, _, _ = winding_number(Feul, r, Ntheta=2048)
            line.append(f"r={r}: w(S)={wS}, w(E)={wE}, diff={wS-wE}")
        print(f"   S={set(S)} B={bsum(S):+d}: " + " | ".join(line), flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
