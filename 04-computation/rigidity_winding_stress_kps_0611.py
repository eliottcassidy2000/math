#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-485 / HYP-2417 -- STRESS TEST of the rigidity on a LARGE family + robustness.

The diagnostic established the truncation-safe annulus: at Kmax>=20 the truncated
Euler product has winding 0 for r <= 0.95, so a winding >= 1 there CERTIFIES a
genuine interior zero (modulus < r < 1) of the true F_S, free of boundary-truncation
artifacts. At r=0.93 all 255 hard-half sets (|S|<=4, k<=10, B(S)>=0) were certified.

Here we push HARDER to look for any counterexample to the rigidity:
  (1) bigger flip support: |S| <= 6 on k <= 12, hard half (B(S) >= 0);
  (2) robustness: certify each set at TWO independent safe radii (0.90, 0.93) and
      with TWO mesh densities, and require AGREEMENT of the integer winding;
  (3) an explicit minimum-|F_S| floor on the safe circle (must stay >> tail bound),
      so the integer count cannot be a near-zero aliasing fluke;
  (4) targeted "adversarial" S: the ones with the LARGEST positive boundary sum
      (deepest in the hard half, farthest from the IVT regime) and the all-even and
      all-odd supports.
A single winding-0 (no interior zero) survivor would REFUTE the rigidity.
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


def winding_number(F, r, Ntheta=2048, refine_depth=26):
    TWO_PI = 2 * math.pi
    THRESH = math.pi / 4
    pts = []
    minmod = float("inf")
    for i in range(Ntheta):
        th = TWO_PI * i / Ntheta
        f = F(r * cmath.exp(1j * th))
        m = abs(f)
        if m < minmod:
            minmod = m
        pts.append(f)

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
    return round(total / TWO_PI), total / TWO_PI, minmod


def bsum(S):
    return sum((1 if k % 2 == 1 else -1) for k in S)


def tail_bound(r, Kmax):
    k0 = Kmax + 1
    gmin = k0 * (3 * k0 - 1) // 2
    return 2.0 * (r ** gmin) / (1.0 - r) if r < 1 else float("inf")


def certify(S, Kmax, radii, meshes):
    """Robust certificate: winding>=1 with integer AGREEMENT across radii & meshes,
    and min|F| floor >> tail. Returns (ok, detail)."""
    F = make_F(S, Kmax)
    windings = []
    minmods = []
    for r in radii:
        for N in meshes:
            w, wf, mm = winding_number(F, r, Ntheta=N)
            windings.append(w)
            minmods.append((mm, tail_bound(r, Kmax)))
    allpos = all(w >= 1 for w in windings)
    floor_ok = all(mm > 1e4 * tb and mm > 1e-9 for mm, tb in minmods)
    return (allpos and floor_ok), (windings, min(mm for mm, _ in minmods))


def main():
    t0 = time.time()
    Kmax = 40
    safe_radii = [0.90, 0.93]
    meshes = [1024, 2048]

    print("=" * 78, flush=True)
    print("STRESS TEST: rigidity on |S|<=6, k<=12, hard half, robust winding cert",
          flush=True)
    print(f"  Kmax={Kmax}, safe radii={safe_radii}, meshes={meshes}", flush=True)
    print(f"  (Euler winding = 0 confirmed for r<=0.95 at Kmax>=20)", flush=True)
    print("=" * 78, flush=True)

    # control re-confirm at this Kmax
    Feul = make_F(frozenset(), Kmax)
    print("\nEuler control (must be winding 0 at safe radii):", flush=True)
    for r in safe_radii:
        w, _, mm = winding_number(Feul, r, Ntheta=2048)
        print(f"   r={r}: winding={w}, min|F|={mm:.3e}, tail<= {tail_bound(r, Kmax):.2e}",
              flush=True)

    KFLIP, MAXSIZE = 12, 6
    hard = []
    for sz in range(1, MAXSIZE + 1):
        for S in itertools.combinations(range(1, KFLIP + 1), sz):
            S = frozenset(S)
            if bsum(S) >= 0:
                hard.append(S)
    print(f"\nHARD half (B(S)>=0) on k<={KFLIP}, |S|<={MAXSIZE}: {len(hard)} sets",
          flush=True)

    certified = 0
    failures = []
    wdist = Counter()
    t_loop = time.time()
    for i, S in enumerate(hard):
        ok, (windings, minmod) = certify(S, Kmax, safe_radii, meshes)
        wdist[min(windings)] += 1
        if ok:
            certified += 1
        else:
            failures.append((set(S), bsum(S), windings, minmod))
        if (i + 1) % 500 == 0:
            print(f"   ... {i+1}/{len(hard)} done "
                  f"({time.time()-t_loop:.0f}s, {certified} certified so far)",
                  flush=True)

    print(f"\nROBUSTLY CERTIFIED (winding>=1 at all radii*meshes, |F| floor ok): "
          f"{certified}/{len(hard)}", flush=True)
    print(f"min-winding distribution: "
          f"{dict(sorted(wdist.items()))}", flush=True)

    if failures:
        print(f"\n*** {len(failures)} NON-CERTIFIED (rigidity counterexample candidates) ***",
              flush=True)
        for S, b, ws, mm in failures[:40]:
            print(f"   S={S} B={b:+d} windings={ws} min|F|={mm:.2e}", flush=True)
    else:
        print("\nNo counterexample on the entire hard half. Rigidity HOLDS here.",
              flush=True)

    # targeted adversarial probes: deepest hard sets (max boundary sum) and parity-pure
    print("\n" + "-" * 78, flush=True)
    print("Adversarial probes (deepest hard half + parity-pure supports):", flush=True)
    print("-" * 78, flush=True)
    probes = [
        frozenset({1, 3, 5, 7, 9, 11}),   # all-odd k -> B = +6, deepest hard
        frozenset({1, 3, 5}),             # B = +3
        frozenset({2, 4, 6, 8, 10, 12}),  # all-even k -> B = -6 (easy, but check)
        frozenset({1, 2, 3, 4, 5, 6}),    # mixed B = 0
        frozenset({11}),                  # single deep odd
        frozenset({1, 3, 5, 7, 9}),       # B = +5
        frozenset({1, 11}),               # B = +2
    ]
    for S in probes:
        ok, (windings, mm) = certify(S, Kmax, safe_radii, meshes)
        print(f"   S={set(S)} B={bsum(S):+d}: certified={ok}, windings={windings}, "
              f"min|F|={mm:.2e}", flush=True)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
