#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""What a DISPROOF candidate looks like across the analogous (esp. geometric) forms of the
LRC, and the cross-form patterns. (mac-mini-2026-06-29-S16)

LRC(n): n-1 speeds, gap 1/n. Equivalent forms (disproof candidate in each):
  - Diophantine: integer (n-1)-set S, M(S)=max_t min_i ||v_i t|| < 1/n.
  - Covering combs: danger combs (width 1/(... )) cover [0,1).
  - VIEW-OBSTRUCTION (geometric): the closed geodesic {t*(v_1..v_{n-1}) mod 1} on the torus
    T^{n-1} AVOIDS the central cube C_n = [1/n, 1-1/n]^{n-1}.
  - BILLIARD: unfolded trajectory trapped in the boundary slab (within 1/n of a wall).
  - LATTICE/approx (Littlewood-sibling): a rational direction whose geodesic stays within 1/n
    of the coordinate-hyperplane arrangement for ALL t.

This script verifies the cross-form PATTERNS:
 P1 equidistribution dichotomy: central cube has volume (1-2/n)^{n-1} -> e^-2 > 0, so irrational
    directions ALWAYS visit it (Weyl) => disproof candidates are RATIONAL/integer only.
 P3 universal extremal: AP {1..n-1}; lonely set = the phi(n) UNITS mod n, in phi(n)/2 antipodal
    pairs (saddle index), touching at multiplicative-inverse runner pairs.
 P4 covering = geodesic 'pinned to walls' at the Farey points 1/q (q<=n); non-covering killed
    by the t=1/q witness.
 P2/P6 the R-antipodal x->1-x (=t->-t=complement) and hyperoctahedral B_{n-1} symmetry.
"""
from __future__ import annotations
import functools, math
print = functools.partial(print, flush=True)


def frac(x): return x - math.floor(x)
def dist(x): f = frac(x); return min(f, 1 - f)
def phi(n): return sum(1 for a in range(1, n) if math.gcd(a, n) == 1)


def lonely_set_AP(n):
    """lonely points of the AP {1..n-1}: t=a/n with min_i ||i a / n|| >= 1/n."""
    S = list(range(1, n)); pts = []
    for a in range(1, n):
        t = a / n
        md = min(dist(v * t) for v in S)
        if md >= 1/n - 1e-12:
            g = math.gcd(a, n)
            touch = [i for i in range(1, n) if abs(dist(i*a/n) - 1/n) < 1e-12]
            pts.append((a, g, touch))
    return pts


def main():
    print("=" * 80)
    print("DISPROOF candidates across the analogous forms of the LRC, and the patterns")
    print("=" * 80)

    # P1: central cube volume -> e^-2
    print("\n[P1] EQUIDISTRIBUTION DICHOTOMY: central cube C_n=[1/n,1-1/n]^{n-1}, vol=(1-2/n)^{n-1}")
    print(f"     {'n':>3} {'vol(C_n)':>10}   (-> e^-2 = {math.exp(-2):.5f})")
    for n in [4, 6, 8, 14, 30, 100, 1000]:
        vol = (1 - 2/n) ** (n - 1)
        print(f"     {n:>3} {vol:>10.5f}")
    print("     vol>0 always => irrational directions EQUIDISTRIBUTE (Weyl) & VISIT C_n => not")
    print("     disproofs. DISPROOF CANDIDATES = rational/integer directions ONLY (closed geodesics).")
    print("     [Littlewood-sibling: a disproof = a direction badly-distributed toward the walls.]")

    # P3: the universal extremal -- AP lonely set = units mod n
    print("\n[P3] UNIVERSAL EXTREMAL = AP {1..n-1}; lonely set = UNITS mod n, phi(n)/2 antipodal pairs:")
    print(f"     {'n':>3} {'phi(n)':>6} {'#pairs':>6} {'lonely a (units)':>22} {'saddle idx phi/2':>16}")
    for n in range(4, 15):
        pts = lonely_set_AP(n)
        units = [a for (a, g, t) in pts]
        allunit = all(g == 1 for (a, g, t) in pts)
        npairs = len(units) // 2
        print(f"     {n:>3} {phi(n):>6} {npairs:>6} {str(units):>22} {phi(n)//2:>16}"
              + ("" if allunit and len(units) == phi(n) else "  <-- mismatch!"))
    print("     PATTERN: lonely set of the AP = EXACTLY the phi(n) units mod n (proof: t=a/n is")
    print("     lonely iff gcd(a,n)=1, since gcd=g>1 gives runner i=n/g hitting 0). Touch at i=a^{-1}.")
    print("     saddle index = phi(n)/2 = #antipodal unit pairs. n=14: 3 = (7-1)/2 (apex 7).")

    # the inverse-touch structure at n=14
    print("\n     n=14 inverse-touch detail (a, touching runners = +-a^{-1} mod 14):")
    for (a, g, touch) in lonely_set_AP(14):
        print(f"       t={a}/14: runners {touch} touch (a^{{-1}} mod 14 = {pow(a, -1, 14)})")

    # P4: covering = pinned at Farey points 1/q; non-covering witness
    print("\n[P4] COVERING = geodesic PINNED to walls at Farey points 1/q (q<=n):")
    print("     non-covering S (omits all multiples of some q<=n) is killed by the t=1/q WITNESS:")
    print("     at t=1/q, every v in S has v/q non-integer => ||v/q|| >= 1/q >= 1/n. Verify:")
    import random
    rng = random.Random(16)
    ok = 0; tot = 0
    for _ in range(2000):
        n = 14
        S = sorted(set(rng.randint(1, 80) for _ in range(13)))[:13]
        if len(S) != 13: continue
        # find a q<=14 with no multiple in S (non-covering)
        q = next((q for q in range(2, 15) if all(v % q for v in S)), None)
        if q is None: continue
        tot += 1
        Mwit = min(dist(v * (1/q)) for v in S)   # min dist at the witness t=1/q
        if Mwit >= 1/14 - 1e-12: ok += 1
    print(f"     non-covering sets tested: {tot}; witness t=1/q gives min-dist >= 1/14: {ok}/{tot} "
          f"({'CONFIRMED' if ok==tot else 'FAIL'})")
    print("     => a DISPROOF geodesic must pass through a wall at EVERY Farey point 1/2,1/3,...,1/n")
    print("     (covering) AND stay in the wall-tube between them (M<1/n). Pinned + trapped.")

    # P2/P6: symmetry
    print("\n[P2/P6] SYMMETRY of the central cube C_n (and the disproof-candidate set):")
    print("     hyperoctahedral B_{n-1} = (Z/2)^{n-1} >| S_{n-1}, order 2^{n-1}(n-1)!:")
    print("       - S_{n-1}: permute coordinates = relabel speeds.")
    print("       - (Z/2)^{n-1}: reflect x_i->1-x_i = flip speed sign v_i->-v_i.")
    print("     The GLOBAL reflection x->1-x (all coords) = t->-t = the COMPLEMENT/REVERSAL R.")
    print("     => the lonely set is R-symmetric (verified S15); the central-cube antipodal map IS")
    print("     the tournament complement (klein THM-584). Disproof candidates = B_{n-1}-orbits of")
    print("     rational directions, R-antipodally paired -- the SAME signed-permutation theme as")
    print("     the tournament metagraph signed cycle index (klein THM-587).")

    print("\n" + "=" * 80)
    print("CROSS-FORM PATTERNS: disproof candidates are (P1) rational/integer directions only")
    print("(closed geodesics; irrationals equidistribute into the vol->e^-2 central cube); (P3)")
    print("extremal = AP/staircase, lonely set = phi(n) units in phi(n)/2 antipodal pairs; (P4)")
    print("covering = pinned to walls at Farey points 1/q + trapped in the wall-tube; (P2/P6)")
    print("organized by hyperoctahedral B_{n-1} with the R-antipodal (=complement) involution.")
    print("=" * 80)


if __name__ == "__main__":
    main()
