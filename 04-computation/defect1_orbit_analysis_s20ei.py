#!/usr/bin/env python3
"""
defect1_orbit_analysis_s20ei.py — Detailed orbit analysis for defect-1 Burnside
kind-pasteur-2026-03-23-S20ei

For each cycle type with at least one even cycle, analyze the ARC ORBIT structure:
- How many orbits?
- What are their sizes?
- How many reversals per orbit?
- Which orbits are "inconsistent" (odd reversals)?

The defect-1 condition: sigma(T flip e) = T means exactly:
  sigma(T) differs from T at exactly the arc sigma(e),
  and this difference is compensated by the flip at e.

Equivalently: there exists exactly ONE arc orbit where the T-assignment
is "inconsistent" with sigma, and that orbit is the one containing e,
and flipping e within that orbit makes it consistent.

For an orbit of size d with r reversals:
- If r is even: orbit CAN be consistent (2 choices) or fully inconsistent (2^d - 2 choices)
  When consistent: 0 defect. When inconsistent: ALL d positions have defect.
- If r is odd: orbit is ALWAYS "half-inconsistent" — exactly d/2 positions agree, d/2 disagree?
  No: with odd reversals in a cycle of length d, the chain constraint has a parity mismatch.
  Going around the cycle: T(a_d) = (-1)^r * T(a_0) = -T(a_0). Contradiction if d is odd.
  For ANY starting value, exactly ceil(d/2) or floor(d/2) positions disagree.

Actually let me think more carefully about the orbit defect structure.
"""

import sys
from math import factorial, gcd
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  ORBIT ANALYSIS FOR DEFECT-1")
print("  kind-pasteur-2026-03-23-S20ei")
print("=" * 80)

def build_sigma(ct, n):
    sigma = [0] * n
    pos = 0
    for c in ct:
        for i in range(c - 1):
            sigma[pos + i] = pos + i + 1
        sigma[pos + c - 1] = pos
        pos += c
    return sigma

def analyze_orbits(ct, n):
    """Analyze arc orbits under sigma with cycle type ct."""
    sigma = build_sigma(ct, n)
    sigma_inv = [0] * n
    for i in range(n): sigma_inv[sigma[i]] = i

    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    arc_idx = {a: i for i, a in enumerate(ALL_ARCS)}

    visited = [False] * m
    orbits = []

    for start in range(m):
        if visited[start]: continue
        orbit = []
        k = start
        while not visited[k]:
            visited[k] = True
            u, v = ALL_ARCS[k]
            su, sv = sigma[u], sigma[v]
            target_arc = (min(su,sv), max(su,sv))
            target_idx = arc_idx[target_arc]
            # Direction: preserved if sigma maintains order, reversed if not
            preserves = (su < sv)
            orbit.append((k, preserves))
            k = target_idx
        orbits.append(orbit)

    results = []
    for orbit in orbits:
        d = len(orbit)
        reversals = sum(1 for _, p in orbit if not p)
        parity = "even" if reversals % 2 == 0 else "odd"
        arcs = [ALL_ARCS[k] for k, _ in orbit]
        results.append({
            'size': d,
            'reversals': reversals,
            'parity': parity,
            'arcs': arcs,
            'directions': [p for _, p in orbit]
        })

    return results

def partitions(n, max_val=None):
    if max_val is None: max_val = n
    if n == 0: yield (); return
    for i in range(min(n, max_val), 0, -1):
        for rest in partitions(n - i, i):
            yield (i,) + rest

for n in range(3, 8):
    print(f"\n{'#'*60}")
    print(f"  n = {n}")
    print(f"{'#'*60}")

    for ct in partitions(n):
        has_even = any(c % 2 == 0 for c in ct)
        if not has_even: continue

        orbits = analyze_orbits(ct, n)
        even_orbits = [o for o in orbits if o['parity'] == 'even']
        odd_orbits = [o for o in orbits if o['parity'] == 'odd']

        print(f"\n  ct = {ct}")
        print(f"    Total orbits: {len(orbits)}")
        print(f"    Even-reversal orbits: {len(even_orbits)} (sizes: {[o['size'] for o in even_orbits]})")
        print(f"    Odd-reversal orbits: {len(odd_orbits)} (sizes: {[o['size'] for o in odd_orbits]})")

        for i, o in enumerate(orbits):
            tag = "EVEN" if o['parity'] == 'even' else "ODD*"
            dirs = ''.join('P' if p else 'R' for p in o['directions'])
            print(f"      Orbit {i}: size={o['size']}, rev={o['reversals']}, [{tag}] dirs={dirs}")
            if o['size'] <= 6:
                print(f"        arcs: {o['arcs']}")

        # The defect-1 analysis:
        # For T to have defect exactly 1 under sigma:
        # - All even-reversal orbits: T must be consistent (2 choices each)
        # - All odd-reversal orbits EXCEPT ONE: T must be... wait.
        #   Odd-reversal orbits are ALWAYS inconsistent (defect = d).
        #   So if there are k odd-reversal orbits, defect >= k*d_min.
        #   For defect = 1: need exactly 1 odd orbit of size 1.
        # But size-1 orbits are FIXED ARCS. They have 0 or 1 reversals.
        # A fixed arc with 1 reversal = odd = sigma reverses it = a TRANSPOSITION arc.
        # So: defect 1 from odd-reversal orbits requires exactly one size-1 odd orbit
        # and no other odd orbits.

        # BUT WAIT: for the flip compensation, the arc e that we flip can be in ANY orbit.
        # The condition sigma(T flip e) = T means sigma acts on (T flip e) to give T.
        # Let's denote T' = T flip e. Then sigma(T') = T.
        # Defect of (sigma, T) at arc a = 1 iff sigma(T)(a) != T(a).
        # sigma(T') = T means sigma(T') agrees with T everywhere.
        # T' = T flip e means T' differs from T only at arc e.
        # sigma(T') differs from sigma(T) only at arc sigma(e).
        # So: sigma(T') = sigma(T) except at sigma(e), where it's flipped.
        # For sigma(T') = T: sigma(T) = T everywhere except at sigma(e),
        # where sigma(T)(sigma(e)) != T(sigma(e)) but sigma(T')(sigma(e)) = T(sigma(e)).
        # So: defect(T, sigma) = 1, at position sigma(e).
        # And the flip at e compensates: sigma(T flip e) = T.

        # So the DEFECT of (T, sigma) must be EXACTLY 1.
        # Defect(T, sigma) = number of arcs a where sigma(T)(a) != T(a).
        # For each orbit: either ALL arcs agree (defect 0 from orbit) or some disagree.
        # The number of disagreements within an orbit depends on T values AND orbit structure.

        # For even-reversal orbit of size d: if T is "consistent" (follows the chain),
        # defect = 0. If T is "inconsistent", defect = d (all positions disagree).
        # So: defect from even orbit is 0 or d.

        # For odd-reversal orbit of size d: T can NEVER be fully consistent.
        # The chain forces T(a_0) = (-1)^r * T(a_0) where r is odd, so T(a_0) = -T(a_0).
        # In binary {0,1}: there's no consistent assignment.
        # But defect is not necessarily d! The chain can have some agreeing and some disagreeing.

        # Let me reconsider: in the orbit [a_0, a_1, ..., a_{d-1}],
        # sigma maps a_k to a_{k+1 mod d}.
        # sigma(T)(a_{k+1}) = T(a_k) if preserving, = 1-T(a_k) if reversing.
        # Defect at a_{k+1}: sigma(T)(a_{k+1}) != T(a_{k+1}).
        # So defect at a_{k+1} iff T(a_k) != T(a_{k+1}) (if preserving)
        #                      or T(a_k) = T(a_{k+1}) (if reversing).

        # Let v_k = T(a_k) in {0,1}. The defect at position k+1 depends on (v_k, v_{k+1}).
        # Going around the orbit, the TOTAL reversals r is odd,
        # so the product of (-1)^{reversal at each position} = (-1)^r = -1.
        # This means: the number of "direction changes" has odd parity mismatch.
        # So: an odd number of positions have defect. But HOW MANY?

        # For a specific T: assign v_0 = 0 or 1. Then v_1, v_2, ... are determined
        # by the defect=0 condition AT EACH LINK, EXCEPT one link must break.
        # With r odd: going around forces v_0 = 1-v_0. So exactly 1 link must break.
        # NO: "forced v_0 = 1-v_0" means the chain is inconsistent with period d.
        # For any assignment of v_0...v_{d-1}: the number of defects (disagreeing links)
        # is always ODD. But it could be 1, 3, 5, ...

        num_odd = len(odd_orbits)
        num_even = len(even_orbits)
        if num_odd == 0:
            print(f"    -> No odd orbits. Defect from even orbits only (0 or d). Defect=1 impossible unless d=1 even orbit with forced reversal... but d=1 even has 0 reversals. So DEFECT=1 IMPOSSIBLE for this type.")
        elif num_odd == 1 and odd_orbits[0]['size'] == 1:
            print(f"    -> Exactly 1 odd orbit of size 1 (transposition arc). Defect=1 from this. Other orbits must have defect 0.")
            free_bits = num_even  # each even orbit contributes 1 free bit
            print(f"    -> Free bits from even orbits: {free_bits}")
            print(f"    -> R_theory = 2 * 2^{free_bits} = {2 * 2**free_bits}")
            # The factor of 2: the flipped arc e = the transposition arc, and T(e) can be 0 or 1
            # Wait: the flip arc e must satisfy sigma(e) = defect position.
            # For size-1 orbit: sigma(a) = a. So sigma(e) = e. The defect is at e itself.
            # And flipping e: T(e) -> 1-T(e). For sigma(T')(e) = T(e): sigma(T)(e) = 1-T(e) != T(e). Check.
            # sigma(T)(e) = 1-T(e) because sigma reverses the arc.
            # So defect at e: sigma(T)(e) = 1-T(e) != T(e). Yes, defect=1 at e.
            # Flip at e: T' has T'(e) = 1-T(e). sigma(T')(e) = 1-T'(e) = 1-(1-T(e)) = T(e). Match!
            # For all other orbits: sigma(T) must agree with T. Even orbits: T consistent.
            # So R = (# choices for e) * (# consistent T for other orbits) * (2 for T(e))
            # = 1 * 2^{num_even} * 2 = 2^{num_even + 1}
            print(f"    -> R_theory = 2^{num_even + 1} = {2**(num_even+1)}")
        else:
            print(f"    -> {num_odd} odd orbits (sizes: {[o['size'] for o in odd_orbits]})")
            print(f"    -> Complex case: need exactly 1 total defect across all orbits.")
            print(f"    -> Each odd orbit contributes >=1 defect. Multiple odd orbits => defect >= {num_odd}.")
            if num_odd > 1:
                print(f"    -> DEFECT >= {num_odd} > 1. R should be 0 for defect-1.")
            else:
                d = odd_orbits[0]['size']
                print(f"    -> Single odd orbit of size {d}. Defect from this orbit is odd (1,3,5...).")
                print(f"    -> For defect=1: this orbit contributes 1 defect, all even orbits contribute 0.")
                # Need to count T assignments where this orbit has exactly 1 defect

print("\nDONE.")
print("=" * 80)
