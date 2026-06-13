#!/usr/bin/env python3
"""
sl_mine_efficient_s20eg.py — Efficient SL_mine via orbit direction analysis
kind-pasteur-2026-03-23-S20eg

THE TRICK: For sigma with cycle type ct, each arc orbit under sigma has
a DIRECTION SIGNATURE: at each position in the orbit, sigma either
PRESERVES or REVERSES the arc direction. A "mixed" orbit (some preserve,
some reverse) contributes defect to sigma(T) vs T.

For defect = 1: exactly one position in one orbit is "always-reversing",
and the arc flip compensates this. The count R(ct) is computable from
the cycle type alone by counting mixed orbits.

For a PRACTICAL formula:
  R(ct) = (number of "compensable" arc positions) * 2^{a(ct) - (correction)}
where "compensable" means the flip at that position fixes a mixed orbit.

We compute R(ct) for each partition ct with exactly one even part,
summing count(ct) * R(ct) / n! to get SL_mine.
"""

import sys
from math import factorial, comb, gcd
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EFFICIENT SL_mine COMPUTATION")
print("  kind-pasteur-2026-03-23-S20eg")
print("=" * 80)

def partitions(n, max_val=None):
    if max_val is None: max_val = n
    if n == 0: yield (); return
    for i in range(min(n, max_val), 0, -1):
        for rest in partitions(n - i, i):
            yield (i,) + rest

def count_perms(ct, n):
    mult = Counter(ct)
    d = 1
    for c, m in mult.items(): d *= (c**m) * factorial(m)
    return factorial(n) // d

def compute_R(ct, n):
    """Compute R(ct) = #{(T, e) : sigma(T flip e) = T} for sigma with cycle type ct.

    Strategy: build a representative sigma, analyze its arc orbits,
    determine which are "mixed" (have both preserving and reversing positions),
    and count the free parameters.
    """
    # Build representative sigma
    cycles_list = []
    pos = 0
    for c in ct:
        cycle = list(range(pos, pos + c))
        cycles_list.append(cycle)
        pos += c

    sigma = [0] * n
    for cycle in cycles_list:
        for i in range(len(cycle) - 1):
            sigma[cycle[i]] = cycle[i+1]
        sigma[cycle[-1]] = cycle[0]

    sigma_inv = [0] * n
    for i in range(n): sigma_inv[sigma[i]] = i

    # Enumerate arc positions and build orbits
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    arc_idx = {a: i for i, a in enumerate(ALL_ARCS)}

    visited = [False] * m
    orbits = []  # each orbit: list of (arc_index, preserves_direction)

    for start in range(m):
        if visited[start]: continue
        orbit = []
        k = start
        while not visited[k]:
            visited[k] = True
            u, v = ALL_ARCS[k]
            su, sv = sigma[u], sigma[v]
            # sigma maps arc (u,v) to arc (su, sv)
            # In our encoding: (min, max)
            target_arc = (min(su,sv), max(su,sv))
            target_idx = arc_idx[target_arc]
            # Direction: preserved if su < sv (same order), reversed if su > sv
            preserves = (su < sv)
            orbit.append((k, preserves))
            k = target_idx
        orbits.append(orbit)

    # Analyze each orbit
    # Standard orbit: all positions preserve direction -> sigma-invariant iff all same value
    # Reversed orbit: all positions reverse -> forced inconsistency (defect = orbit_size)
    # Mixed orbit: some preserve, some reverse -> partial defect

    # For defect-1 with flip:
    # We need sigma(T flip e) = T. This means: sigma(T) differs from T at exactly
    # the arc sigma(e), and the flip at e compensates.
    # So we need: defect(T, sigma) = 1, and the single defect is at sigma(e).
    # The flip at e changes T at position e, and sigma maps e to sigma(e).
    # After the flip: sigma(T flip e) at sigma(e) = (T flip e)(e) = 1-T(e).
    # For this to equal T(sigma(e)): 1-T(e) = T(sigma(e)).
    # And sigma(T)(sigma(e)) = T(e) [from the orbit structure].
    # So the defect at sigma(e): sigma(T)(sigma(e)) = T(e) != T(sigma(e)) iff T(e) != T(sigma(e)).
    # The flip condition: 1-T(e) = T(sigma(e)), i.e., T(e) != T(sigma(e)).
    # So the defect at sigma(e) = 1, and the flip compensates it.

    # For ALL other positions: defect = 0 (sigma(T) = T).

    # This means: we need an orbit where ONLY the position sigma(e) has defect.
    # In the orbit containing e: sigma maps e to some next position.
    # The orbit traces: e -> sigma(e) -> sigma^2(e) -> ...
    # Defect at each position p in the orbit: sigma(T)(p) != T(p).
    # sigma(T)(p) = T(sigma^{-1}(p)). So defect at p iff T(sigma^{-1}(p)) != T(p).
    # But we need to account for direction: sigma may reverse arcs.

    # For the orbit of arc e: let the orbit be [a0, a1, ..., a_{d-1}]
    # where sigma maps a_k to a_{k+1} (indices mod d).
    # At a_k: sigma(T)(a_k) = T(sigma^{-1}(a_k)) adjusted for direction.
    # If sigma preserves direction at a_k: sigma(T)(a_k) = T(a_{k-1}).
    # If sigma reverses at a_k: sigma(T)(a_k) = 1-T(a_{k-1}).

    # For defect = 0 at a_k (k != the special one):
    # T(a_k) = sigma(T)(a_k) = T(a_{k-1}) or 1-T(a_{k-1}) depending on direction.

    # This creates a CHAIN of constraints around the orbit.
    # For the orbit to have defect at EXACTLY one position:
    # the chain must be consistent everywhere except one link.

    # The number of "reversals" in the orbit's direction signature:
    # r = number of positions where sigma reverses direction.
    # If r is EVEN: the orbit is self-consistent (defect 0) for SOME T values.
    # If r is ODD: the orbit is ALWAYS inconsistent (defect >= 1).

    # For defect = 1 at the special position sigma(e):
    # We need r to be ODD in the orbit of e, and the flip at e
    # to make the chain consistent.

    # Count per orbit: reversals
    total_R = 0

    for orbit in orbits:
        d = len(orbit)
        reversals = sum(1 for _, preserves in orbit if not preserves)

        if reversals % 2 == 0:
            # Even reversals: orbit can be fully consistent.
            # No defect contribution. Cannot be used for flip-compensation.
            continue

        # Odd reversals: orbit ALWAYS has defect >= 1.
        # The defect is exactly 1 at one position (the "mismatch point").
        # Flipping at the RIGHT arc compensates this.
        # Number of valid flip positions in this orbit: d (any position works,
        # as the chain can be compensated by flipping at any link).
        # Actually: only specific positions work. Let me count.

        # For each position a_k in the orbit: flipping at a_{k-1} (the previous)
        # might compensate the defect at a_k.
        # The valid flip arc: e such that sigma(e) = the defect position.
        # There's exactly 1 such e per defect position.
        # And the defect position is the ONE position where the chain breaks.
        # With odd reversals: exactly one chain break, at a specific position.
        # So: 1 valid flip arc per odd-reversal orbit.

        # Wait: with d positions and odd reversals, the chain constraints
        # propagate around the cycle. The break point depends on the
        # specific arrangement of preserving/reversing links.
        # But for COUNTING: each orbit with odd reversals contributes
        # d possible (T, e) pairs where e is any arc in the orbit.
        # No: the flip must be at the specific arc whose sigma-image
        # is the defect point.

        # Actually from the brute force: for 4-cycle (0,1,2,3),
        # orbit {{0,2},{1,3}} has 2 arcs and 1 reversal (odd).
        # Both arcs {0,2} and {1,3} gave 8 solutions each.
        # So BOTH positions in the orbit are valid flip arcs.
        # That's d = 2 valid arcs.

        # So: each odd-reversal orbit of size d contributes d valid flip arcs.
        # For each valid flip arc: the number of T satisfying the constraints
        # is 2^{(free bits from other orbits) + (free bits from this orbit)}.
        # The orbit itself: d positions with d-1 independent constraints
        # (chain relations), plus the flip gives 1 free bit.
        # So: orbit contributes 1 free bit.
        # Other orbits: each contributes their standard free bits.
        # Even-reversal orbits contribute 1 free bit (all-same or all-different).

        # Total free bits = 1 (this orbit) + sum over other orbits of their free bits.
        total_R += d  # d valid flip positions per odd-reversal orbit

    # Total free bits for T: count standard free bits across all orbits
    free_bits_standard = 0
    for orbit in orbits:
        d = len(orbit)
        reversals = sum(1 for _, preserves in orbit if not preserves)
        if reversals % 2 == 0:
            free_bits_standard += 1  # orbit is consistent: 1 free bit

    # For each odd-reversal orbit that we use for the flip:
    # that orbit contributes 1 free bit (the starting value, rest determined)
    # Total free bits = free_bits_standard + 1 (from the odd orbit)
    total_free = free_bits_standard + 1

    # R = (total valid flip positions across all odd orbits) * 2^{total_free}
    R = total_R * (2 ** total_free)

    return R


# Compute SL_mine for each n
print()
for n in range(3, 16):
    t0 = time.time()
    total = 0

    for ct in partitions(n):
        # Check: exactly one even cycle
        even_count = sum(1 for c in ct if c % 2 == 0)
        if even_count != 1: continue

        nperms = count_perms(ct, n)
        R = compute_R(ct, n)
        total += nperms * R

    sl_mine = total // factorial(n)

    # Known values
    T_known = {3:4,4:16,5:88,6:704,7:8912,8:188288,9:6847200}
    E_known = {3:1,4:5,5:30,6:290,7:4086,8:91161,9:3380751}
    sl_expected = {n: T_known[n]-2*E_known[n] for n in T_known if n in E_known}

    expected = sl_expected.get(n)
    check = f"expected={expected} {'EXACT!' if sl_mine==expected else f'off by {sl_mine-expected}'}" if expected else ""
    print(f"  n={n:2d}: SL_mine={sl_mine:>14d}  {check}  ({time.time()-t0:.3f}s)")

print()
print("DONE.")
print("=" * 80)
