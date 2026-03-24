#!/usr/bin/env python3
"""
sl_mine_exact_s20ef.py — Compute SL_mine EXACTLY via defect-1 Burnside
kind-pasteur-2026-03-23-S20ef

SL_mine(n) = (1/n!) * sum_{all sigma in S_n} sum_{all labeled T} 1[defect(T,sigma)=1]

For a fixed sigma with cycle type ct (all odd cycles, otherwise 0):
  sum_T 1[defect(T,sigma)=1] = #{T: exactly 1 arc position disagrees with sigma(T)}

The arc positions split into orbits under sigma.
sigma(T) = T at an orbit iff all arcs in the orbit have consistent values.
sigma(T) != T at an orbit iff the orbit has an INCONSISTENCY.

Defect = 1 means: exactly 1 arc position disagrees.
This means: exactly 1 arc ORBIT has an inconsistency, and that orbit has size 1
(a fixed arc), and the fixed arc has the "wrong" value.

Wait: if an orbit has size k > 1, and one position disagrees, then sigma
cycles the positions within the orbit. If position p disagrees, then
sigma maps p to sigma(p), and sigma(T)(sigma(p)) = T(p) != T(sigma(p))
(since T is NOT sigma-invariant at p). So sigma(p) also disagrees.
Continuing: ALL positions in the orbit disagree. So defect >= orbit_size.

THEREFORE: defect = 1 is possible ONLY if the disagreeing orbit has size 1.
An orbit of size 1 = a FIXED arc (both endpoints fixed by sigma, or
the endpoints are swapped by sigma).

For TOURNAMENT arc {u,v}: sigma fixes {u,v} as a SET iff
either sigma(u)=u,sigma(v)=v (both fixed) or sigma(u)=v,sigma(v)=u (swapped).

If both fixed: sigma(T)(u,v) = T(u,v). Agreement always. Cannot have defect here.

If swapped: sigma(T)(u,v) = T(v,u) = 1-T(u,v). Defect = 1 at this arc.
THIS ALWAYS DISAGREES for tournaments! So a "swapped" fixed arc ALWAYS has defect 1.

THEREFORE: defect(T, sigma) = (number of swapped arc pairs in sigma).

A "swapped pair" = a 2-element subset {u,v} where sigma(u)=v and sigma(v)=u.
This happens iff u and v are in the SAME 2-cycle of sigma.
Since sigma has all odd cycles (for Fix > 0), there are NO 2-cycles.

Wait: if sigma has all ODD cycles, there are no 2-cycles! So no swapped pairs!
So defect = 0 always for all-odd sigma? That can't be right...

The issue: for all-odd-cycle sigma, sigma(T) = T means T is sigma-invariant.
Each arc orbit under sigma has all-equal values.
The number of T with sigma(T) = T is 2^{a(sigma)} (one free bit per orbit).
But for SOME T: sigma(T) != T. The defect is the number of disagreeing positions.

For sigma with all-odd cycles: an arc {u,v} where u,v are in the SAME cycle.
sigma maps {u,v} to {sigma(u), sigma(v)}, which is a DIFFERENT pair (since the
cycle has length >= 3, sigma doesn't fix or swap u,v).
So this arc is in an orbit of size > 1. If ONE position in the orbit disagrees,
ALL positions disagree. Defect contribution = orbit_size.

For an arc between u (in cycle C1) and v (in cycle C2, C1 != C2):
sigma maps to {sigma(u), sigma(v)}. Orbit size = lcm(|C1|, |C2|) / gcd...
It's gcd(|C1|,|C2|). Actually orbit size = lcm(|C1|,|C2|).

For sigma with EVEN cycles: sigma has a 2-cycle (u,v).
sigma maps arc {u,v} to {v,u} = {u,v}. Fixed as a set!
But sigma(T)(u,v) = T(v,u) = 1-T(u,v). ALWAYS disagrees. Defect += 1.

So: for sigma with a 2-cycle (u,v): the arc {u,v} ALWAYS has defect 1,
regardless of T. And other arcs may have additional defect.

For defect EXACTLY 1: we need the arc {u,v} (from the 2-cycle) to disagree,
AND all other arcs to agree.

This means: sigma has exactly ONE 2-cycle, and all other cycles are odd.
AND: T is sigma-invariant on all arc orbits except the one containing {u,v}.
Since {u,v} is a fixed arc (swapped pair), it's an orbit of size 1.
So: T must be sigma-invariant on all OTHER orbits, and {u,v} auto-disagrees.

For T to be sigma-invariant on all other orbits: 2^{a(sigma) - 1} choices
(all orbits free except the fixed arc, which is forced).

Wait: the "fixed arc" {u,v} has orbit size 1 (sigma swaps u and v, mapping
{u,v} to {v,u} = {u,v}). For this orbit: sigma(T)(u,v) = 1-T(u,v) ALWAYS.
So this orbit always disagrees. There are 0 free bits for this orbit
(T(u,v) can be anything, but sigma always flips it).

For the other a(sigma)-1 orbits: each has a free bit (2 choices), and
we need sigma(T) = T on these orbits. The number of T satisfying this:
2^{a(sigma)-1} * 2 = 2^{a(sigma)} (the {u,v} bit is free, just always flipped).

Wait: for defect EXACTLY 1: we need sigma(T) to agree with T on ALL orbits
except {u,v}. The number of such T = 2^{a(sigma)-1} * 2 = 2^{a(sigma)}.
Hmm, that's the FULL Fix(sigma) count, which should be 0 for even-cycle sigma.

I'm confusing myself. Let me restart.

For sigma with cycle type including a 2-cycle:
Fix(sigma) = 0 (no tournament is fixed by sigma, because the 2-cycle
forces T(u,v) = 1-T(u,v) which is impossible).

For the DEFECT computation: defect(T, sigma) counts arc positions where
sigma(T) != T. For the arc {u,v} (swapped pair): ALWAYS disagrees (defect +1).
For other arcs in sigma orbits of odd size: defect = 0 iff orbit is consistent,
= orbit_size otherwise.

So: defect(T, sigma) >= 1 for ALL T (the swapped arc always contributes).
Defect = exactly 1 iff ALL other orbits are consistent.
Number of such T: 2^{a(sigma) - 1} * 2 (free choices for consistent orbits,
plus 2 choices for the {u,v} arc value).
= 2^{a(sigma)}.

But wait: I said Fix(sigma) = 0 and now defect-1 count = 2^{a(sigma)} > 0?

Yes! Fix(sigma) counts T with defect = 0. That's 0 (impossible for even cycles).
But defect = 1 is POSSIBLE: the swapped arc always disagrees, and all others agree.

So: for sigma with EXACTLY ONE 2-cycle and all other cycles odd:
  #{T: defect(T,sigma) = 1} = 2^{a(sigma)}
  where a(sigma) = total arc orbits including the {u,v} orbit.

For sigma with MORE THAN ONE 2-cycle: each 2-cycle adds a forced disagreement.
So defect >= #{2-cycles}. Defect = #{2-cycles} iff all other orbits agree.
For defect = 1: need exactly one 2-cycle.

FORMULA:
  total_d1 = sum over sigma with exactly one 2-cycle and all other cycles odd
             of 2^{a(sigma)}

This is BURNSIDE-COMPUTABLE! We sum over cycle types (2, odd, odd, ..., odd).
"""

import sys
from math import factorial, comb, gcd
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EXACT SL_mine VIA DEFECT-1 BURNSIDE")
print("  kind-pasteur-2026-03-23-S20ef")
print("=" * 80)

def partitions(n, max_val=None):
    if max_val is None: max_val = n
    if n == 0: yield (); return
    for i in range(min(n, max_val), 0, -1):
        for rest in partitions(n - i, i):
            yield (i,) + rest

def count_perms(ct, n):
    mult = Counter(ct)
    denom = 1
    for c, m in mult.items():
        denom *= (c ** m) * factorial(m)
    return factorial(n) // denom

def arc_orbits(ct, n):
    """Number of arc orbits under perm with cycle type ct."""
    cycles = list(ct)
    total = sum(c // 2 for c in cycles)
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            total += gcd(cycles[i], cycles[j])
    return total

print("\nFORMULA: SL_mine(n) = (1/n!) * sum over cycle types (2, odd^*)")
print("         of count(type) * 2^{a(type)}")
print()

# The cycle types with exactly one 2-cycle and all other cycles odd:
# (2, c1, c2, ...) where each ci is odd and 2 + sum(ci) = n.

for n in range(3, 16):
    t0 = time.time()
    total = 0

    for ct in partitions(n):
        # Check: exactly one 2-cycle, all others odd
        num_2 = sum(1 for c in ct if c == 2)
        all_others_odd = all(c % 2 == 1 for c in ct if c != 2)

        if num_2 == 1 and all_others_odd:
            nperms = count_perms(ct, n)
            a = arc_orbits(ct, n)
            contrib = nperms * (2 ** a)
            total += contrib

    sl_mine = total // factorial(n)

    # Known values for verification
    T_known = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
    E_known = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

    T = T_known.get(n)
    E = E_known.get(n)
    sl_expected = T - 2*E if T and E else None

    check = f"expected={sl_expected} {'MATCH!' if sl_mine == sl_expected else 'MISMATCH'}" if sl_expected else ""
    print(f"  n={n:2d}: SL_mine = {sl_mine:>12d}  {check}  ({time.time()-t0:.3f}s)")

print()
print("E(G_n) = (T_n - SL_mine) / 2 where BOTH are Burnside-computable!")
print()

# Now compute E exactly
print("EXACT E(G_n):")
for n in range(3, 16):
    # T_n via standard Burnside
    T = 0
    for ct in partitions(n):
        if any(c % 2 == 0 for c in ct): continue
        nperms = count_perms(ct, n)
        a = arc_orbits(ct, n)
        fixed_arcs = comb(sum(1 for c in ct if c == 1), 2)
        T += nperms * (2 ** a) * fixed_arcs
    T //= factorial(n)

    # SL_mine via defect-1
    SL = 0
    for ct in partitions(n):
        num_2 = sum(1 for c in ct if c == 2)
        all_others_odd = all(c % 2 == 1 for c in ct if c != 2)
        if num_2 == 1 and all_others_odd:
            SL += count_perms(ct, n) * (2 ** arc_orbits(ct, n))
    SL //= factorial(n)

    E = (T - SL) // 2
    E_known_val = E_known.get(n)
    check = f"known={E_known_val} {'EXACT!' if E == E_known_val else f'off by {E - E_known_val}'}" if E_known_val else ""

    print(f"  n={n:2d}: T={T:>16d}  SL={SL:>12d}  E=(T-SL)/2={E:>14d}  {check}")

print()
print("DONE.")
print("=" * 80)
