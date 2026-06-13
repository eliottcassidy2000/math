#!/usr/bin/env python3
"""
overnight_round1_s290.py — Round 1: D(n) formula for SL_mine
opus-2026-03-24-S290

kind-pasteur S20ef/S20eg: SL_mine computable via Burnside on cycle types
with EXACTLY ONE even cycle. The formula:

D(n) = (1/n!) Σ_{ct with exactly 1 even cycle 2k} count(ct) × k × 2^{a(ct)}

where a(ct) = arc orbits for the REMAINING odd cycles.

If D(n) ≈ SL_mine, then E(G_n) ≈ (T - D(n))/2 is exact up to
the multi-edge surplus.

Let me compute D(n) at n=3..12 and compare with known values.
"""

import sys
from math import comb, factorial, gcd
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

def gen_all_parts(n):
    result = []
    def _g(rem, mx, cur):
        if rem == 0: result.append(tuple(cur)); return
        for p in range(min(rem, mx), 0, -1): _g(rem - p, p, cur + [p])
    _g(n, n, [])
    return result

def ccs(n, ct):
    c = Counter(ct); r = factorial(n)
    for l, k in c.items(): r //= pow(l, k) * factorial(k)
    return r

def arc_orbits_odd(ct):
    """Arc orbits for odd-cycle-only part of ct."""
    s = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)): s += gcd(ct[i], ct[j])
    return s

def compute_D(n):
    """D(n) = (1/n!) Σ_{ct with exactly 1 even cycle} ccs × k × 2^{a}"""
    nf = factorial(n)
    D_sum = 0

    for ct in gen_all_parts(n):
        ct_list = list(ct)
        # Find cycle types with EXACTLY one even cycle
        even_cycles = [c for c in ct_list if c % 2 == 0]
        odd_cycles = [c for c in ct_list if c % 2 == 1]

        if len(even_cycles) != 1:
            continue

        k = even_cycles[0] // 2  # half the even cycle length
        cc = ccs(n, ct_list)

        # Arc orbits for the odd cycles
        if odd_cycles:
            ao = arc_orbits_odd(odd_cycles)
        else:
            ao = 0

        # Also need arc orbits BETWEEN the even cycle and odd cycles
        # For an even cycle of length 2k and odd cycles:
        # Arcs between them have gcd(2k, c_i) orbits each
        for c in odd_cycles:
            ao += gcd(even_cycles[0], c)

        # Self-pairing within the even cycle: k orbits (floor(2k/2) = k)
        # But some are self-conjugate...
        # Actually, the even cycle has 2k vertices. Arc orbits within it:
        # floor(2k/2) = k for intra-cycle
        ao += k  # intra-cycle contribution

        # The formula: contribution = ccs × k × 2^{ao}
        # Wait — kind-pasteur's formula is:
        # D(n) = (1/n!) Σ count(ct) × k × 2^{a(ct)}
        # where a(ct) includes ALL arc orbits of the combined permutation
        # MINUS the self-paired orbit from the even cycle that causes
        # the neutral arc.

        # Actually, I think the formula is simpler.
        # For a permutation with exactly one even cycle of length 2k:
        # It fixes tournaments only if the even cycle is "consistent."
        # The neutral arc is the one paired by the half-rotation of the 2k-cycle.
        # The other arcs contribute 2^{total_arc_orbits - 1} fixed tournaments.
        # And the neutral arc can be flipped → 1 self-loop.

        # So: SL contribution = ccs × 2^{total_arc_orbits - 1} × 1
        # where total_arc_orbits counts ALL orbits including the even cycle.

        # Let me just compute the full arc orbits for the entire cycle type
        full_ao = 0
        for i in range(len(ct_list)):
            full_ao += ct_list[i] // 2
            for j in range(i+1, len(ct_list)):
                full_ao += gcd(ct_list[i], ct_list[j])

        # But for even cycle types, Fix_tournament = 0 in general.
        # The key: we're not counting Fix_tournament. We're counting
        # how many (tournament, arc) pairs have σ(T) = T AND σ(arc flip) gives same class.

        # Actually, kind-pasteur's D(n) formula counts something specific.
        # Let me just use the values they computed:
        # D(3..12) = 2, 6, 16, 60, 328, 3160, 54928, 1722992, 97323552, 9941203552

        D_sum += cc * k * pow(2, full_ao - 1)  # tentative formula

    return D_sum // nf

# Known values from kind-pasteur
known_D = {3:2, 4:6, 5:16, 6:60, 7:328, 8:3160, 9:54928, 10:1722992}
known_T = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
known_E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}
known_twin = {3:2, 4:4, 5:12, 6:48, 7:296, 8:3040, 9:54256}

print("=" * 72)
print("  ROUND 1: D(n) AND E(G_n)")
print("=" * 72)

print(f"\n  {'n':>3} {'T':>12} {'twin_SL':>10} {'D(n)_known':>12} {'D(n)_comp':>12} "
      f"{'E_known':>10} {'(T-D)/2':>10} {'error':>8}")

for n in range(3, 11):
    # Compute D(n) with my formula
    D_comp = compute_D(n)

    T = known_T.get(n, 0)
    E = known_E.get(n)
    D_known = known_D.get(n)
    twin = known_twin.get(n, 0)

    E_from_D = (T - D_comp) // 2 if T > 0 and D_comp > 0 else 0
    err = abs(E_from_D - E) if E else '?'

    print(f"  {n:>3} {T:>12} {twin:>10} {str(D_known or '?'):>12} {D_comp:>12} "
          f"{str(E or '?'):>10} {E_from_D:>10} {str(err):>8}")

# The multi-edge surplus
print(f"\n  MULTI-EDGE SURPLUS (from kind-pasteur):")
print(f"  n=3: 0, n=4: 0, n=5: 12, n=6: 66, n=7: 416")
print(f"  T - 2E - D: ", end="")
for n in range(3, 10):
    T = known_T.get(n, 0); E = known_E.get(n, 0)
    D = known_D.get(n, 0)
    surplus = T - 2*E - D if T and E and D else '?'
    print(f"n={n}:{surplus}, ", end="")
print()

# The residual R/2 = (T - twin_SL - 2E) / 2
print(f"\n  RESIDUAL R/2:")
for n in range(3, 10):
    T = known_T.get(n, 0); E = known_E.get(n, 0)
    twin = known_twin.get(n, 0)
    R = T - twin - 2*E
    print(f"  n={n}: R={R}, R/2={R//2}")

print(f"\n  R/2 sequence: 0, 1, 8, 38, 222, 1463, 15721")
print(f"  D - twin_SL: ", end="")
for n in range(3, 10):
    D = known_D.get(n, 0); twin = known_twin.get(n, 0)
    print(f"{D - twin}, ", end="")
print()

# D - twin_SL should = complex_SL (the non-twin neutral arcs)
print(f"\n  D - twin_SL = complex_SL:")
for n in range(3, 10):
    D = known_D.get(n, 0); twin = known_twin.get(n, 0)
    cs = D - twin
    print(f"  n={n}: complex_SL = {cs}")

print("\nDONE.")
