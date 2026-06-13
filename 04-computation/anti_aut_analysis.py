#!/usr/bin/env python3
"""
Anti-automorphism structure analysis for SC tournaments.
kind-pasteur-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, opposite_tournament, find_odd_cycles, conflict_graph
from itertools import permutations

def find_all_anti_auts(T):
    n = len(T)
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j and T[perm[i]][perm[j]] != (1 - T[i][j]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            auts.append(perm)
    return auts

def is_involution(sigma):
    return all(sigma[sigma[i]] == i for i in range(len(sigma)))

def order_of(sigma):
    n = len(sigma)
    x = list(range(n))
    for k in range(1, n + 2):
        x = [sigma[xi] for xi in x]
        if all(x[i] == i for i in range(n)):
            return k
    return None

# n=6, the H=45 SC tournament
print("=" * 70)
print("ANTI-AUTOMORPHISM STRUCTURE")
print("=" * 70)

n = 6
T = tournament_from_bits(n, 408)
auts = find_all_anti_auts(T)
print(f"\nn=6, bits=408 (H=45): {len(auts)} anti-automorphisms")
for sigma in auts:
    inv = is_involution(sigma)
    ord_s = order_of(sigma)
    fixed = [i for i in range(n) if sigma[i] == i]
    twos = [(i, sigma[i]) for i in range(n) if sigma[i] > i]
    print(f"  sigma={sigma}, order={ord_s}, involution={inv}, fixed={fixed}, 2-cycles={twos}")

# n=5, regular SC tournament
n = 5
T = tournament_from_bits(n, 76)
auts = find_all_anti_auts(T)
print(f"\nn=5, bits=76 (H=15, regular): {len(auts)} anti-automorphisms")
for sigma in auts:
    inv = is_involution(sigma)
    ord_s = order_of(sigma)
    fixed = [i for i in range(n) if sigma[i] == i]
    twos = [(i, sigma[i]) for i in range(n) if sigma[i] > i]
    print(f"  sigma={sigma}, order={ord_s}, involution={inv}, fixed={fixed}, 2-cycles={twos}")

# Key theoretical observation:
print("\n" + "=" * 70)
print("THEORETICAL ANALYSIS")
print("=" * 70)
print("""
For SC tournament T with anti-automorphism sigma:
- sigma reverses ALL arcs: T[sigma(i)][sigma(j)] = 1-T[i][j]
- If C = (v0->v1->...->v_{k-1}->v0) is an odd cycle,
  then sigma(C) = (sigma(v_{k-1})->...->sigma(v0)->sigma(v_{k-1})) is also an odd cycle
- sigma maps cycles to cycles (bijection on Omega(T))
- This induces sigma* on the conflict graph Omega(T)

KEY LEMMA: If sigma is an involution with no fixed points (possible at even n),
and C is a 3-cycle whose vertices come from 3 different sigma-orbits,
then C and sigma(C) are vertex-disjoint.

PROOF: If {v0,v1,v2} picks one from each sigma-orbit {a,sigma(a)},{b,sigma(b)},{c,sigma(c)},
then sigma({v0,v1,v2}) = {sigma(v0),sigma(v1),sigma(v2)} picks the OTHER element
from each orbit. So the vertex sets are disjoint, covering all 6 vertices.

At even n: sigma can be a fixed-point-free involution, creating n/2 orbits.
3-cycles picking one from each of 3 orbits have complementary pairs.
The number of such pairs is C(n/2, 3) * 2^3 / 2 = C(n/2,3) * 4.
At n=6: C(3,3) * 4 = 4 pairs. Matches computation!

At odd n: sigma must have at least one fixed point (since n is odd).
3-cycles through the fixed point are SELF-PAIRED (sigma(C) involves
same vertices modulo fixed point), reducing disjoint pair count.
""")

# Verify the count at n=6
from itertools import combinations
n = 6
sigma = (1, 0, 5, 4, 3, 2)
orbits = [(i, sigma[i]) for i in range(n) if i < sigma[i]]
print(f"n=6, sigma={sigma}")
print(f"Orbits: {orbits}")

# 3-cycles picking one from each orbit
cycles_from_orbits = []
for choice in range(8):  # 2^3 choices
    verts = []
    for k, (a, b) in enumerate(orbits):
        if choice & (1 << k):
            verts.append(b)
        else:
            verts.append(a)
    cycles_from_orbits.append(tuple(sorted(verts)))

print(f"Vertex sets from orbit choices: {cycles_from_orbits}")

# How many are complementary pairs?
pairs = []
for i in range(len(cycles_from_orbits)):
    comp = tuple(sorted([sigma[v] for v in cycles_from_orbits[i]]))
    for j in range(i + 1, len(cycles_from_orbits)):
        if cycles_from_orbits[j] == comp:
            pairs.append((cycles_from_orbits[i], cycles_from_orbits[j]))

print(f"Complementary disjoint pairs: {pairs}")
print(f"Count: {len(pairs)}")

# Now: of these 8 vertex sets, which are ACTUAL 3-cycles in the tournament?
T = tournament_from_bits(6, 408)
actual_cycles = find_odd_cycles(T)
c3_sets = [frozenset(c) for c in actual_cycles if len(c) == 3]
print(f"\nActual 3-cycles: {[sorted(c) for c in c3_sets]}")

orbit_sets = [frozenset(c) for c in cycles_from_orbits]
print(f"Orbit vertex sets that are actual 3-cycles:")
for vs in orbit_sets:
    is_cycle = vs in c3_sets
    print(f"  {sorted(vs)}: {'YES' if is_cycle else 'NO'}")

print("\nDone.")
