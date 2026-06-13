#!/usr/bin/env python3
"""
How many directed 5-cycles can a tournament on 5 vertices have?
opus-2026-03-14-S71f

A directed 5-cycle on {0,1,2,3,4} is a cyclic ordering (a,b,c,d,e)
where a→b→c→d→e→a in the tournament.
Two orderings are the same cycle iff they're cyclic rotations.
So there are 4!/1 = 24 possible directed cycles (fix start vertex,
then 4! orderings, divide by... wait, fixing start vertex already
removes rotations).

Actually: fix vertex 0 as start. Then try all 4! = 24 orderings of
{1,2,3,4}. Each gives a potential cycle 0→p₁→p₂→p₃→p₄→0.
The cycle 0→1→2→3→4→0 is DIFFERENT from 0→4→3→2→1→0 (reverse).
"""

from itertools import permutations
from collections import Counter

def count_directed_5cycles_on_5_verts(A):
    """Count directed 5-cycles on {0,1,2,3,4}."""
    count = 0
    for p in permutations([1,2,3,4]):
        order = [0] + list(p)
        ok = True
        for i in range(5):
            if A[order[i]][order[(i+1) % 5]] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

# Exhaustive over all 2^10 = 1024 tournaments on 5 vertices
dist = Counter()
for mask in range(1024):
    A = [[0]*5 for _ in range(5)]
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if mask & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    c = count_directed_5cycles_on_5_verts(A)
    dist[c] += 1

print("Directed 5-cycle count on n=5 tournaments:")
for c in sorted(dist):
    print(f"  {c} directed 5-cycles: {dist[c]} tournaments ({dist[c]/1024*100:.1f}%)")

print("\nNote: each cyclic rotation gives the same cycle, and we fixed")
print("start vertex=0, so these are DISTINCT directed cycles.")
print(f"\nTotal with ≥2: {sum(v for k,v in dist.items() if k >= 2)}")

# Also check: at n=6, a subtournament on 5 vertices can have how many?
# It's the same distribution since the subtournament IS a tournament on 5 vertices.
# Key question: how does multiplicity affect OCF?

# Let's verify: for those with multiple 5-cycles, the cycles all share all vertices,
# so they're all adjacent in Ω. They form a clique of size (count).
# I(K_c, 2) = 1 + 2c (independent sets: empty + c singletons)
# But we also have 3-cycles that share vertices with these 5-cycles.

# For a VERTEX SET approach (counting vertex sets, not directed cycles):
# I_vs(Ω, 2) would use α₁ = # vertex sets with ANY directed cycle
# For DIRECTED CYCLE approach (OCF uses this):
# I_dc(Ω, 2) would use α₁ = # directed cycles (proper OCF)

# The difference: a vertex set with c directed 5-cycles contributes
#   c to α₁(dc) but only 1 to α₁(vs)
# Since all c cycles are mutually adjacent (same vertices), they form K_c in Ω
# I(K_c, 2) = 1 + 2c = 1 + 2·α₁_from_this_set

print(f"\n{'='*60}")
print("Impact on OCF:")
print(f"{'='*60}")
print("If a 5-element set has c directed 5-cycles:")
print("  Vertex-set approach: contributes 1 to α₁ (wrong)")
print("  Directed-cycle approach: contributes c to α₁ (correct)")
print("  These c cycles are all in a clique, so they contribute")
print("  I(K_c, 2) - 1 = 2c to H (not 2·1 = 2)")
print()
for c in sorted(dist):
    if c > 0:
        print(f"  c={c}: correct contribution = 2·{c} = {2*c}")
        print(f"         vs-approach contribution = 2·1 = 2  {'← ERROR!' if c > 1 else ''}")
