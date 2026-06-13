#!/usr/bin/env python3
"""
Proving the order-2 coefficient formula for shared-vertex tile pairs.
opus-2026-04-03-S28

CLAIM: For tiles (x₁,y₁) and (x₂,y₂) sharing vertex v in the transitive tournament:
  Case 1 (same-end): v = y₁ = y₂ (both lower) or v = x₁ = x₂ (both upper)
    → c₂ = -2

  Case 2 (cross-end): v = x₁ = y₂ or v = y₁ = x₂
    → c₂ = +2^(s₁ + s₂ - 2)  where s₁ = x₁-y₁, s₂ = x₂-y₂

PROOF STRATEGY:
We need to count H(flip both) - H(flip 1) - H(flip 2) + H(none).

H(none) = 1 (transitive)
H(flip i) = 1 + 2^(sᵢ-1)

For H(flip both), we need to count paths in the tournament where both arcs
y₁→x₁ and y₂→x₂ exist (instead of x₁→y₁ and x₂→y₂).

Case 1 (same lower vertex, v = y₁ = y₂, x₁ < x₂ WLOG):
  Both arcs go OUT of v: v→x₁ and v→x₂.
  A path can use AT MOST ONE of these arcs (since v has only one position in any HP).
  So any HP uses v→x₁ or v→x₂ or neither.

  The "neither" paths: same as transitive, count = 1.
  The "v→x₁ only" paths: must use v→x₁ but not v→x₂.
  The "v→x₂ only" paths: must use v→x₂ but not v→x₁.
  The "both" paths: IMPOSSIBLE (v visited only once).

  So H(flip both) = 1 + |A-only| + |B-only|

  But how does flipping the SECOND arc affect the |A-only| count?
  When only arc (x₁,y₁) is flipped, the A-type paths number 2^(s₁-1).
  When BOTH arcs are flipped, some A-type paths may be DESTROYED because
  the second arc (now y₂→x₂ instead of x₂→y₂) might block some paths.

  Actually NO — flipping x₂→y₂ to y₂→x₂ doesn't remove any arcs that
  the A-type paths use. The A-type paths go through v→x₁, and the arc
  between x₂ and y₂ (=v) is now v→x₂ instead of x₂→v.

  Wait: in the transitive tournament, x₂→v (since x₂ > v).
  After flipping: v→x₂.

  The A-type paths descend to v, jump to x₁, then continue descending.
  In the transitive+flip1 tournament, x₂→v is still available.
  But in the transitive+flip_both tournament, x₂→v is REMOVED (now v→x₂).
  So any A-type path that had x₂ visiting v via x₂→v is broken!

  But A-type paths go THROUGH v (they must reach v and then jump to x₁).
  The predecessor of v in an A-type path must be some vertex z > v.
  If z = x₂, then the path used arc x₂→v, which is now broken.
  If z ≠ x₂, the path is unaffected.

  How many A-type paths have x₂ as v's predecessor?
  In the "between-vertex partition" for arc v→x₁:
    Between-vertices = {v+1, ..., x₁-1}
    x₂ is one of the "above v" vertices (since x₂ > v = y₁)

  This is getting complex. Let me just verify the formulas computationally
  and look for patterns.
"""

def tournament_matrix(n, flips=None):
    T = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i > j:
                T[i][j] = 1
    if flips:
        for (x, y) in flips:
            T[x][y] = 0
            T[y][x] = 1
    return T

def count_hp(T, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0:
                continue
            for u in range(n):
                if not (mask & (1 << u)) and T[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[full])

from itertools import permutations

def list_hps(T, n):
    paths = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            paths.append(perm)
    return paths

# Systematic test: same-end pairs with various skip values
print("="*70)
print("SAME-END PAIRS (both lower vertex = v)")
print("="*70)

for n in range(4, 10):
    for v in range(n):
        # All pairs of tiles (x₁, v) and (x₂, v) where x₁-v >= 2, x₂-v >= 2
        upper_vertices = [x for x in range(v+2, n)]
        for i, x1 in enumerate(upper_vertices):
            for x2 in upper_vertices[i+1:]:
                s1, s2 = x1-v, x2-v

                H0 = 1
                H1 = count_hp(tournament_matrix(n, [(x1, v)]), n)
                H2 = count_hp(tournament_matrix(n, [(x2, v)]), n)
                H12 = count_hp(tournament_matrix(n, [(x1, v), (x2, v)]), n)
                c2 = H12 - H1 - H2 + H0

                if c2 != -2:
                    print(f"  n={n} v={v+1}: ({x1+1},{v+1})×({x2+1},{v+1}) "
                          f"skip=({s1},{s2}): c₂={c2} ≠ -2!")
                elif n <= 7:
                    print(f"  n={n} v={v+1}: ({x1+1},{v+1})×({x2+1},{v+1}) "
                          f"skip=({s1},{s2}): c₂={c2} ✓")

print(f"\n{'='*70}")
print("SAME-END PAIRS (both upper vertex = v)")
print("="*70)

for n in range(4, 10):
    for v in range(n):
        # All pairs of tiles (v, y₁) and (v, y₂) where v-y₁ >= 2, v-y₂ >= 2
        lower_vertices = [y for y in range(0, v-1)]
        for i, y1 in enumerate(lower_vertices):
            for y2 in lower_vertices[i+1:]:
                s1, s2 = v-y1, v-y2

                H0 = 1
                H1 = count_hp(tournament_matrix(n, [(v, y1)]), n)
                H2 = count_hp(tournament_matrix(n, [(v, y2)]), n)
                H12 = count_hp(tournament_matrix(n, [(v, y1), (v, y2)]), n)
                c2 = H12 - H1 - H2 + H0

                if c2 != -2:
                    print(f"  n={n} v={v+1}: ({v+1},{y1+1})×({v+1},{y2+1}) "
                          f"skip=({s1},{s2}): c₂={c2} ≠ -2!")
                elif n <= 7:
                    print(f"  n={n} v={v+1}: ({v+1},{y1+1})×({v+1},{y2+1}) "
                          f"skip=({s1},{s2}): c₂={c2} ✓")

print(f"\n{'='*70}")
print("CROSS-END PAIRS (v = x₁ = y₂, i.e., 'through' vertex)")
print("="*70)

for n in range(4, 10):
    for v in range(n):
        # Tile 1: (v, y₁) with v-y₁ >= 2 => y₁ <= v-2
        # Tile 2: (x₂, v) with x₂-v >= 2 => x₂ >= v+2
        for y1 in range(0, v-1):
            for x2 in range(v+2, n):
                s1 = v - y1  # skip of tile 1
                s2 = x2 - v  # skip of tile 2

                H0 = 1
                H1 = count_hp(tournament_matrix(n, [(v, y1)]), n)
                H2 = count_hp(tournament_matrix(n, [(x2, v)]), n)
                H12 = count_hp(tournament_matrix(n, [(v, y1), (x2, v)]), n)
                c2 = H12 - H1 - H2 + H0

                predicted = 2**(s1 + s2 - 2)
                match = "✓" if c2 == predicted else "✗"

                if n <= 7 or c2 != predicted:
                    print(f"  n={n} v={v+1}: ({v+1},{y1+1})×({x2+1},{v+1}) "
                          f"skip=({s1},{s2}): c₂={c2}, predicted={predicted} {match}")
