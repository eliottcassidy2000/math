#!/usr/bin/env python3
"""
Recursive structure: new HPs from single flip ↔ binary strings ↔ smaller staircases.
opus-2026-04-03-S28

When we flip tile (x,y) from the transitive tournament:
- The 2^(s-1) new HPs correspond to binary strings in {0,1}^{s-1}
  where bit k = 0 means "visit between-vertex y+k before the jump"
  and bit k = 1 means "visit between-vertex y+k after the jump"
- This binary string IS a tiling of a smaller interval [y, x]!

QUESTION: When we flip TWO tiles simultaneously, the new HPs that use
BOTH reversed arcs are counted by c₂. Is c₂ related to counting
binary strings that satisfy BOTH tile constraints simultaneously?

This would give the multilinear polynomial a RECURSIVE structure:
H at level n involves products of H-like counts at smaller levels.

Let me investigate this.
"""

from itertools import permutations, product
from collections import defaultdict

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

def list_hps(T, n):
    paths = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            paths.append(perm)
    return paths

# For n=5, study the HPs from flipping BOTH the apex (5,1) and a skip-2 tile

n = 5
print("="*70)
print(f"RECURSIVE TILING STRUCTURE: n={n}")
print("="*70)

# Apex flip A=(5,1) creates 8 new paths
# Each is a binary string of length 3 (between-vertices: 2, 3, 4)
# Skip-2 tile E=(4,2) creates 2 new paths
# Between-vertex: 3

# Double flip AE: how do the binary strings interact?
T_ae = tournament_matrix(n, flips=[(4, 0), (3, 1)])  # (5,1) and (4,2) 0-indexed
hps = list_hps(T_ae, n)

print(f"\nDouble flip A=(5,1) + E=(4,2): H = {len(hps)}")
print(f"Expected: 1 + 2^3 + 2^1 + c₂ = 1 + 8 + 2 + 4 = 15")

# Categorize each HP
T_trans = tournament_matrix(n)
for p in hps:
    # Does it use arc 1->5 (i.e., 0->4 in 0-indexed)?
    uses_A = any(p[k] == 0 and p[k+1] == 4 for k in range(n-1))
    # Does it use arc 2->4 (i.e., 1->3 in 0-indexed)?
    uses_E = any(p[k] == 1 and p[k+1] == 3 for k in range(n-1))

    category = ""
    if uses_A and uses_E: category = "BOTH"
    elif uses_A: category = "A only"
    elif uses_E: category = "E only"
    else: category = "neither (original)"

    path_str = ' -> '.join(str(v+1) for v in p)
    print(f"  {path_str}  [{category}]")

# The "BOTH" paths should number c₂ = 4
both_count = sum(1 for p in hps
                 if any(p[k]==0 and p[k+1]==4 for k in range(n-1))
                 and any(p[k]==1 and p[k+1]==3 for k in range(n-1)))
print(f"\nPaths using BOTH reversed arcs: {both_count}")

# The STRUCTURE of "BOTH" paths:
# They must contain ... -> 1 -> 5 -> ... and ... -> 2 -> 4 -> ...
# In vertices 1-5, the path visits all 5 vertices using these two jumps.
# The "between" vertices for arc 1->5 are {2,3,4}.
# The arc 2->4 is INSIDE this between-set!
# So the binary string for the 1->5 jump must have vertex 2 "before" and vertex 4 "after"
# (since the HP goes ...->2->... then ...->1->5->...->4->... )
# Wait, that's not right. Let me trace the paths.

print(f"\nDetailed BOTH paths:")
for p in hps:
    uses_A = any(p[k]==0 and p[k+1]==4 for k in range(n-1))
    uses_E = any(p[k]==1 and p[k+1]==3 for k in range(n-1))
    if uses_A and uses_E:
        path_str = ' -> '.join(str(v+1) for v in p)
        # Find positions
        a_pos = next(k for k in range(n-1) if p[k]==0 and p[k+1]==4)
        e_pos = next(k for k in range(n-1) if p[k]==1 and p[k+1]==3)
        print(f"  {path_str}  (1->5 at pos {a_pos}, 2->4 at pos {e_pos})")

# Let me now study the GENERAL pattern: how many "BOTH" paths for all pairs?
print(f"\n{'='*70}")
print(f"BOTH-ARC PATH COUNTS FOR ALL PAIRS")
print(f"{'='*70}")

tiles_1idx = [(5,1), (4,1), (3,1), (5,2), (4,2), (5,3)]
m = len(tiles_1idx)

for ti in range(m):
    for tj in range(ti+1, m):
        x1, y1 = tiles_1idx[ti]
        x2, y2 = tiles_1idx[tj]

        T = tournament_matrix(n, flips=[(x1-1, y1-1), (x2-1, y2-1)])
        hps = list_hps(T, n)

        both = sum(1 for p in hps
                   if any(p[k]==y1-1 and p[k+1]==x1-1 for k in range(n-1))
                   and any(p[k]==y2-1 and p[k+1]==x2-1 for k in range(n-1)))

        a_only = sum(1 for p in hps
                    if any(p[k]==y1-1 and p[k+1]==x1-1 for k in range(n-1))
                    and not any(p[k]==y2-1 and p[k+1]==x2-1 for k in range(n-1)))

        b_only = sum(1 for p in hps
                    if not any(p[k]==y1-1 and p[k+1]==x1-1 for k in range(n-1))
                    and any(p[k]==y2-1 and p[k+1]==x2-1 for k in range(n-1)))

        neither = sum(1 for p in hps
                     if not any(p[k]==y1-1 and p[k+1]==x1-1 for k in range(n-1))
                     and not any(p[k]==y2-1 and p[k+1]==x2-1 for k in range(n-1)))

        H_i = 1 + 2**(x1-y1-1)
        H_j = 1 + 2**(x2-y2-1)
        c2 = len(hps) - H_i - H_j + 1

        shared = set(tiles_1idx[ti]) & set(tiles_1idx[tj])

        label = chr(65+ti) + chr(65+tj)
        print(f"  {label} ({x1},{y1})×({x2},{y2}): "
              f"H={len(hps)} = 1 + {a_only}(A) + {b_only}(B) + {both}(AB) + {neither-1}(extra)")
        print(f"    c₂={c2}  |  a_only={a_only}  b_only={b_only}  both={both}  original={neither}")

        # KEY: c₂ = both - (a_only - (2^s1 - 1)) - (b_only - (2^s2 - 1))???
        # Actually: H = neither + a_only + b_only + both
        # H(0,0) = 1, H(1,0) = 1 + 2^(s1-1), H(0,1) = 1 + 2^(s2-1)
        # c₂ = H(1,1) - H(1,0) - H(0,1) + H(0,0)
        # = (neither + a_only + b_only + both) - (1 + 2^(s1-1)) - (1 + 2^(s2-1)) + 1
        # = neither + a_only + b_only + both - 1 - 2^(s1-1) - 2^(s2-1)
        # But neither >= 1 (original path always there)
        # a_only <= 2^(s1-1) (paths from flip 1 only)
        # b_only <= 2^(s2-1) (paths from flip 2 only)
        # So c₂ = (neither - 1) + (a_only - 2^(s1-1)) + (b_only - 2^(s2-1)) + both

        excess_neither = neither - 1
        excess_a = a_only - (2**(x1-y1-1))
        excess_b = b_only - (2**(x2-y2-1))
        print(f"    c₂ = {excess_neither}(extra_orig) + {excess_a}(extra_A) + {excess_b}(extra_B) + {both}(AB)")
        check = excess_neither + excess_a + excess_b + both
        print(f"    Check: {check} == {c2}? {check == c2}")
