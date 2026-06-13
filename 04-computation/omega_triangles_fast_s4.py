"""
Fast targeted check: does Omega(T) ever contain K3 (triangle subgraph)?
Session: opus-2026-05-28-S4

Key question: among all 3 pairwise-conflicting odd cycles in Omega(T),
does there always exist a 4th cycle? (i.e., is Omega K3-free?)

This is STRONGER than HYP-1748 (K3 as complete graph) — this asks
whether K3 can appear as an INDUCED SUBGRAPH when alpha_1 > 3.
"""

from collections import defaultdict
import itertools

def get_tiles(n):
    return [(x, y) for y in range(n-1) for x in range(y+2, n)]

def tiling_to_adj(n, tiles, bits):
    adj = [[0]*n for _ in range(n)]
    for k in range(1, n):
        adj[k][k-1] = 1
    for idx, (x, y) in enumerate(tiles):
        if (bits >> idx) & 1:
            adj[x][y] = 1
        else:
            adj[y][x] = 1
    return adj

def normalize_cycle(path):
    min_idx = path.index(min(path))
    return tuple(path[min_idx:] + path[:min_idx])

def find_all_odd_cycles(adj, n):
    odd_cycles = set()
    def dfs(start, current, path, in_path):
        for nxt in range(n):
            if not adj[current][nxt]:
                continue
            if nxt == start and len(path) >= 3:
                if len(path) % 2 == 1:
                    odd_cycles.add(normalize_cycle(list(path)))
            elif nxt not in in_path and nxt > start:
                in_path.add(nxt)
                path.append(nxt)
                dfs(start, nxt, path, in_path)
                path.pop()
                in_path.remove(nxt)
    for start in range(n):
        dfs(start, start, [start], {start})
    return odd_cycles

def has_triangle_in_conflict(cycles):
    """Check if any 3 cycles are mutually conflicting (K3 in Omega)."""
    cyc_list = list(cycles)
    vsets = [set(c) for c in cyc_list]
    nc = len(cyc_list)
    for i in range(nc):
        for j in range(i+1, nc):
            if not (vsets[i] & vsets[j]):
                continue
            # i and j conflict, find k that conflicts with both
            for k in range(j+1, nc):
                if (vsets[i] & vsets[k]) and (vsets[j] & vsets[k]):
                    return True, (cyc_list[i], cyc_list[j], cyc_list[k])
    return False, None

# ============================================================
# MAIN QUESTION: Is Omega EVER triangle-free? Or does it always contain K3?
# ============================================================
print("=" * 60)
print("TRIANGLE CHECK: Does Omega(T) contain K3 subgraph?")
print("=" * 60)

first_triangle = None
for n in [5, 6, 7]:
    tiles = get_tiles(n)
    m = len(tiles)
    tri_count = 0

    for bits in range(1 << m):
        adj = tiling_to_adj(n, tiles, bits)
        cycles = find_all_odd_cycles(adj, n)
        if len(cycles) < 3:
            continue
        has_tri, tri_ex = has_triangle_in_conflict(cycles)
        if has_tri:
            tri_count += 1
            if first_triangle is None:
                first_triangle = (n, bits, tri_ex, len(cycles))

    print(f"n={n}: tournaments with K3 in Omega: {tri_count} / {1<<m}")

if first_triangle is None:
    print("\n*** OMEGA(T) IS ALWAYS TRIANGLE-FREE (K3-free)! ***")
    print("This is a MAJOR structural result for TRRT!")
else:
    n0, bits0, (C1, C2, C3), nc0 = first_triangle
    print(f"\nFirst K3 found: n={n0}, alpha_1={nc0}")
    print(f"  C1={C1}")
    print(f"  C2={C2}")
    print(f"  C3={C3}")
    print(f"  Shared vertices: C1∩C2={set(C1)&set(C2)}, C1∩C3={set(C1)&set(C3)}, C2∩C3={set(C2)&set(C3)}")

# ============================================================
# KEY INSIGHT TEST: The proof shows 3 pairwise-conflicting cycles
# ALWAYS force a 4th. So K3 in Omega always has more cycles nearby.
# The question is whether K3 can exist as a SUBGRAPH (not the full Omega).
# ============================================================
print("\n" + "=" * 60)
print("SUBGRAPH vs FULL GRAPH DISTINCTION")
print("=" * 60)

print("""
THM-343 shows: Omega(T) cannot EQUAL K3 (cannot have exactly 3 mutually
conflicting cycles with no others).

But Omega(T) CAN contain K3 as a SUBGRAPH when alpha_1 >= 4.

Example: if 4 cycles C1,C2,C3 (mutually conflicting) and C4 (disjoint from all),
then Omega has vertices {C1,C2,C3,C4}, edges {C1-C2, C1-C3, C2-C3}.
This IS K3 + isolated vertex.

The question above checks for this general case (K3 as subgraph).
""")

# ============================================================
# PROOF OF ALGEBRAIC LEMMA: Three pairwise-intersecting 3-cycles
# Let's verify the specific case in the proof of THM-343.
# ============================================================
print("=" * 60)
print("ALGEBRAIC PROOF VERIFICATION")
print("=" * 60)

print("""
Proof of Key Lemma (Type 2 case):
C1=(a,b,c), C2=(a,d,e), C3=(b,d,f) with {a,b,c,d,e,f} distinct.

Arcs: a→b, b→c, c→a (C1); a→d, d→e, e→a (C2); b→d, d→f, f→b (C3).

Claim: the arc between c and d forces a 4th odd cycle.
  c→d: 5-cycle (c,d,e,a,b) formed.
  d→c: 3-cycle (d,c,a) formed.

Let's verify this with the tiling model at n=7:
Find all tournaments where {C1=(0,2,1), C2=(0,3,4), C3=(2,3,5)} all exist.
(a=0, b=2, c=1, d=3, e=4, f=5)
""")

n = 7
tiles = get_tiles(n)
m = len(tiles)

# Find tournaments where C1=(0,2,1), C2=(0,3,4), C3=(2,3,5) all exist
# (0,2,1) means cycle 0→2→1→0: arcs 0→2, 2→1, 1→0 (or normalize: 0→2,2→1,1→0)
# Actually normalizeation says (0,2,1) = starts at 0 (min), then 2, then 1.
# So the directed cycle is 0→2→1→0.

target_cycles = {
    (0, 2, 1),  # C1: 0→2→1→0
    (0, 3, 4),  # C2: 0→3→4→0? Wait: (0,3,4) normalized means 0→3→4→0
    (2, 3, 5),  # C3: 2→3→5→2? Normalized (2,3,5) means 2→3→5→2
}

# Actually, the cycle (a,b,c) means a→b→c→a.
# (0,2,1): 0→2→1→0
# (0,3,4): 0→3→4→0? But wait, normalized means minimum vertex first.
# (0,3,4): 0 is min, then 3, then 4. So 0→3→4→0.
# But for this to be in a tournament: 0→3, 3→4, 4→0.
# Hmm, but 4>0 so arc 4→0 means 4 beats 0, but 0 beats 3 (score 0 has high out-degree???)
# Let me re-examine. In our tiling model, vertices 0..n-1.
# The BASE PATH is n-1→n-2→...→1→0.
# So arc k→k-1 is always present in our model.
# "0→3" would be an arc from vertex 0 to vertex 3, which is a "tile" in our model
# (since 3 ≥ 0+2). This tile (3,0) has upward direction = 3→0 or downward = 0→3.
# Wait: tile (x,y) means x > y+1. Direction: upward = x→y, downward = y→x.
# So tile (3,0): upward = 3→0, downward = 0→3.
# For "0→3": we need tile (3,0) in downward direction = 0→3? No wait.
# In our model: tile (x,y) with x > y+1. The upward direction = x→y (going "up" the staircase),
# downward = y→x (going "down").
# So tile (3,0): upward means 3→0, downward means 0→3.
# For the cycle 0→2→1→0 (normalized as (0,2,1)):
#   0→2: from tile (2,0), need it to be in downward direction (0→2). Wait:
#         tile (2,0) upward = 2→0, downward = 0→2. So downward gives 0→2. bit=0.
#   2→1: from base path (2→1 is the arc 2→2-1=1). This is the base path arc, always present!
#   1→0: from base path (1→0). Always present.
# So cycle (0,2,1) = 0→2→1→0 requires tile (2,0) to be in downward position (bit for tile (2,0) = 0).

# For cycle (0,3,4):
# Normalized: 0→3→4→0. Wait, 4→0 is the arc from 4 to 0. But the base path has 4→3 and 3→2 and 2→1 and 1→0.
# 4→0 would require tile (4,0) upward direction? tile (4,0): upward = 4→0, downward = 0→4.
# Yes! tile (4,0) upward = 4→0. So arc 4→0 means tile (4,0) is in upward position.
# 0→3: tile (3,0) downward = 0→3.
# 3→4: tile (4,3) = WRONG! Wait, tiles are (x,y) with x≥y+2. So (4,3) has 4≥3+1=4, which means 4≥4 ✓.
#       But x=4, y=3: x-y=1 which is exactly 1, not ≥2. So (4,3) is NOT a tile! It's a base path arc!
# 3→4: base path has 4→3 (n-1→n-2→...→1→0), so 4→3 is always there. NOT 3→4.
# Hmm, so 3→4 (3 beats 4) is impossible in our tiling model since the base path has 4→3!
# Wait, I'm confused about the tiling model direction.

# Base path direction: n-1 → n-2 → ... → 1 → 0.
# So at n=7: 6→5→4→3→2→1→0.
# Arc k→k-1 is ALWAYS present (k beats k-1).
# So 4→3 is always present. 3→4 would require reversing this arc, but the base path arcs are fixed.

# This means the cycle (0,3,4,0) = 0→3→4→0 requires arc 3→4. But 3→4 is IMPOSSIBLE in our model!
# So cycle (0,3,4) cannot exist in our tiling model.

# I need to recheck which cycles can exist given the base path structure.

print("Checking which canonical 3-cycles can appear in n=7 tilings...")
n = 7
tiles_7 = get_tiles(n)
m7 = len(tiles_7)

# Sample some tilings to find 3-cycles
from collections import Counter
three_cycle_freq = Counter()
for bits in range(0, 1 << m7, 256):  # sample every 256th
    adj = tiling_to_adj(n, tiles_7, bits)
    cycles = find_all_odd_cycles(adj, n)
    for c in cycles:
        if len(c) == 3:
            three_cycle_freq[c] += 1

print("Most common 3-cycles in n=7 tilings:")
for cyc, cnt in three_cycle_freq.most_common(10):
    print(f"  {cyc}: {cnt} tilings")

# ============================================================
# QUICK SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("KEY RESULTS SUMMARY")
print("=" * 60)
print("""
1. H=7 IMPOSSIBLE (THM-343): For ALL n tested (n≤7, sampled n=8).
   - Algebraic proof for n=5,6 via score sequences.
   - At n≥7: alpha_1=3 exists but always has alpha_2≥1 (K3 impossible as full graph).
   - Core lemma: 3 pairwise-conflicting 3-cycles force a 4th (proved case-by-case).

2. Alpha_1=3 structure at n=7: Always H=15, always K1⊔K2 conflict graph.

3. HYP-1751 REFUTED for p=7: H mod 7 achieves all residues at n=7.
   The mod-5 gap at n=5 is specific to n=p=5 and the universal H=7 impossibility.

4. H mod p results:
   - All n: H ≡ 1 (mod 2) (H always odd — proved by Rédei)
   - n=5: H mod 5 misses {2} (because H=7 impossible AND H_max=15)
   - n=6: H mod 5 achieves all residues (H=17≡2 mod 5 is achievable)
   - n=7: H mod 7 achieves ALL residues {0,1,2,3,4,5,6}

5. Omega(T) triangle check: PENDING (too slow for n=7 exhaustive).
   Key theoretical question: is Omega always K3-free (as a subgraph)?
   Our proof of THM-343 only shows K3 can't be the ENTIRE Omega.
""")
