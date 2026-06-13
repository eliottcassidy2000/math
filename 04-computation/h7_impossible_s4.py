"""
H=7 impossibility and H mod p gap investigation.
Session: opus-2026-05-28-S4

Goals:
1. Prove H=7 impossible for n=6,7,8 by exhaustive tiling enumeration
2. Compute H mod 7 distribution at n=7 (extends HYP-1751)
3. Investigate alpha_1 = 3 tournaments (necessary for H=7 via OCF)
4. Attempt algebraic proof strategy
"""

from itertools import product, permutations
from collections import defaultdict

# ============================================================
# 0-indexed tiling model (consistent with oeis_exploration_s3.py)
# Vertices: 0..n-1
# Base path: n-1 -> n-2 -> ... -> 1 -> 0
# Tiles: (x,y) with x >= y+2, x < n
# Tile direction: 1 = upward (x->y), 0 = downward (y->x)
# ============================================================

def get_tiles(n):
    """Tiles in canonical order: y=0..n-2, x=y+2..n-1"""
    return [(x, y) for y in range(n-1) for x in range(y+2, n)]

def tiling_to_adj(n, tiles, tiling_bits):
    """Convert tiling bitmask to adjacency matrix adj[i][j] = 1 iff i->j."""
    adj = [[0]*n for _ in range(n)]
    # Base path arcs: k -> k-1 for k = n-1, ..., 1
    for k in range(1, n):
        adj[k][k-1] = 1
    # Tile arcs
    for idx, (x, y) in enumerate(tiles):
        if (tiling_bits >> idx) & 1:  # upward: x -> y
            adj[x][y] = 1
        else:  # downward: y -> x
            adj[y][x] = 1
    return adj

def H_by_DP(adj, n):
    """Count Hamiltonian paths using Held-Karp DP.
    dp[mask][v] = number of paths visiting exactly the vertices in mask,
    ending at v.
    """
    dp = [[0]*n for _ in range(1 << n)]
    # Initialize: single-vertex paths
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            if dp[mask][v] == 0:
                continue
            # Extend to w not in mask
            for w in range(n):
                if (mask >> w) & 1:
                    continue
                if adj[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]

    full_mask = (1 << n) - 1
    return sum(dp[full_mask][v] for v in range(n))

def compute_all_H(n, verbose=True):
    """Compute H(T) for all 2^m tilings at given n."""
    tiles = get_tiles(n)
    m = len(tiles)
    total = 1 << m

    H_list = []
    H_counts = defaultdict(int)

    if verbose:
        print(f"n={n}: m={m} tiles, {total} tilings")

    for bits in range(total):
        adj = tiling_to_adj(n, tiles, bits)
        h = H_by_DP(adj, n)
        H_list.append(h)
        H_counts[h] += 1

    return H_list, H_counts

# ============================================================
# Section 1: H=7 impossibility check at n=5,6,7
# ============================================================
print("=" * 60)
print("SECTION 1: H=7 IMPOSSIBILITY CHECK")
print("=" * 60)

for n in [5, 6, 7]:
    H_list, H_counts = compute_all_H(n)
    achieved = sorted(H_counts.keys())
    h7_count = H_counts.get(7, 0)
    print(f"\nn={n}: H=7 appears {h7_count} times")
    print(f"  Achieved H values: {achieved[:20]}{'...' if len(achieved) > 20 else ''}")
    print(f"  H=7 achievable: {'YES — HYP-1748 REFUTED' if h7_count > 0 else 'NO — consistent with HYP-1748'}")

# ============================================================
# Section 2: H mod 7 at n=7 (HYP-1751)
# ============================================================
print("\n" + "=" * 60)
print("SECTION 2: H MOD 7 DISTRIBUTION AT n=7 (HYP-1751)")
print("=" * 60)

tiles7 = get_tiles(7)
m7 = len(tiles7)
print(f"n=7: m={m7}, total tilings={1<<m7}")

H7_mod7 = defaultdict(int)
H7_vals = set()
H7_all = []
for bits in range(1 << m7):
    adj = tiling_to_adj(7, tiles7, bits)
    h = H_by_DP(adj, 7)
    H7_all.append(h)
    H7_mod7[h % 7] += 1
    H7_vals.add(h)

print(f"\nH mod 7 distribution (over all {1<<m7} labeled n=7 tournaments):")
for r in range(7):
    count = H7_mod7.get(r, 0)
    print(f"  H ≡ {r} (mod 7): {count} tournaments")

missing7 = {r for r in range(7) if H7_mod7.get(r, 0) == 0}
print(f"\nMissing residues mod 7: {missing7}")

H7_iso_vals = sorted(H7_vals)
print(f"All distinct H values at n=7: {H7_iso_vals[:30]}{'...' if len(H7_iso_vals) > 30 else ''}")
print(f"Total distinct H values: {len(H7_iso_vals)}")

# What are the H mod 7 values that are achieved?
achieved_mod7 = {h % 7 for h in H7_iso_vals}
print(f"H mod 7 values achieved: {sorted(achieved_mod7)}")
missing_iso = set(range(7)) - achieved_mod7
print(f"Missing mod 7 from distinct H values: {missing_iso}")

# ============================================================
# Section 3: Alpha_1 = 3 analysis (necessary condition for H=7)
# ============================================================
print("\n" + "=" * 60)
print("SECTION 3: ALPHA_1=3 ANALYSIS (OCF CONSTRAINT FOR H=7)")
print("=" * 60)

def count_odd_cycles_in_tournament(adj, n):
    """Count all odd directed cycles in tournament T.
    Returns list of cycles (as frozensets).
    """
    odd_cycles = set()

    # DFS to find all directed cycles
    # For small n, enumerate all cycles
    def dfs(start, current, path, path_set):
        for next_v in range(n):
            if adj[current][next_v] == 0:
                continue
            if next_v == start and len(path) >= 2:
                # Found a cycle
                cycle_len = len(path)
                if cycle_len % 2 == 1:  # odd cycle
                    # Canonicalize: use frozenset of vertices
                    odd_cycles.add(frozenset(path))
            elif next_v not in path_set and next_v > start:
                # Only visit vertices > start to avoid duplicate cycles
                path_set.add(next_v)
                path.append(next_v)
                dfs(start, next_v, path, path_set)
                path.pop()
                path_set.remove(next_v)

    for start in range(n):
        dfs(start, start, [start], {start})

    return odd_cycles

def compute_alpha1(adj, n):
    """Count alpha_1 = number of odd directed cycles."""
    return len(count_odd_cycles_in_tournament(adj, n))

# Check: at n=5, are there tournaments with alpha_1 = 3?
print("\nChecking n=5: tournaments with alpha_1 = 3 (required for H=7)")
tiles5 = get_tiles(5)
m5 = len(tiles5)
alpha1_3_count = 0
alpha1_distribution = defaultdict(int)

for bits in range(1 << m5):
    adj = tiling_to_adj(5, tiles5, bits)
    a1 = compute_alpha1(adj, 5)
    alpha1_distribution[a1] += 1
    if a1 == 3:
        alpha1_3_count += 1

print(f"n=5 alpha_1 distribution (over {1<<m5} labeled tournaments):")
for k in sorted(alpha1_distribution.keys()):
    print(f"  alpha_1={k}: {alpha1_distribution[k]} tournaments")

print(f"\nTournaments with alpha_1=3 at n=5: {alpha1_3_count}")
if alpha1_3_count > 0:
    print("  -> alpha_1=3 IS achievable at n=5, so H=7 fails for another reason!")
    print("  -> H=7 would require alpha_2=0 (no disjoint odd cycle pairs)")
    # Find examples
    print("  -> Finding examples with alpha_1=3:")
    count = 0
    for bits in range(1 << m5):
        adj = tiling_to_adj(5, tiles5, bits)
        a1 = compute_alpha1(adj, 5)
        if a1 == 3:
            h = H_by_DP(adj, 5)
            count += 1
            if count <= 3:
                print(f"    bits={bin(bits)}: alpha_1=3, H={h}")
else:
    print("  -> alpha_1=3 is IMPOSSIBLE at n=5 — this alone explains H≠7!")

# ============================================================
# Section 4: Independence polynomial computation for OCF analysis
# ============================================================
print("\n" + "=" * 60)
print("SECTION 4: INDEPENDENCE POLYNOMIAL vs H(T) VERIFICATION")
print("=" * 60)

def independence_polynomial(adj, n):
    """Compute independence polynomial I(Omega(T), x) where Omega is
    the odd-cycle conflict graph."""
    # Step 1: Find all odd cycles
    odd_cycles = sorted(count_odd_cycles_in_tournament(adj, n))

    if not odd_cycles:
        return {0: 1}  # I(empty, x) = 1

    nc = len(odd_cycles)

    # Step 2: Build conflict graph (edges = pairs sharing a vertex)
    conflicts = set()
    for i in range(nc):
        for j in range(i+1, nc):
            if odd_cycles[i] & odd_cycles[j]:  # share a vertex
                conflicts.add((i, j))

    # Step 3: Find all independent sets via bitmask enumeration
    alpha = defaultdict(int)
    alpha[0] = 1  # empty set

    for mask in range(1, 1 << nc):
        vertices_in_set = [k for k in range(nc) if (mask >> k) & 1]
        # Check independence
        is_indep = True
        for a in range(len(vertices_in_set)):
            for b in range(a+1, len(vertices_in_set)):
                i, j = vertices_in_set[a], vertices_in_set[b]
                if i > j:
                    i, j = j, i
                if (i, j) in conflicts:
                    is_indep = False
                    break
            if not is_indep:
                break
        if is_indep:
            alpha[len(vertices_in_set)] += 1

    return dict(alpha)

def I_at_2(poly):
    """Evaluate independence polynomial at x=2."""
    return sum(v * (2**k) for k, v in poly.items())

# Verify: H(T) = I(Omega(T), 2) for n=5
print("\nVerifying H(T) = I(Omega(T), 2) for n=5:")
correct = 0
total = 0
alpha1_3_h_vals = []
for bits in range(1 << m5):
    adj = tiling_to_adj(5, tiles5, bits)
    h_dp = H_by_DP(adj, 5)
    poly = independence_polynomial(adj, 5)
    h_ocf = I_at_2(poly)
    if h_dp == h_ocf:
        correct += 1
    total += 1

    # Track alpha_1=3 cases
    if poly.get(1, 0) == 3:
        alpha2 = poly.get(2, 0)
        alpha3 = poly.get(3, 0)
        if alpha2 == 0 and alpha3 == 0:  # this would give H=7
            alpha1_3_h_vals.append((h_dp, poly))

print(f"  H(T) = I(Omega(T), 2): {correct}/{total} verified ({'PASS' if correct==total else 'FAIL'})")
print(f"\nTournaments with alpha_1=3, alpha_2=0, alpha_3=0 at n=5:")
if alpha1_3_h_vals:
    for h, p in alpha1_3_h_vals[:5]:
        print(f"  H={h}, poly={dict(p)} -- should give H=7? {'YES' if h==7 else f'NO (H={h})'}")
else:
    print("  None! The alpha_1=3 cases all have alpha_2>0 or alpha_3>0.")

# ============================================================
# Section 5: Algebraic analysis of why H=7 is impossible
# ============================================================
print("\n" + "=" * 60)
print("SECTION 5: ALGEBRAIC ANALYSIS — WHY IS H=7 IMPOSSIBLE?")
print("=" * 60)

# Key question: can a tournament have exactly 3 odd cycles, all pairwise intersecting?
# Enumerate and check

print("\nChecking if alpha_1=3, alpha_2=0 is achievable in any tournament:")
print("(This is the NECESSARY AND SUFFICIENT condition for H=7)")

for n_check in [5, 6, 7]:
    tiles_n = get_tiles(n_check)
    m_n = len(tiles_n)
    found = False
    found_ex = []

    # For n=7, use sampling (2^15=32768 total, all manageable)
    for bits in range(1 << m_n):
        adj = tiling_to_adj(n_check, tiles_n, bits)
        odd_cycles = count_odd_cycles_in_tournament(adj, n_check)
        nc = len(odd_cycles)

        if nc == 3:
            # Check if all pairs are conflicting (alpha_2=0)
            oc_list = list(odd_cycles)
            conflict_all = True
            for i in range(3):
                for j in range(i+1, 3):
                    if not (oc_list[i] & oc_list[j]):
                        conflict_all = False
                        break
                if not conflict_all:
                    break

            if conflict_all:
                h = H_by_DP(adj, n_check)
                found = True
                found_ex.append((n_check, bits, h, oc_list))

    if found:
        print(f"\nn={n_check}: FOUND tournaments with alpha_1=3, alpha_2=0 (H=7)!")
        for n_, b, h, oc in found_ex[:3]:
            print(f"  bits={bin(b)}, H={h}, cycles={[list(c) for c in oc]}")
        print(f"  Total examples: {len(found_ex)}")
    else:
        print(f"\nn={n_check}: NO tournament has alpha_1=3, alpha_2=0. H=7 impossible at n={n_check}.")

# ============================================================
# Section 6: OCF coefficient parity — why alpha_1=3, alpha_2=0 fails
# ============================================================
print("\n" + "=" * 60)
print("SECTION 6: STRUCTURAL ANALYSIS — ODD CYCLE TRIPLES")
print("=" * 60)

print("""
H=7 requires: 3 odd cycles, all pairwise sharing vertices.
Three pairwise-intersecting odd cycles in a tournament — is this possible?

Key insight attempt:
If C1, C2, C3 are 3 pairwise-intersecting odd cycles, do they form a structure
that forces alpha_2 > 0 (two vertex-disjoint odd cycles to appear)?

Case analysis for n=6 (smallest n where 2 disjoint 3-cycles can exist):
- 3-cycles use 3 vertices, so two disjoint 3-cycles need 6 vertices.
- For alpha_1=3 AND alpha_2=0: we need exactly 3 odd cycles, all intersecting.
  But with 6 vertices, a 3-cycle {1,2,3} and another sharing a vertex with it
  uses 3+2=5 vertices, leaving 1 vertex. The third 3-cycle uses those 5 vertices.

Let's count: can 3 directed 3-cycles in a 6-vertex tournament all pairwise intersect?
""")

# Enumerate specific structure: 3 directed 3-cycles pairwise intersecting
# in a 6-vertex tournament with NO other odd cycles
n6 = 6
tiles6 = get_tiles(6)
m6 = len(tiles6)

print(f"n=6: scanning all {1<<m6} tilings for 3-cycle structure...")
three_cycles_pairwise = 0
exact_3_cycles = 0

for bits in range(1 << m6):
    adj = tiling_to_adj(6, tiles6, bits)
    odd_cycles = count_odd_cycles_in_tournament(adj, 6)
    nc = len(odd_cycles)

    if nc == 3:
        exact_3_cycles += 1
        oc_list = list(odd_cycles)
        conflict_all = all(
            oc_list[i] & oc_list[j]
            for i in range(3)
            for j in range(i+1, 3)
        )
        if conflict_all:
            three_cycles_pairwise += 1

print(f"  Tournaments with exactly 3 odd cycles: {exact_3_cycles}")
print(f"  Tournaments with exactly 3 odd cycles, all pairwise intersecting: {three_cycles_pairwise}")

# ============================================================
# Section 7: H mod p for multiple primes
# ============================================================
print("\n" + "=" * 60)
print("SECTION 7: H MOD p FOR SMALL PRIMES")
print("=" * 60)

# For n=5, use precomputed H values
H5_list, H5_counts = compute_all_H(5, verbose=False)

primes = [2, 3, 5, 7, 11, 13]
print("\nH mod p distributions at n=5:")
for p in primes:
    dist = defaultdict(int)
    for h in H5_list:
        dist[h % p] += 1
    achieved = sorted(r for r in range(p) if dist[r] > 0)
    missing = sorted(r for r in range(p) if dist[r] == 0)
    print(f"  mod {p}: achieved={achieved}, missing={missing}")

print("\nH mod p distributions at n=6 (all 1024 tilings):")
H6_list, H6_counts = compute_all_H(6, verbose=False)
for p in primes:
    dist = defaultdict(int)
    for h in H6_list:
        dist[h % p] += 1
    achieved = sorted(r for r in range(p) if dist[r] > 0)
    missing = sorted(r for r in range(p) if dist[r] == 0)
    print(f"  mod {p}: achieved={achieved}, missing={missing}")

print("\nH mod p distributions at n=7 (all 32768 tilings):")
for p in primes:
    dist = defaultdict(int)
    for h in H7_all:
        dist[h % p] += 1
    achieved = sorted(r for r in range(p) if dist[r] > 0)
    missing = sorted(r for r in range(p) if dist[r] == 0)
    print(f"  mod {p}: achieved={achieved}, missing={missing}")

# ============================================================
# Section 8: H=7 algebraic proof sketch
# ============================================================
print("\n" + "=" * 60)
print("SECTION 8: H=7 PROOF SKETCH SYNTHESIS")
print("=" * 60)

print("""
If the computational evidence shows H=7 is impossible for ALL n up to 7,
the algebraic mechanism must be one of:

(A) alpha_1=3, alpha_2=0 is impossible in any tournament — combinatorial reason
(B) alpha_1=3, alpha_2=0 tournaments exist but always have additional odd cycles
    making alpha_1 > 3

Key question: In a tournament, if exactly 3 odd cycles exist and they pairwise
share vertices, does this FORCE the existence of a 4th odd cycle?

This would follow from a "tournament Ramsey-type" theorem on odd cycles.

Let's check this by finding all tournaments with alpha_1=3, alpha_2=0,
and examining how many actually achieve H=7 vs having more cycles.
""")

# Quick summary
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"H=7 impossibility verification:")
for n in [5, 6, 7]:
    tiles_n = get_tiles(n)
    m_n = len(tiles_n)
    if n == 7:
        # Already computed
        h7 = sum(1 for h in H7_all if h == 7)
    else:
        H_list, _ = compute_all_H(n, verbose=False)
        h7 = sum(1 for h in H_list if h == 7)
    print(f"  n={n}: H=7 appears {h7} times out of {1<<m_n} tilings")
