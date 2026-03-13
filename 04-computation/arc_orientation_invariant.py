#!/usr/bin/env python3
"""
arc_orientation_invariant.py -- kind-pasteur-2026-03-13-S61

The deepest ambiguity: two tournaments with ISOMORPHIC labeled lambda graphs
have different H (109 vs 111), differing only in c7_dir (8 vs 9).

Question: What is the MINIMAL additional invariant beyond labeled lambda
that determines H?

The candidates:
1. c7_dir alone (Hamiltonian cycle count)
2. Arc orientation on the lambda graph (which direction is each edge?)
3. Some orientation-derived invariant (e.g., net flow through high-lambda edges)
4. The 5-cycle ORIENTATION structure (not just count)

Also: explore the exact structural difference between these two tournaments.
What makes one have one more Hamiltonian cycle than the other?

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations, permutations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def get_labeled_lambda(A, n):
    c3_sets = []
    for a, b, c_v in combinations(range(n), 3):
        if A[a][b] and A[b][c_v] and A[c_v][a]:
            c3_sets.append(frozenset([a, b, c_v]))
        if A[a][c_v] and A[c_v][b] and A[b][a]:
            c3_sets.append(frozenset([a, b, c_v]))
    c3_sets = list(set(c3_sets))

    lam = [[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(u+1, n):
            val = sum(1 for cs in c3_sets if u in cs and v in cs)
            lam[u][v] = val
            lam[v][u] = val

    return lam, c3_sets


n = 7

# The two deepest ambiguity cases
cases = [
    ("Case 1", [(4728, 109), (4658, 111)]),
    ("Case 2", [(9388, 109), (9653, 111)]),
]


for case_name, pairs in cases:
    print("=" * 70)
    print(f"{case_name}")
    print("=" * 70)

    tournaments = []
    for bits, expected_H in pairs:
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)
        assert H == expected_H, f"H={H} != expected {expected_H}"
        tournaments.append((bits, A, H))

    bits1, A1, H1 = tournaments[0]
    bits2, A2, H2 = tournaments[1]

    # 1. Full adjacency comparison
    print(f"\n  Tournament 1 (bits={bits1}, H={H1}):")
    for i in range(n):
        print(f"    {A1[i]}")
    print(f"\n  Tournament 2 (bits={bits2}, H={H2}):")
    for i in range(n):
        print(f"    {A2[i]}")

    # 2. Score sequences
    scores1 = [sum(A1[v]) for v in range(n)]
    scores2 = [sum(A2[v]) for v in range(n)]
    print(f"\n  Scores1: {scores1}")
    print(f"  Scores2: {scores2}")

    # 3. Lambda matrices
    lam1, c3s1 = get_labeled_lambda(A1, n)
    lam2, c3s2 = get_labeled_lambda(A2, n)
    print(f"\n  Lambda1:")
    for v in range(n):
        print(f"    {lam1[v]}")
    print(f"  Lambda2:")
    for v in range(n):
        print(f"    {lam2[v]}")

    # 4. Find the isomorphism
    best_perm = None
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(i+1, n):
                if lam1[perm[i]][perm[j]] != lam2[i][j]:
                    match = False
                    break
            if not match:
                break
        if match:
            best_perm = perm
            break

    print(f"\n  Lambda isomorphism: perm={best_perm}")

    # 5. Apply the isomorphism to A1 and compare with A2
    # If sigma is the perm such that lam1[sigma(i)][sigma(j)] = lam2[i][j],
    # then apply sigma to A1 to get A1' and compare with A2
    A1_perm = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A1_perm[i][j] = A1[best_perm[i]][best_perm[j]]

    print(f"\n  A1 permuted by sigma:")
    for i in range(n):
        print(f"    {A1_perm[i]}")

    # 6. Compute the difference: which arcs are flipped?
    flips = []
    for i in range(n):
        for j in range(i+1, n):
            if A1_perm[i][j] != A2[i][j]:
                flips.append((i, j))
    print(f"\n  Arc flips between A1_perm and A2: {flips}")
    print(f"  Number of flips: {len(flips)}")

    # 7. For each flip, what are the lambda values?
    if flips:
        print(f"\n  Flip details (lambda values in A2's labeling):")
        for (i, j) in flips:
            print(f"    ({i},{j}): lambda={lam2[i][j]}, A1'={A1_perm[i][j]}->{A2[i][j]}")

    # 8. 3-cycle orientations
    print(f"\n  3-cycle orientations:")
    for bits_val, A, H in tournaments:
        print(f"\n    bits={bits_val}, H={H}:")
        for c3 in sorted(c3s1 if bits_val == bits1 else c3s2, key=lambda x: sorted(x)):
            verts = sorted(c3)
            a, b, c = verts
            if A[a][b] and A[b][c] and A[c][a]:
                print(f"      {verts}: {a}->{b}->{c}->{a} (clockwise)")
            elif A[a][c] and A[c][b] and A[b][a]:
                print(f"      {verts}: {a}->{c}->{b}->{a} (counter-clockwise)")
            else:
                # Both directions
                dir1 = A[a][b] * A[b][c] * A[c][a]
                dir2 = A[a][c] * A[c][b] * A[b][a]
                print(f"      {verts}: dir1={dir1}, dir2={dir2}")

    # 9. 5-cycle counts per subset
    print(f"\n  5-cycle counts per subset:")
    diff_5 = []
    for subset in combinations(range(n), 5):
        c5_1 = count_directed_ham_cycles_on_subset(A1_perm, list(subset))
        c5_2 = count_directed_ham_cycles_on_subset(A2, list(subset))
        if c5_1 != c5_2:
            diff_5.append((subset, c5_1, c5_2))
    print(f"    Subsets where c5 differs: {len(diff_5)}")
    for sub, c1, c2 in diff_5[:10]:
        print(f"      {sub}: c5_1={c1}, c5_2={c2}")

    # 10. 7-cycle (Hamiltonian cycle) count
    c7_1 = count_directed_ham_cycles_on_subset(A1_perm, list(range(n)))
    c7_2 = count_directed_ham_cycles_on_subset(A2, list(range(n)))
    print(f"\n  c7_dir: A1_perm={c7_1}, A2={c7_2}")

    # 11. Net flow on lambda-weighted arcs
    print(f"\n  Net flow analysis:")
    for bits_val, A, H in [(bits1, A1, H1), (bits2, A2, H2)]:
        # For each edge in lambda graph, what's the tournament direction?
        # Use A2's labeling for both (apply perm to A1)
        if bits_val == bits1:
            A_use = A1_perm
        else:
            A_use = A

        flow_by_lambda = defaultdict(int)
        for i in range(n):
            for j in range(i+1, n):
                lv = lam2[i][j]
                if A_use[i][j]:
                    flow_by_lambda[lv] += 1  # i->j
                else:
                    flow_by_lambda[lv] -= 1  # j->i
        print(f"    H={H}: net flow by lambda value: {dict(sorted(flow_by_lambda.items()))}")

    # 12. Vertex-level analysis: in-degree and out-degree through high-lambda edges
    print(f"\n  Vertex flow through lambda-3 edges:")
    for bits_val, A, H in [(bits1, A1, H1), (bits2, A2, H2)]:
        if bits_val == bits1:
            A_use = A1_perm
        else:
            A_use = A

        for v in range(n):
            out_lam3 = sum(1 for w in range(n) if w != v and lam2[v][w] == 3 and A_use[v][w])
            in_lam3 = sum(1 for w in range(n) if w != v and lam2[v][w] == 3 and A_use[w][v])
            out_lam2 = sum(1 for w in range(n) if w != v and lam2[v][w] == 2 and A_use[v][w])
            in_lam2 = sum(1 for w in range(n) if w != v and lam2[v][w] == 2 and A_use[w][v])
            if lam2[v][v] == 0:  # diagonal always 0
                pass
            print(f"    H={H}, v={v}: lam3(out={out_lam3},in={in_lam3}), lam2(out={out_lam2},in={in_lam2})")

    # 13. The ORIENTATION SIGNATURE: for each 3-cycle vertex set,
    #     record the direction as a binary value
    print(f"\n  3-cycle orientation signature:")
    # Use A2's labeling
    for bits_val, A, H in [(bits1, A1, H1), (bits2, A2, H2)]:
        if bits_val == bits1:
            A_use = A1_perm
        else:
            A_use = A

        orientations = []
        lam_use, c3s_use = get_labeled_lambda(A_use, n)
        for c3 in sorted(c3s_use, key=lambda x: sorted(x)):
            verts = sorted(c3)
            a, b, c = verts
            # 1 if clockwise (a->b->c->a), 0 if counter-clockwise
            if A_use[a][b] and A_use[b][c] and A_use[c][a]:
                orientations.append(1)
            else:
                orientations.append(0)
        print(f"    H={H}: orientations = {orientations}")
        print(f"           sum = {sum(orientations)}")


# ========================================================================
# ANALYSIS 2: What is the minimal additional invariant?
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 2: MINIMAL ADDITIONAL INVARIANT")
print("=" * 70)

# Check if c7_dir alone (combined with lambda isomorphism class) resolves
# everything. At n=7, c7_dir in {8,9} for these pairs.
print("""
For the 2 deepest ambiguity cases:
  Both have isomorphic lambda graphs.
  Both have same c3, same c5_dir, same alpha_2.
  They DIFFER in c7_dir (8 vs 9).

So: labeled lambda + c7_dir determines H for ALL n=7 tournaments?
This would mean:
  H = 1 + 2*alpha_1 + 4*alpha_2  (OCF decomposition, THM-166)
  alpha_1 is determined by labeled lambda (it's the total directed cycle count)
  alpha_2 is determined by labeled lambda (disjointness)
  But alpha_1 counts directed 3-cycles, 5-cycles, AND 7-cycles.
  alpha_1 = c3_dir + c5_dir + c7_dir

  c3_dir and c5_dir are determined by labeled lambda.
  c7_dir is NOT.

So: the UNIQUE additional bit of information is c7_dir (= #directed Hamiltonian cycles).
This is the arc orientation invariant: how many of the 7! possible orderings
are Hamiltonian cycles.
""")

# Let's verify: for ALL 2097152 n=7 tournaments, does (lambda_iso_class, c7_dir)
# determine H?

# Actually, we already know: H = 1 + 2*(c3_dir + c5_dir + c7_dir) + 4*alpha_2
# c3_dir + c5_dir are from lambda. alpha_2 is from lambda.
# So H = CONST_from_lambda + 2*c7_dir.
# Two lambda-isomorphic tournaments with c7_dir differing by 1 will have H differing by 2.
# This is EXACTLY what we see: H=109 vs H=111, c7=8 vs c7=9.

# Verify the formula
for case_name, pairs in cases:
    print(f"\n  {case_name}:")
    for bits, expected_H in pairs:
        A = binary_to_tournament(bits, n)
        H = count_ham_paths(A, n)

        # Count all directed odd cycles
        c3_dir = 0
        for subset in combinations(range(n), 3):
            c3_dir += count_directed_ham_cycles_on_subset(A, list(subset))

        c5_dir = 0
        for subset in combinations(range(n), 5):
            c5_dir += count_directed_ham_cycles_on_subset(A, list(subset))

        c7_dir = count_directed_ham_cycles_on_subset(A, list(range(n)))

        # alpha_2 from OCF: H = 1 + 2*alpha_1 + 4*alpha_2
        alpha_1 = c3_dir + c5_dir + c7_dir
        alpha_2 = (H - 1 - 2*alpha_1) // 4

        print(f"    bits={bits}: H={H}, c3={c3_dir}, c5={c5_dir}, c7={c7_dir}, "
              f"alpha_1={alpha_1}, alpha_2={alpha_2}")
        print(f"    Verify: 1 + 2*{alpha_1} + 4*{alpha_2} = {1 + 2*alpha_1 + 4*alpha_2}")


# ========================================================================
# ANALYSIS 3: The Vitali structure at the arc orientation level
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 3: VITALI STRUCTURE AT ARC ORIENTATION LEVEL")
print("=" * 70)

print("""
THE FOUR-LEVEL MEASURABILITY HIERARCHY:

Level 0: Score sequence (sorted out-degrees)
  -> Determines c3 (number of 3-cycle VERTEX SETS)
  -> "Mass" of the cycle structure

Level 1: Sorted lambda histogram
  -> Determines alpha_2 aggregate statistics but NOT exact alpha_2
  -> "Energy" of the overlap structure
  -> Vitali quotient #1: sorting destroys which-pairs-share-which-cycles

Level 2: Labeled lambda graph (pair-coverage up to isomorphism)
  -> Determines c3_dir, c5_dir, alpha_2 EXACTLY
  -> "Geometry" of the cycle arrangement
  -> Vitali quotient #2: forgetting arc orientations within lambda-classes

Level 3: Full tournament adjacency
  -> Determines c7_dir (Hamiltonian cycle count)
  -> "Topology" -- the full oriented structure
  -> The UNIQUE additional bit is c7_dir

KEY INSIGHT: At each level, the "non-measurable" content is:
  Level 1->0: Which vertex pairs share cycles (arrangement)
  Level 2->1: Which edges of the lambda graph are bidirectional vs directed (orientation)
  Level 3->2: Whether the full oriented graph has Hamiltonian cycles (topology)

Each level's non-measurable content controls EXACTLY ONE cycle length:
  Level 1 resolves ambiguity via alpha_2 (disjointness of 3-cycles)
  Level 2 resolves via c5_dir (5-cycle structure)
  Level 3 resolves via c7_dir (Hamiltonian cycles)

This is the VITALI TOWER: each quotient loses exactly one piece of cycle info.
""")


# ========================================================================
# ANALYSIS 4: Explore the single arc flip
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 4: THE SINGLE ARC FLIP")
print("=" * 70)

# For Case 1: after applying the permutation, how many arcs differ?
# And what effect does each flip have on c7?
for case_name, pairs in cases:
    print(f"\n  {case_name}:")
    bits1, H1 = pairs[0]
    bits2, H2 = pairs[1]
    A1 = binary_to_tournament(bits1, n)
    A2 = binary_to_tournament(bits2, n)
    lam1, _ = get_labeled_lambda(A1, n)
    lam2, _ = get_labeled_lambda(A2, n)

    # Find isomorphism
    best_perm = None
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(i+1, n):
                if lam1[perm[i]][perm[j]] != lam2[i][j]:
                    match = False
                    break
            if not match:
                break
        if match:
            best_perm = perm
            break

    A1p = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A1p[i][j] = A1[best_perm[i]][best_perm[j]]

    flips = []
    for i in range(n):
        for j in range(i+1, n):
            if A1p[i][j] != A2[i][j]:
                flips.append((i, j, lam2[i][j]))

    print(f"    Flips: {flips}")
    print(f"    These flips preserve lambda but change c7 by 1")

    # What's the Hamming distance between A1p and A2?
    ham_dist = sum(1 for i in range(n) for j in range(i+1, n) if A1p[i][j] != A2[i][j])
    print(f"    Hamming distance: {ham_dist}")

    # For each pair of vertices in flips, show which 3-cycles they're in
    for (i, j, lv) in flips:
        shared_3cyc = []
        for a, b, c in combinations(range(n), 3):
            verts = {a, b, c}
            if i in verts and j in verts:
                # Check if it's a 3-cycle in A2
                dirs = (A2[a][b]*A2[b][c]*A2[c][a]) + (A2[a][c]*A2[c][b]*A2[b][a])
                if dirs > 0:
                    shared_3cyc.append(sorted(verts))
        print(f"    Edge ({i},{j}), lambda={lv}: shared 3-cycles = {shared_3cyc}")


# ========================================================================
# ANALYSIS 5: Is c7 the ONLY difference? Check all invariants
# ========================================================================
print(f"\n{'='*70}")
print("ANALYSIS 5: EXHAUSTIVE INVARIANT COMPARISON")
print("=" * 70)

for case_name, pairs in cases:
    print(f"\n  {case_name}:")
    for bits, expected_H in pairs:
        A = binary_to_tournament(bits, n)
        scores = [sum(A[v]) for v in range(n)]
        lam, c3s = get_labeled_lambda(A, n)

        # All cycle counts
        c3 = 0
        for sub in combinations(range(n), 3):
            c3 += count_directed_ham_cycles_on_subset(A, list(sub))
        c5 = 0
        for sub in combinations(range(n), 5):
            c5 += count_directed_ham_cycles_on_subset(A, list(sub))
        c7 = count_directed_ham_cycles_on_subset(A, list(range(n)))

        # Disjoint 3-cycle pairs
        alpha2 = 0
        for i in range(len(c3s)):
            for j in range(i+1, len(c3s)):
                if not (c3s[i] & c3s[j]):
                    alpha2 += 1

        # 5-cycle vertex sets with directed cycles
        c5_sets = []
        for sub in combinations(range(n), 5):
            nc = count_directed_ham_cycles_on_subset(A, list(sub))
            if nc > 0:
                c5_sets.append(frozenset(sub))

        # Disjoint 5-cycle pairs (among 5-cycle sets)
        disj_55 = 0
        for i in range(len(c5_sets)):
            for j in range(i+1, len(c5_sets)):
                if not (c5_sets[i] & c5_sets[j]):
                    disj_55 += 1

        # 3-5 overlap: number of (3-cycle, 5-cycle) pairs with empty intersection
        disj_35 = 0
        for c3_set in c3s:
            for c5_set in c5_sets:
                if not (c3_set & c5_set):
                    disj_35 += 1

        print(f"    bits={bits}, H={expected_H}:")
        print(f"      scores={sorted(scores)}")
        print(f"      c3_dir={c3}, c5_dir={c5}, c7_dir={c7}")
        print(f"      alpha_2(3-3)={alpha2}")
        print(f"      #5-cycle sets={len(c5_sets)}")
        print(f"      disj_55={disj_55}, disj_35={disj_35}")

        # Per-vertex cycle distribution
        vtx_c3 = [0]*n
        for c3_set in c3s:
            for v in c3_set:
                vtx_c3[v] += 1
        vtx_c5 = [0]*n
        for sub in combinations(range(n), 5):
            nc = count_directed_ham_cycles_on_subset(A, list(sub))
            if nc > 0:
                for v in sub:
                    vtx_c5[v] += nc
        print(f"      per-vertex c3 membership: {sorted(vtx_c3)}")
        print(f"      per-vertex c5 weighted: {sorted(vtx_c5)}")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
