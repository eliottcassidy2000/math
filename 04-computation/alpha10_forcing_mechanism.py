"""
alpha10_forcing_mechanism.py -- kind-pasteur-2026-03-14-S66

WHY does alpha_1=10 force alpha_2>=2?

Strategy: At n=6 (exhaustive), find ALL tournaments with alpha_1=10.
For each, analyze the cycle structure:
- How many vertex sets have cycles?
- What sizes are the cycles?
- Which vertex sets are they on?
- Why can't we reduce alpha_2 to 0 or 1?

KEY INSIGHT TO INVESTIGATE: alpha_1=10 means exactly 10 directed odd cycles.
At n=6, the only odd cycle sizes are 3 and 5.
If c3 = number of 3-cycle vertex sets, c5 = number of 5-cycle vertex sets,
then directed cycles come from these vertex sets.
Each 3-cycle vertex set has exactly 1 directed cycle.
Each 5-cycle vertex set has 0, 1, 2, or 3 directed Hamiltonian cycles.

So alpha_1 = c3_directed + sum(dc5_i) where dc5_i is the number of
directed Hamiltonian cycles in the i-th 5-cycle vertex set.

The constraint alpha_2=0 means ALL pairs of directed cycles share a vertex.
For 3-cycles on 3 vertices in a 6-vertex tournament, two 3-cycles are
disjoint iff their vertex sets are disjoint — and 3+3=6 exactly covers all vertices.

This is the KEY: at n=6, any two vertex-disjoint 3-cycles partition all 6 vertices.
"""

import numpy as np
from itertools import combinations
from math import comb

def tournament_from_bits(n, bits):
    A = np.zeros((n, n), dtype=np.int8)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[j][i] = 1
            else:
                A[i][j] = 1
            idx += 1
    return A

def count_directed_hamcycles(A, vertices):
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0
    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = int(A[vlist[i]][vlist[j]])
    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def analyze_tournament(A, n):
    """Full cycle analysis of a tournament."""
    cycles_3 = []  # (vertex_set, directed_count)
    cycles_5 = []

    for subset in combinations(range(n), 3):
        cnt = count_directed_hamcycles(A, list(subset))
        if cnt > 0:
            cycles_3.append((frozenset(subset), cnt))

    for subset in combinations(range(n), 5):
        cnt = count_directed_hamcycles(A, list(subset))
        if cnt > 0:
            cycles_5.append((frozenset(subset), cnt))

    # Compute alpha_1 and alpha_2
    all_cycles = cycles_3 + cycles_5
    alpha_1 = sum(cnt for _, cnt in all_cycles)

    # Alpha_2: count disjoint DIRECTED cycle pairs
    alpha_2 = 0
    disjoint_pairs = []
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            vs_i, cnt_i = all_cycles[i]
            vs_j, cnt_j = all_cycles[j]
            if len(vs_i & vs_j) == 0:
                alpha_2 += cnt_i * cnt_j
                disjoint_pairs.append((i, j, cnt_i * cnt_j))

    return {
        'cycles_3': cycles_3,
        'cycles_5': cycles_5,
        'alpha_1': alpha_1,
        'alpha_2': alpha_2,
        'all_cycles': all_cycles,
        'disjoint_pairs': disjoint_pairs
    }

def main():
    print("=" * 70)
    print("WHY DOES alpha_1=10 FORCE alpha_2>=2?")
    print("=" * 70)

    # Part 1: At n=6, find all tournaments with alpha_1=10
    print("\n--- Part 1: All n=6 tournaments with alpha_1=10 ---")
    n = 6
    num_edges = n*(n-1)//2  # 15
    total = 1 << num_edges   # 32768

    alpha1_10_tours = []
    alpha1_counts = {}

    for bits in range(total):
        A = tournament_from_bits(n, bits)
        info = analyze_tournament(A, n)
        a1 = info['alpha_1']
        a2 = info['alpha_2']

        key = (a1, a2)
        alpha1_counts[key] = alpha1_counts.get(key, 0) + 1

        if a1 == 10:
            alpha1_10_tours.append((bits, info))

    print(f"Total tournaments: {total}")
    print(f"Tournaments with alpha_1=10: {len(alpha1_10_tours)}")

    # Show alpha_2 distribution for alpha_1=10
    a2_dist = {}
    for bits, info in alpha1_10_tours:
        a2 = info['alpha_2']
        a2_dist[a2] = a2_dist.get(a2, 0) + 1

    print(f"\nalpha_2 distribution when alpha_1=10:")
    for a2 in sorted(a2_dist.keys()):
        print(f"  alpha_2={a2}: {a2_dist[a2]} tournaments (H={1+20+4*a2})")

    # Part 2: Detailed analysis of alpha_1=10 tournaments
    print(f"\n--- Part 2: Cycle structure of alpha_1=10 tournaments ---")

    # Classify by (c3_sets, c5_sets, c3_directed, c5_directed)
    structure_types = {}
    for bits, info in alpha1_10_tours:
        c3_sets = len(info['cycles_3'])
        c5_sets = len(info['cycles_5'])
        c3_dir = sum(cnt for _, cnt in info['cycles_3'])
        c5_dir = sum(cnt for _, cnt in info['cycles_5'])
        key = (c3_sets, c5_sets, c3_dir, c5_dir)
        if key not in structure_types:
            structure_types[key] = []
        structure_types[key].append((bits, info))

    for key in sorted(structure_types.keys()):
        c3s, c5s, c3d, c5d = key
        examples = structure_types[key]
        print(f"\n  Type ({c3s} 3-sets, {c5s} 5-sets, {c3d} dir-3, {c5d} dir-5): "
              f"{len(examples)} tournaments")

        # Show first example in detail
        bits, info = examples[0]
        print(f"    Example: bits={bits}")
        print(f"    3-cycles: ", end="")
        for vs, cnt in info['cycles_3']:
            print(f"{set(vs)}({cnt}d) ", end="")
        print()
        if info['cycles_5']:
            print(f"    5-cycles: ", end="")
            for vs, cnt in info['cycles_5']:
                print(f"{set(vs)}({cnt}d) ", end="")
            print()

        # Show disjoint pair structure
        print(f"    alpha_2 = {info['alpha_2']}")
        if info['disjoint_pairs']:
            print(f"    Disjoint pairs ({len(info['disjoint_pairs'])}):")
            for i, j, product in info['disjoint_pairs']:
                vs_i = info['all_cycles'][i][0]
                vs_j = info['all_cycles'][j][0]
                print(f"      {set(vs_i)} x {set(vs_j)} -> {product} directed pairs")

    # Part 3: The key constraint — complementary 3-cycles
    print(f"\n{'='*70}")
    print("THE COMPLEMENTARY 3-CYCLE MECHANISM")
    print(f"{'='*70}")

    # At n=6, any two disjoint 3-vertex sets partition {0,1,2,3,4,5}
    # So a disjoint 3-cycle pair = complementary partition
    print("\nAt n=6: two disjoint 3-vertex sets partition all 6 vertices.")
    print("There are C(6,3)/2 = 10 complementary pairs of 3-vertex sets.")
    print()

    # For alpha_1=10 tournaments, how many complementary pairs are both cycles?
    comp_pair_data = {}
    for bits, info in alpha1_10_tours:
        # Get the 3-cycle vertex sets
        c3_sets = set(vs for vs, _ in info['cycles_3'])

        # Count complementary pairs where both are 3-cycles
        all_verts = frozenset(range(6))
        comp_both = 0
        for vs in c3_sets:
            comp = all_verts - vs
            if comp in c3_sets and vs < comp:  # avoid double count
                comp_both += 1

        if comp_both not in comp_pair_data:
            comp_pair_data[comp_both] = {'count': 0, 'a2_vals': []}
        comp_pair_data[comp_both]['count'] += 1
        comp_pair_data[comp_both]['a2_vals'].append(info['alpha_2'])

    print("Complementary 3-cycle pairs in alpha_1=10 tournaments:")
    for cp in sorted(comp_pair_data.keys()):
        a2_vals = comp_pair_data[cp]['a2_vals']
        print(f"  {cp} complementary cycle pairs: {comp_pair_data[cp]['count']} tournaments, "
              f"alpha_2 range [{min(a2_vals)}, {max(a2_vals)}]")

    # Part 4: Why 10 specifically?
    print(f"\n{'='*70}")
    print("WHY 10 = C(5,2) SPECIFICALLY?")
    print(f"{'='*70}")

    # At n=6: C(6,3) = 20 possible 3-vertex subsets,
    # C(6,5) = 6 possible 5-vertex subsets
    # Each 3-subset gives 0 or 1 directed 3-cycle
    # Each 5-subset gives 0, 1, 2, or 3 directed 5-cycles

    # With alpha_1=10, what are the possible decompositions?
    print("\nalpha_1 = 10 decompositions (c3_directed + c5_directed = 10):")
    for c3d in range(min(20, 10), -1, -1):
        c5d = 10 - c3d
        if c5d < 0:
            continue
        # c3d directed 3-cycles means c3d 3-vertex sets are 3-cycles
        # c5d directed 5-cycles can come from at most 6 5-vertex sets
        # Each 5-vertex set gives at most 3 directed 5-cycles (for n=5 regular)
        max_5sets = min(6, c5d)  # need at least c5d/3 5-sets
        min_5sets = (c5d + 2) // 3  # ceiling division
        if c5d > 0 and c5d <= 6*3:
            print(f"  c3_dir={c3d}, c5_dir={c5d}: valid if {min_5sets}-{max_5sets} 5-vertex sets")
        elif c5d == 0:
            print(f"  c3_dir={c3d}, c5_dir=0: needs exactly {c3d} 3-cycles")

    # What actually occurs?
    print("\nActual decompositions observed:")
    for key in sorted(structure_types.keys()):
        c3s, c5s, c3d, c5d = key
        print(f"  c3_dir={c3d}, c5_dir={c5d} ({c3s} 3-sets, {c5s} 5-sets): "
              f"{len(structure_types[key])} tournaments")

    # Part 5: The 10 = C(5,2) connection
    print(f"\n{'='*70}")
    print("THE PIGEONHOLE ARGUMENT")
    print(f"{'='*70}")
    print()
    print("At n=6, there are C(6,3)/2 = 10 complementary pairs of 3-sets.")
    print("CLAIM: If a tournament has c3 >= 8 3-cycle vertex sets,")
    print("then some complementary pair are both 3-cycles, giving alpha_2 >= 1.")
    print()

    # Verify: for each c3 value, what's the min number of complementary cycle pairs?
    print("Min complementary cycle pairs by c3 (3-cycle vertex set count), n=6:")
    c3_to_min_comp = {}
    c3_to_a2_range = {}

    for bits in range(total):
        A = tournament_from_bits(n, bits)
        # Fast c3 count
        c3_sets = set()
        for a, b, c in combinations(range(n), 3):
            if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                c3_sets.add(frozenset([a, b, c]))

        c3 = len(c3_sets)

        # Count complementary pairs
        all_v = frozenset(range(6))
        comp = sum(1 for vs in c3_sets if (all_v - vs) in c3_sets and vs < (all_v - vs))

        if c3 not in c3_to_min_comp or comp < c3_to_min_comp[c3]:
            c3_to_min_comp[c3] = comp
        if c3 not in c3_to_a2_range:
            c3_to_a2_range[c3] = [comp, comp]
        else:
            c3_to_a2_range[c3][0] = min(c3_to_a2_range[c3][0], comp)
            c3_to_a2_range[c3][1] = max(c3_to_a2_range[c3][1], comp)

    for c3 in sorted(c3_to_min_comp.keys()):
        lo, hi = c3_to_a2_range[c3]
        print(f"  c3={c3:2d}: min comp pairs = {c3_to_min_comp[c3]}, range [{lo}, {hi}]")

    # Part 6: The splicing argument revisited
    print(f"\n{'='*70}")
    print("THE VERTEX-COVER ARGUMENT FOR alpha_1=10")
    print(f"{'='*70}")
    print()
    print("Each 3-cycle uses 3 of 6 vertices. Each vertex appears in C(5,2)=10")
    print("of the C(6,3)=20 possible 3-vertex subsets.")
    print()
    print("If c3 3-vertex subsets are 3-cycles, then each of the 10 complementary")
    print("pairs (S, S^c) has probability roughly (c3/20)^2 of both being cycles.")
    print("Expected complementary cycle pairs ~ 10 * (c3/20)^2.")
    print()
    print("For alpha_1=10 with only 3-cycles (c5=0):")
    print("  c3=10 directed 3-cycles means c3=10 vertex sets are 3-cycles.")
    print("  Expected comp pairs ~ 10 * (10/20)^2 = 10/4 = 2.5")
    print("  MINIMUM comp pairs = ?")
    print()

    # Check: at c3=10, what's the minimum complementary pairs?
    if 10 in c3_to_min_comp:
        print(f"  VERIFIED: c3=10 has min comp pairs = {c3_to_min_comp[10]}")

    # Part 7: Graph coloring perspective
    print(f"\n{'='*70}")
    print("COMBINATORIAL LOWER BOUND ON COMPLEMENTARY PAIRS")
    print(f"{'='*70}")
    print()
    print("Model: 20 3-vertex subsets of {0,...,5}, paired into 10 complementary pairs.")
    print("Choose c3 of the 20 subsets. How few complementary pairs can be hit?")
    print()

    # This is equivalent to: given 10 edges (complementary pairs) in a graph on 20 vertices,
    # choose c3 vertices. Min number of edges with both endpoints chosen.
    # By inclusion-exclusion: if we choose c3 subsets, the complementary graph has 10 edges.
    # An independent set has no complementary pairs. Max independent set = 10 (one from each pair).
    # So for c3 > 10, at least c3-10 complementary pairs must be hit.

    for c3 in range(0, 21):
        # Lower bound: max(0, c3 - 10) pairs minimum
        # Because each pair can contribute at most 1 to the independent part
        lb = max(0, c3 - 10)
        actual = c3_to_min_comp.get(c3, '?')
        print(f"  c3={c3:2d}: lower bound = max(0, {c3}-10) = {lb}, actual min = {actual}")

    print()
    print("KEY INSIGHT: the lower bound max(0, c3-10) comes from")
    print("the 10 complementary pairs being a PERFECT MATCHING on 20 items.")
    print("Any selection of c3 > 10 items from 20 must include both")
    print("endpoints of at least c3-10 matched pairs.")
    print()
    print("For alpha_1=10 with ONLY 3-cycles: c3=10 of the 20 subsets.")
    print("  => at least max(0, 10-10)=0 complementary pairs FORCED.")
    print("  But the ACTUAL minimum is higher because tournament structure constrains")
    print("  which 10 of the 20 subsets can be 3-cycles.")
    print()
    print("CRITICAL: Tournament structure prevents choosing an independent set!")
    print("The achievable 3-cycle configurations are NOT arbitrary subsets.")
    print("Regular tournament at n=5 (unique): c3=2, so 3-cycle vertex sets")
    print("are highly constrained by the tournament orientation.")

    # Part 8: Score sequence constraint
    print(f"\n{'='*70}")
    print("SCORE SEQUENCE AND 3-CYCLE COUNT")
    print(f"{'='*70}")

    # At n=6, score sequences and c3 values
    score_c3 = {}
    for bits in range(total):
        A = tournament_from_bits(n, bits)
        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        c3 = 0
        for a, b, c_v in combinations(range(n), 3):
            if (A[a][b] and A[b][c_v] and A[c_v][a]) or \
               (A[a][c_v] and A[c_v][b] and A[b][a]):
                c3 += 1
        if scores not in score_c3:
            score_c3[scores] = set()
        score_c3[scores].add(c3)

    print("\nScore sequences and possible c3 values at n=6:")
    for scores in sorted(score_c3.keys()):
        c3_vals = sorted(score_c3[scores])
        # c3 formula: c3 = C(n,3) - sum C(s_i, 2) = 20 - sum s_i*(s_i-1)/2
        c3_formula = 20 - sum(s*(s-1)//2 for s in scores)
        print(f"  scores={scores}: c3={c3_formula} (formula), observed {c3_vals}")

    print()
    print("OBSERVATION: c3 is UNIQUELY determined by score sequence!")
    print("c3 = C(n,3) - sum C(s_i, 2)  (classical result)")
    print()
    print("For c3=10 at n=6: need sum C(s_i,2) = 20-10 = 10")
    print("Score sequences with sum C(s_i,2) = 10:")
    for scores in sorted(score_c3.keys()):
        c3_val = 20 - sum(s*(s-1)//2 for s in scores)
        if c3_val == 10:
            print(f"  {scores}: c3={c3_val}")

    # Part 9: The definitive mechanism
    print(f"\n{'='*70}")
    print("DEFINITIVE MECHANISM: WHY alpha_1=10 FORCES alpha_2>=2")
    print(f"{'='*70}")

    # For alpha_1=10 tournaments at n=6, show the full accounting
    if alpha1_10_tours:
        # Group by score sequence
        score_groups = {}
        for bits, info in alpha1_10_tours:
            A = tournament_from_bits(n, bits)
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            if scores not in score_groups:
                score_groups[scores] = []
            score_groups[scores].append((bits, info))

        for scores, group in sorted(score_groups.items()):
            c3 = 20 - sum(s*(s-1)//2 for s in scores)
            print(f"\n  Score {scores} (c3={c3}): {len(group)} tournaments with alpha_1=10")

            # For this score, what's the directed 5-cycle contribution?
            c5_data = {}
            for bits, info in group:
                c5_dir = sum(cnt for _, cnt in info['cycles_5'])
                needed_from_3 = 10 - c5_dir
                key = (c3, c5_dir, needed_from_3)
                if key not in c5_data:
                    c5_data[key] = 0
                c5_data[key] += 1

            for (c3v, c5d, need3), count in sorted(c5_data.items()):
                print(f"    c3={c3v} vertex sets, {c5d} directed 5-cycles, "
                      f"need {need3} directed 3-cycles: {count} tournaments")
                print(f"    But c3 3-cycles means {c3v} directed 3-cycles (each vertex set gives exactly 1).")
                print(f"    So alpha_1 = {c3v} + {c5d} = {c3v+c5d}")

            # Show WHY alpha_2>=2 for this score
            print(f"\n    COMPLEMENTARY PAIR ANALYSIS:")
            all_v = frozenset(range(6))
            for bits, info in group[:3]:  # First 3 examples
                c3_sets = set(vs for vs, _ in info['cycles_3'])
                comp_pairs = [(vs, all_v-vs) for vs in c3_sets
                              if (all_v-vs) in c3_sets and vs < (all_v-vs)]
                print(f"    bits={bits}: {len(c3_sets)} 3-cycle sets, "
                      f"{len(comp_pairs)} complementary pairs -> alpha_2={info['alpha_2']}")
                for vs, comp in comp_pairs:
                    print(f"      {set(vs)} <-> {set(comp)}")

if __name__ == "__main__":
    main()
