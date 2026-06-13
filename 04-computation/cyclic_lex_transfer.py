"""
cyclic_lex_transfer.py -- kind-pasteur-2026-03-14-S82
Understanding WHY cyclic lex products give huge H values.

KEY QUESTION: In T_3 lex T_2 (6 vertices), why H=45 instead of 3?

The 3-cycle lex transitive-2:
Vertices: (0,0),(0,1), (1,0),(1,1), (2,0),(2,1)
Arcs:
  Within group 0: (0,0)->(0,1) [from T2]
  Within group 1: (1,0)->(1,1) [from T2]
  Within group 2: (2,0)->(2,1) [from T2]
  Between groups: T1 = 3-cycle: 0->1, 1->2, 2->0
    So: (0,j)->(1,j'), (1,j)->(2,j'), (2,j)->(0,j') for ALL j,j'

The lex product makes ALL cross-group arcs follow T1.
So the tournament looks like:
  Group 0 beats Group 1 (all 4 arcs: (0,0)->(1,0), (0,0)->(1,1), (0,1)->(1,0), (0,1)->(1,1))
  Group 1 beats Group 2 (all 4 arcs)
  Group 2 beats Group 0 (all 4 arcs)
  Within each group: one directed arc.

A Hamiltonian path visits all 6 vertices. It can:
- Stay within a group (using the internal arc)
- Jump between groups (using cross arcs)

The 45 paths arise from the ENORMOUS number of valid orderings.

TRANSFER MATRIX APPROACH:
Define the "state" as (current group, set of visited vertices in current group).
The transfer matrix M captures transitions between groups.

Actually simpler: enumerate ALL 45 paths and classify them.
"""

import numpy as np
from itertools import permutations
from collections import Counter, defaultdict
import sys

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("CYCLIC LEX PRODUCT — PATH ANATOMY")
    print("kind-pasteur-2026-03-14-S82")
    print("=" * 70)

    # Build T_3 lex T_2
    # Vertices: 0=(0,0), 1=(0,1), 2=(1,0), 3=(1,1), 4=(2,0), 5=(2,1)
    n = 6
    # Group membership
    group = [0, 0, 1, 1, 2, 2]
    # Within group: T2 = transitive (0->1)
    # (0,0)->(0,1): vertex 0->1
    # (1,0)->(1,1): vertex 2->3
    # (2,0)->(2,1): vertex 4->5
    # Between groups: 3-cycle 0->1->2->0
    # Group 0 beats Group 1: 0->2, 0->3, 1->2, 1->3
    # Group 1 beats Group 2: 2->4, 2->5, 3->4, 3->5
    # Group 2 beats Group 0: 4->0, 4->1, 5->0, 5->1

    A = np.zeros((n, n), dtype=int)
    # Internal
    A[0][1] = 1  # (0,0)->(0,1)
    A[2][3] = 1  # (1,0)->(1,1)
    A[4][5] = 1  # (2,0)->(2,1)
    # Cross: 0->1
    A[0][2] = A[0][3] = A[1][2] = A[1][3] = 1
    # Cross: 1->2
    A[2][4] = A[2][5] = A[3][4] = A[3][5] = 1
    # Cross: 2->0
    A[4][0] = A[4][1] = A[5][0] = A[5][1] = 1

    print(f"\n  Adjacency matrix:")
    for i in range(n):
        print(f"    {list(A[i])}")

    scores = A.sum(axis=1)
    print(f"  Scores: {list(scores.astype(int))}")

    # Enumerate ALL Hamiltonian paths
    paths = []
    for perm in permutations(range(n)):
        valid = all(A[perm[i]][perm[i+1]] for i in range(n-1))
        if valid:
            paths.append(perm)

    print(f"\n  Total Ham paths: {len(paths)}")

    # Classify paths by their group sequence
    group_seq_count = Counter()
    for path in paths:
        gseq = tuple(group[v] for v in path)
        group_seq_count[gseq] += 1

    print(f"\n  PATH CLASSIFICATION BY GROUP SEQUENCE:")
    print(f"  (The group each vertex belongs to, in path order)")
    for gseq, count in sorted(group_seq_count.items(), key=lambda x: -x[1]):
        # Check: how many group switches?
        switches = sum(1 for i in range(5) if gseq[i] != gseq[i+1])
        # Is it "block" (all of one group, then all of next)?
        # blocks = set of contiguous group runs
        runs = []
        for i in range(6):
            if i == 0 or gseq[i] != gseq[i-1]:
                runs.append(gseq[i])
        print(f"    gseq={gseq}: {count} paths, {switches} switches, group runs={runs}")

    # How many paths stay within groups (block paths)?
    block_paths = sum(count for gseq, count in group_seq_count.items()
                      if all(gseq[i] <= gseq[i+1] or gseq[i] >= gseq[i+1] for i in range(5))
                      and len(set(gseq[i] for i in range(6))) == 3)

    # Actually "block" = groups visited contiguously:
    # group sequence like (0,0,1,1,2,2) or (2,2,0,0,1,1) etc.
    # This means each group appears as a contiguous block
    def is_block(gseq):
        seen = set()
        prev = -1
        for g in gseq:
            if g != prev:
                if g in seen:
                    return False
                seen.add(g)
                prev = g
        return True

    block_count = sum(count for gseq, count in group_seq_count.items() if is_block(gseq))
    non_block_count = len(paths) - block_count

    print(f"\n  Block paths (each group contiguous): {block_count}")
    print(f"  Non-block paths (groups interleaved): {non_block_count}")

    # For block paths: # = H(T1) * prod H(T2_copy) = 3 * 1^3 = 3?
    # Wait: block paths need groups in an ORDER consistent with T1.
    # For the 3-cycle: valid group orderings are (0,1,2), (1,2,0), (2,0,1)
    # For each: 2 internal orderings per group (but only 1 for transitive T2)
    # So block paths = 3 group orderings * 1 * 1 * 1 = 3. Let me check.

    block_seqs = [gseq for gseq, count in group_seq_count.items() if is_block(gseq)]
    print(f"\n  Block group sequences: {block_seqs}")
    for gseq in block_seqs:
        print(f"    {gseq}: {group_seq_count[gseq]} paths")

    # The REMAINING 42 paths are non-block — they interleave groups!
    # These are paths like (0,0, 1,0, 0,1, 1,1, 2,0, 2,1) where
    # the path goes back to group 0 after visiting group 1.

    # This is the KEY: the cyclic structure of T1 allows BACKTRACKING
    # between groups, creating many more valid orderings.

    print(f"\n  *** KEY INSIGHT ***")
    print(f"  The simple formula H(T1)*H(T2)^|V1| counts only BLOCK paths")
    print(f"  (where each group is visited contiguously).")
    print(f"  But the cyclic T1 = 3-cycle allows interleaving: visiting")
    print(f"  part of group A, then group B, then back to group A.")
    print(f"  This creates {non_block_count} extra paths beyond the {block_count} block paths!")
    print(f"  Total: {block_count} + {non_block_count} = {len(paths)} = 45.")

    # ========================================
    # PART 2: Path interleaving patterns
    # ========================================
    print(f"\n{'='*70}")
    print("PART 2: INTERLEAVING PATTERNS")
    print("  Non-block paths interleave groups. How?")
    print(f"{'='*70}")

    # Group the non-block paths by their "interleaving signature"
    # = the sequence of group RUNS
    interleave_patterns = Counter()
    for path in paths:
        gseq = tuple(group[v] for v in path)
        runs = []
        for i in range(6):
            if i == 0 or gseq[i] != gseq[i-1]:
                runs.append(gseq[i])
        interleave_patterns[tuple(runs)] += 1

    print(f"\n  Interleaving patterns (group run sequences):")
    for pattern, count in sorted(interleave_patterns.items(), key=lambda x: (-len(x[0]), x[0])):
        is_blk = len(pattern) == 3  # exactly 3 runs = block
        bstr = "BLOCK" if is_blk else "INTERLEAVED"
        print(f"    {pattern}: {count} paths [{bstr}]")

    # ========================================
    # PART 3: Transfer matrix approach
    # ========================================
    print(f"\n{'='*70}")
    print("PART 3: TRANSFER MATRIX FOR LEX PRODUCT")
    print("  Define states as (group, position_within_group)")
    print("  The transfer matrix M[s1, s2] = 1 iff transition s1->s2 is an arc")
    print("  H = #{Ham paths in M} = same as our tournament")
    print(f"{'='*70}")

    # States: (group_id, vertex_within_group)
    # The tournament on 6 vertices IS the transfer system.
    # H = 45 is already computed.

    # More useful: BLOCK the transfer matrix by groups
    # M_ij = 2x2 block describing arcs from group i to group j
    print(f"\n  Block structure of adjacency matrix:")
    for gi in range(3):
        for gj in range(3):
            block = [[int(A[gi*2+a][gj*2+b]) for b in range(2)] for a in range(2)]
            print(f"    M[group {gi} -> group {gj}] = {block}")

    # Internal blocks: M[0->0] = [[0,1],[0,0]] (transitive)
    # Cross blocks: M[0->1] = [[1,1],[1,1]] (all ones)
    # Cross blocks: M[2->0] = [[1,1],[1,1]] (all ones)
    # Cross blocks: M[1->0] = [[0,0],[0,0]] (all zeros)

    # The "all ones" cross blocks mean: ANY vertex in group i can reach
    # ANY vertex in group j (in the direction T1[i][j]=1).
    # This is why interleaving works: you can leave group i at any point
    # and enter group j at any point.

    print(f"\n  Cross-group blocks are ALL-ONES for winning direction,")
    print(f"  ALL-ZEROS for losing direction. This creates maximum flexibility")
    print(f"  for path construction, explaining the high H value.")

    # ========================================
    # PART 4: What if T2 is NOT transitive?
    # ========================================
    print(f"\n{'='*70}")
    print("PART 4: T_3 lex T_3 — CYCLIC LEX CYCLIC")
    print("  T_3 lex T_3 gives n=9, H=3159")
    print("  The internal blocks are 3-cycles, not transitive")
    print(f"{'='*70}")

    # Build T_3 lex T_3
    cycle3 = np.array([[0,1,0],[0,0,1],[1,0,0]])
    n9 = 9
    A9 = np.zeros((n9, n9), dtype=int)
    for i1 in range(3):
        for j1 in range(3):
            for i2 in range(3):
                for j2 in range(3):
                    v1 = i1 * 3 + j1
                    v2 = i2 * 3 + j2
                    if v1 == v2: continue
                    if cycle3[i1][i2] == 1:
                        A9[v1][v2] = 1
                    elif i1 == i2 and cycle3[j1][j2] == 1:
                        A9[v1][v2] = 1

    scores9 = A9.sum(axis=1)
    print(f"  Scores: {list(scores9.astype(int))}")
    print(f"  All scores = 4? {all(s == 4 for s in scores9)} (REGULAR)")

    # Compute H
    from itertools import permutations as perm9
    # This would be too slow with brute force for n=9 (9! = 362880)
    # Use DP instead
    def compute_H_dp(A, n):
        dp = {}
        for v in range(n): dp[(1 << v, v)] = 1
        for ms in range(2, n+1):
            for mask in range(1 << n):
                if bin(mask).count('1') != ms: continue
                for v in range(n):
                    if not (mask & (1 << v)): continue
                    pm = mask ^ (1 << v)
                    t = sum(dp.get((pm, u), 0) for u in range(n) if (pm & (1 << u)) and A[u][v])
                    if t: dp[(mask, v)] = t
        return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

    H9 = compute_H_dp(A9.tolist(), n9)
    print(f"  H(T_3 lex T_3) = {H9}")
    print(f"  max_H(9) = 3357")
    print(f"  Ratio: {H9/3357:.4f}")
    print(f"  T_3 lex T_3 achieves {100*H9/3357:.1f}% of maximum!")

    # Compare: is T_3 lex T_3 the REGULAR tournament at n=9 closest to max?
    # Regular tournaments at n=9 have scores (4,4,4,4,4,4,4,4,4).
    # There are many such tournaments. T_3 lex T_3 is one.

    print(f"\n  NOTE: T_3 lex T_3 is REGULAR (all scores 4).")
    print(f"  The Paley tournament at p=7 has n=7 vertices (wrong size).")
    print(f"  At n=9, the maximizer likely has a different structure.")

    # ========================================
    # PART 5: General pattern — when does lex product maximize?
    # ========================================
    print(f"\n{'='*70}")
    print("PART 5: WHEN DOES LEX PRODUCT ACHIEVE max_H?")
    print(f"{'='*70}")

    results = {
        3: ("T_3 (trivially)", 3, 3, True),
        4: ("no lex product", None, 5, False),
        6: ("T_3 lex T_2", 45, 45, True),
        9: ("T_3 lex T_3", 3159, 3357, False),
    }

    print(f"  {'n':>3} {'Construction':>20} {'H_lex':>8} {'max_H':>8} {'Achieves max?':>15}")
    for n_val, (desc, h_lex, h_max, achieves) in sorted(results.items()):
        h_str = str(h_lex) if h_lex else 'N/A'
        print(f"  {n_val:3d} {desc:>20} {h_str:>8} {h_max:>8} {'YES' if achieves else 'NO':>15}")

    print(f"\n  CONCLUSION: Lex products give maximizers at n=3,6 (=3*1, 3*2)")
    print(f"  but NOT at n=4,9. The pattern is NOT simply 'blow up T_3'.")
    print(f"  At odd prime n=7,11: Paley tournaments achieve max, not lex products.")
    print(f"  At even n: lex products of T_3 with smaller tournaments can help.")

    print(f"\n{'='*70}")
    print("DONE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
