#!/usr/bin/env python3
"""
test_perpath_n7.py -- Test OPEN-Q-010: per-path formula at n=7 with 3+5 cycles.

CONTEXT (from THM-008):
  At n=7, for any 5-cycle C = (v, a1, a2, a3, a4) through v:
    V \ {v, a1, a2, a3, a4} has exactly 2 vertices.
    No odd cycle fits on 2 vertices.
    => mu(C) = 1 for ALL 5-cycles at n=7.

So at n=7, only 3-cycles can have mu > 1. But the per-path identity only sums
over 3-cycles (by THM-005), and at n=7, 5-cycles also contribute to Claim A's
RHS with mu=1. The per-path identity counts 3-cycles only.

The question (OPEN-Q-010): if we extend the per-path sum to include both
  - 3-cycle embeddings (each with their mu weight)
  - 5-cycle embeddings (each with mu=1, since they're trivial at n=7)

does the resulting formula equal (inshat-1)/2?

Per-path formula: for P' a Ham path of T-v,
  (inshat(v, P') - 1) / 2 = sum over TypeII positions j: mu(v, P'[j], P'[j+1])
  (this is the existing per-path identity using only 3-cycles at TypeII positions)

Extended formula attempt: sum over ALL consecutive pairs (a,b) in P' such that
  (v,a,b) forms a directed 3-cycle OR there's a 5-cycle embedding starting with v
  at a "relevant" position...

Wait, the per-path identity (THM-005) is:
  (inshat-1)/2 = #{TypeII positions} = #{3-cycle embeddings where (a,b) consecutive in P'}

At n=7, the 5-cycles also have mu=1. Claim A says:
  H(T) - H(T-v) = 2 * sum_C mu(C)  [over all odd cycles C through v]
              = 2 * (sum_{3-cycles} 1 + sum_{5-cycles} 1)  [at n=7, all mu=1]

But per-path identity gives per path P':
  (inshat-1)/2 = #{3-cycles (v,a,b) with (a,b) consecutive in P'}

Summing over all P' in Ham(T-v):
  H(T) - H(T-v) = sum_{P'} (inshat-1)/2 is NOT what Claim A says.

Actually, let me re-read. Claim A says:
  H(T) - H(T-v) = 2 * Sigma_C mu(C)

This is a SUM over cycles, not a per-path formula. The per-path identity is a
different statement. Let me think about what OPEN-Q-010 is actually asking.

From opus session log: "A per-path formula summing over both 3-cycle and 5-cycle
embeddings (each with their mu weights) might work at n=7."

So the question is: at n=7, is there a per-path formula
  (inshat(v,P') - 1) / 2 = f(v, P', T)
where f sums over both 3-cycle and 5-cycle embeddings?

For 5-cycles: the Type-II count for 5-cycle embeddings would need to be defined
somehow. A 5-cycle (v, a, b, c, d) is "embedded" at position j in P' if what?
Perhaps: there's a window of 4 consecutive vertices in P' that together with v
form a directed 5-cycle. The Type-II contribution from a 5-cycle window...

Actually, let's think from first principles. The per-path identity:
  (inshat-1)/2 = #{TypeII positions}

TypeII position j: s[j]=1, s[j+1]=0 (v->P'[j], P'[j+1]->v).
Each TypeII position (j) at a pair (a=P'[j], b=P'[j+1]) with 3-cycle (v,a,b)
contributes 1 to the RHS.

For a 5-cycle window of length 4 at positions (j, j+1, j+2, j+3), the
contribution from THM-007 is sum over k TypeII positions in window = C(2, 2k-1)
patterns. The mean TypeII count per 5-cycle window is (5-4)/4 = 0.25.

But we want a formula where each 5-cycle embedding contributes a fixed amount
(perhaps its mu=1), matching what Claim A needs.

Let me just test empirically: compute for each (T, v, P') at n=7:
  LHS = (inshat(v,P') - 1) / 2
  RHS_3 = #{3-cycle (v,a,b) with (a,b) consecutive in P'}  [existing formula]
  RHS_35 = #{3-cycle embeddings} + #{5-cycle embeddings at P'}

where "5-cycle embedding" means: there is a directed 5-cycle (v, P'[j], P'[j+1],
P'[j+2], P'[j+3]) for some j. (Or all directed 5-cycles through v that can be
embedded contiguously in P'.)

Actually let me just compute LHS - RHS_3 (the "gap" from 3-cycles only) and see
if this gap equals the 5-cycle count, or some other function of 5-cycles.
"""

import sys
sys.path.insert(0, '.')
from tournament_lib import (
    random_tournament, all_tournaments, delete_vertex,
    find_odd_cycles, hamiltonian_path_count, mu
)
from itertools import permutations
import random


def inshat(T, v, path):
    """Compute inshat(v, P') = signed count of insertions of v into path.
    path is a list of vertices (a Ham path of T-v).
    Returns (inshat_value, signature, typeii_positions).

    Per definitions.md:
      s[j] = 1 if v->P'[j], 0 if P'[j]->v
      TypeII at j: s[j]=1 and s[j+1]=0
      b = s[0] + (1 - s[n-2])
      inshat = b + #{TypeII} + #{TypeI}
        where TypeI at j: s[j]=0 and s[j+1]=1

    Actually inshat = b + #{valid insertion positions}
    where a position k is valid if P'[k-1]->v->P'[k] is consistent.
    A position between P'[j] and P'[j+1] is valid if:
      either (boundary: before first or after last) ...
    Let me use the formula: inshat = sum over positions, each counted +1 or 0.

    From the definitions: inshat = b + #{TypeI} + #{TypeII}
    Note TypeI: s[j]=0, s[j+1]=1 means P'[j]->v and v->P'[j+1]
    TypeII: s[j]=1, s[j+1]=0 means v->P'[j] and P'[j+1]->v

    But the actual count of VALID insertions is different. Let me use the
    fact that (inshat-1)/2 = #{TypeII} by THM-004.
    So I only need to count TypeII positions.
    """
    n = len(T)
    # Build signature
    sig = [1 if T[v][path[j]] else 0 for j in range(len(path))]
    # TypeII positions: sig[j]=1, sig[j+1]=0
    typeii = [(j, path[j], path[j+1])
              for j in range(len(path)-1)
              if sig[j] == 1 and sig[j+1] == 0]
    # TypeI: sig[j]=0, sig[j+1]=1
    typei = sum(1 for j in range(len(path)-1) if sig[j]==0 and sig[j+1]==1)
    # Boundary
    b = sig[0] + (1 - sig[-1])
    inshat_val = b + typei + len(typeii)
    return inshat_val, sig, typeii


def ham_paths_tv(T, v):
    """All Hamiltonian paths of T-v."""
    Tv, labels = delete_vertex(T, v)
    n = len(T)
    others = [u for u in range(n) if u != v]
    paths = []
    for perm in permutations(others):
        # Check if it's a valid directed Ham path in T
        if all(T[perm[i]][perm[i+1]] for i in range(len(perm)-1)):
            paths.append(list(perm))
    return paths


def find_cycles_through_v_of_length(T, v, length):
    """Find all directed cycles of given length through v in T."""
    n = len(T)
    others = [u for u in range(n) if u != v]
    cycles = []
    # Try all ordered sequences of (length-1) other vertices
    for perm in permutations(others, length - 1):
        seq = [v] + list(perm)
        # Check cycle: seq[0]->seq[1]->...->seq[-1]->seq[0]
        if all(T[seq[i]][seq[(i+1) % length]] for i in range(length)):
            cycles.append(tuple(seq))
    return cycles


def test_perpath_n7(n_samples=200, seed=42):
    """Test the per-path formula at n=7 for random tournaments."""
    rng = random.Random(seed)
    n = 7

    total = 0
    matches_3only = 0
    gap_stats = {}  # gap -> count

    print(f"Testing per-path formula at n={n} ({n_samples} random tournaments)")
    print(f"For each (T, v, P'):")
    print(f"  LHS = (inshat-1)/2")
    print("  RHS_3 = #3-cycle embeddings at TypeII positions")
    print(f"  gap = LHS - RHS_3")
    print()

    for _ in range(n_samples):
        T = random_tournament(n, rng)

        for v in range(n):
            paths = ham_paths_tv(T, v)
            # Find all 5-cycles through v
            cycles5 = find_cycles_through_v_of_length(T, v, 5)
            # Set of 5-cycle vertex sets through v (as frozensets)
            cycle5_sets = [frozenset(c) for c in cycles5]

            for P in paths:
                val, sig, typeii_pos = inshat(T, v, P)
                lhs = (val - 1) // 2

                # RHS_3: count TypeII positions that have 3-cycle (v, P[j], P[j+1])
                rhs3 = sum(1 for (j, a, b) in typeii_pos if T[v][a] and T[a][b] and T[b][v])
                # (TypeII already means v->a and b->v, so we need a->b for a 3-cycle v->a->b->v)
                # Actually: 3-cycle (v, a, b) means v->a, a->b, b->v. TypeII pos j means v->P[j] and P[j+1]->v.
                # So a=P[j], b=P[j+1]. We need a->b = T[P[j]][P[j+1]] = T[a][b]. This is already given
                # since P is a valid Ham path in T-v, so T[P[j]][P[j+1]] must be 1. So ALL TypeII positions
                # have T[a][b] = 1 (as P is a valid directed path).
                # Wait... a Ham path means each consecutive pair has an arc. So T[P[j]][P[j+1]] = 1 always.
                # That means every TypeII position is automatically a 3-cycle (v, a, b)! So rhs3 = len(typeii_pos).

                rhs3 = len(typeii_pos)
                gap = lhs - rhs3

                total += 1
                if gap == 0:
                    matches_3only += 1
                gap_stats[gap] = gap_stats.get(gap, 0) + 1

    print(f"Total (T, v, P') triples: {total}")
    print(f"Gap = 0 (3-cycle formula works): {matches_3only} ({100*matches_3only/total:.1f}%)")
    print(f"Gap distribution: {sorted(gap_stats.items())}")
    print()

    if matches_3only < total:
        print("Gap != 0 means 3-cycle formula FAILS for some paths.")
        print("Now investigating: does the gap equal the 5-cycle contribution?")
        print()
        investigate_gap_vs_5cycles(n_samples // 2, seed + 1)


def investigate_gap_vs_5cycles(n_samples, seed):
    """For paths where gap != 0, check if gap = number of 5-cycle embeddings."""
    rng = random.Random(seed)
    n = 7

    total_gap_nonzero = 0
    gap_eq_5cycle_count = 0
    discrepancies = []

    for _ in range(n_samples):
        T = random_tournament(n, rng)

        for v in range(n):
            paths = ham_paths_tv(T, v)
            cycles5 = find_cycles_through_v_of_length(T, v, 5)
            # 5-cycle vertices (as tuples for fast lookup)
            cycle5_vsets = set(frozenset(c) for c in cycles5)

            for P in paths:
                val, sig, typeii_pos = inshat(T, v, P)
                lhs = (val - 1) // 2
                rhs3 = len(typeii_pos)
                gap = lhs - rhs3

                if gap == 0:
                    continue

                total_gap_nonzero += 1

                # Count 5-cycle embeddings: windows of 4 consecutive vertices in P
                # that form a 5-cycle (v, P[j], P[j+1], P[j+2], P[j+3]) for some j
                count5 = 0
                for j in range(len(P) - 3):
                    a, b, c, d = P[j], P[j+1], P[j+2], P[j+3]
                    if T[v][a] and T[a][b] and T[b][c] and T[c][d] and T[d][v]:
                        count5 += 1

                if gap == count5:
                    gap_eq_5cycle_count += 1
                else:
                    discrepancies.append((gap, count5))

    if total_gap_nonzero == 0:
        print("No gaps found (3-cycle formula is exact at n=7).")
        return

    print(f"Triples with gap != 0: {total_gap_nonzero}")
    pct = 100 * gap_eq_5cycle_count / total_gap_nonzero
    print(f"Gap == #5-cycle-embeddings: {gap_eq_5cycle_count} ({pct:.1f}%)")

    if discrepancies:
        from collections import Counter
        disc_counts = Counter(discrepancies)
        print(f"Discrepancies (gap, count5): {dict(list(disc_counts.most_common(10)))}")
    else:
        print("PERFECT MATCH: gap == #contiguous-5-cycle-embeddings in ALL failing cases!")
        print()
        print("CONCLUSION: At n=7, the extended per-path formula works:")
        print("  (inshat-1)/2 = #{TypeII positions} + #{contiguous 5-cycle windows in P'}")
        print("  i.e., adding 5-cycle embeddings exactly closes the gap!")


if __name__ == "__main__":
    test_perpath_n7(n_samples=300)
