#!/usr/bin/env python3
"""
EXHAUSTIVE n=7 VERIFICATION of the Independence Polynomial Formula
opus-2026-03-14-S89b

CONJECTURE (THM-209):
  H(T) = IP(OddCycleDisjointnessGraph(T), 2)
       = Σ_{k≥0} 2^k · i_k(T)

where i_k(T) = number of collections of k pairwise vertex-disjoint
directed odd cycles in T.

At n=7, the possible disjoint collections are:
  k=0: empty set (always 1)
  k=1: any single odd cycle (t₃ + t₅ + t₇)
  k=2: any pair of disjoint odd cycles
    - (3,3) pairs need 6 ≤ 7 vertices ✓
    - (3,5) pairs need 8 > 7 vertices ✗
    - others need even more ✗
  k≥3: impossible (need ≥ 9 vertices)

So: H = 1 + 2(t₃ + t₅ + t₇) + 4·d₃₃

Full enumeration: 2^21 = 2,097,152 tournaments. Feasible with optimized code.
"""

from itertools import combinations, permutations
import sys

def compute_H(n, adj_matrix):
    """DP bitmask Hamiltonian path count using adjacency matrix."""
    dp = [0] * (n * (1 << n))
    for v in range(n):
        dp[(1 << v) * n + v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp[mask * n + v]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj_matrix[v][u]:
                    dp[(mask | (1 << u)) * n + u] += val

    full_mask = (1 << n) - 1
    return sum(dp[full_mask * n + v] for v in range(n))

def count_3cycles(n, adj):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    count += 1
                elif adj[i][k] and adj[k][j] and adj[j][i]:
                    count += 1
    return count

def get_3cycle_triples(n, adj):
    """Return list of frozensets of 3-cycle vertex sets."""
    triples = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    triples.append(frozenset([i,j,k]))
    return triples

def count_5cycles(n, adj):
    """Count directed 5-cycles using matrix method."""
    count = 0
    for combo in combinations(range(n), 5):
        a, b, c, d, e = combo
        # Check all 12 possible directed 5-cycles on these vertices
        # A 5-cycle visits all 5 in some order and returns
        for perm in permutations(combo):
            ok = True
            for idx in range(5):
                if not adj[perm[idx]][perm[(idx+1)%5]]:
                    ok = False
                    break
            if ok:
                count += 1
    return count // 5

def count_7cycles(n, adj):
    """Count directed 7-cycles. Only for n=7 (single combo)."""
    if n < 7:
        return 0
    count = 0
    for combo in combinations(range(n), 7):
        for perm in permutations(combo):
            ok = True
            for idx in range(7):
                if not adj[perm[idx]][perm[(idx+1)%7]]:
                    ok = False
                    break
            if ok:
                count += 1
    return count // 7

def count_disjoint_pairs(triples):
    """Count disjoint pairs of 3-cycles."""
    count = 0
    for i in range(len(triples)):
        for j in range(i+1, len(triples)):
            if len(triples[i] & triples[j]) == 0:
                count += 1
    return count

n = 7
m = n*(n-1)//2  # 21
total = 1 << m   # 2097152

print(f"EXHAUSTIVE n={n} VERIFICATION")
print(f"m={m}, total={total} tournaments")
print(f"Testing: H = 1 + 2(t₃+t₅+t₇) + 4·d₃₃")
print()

# For n=7, counting 7-cycles via permutations of 7 elements (5040) for each
# of C(7,7)=1 combo is manageable but slow. Let's optimize.

# Actually for n=7 there's exactly 1 combo for 7-cycles: all 7 vertices.
# We need to count directed 7-cycles = (7-1)!/2 possible directed cycles,
# but only those where all arcs go forward.

# Optimization: precompute adjacency matrix, use arrays

errors = 0
max_error = 0
progress_step = total // 20

# For 7-cycles on all 7 vertices, we can precompute more efficiently
# Number of directed Hamiltonian cycles = number of cyclic orderings
# Fix vertex 0, permute rest, check if 0->p[0]->p[1]->...->p[5]->0

from itertools import permutations as perms

for bits in range(total):
    if bits % progress_step == 0:
        pct = 100 * bits // total
        print(f"  {pct}% done... (errors so far: {errors})", file=sys.stderr, flush=True)

    # Build adjacency matrix
    adj = [[False]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            idx += 1

    # H via DP
    H = compute_H(n, adj)

    # t₃
    t3 = count_3cycles(n, adj)

    # t₅
    t5 = count_5cycles(n, adj)

    # t₇: count Hamiltonian directed cycles
    # Fix vertex 0, permute 1..6
    t7 = 0
    for perm in perms(range(1, n)):
        ok = adj[0][perm[0]]
        if not ok:
            continue
        for idx2 in range(len(perm)-1):
            if not adj[perm[idx2]][perm[idx2+1]]:
                ok = False
                break
        if ok and adj[perm[-1]][0]:
            t7 += 1
    t7 //= 2  # each directed cycle counted twice (cw and ccw)...
    # Wait: in a tournament, a directed Hamiltonian cycle has a UNIQUE direction.
    # Fix v0, then the (n-1)! orderings include each directed cycle n times
    # (once for each starting vertex), but we fixed v0, so each cycle appears once.
    # No wait — fixing v0 means each directed cycle that includes v0 appears
    # exactly once as 0->perm[0]->...->perm[5]->0.
    # But a directed Hamiltonian cycle visits 0 exactly once, so fixing the start
    # at 0 gives each directed cycle exactly once.
    # However, the "reverse" cycle is a DIFFERENT directed cycle in a tournament.
    # So t7 should NOT be divided by 2.
    t7 *= 2  # undo the division
    # Actually let me think again:
    # A directed 7-cycle: v₀→v₁→...→v₆→v₀
    # Fixing v₀=0: we get each directed cycle exactly once
    # The reverse 0→v₆→v₅→...→v₁→0 is a DIFFERENT cycle (different direction)
    # But in a tournament, only one of these can exist (since arc direction is fixed)
    # So count = number of perms of [1..6] such that 0->p1->p2->...->p6->0
    # This directly counts directed Hamiltonian cycles. No division needed.
    # Let me recount properly:

    t7_correct = 0
    for perm in perms(range(1, n)):
        ok = adj[0][perm[0]]
        if not ok:
            continue
        for idx2 in range(len(perm)-1):
            if not adj[perm[idx2]][perm[idx2+1]]:
                ok = False
                break
        if ok and adj[perm[-1]][0]:
            t7_correct += 1
    # This counts each directed Hamiltonian cycle once (start fixed at 0).
    # But wait — the "standard" definition of t_k counts the number of
    # directed k-cycles as SUBGRAPHS, where each cycle on a vertex set
    # is counted by its edge set. A directed k-cycle on {v0,...,v_{k-1}}
    # has k rotations that give the same subgraph, so if we sum over
    # all labeled orderings we get k times the subgraph count.
    # But by fixing v0=0, we avoid the rotation issue — each cycle appears once.
    # So t7_correct IS the number of directed Hamiltonian 7-cycles.
    t7 = t7_correct

    # d₃₃: disjoint 3-cycle pairs
    triples = get_3cycle_triples(n, adj)
    d33 = count_disjoint_pairs(triples)

    predicted = 1 + 2*(t3 + t5 + t7) + 4*d33

    if predicted != H:
        errors += 1
        err = abs(H - predicted)
        if err > max_error:
            max_error = err
        if errors <= 10:
            print(f"  MISMATCH #{errors}: bits={bits}, H={H}, pred={predicted}, "
                  f"t3={t3}, t5={t5}, t7={t7}, d33={d33}, diff={H-predicted}")

print(f"\n{'='*60}")
if errors == 0:
    print(f"*** PERFECT: H = 1 + 2(t₃+t₅+t₇) + 4·d₃₃ EXACT for ALL {total} n=7 tournaments! ***")
    print(f"*** THM-209 VERIFIED through n=7! ***")
else:
    print(f"FAILED: {errors}/{total} mismatches, max error = {max_error}")
print(f"{'='*60}")
