#!/usr/bin/env python3
"""
Connect the quadratic H formula H = disj^2 - 23*disj + 301
for regular n=7 tournaments to the OCF decomposition.

From kind-pasteur-S61: For ALL 2640 regular n=7 tournaments,
  H = disj_33^2 - 23*disj_33 + 301
where disj_33 = number of vertex-disjoint 3-cycle pairs.

From OCF: H = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
where alpha_k = # independent k-sets in the conflict graph Omega(T).

QUESTION: Can we derive the quadratic formula from OCF?

For regular n=7: every vertex has out-degree 3, in-degree 3.
Number of 3-cycles: c3 = 14 for ALL regular n=7 (why? because
  c3 = n(n-1)(n-2)/6 - n*(s_i choose 2) adjustment... actually,
  c3 is determined by the score sequence for regular tournaments).

So alpha_1 = c3 = 14 always. H = 1 + 28 + 4*alpha_2 + 8*alpha_3 + ...

disj_33 counts pairs of 3-cycles sharing NO vertices.
In the conflict graph Omega, two 3-cycles are ADJACENT iff they share a vertex.
So disj_33 = # NON-EDGES among 3-cycle vertices in Omega.

Wait — disj_33 = # pairs of 3-cycles that DON'T share a vertex.
Total 3-cycle pairs: C(14,2) = 91.
disj_33 = 91 - (# pairs sharing a vertex) = 91 - |E(Omega restricted to 3-cycles)|.

But alpha_2 counts independent SETS in the FULL conflict graph (including 5-cycles, 7-cycles).
An independent pair means two odd cycles sharing no vertex.

Key: at n=7, the odd cycles are 3-cycles, 5-cycles, and 7-cycles.
An independent pair in Omega can be:
  (3,3): two disjoint 3-cycles = disj_33
  (3,5): impossible! 3+5=8 > 7 vertices total, so they MUST share at least 1 vertex
  (3,7): impossible! 3+7=10 > 7
  (5,5): impossible! 5+5=10 > 7
  (5,7): impossible! 5+7=12 > 7
  (7,7): impossible! only 7 vertices, can't have two disjoint 7-cycles

So alpha_2 = disj_33! The only independent pairs are disjoint 3-cycle pairs!

Similarly, alpha_3 = # independent triples. Three pairwise disjoint 3-cycles
need 9 vertices, but n=7. So alpha_3 = 0.
And alpha_k = 0 for k >= 3.

Therefore: H = 1 + 2*14 + 4*disj_33 = 29 + 4*disj_33

But wait, the actual values are H=189 (disj=7), H=171 (disj=10), H=175 (disj=14).
  29 + 4*7 = 57 ≠ 189.

Something is wrong. Let me reconsider.

opus-2026-03-13-S71b
"""

import itertools
import numpy as np
from collections import Counter

def all_tournaments(n):
    """Generate all tournaments on n vertices."""
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for bits in range(2**m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

def is_regular(adj, n):
    """Check if tournament is regular (all out-degrees equal)."""
    for i in range(n):
        if sum(adj[i]) != (n-1)//2:
            return False
    return True

def count_hp(adj, n):
    """Count Hamiltonian paths by DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and adj[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def find_odd_cycles(adj, n):
    """Find all directed odd cycles."""
    cycles = set()
    for length in range(3, n+1, 2):
        for combo in itertools.combinations(range(n), length):
            for perm in itertools.permutations(combo):
                # Check if this is a directed cycle
                valid = True
                for i in range(length):
                    u, v = perm[i], perm[(i+1) % length]
                    if not adj[u][v]:
                        valid = False
                        break
                if valid:
                    # Normalize: start from smallest vertex, choose direction
                    min_idx = perm.index(min(perm))
                    canonical = perm[min_idx:] + perm[:min_idx]
                    cycles.add(canonical)
    return list(cycles)

def build_conflict_graph(cycles, n):
    """Build conflict graph: edges between cycles sharing a vertex."""
    k = len(cycles)
    vertex_sets = [set(c) for c in cycles]
    adj_conflict = [[False]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1, k):
            if vertex_sets[i] & vertex_sets[j]:
                adj_conflict[i][j] = True
                adj_conflict[j][i] = True
    return adj_conflict

def count_independent_sets(adj, k):
    """Count independent sets by size."""
    n = len(adj)
    counts = Counter()
    counts[0] = 1  # empty set
    for size in range(1, n+1):
        count = 0
        for combo in itertools.combinations(range(n), size):
            indep = True
            for i in range(len(combo)):
                for j in range(i+1, len(combo)):
                    if adj[combo[i]][combo[j]]:
                        indep = False
                        break
                if not indep:
                    break
            if indep:
                count += 1
        counts[size] = count
        if count == 0:
            break
    return counts

def independence_poly_at_2(counts):
    """Compute I(G, 2) = sum alpha_k * 2^k."""
    return sum(counts[k] * (2**k) for k in counts)

# ============================================================
n = 7
print(f"Analyzing regular tournaments at n={n}")
print(f"OCF: H = I(Omega, 2) = sum alpha_k * 2^k\n")

# Sample regular tournaments to analyze
# We'll gather the three H-classes
results_by_H = {}
count = 0

for adj in all_tournaments(n):
    if not is_regular(adj, n):
        continue
    count += 1
    if count > 100:  # Sample first 100
        break

    H = count_hp(adj, n)
    if H not in results_by_H:
        cycles = find_odd_cycles(adj, n)

        # Classify cycles by length
        by_length = Counter()
        for c in cycles:
            by_length[len(c)] += 1

        # Count disjoint 3-cycle pairs
        three_cycles = [c for c in cycles if len(c) == 3]
        n3 = len(three_cycles)
        disj_33 = 0
        for i in range(n3):
            for j in range(i+1, n3):
                if not (set(three_cycles[i]) & set(three_cycles[j])):
                    disj_33 += 1

        # Build conflict graph and compute independence polynomial
        conflict_adj = build_conflict_graph(cycles, n)
        alpha = count_independent_sets(conflict_adj, len(cycles))
        I_at_2 = independence_poly_at_2(alpha)

        results_by_H[H] = {
            'cycles_by_length': dict(by_length),
            'disj_33': disj_33,
            'alpha': dict(alpha),
            'I_at_2': I_at_2,
            'n_cycles': len(cycles)
        }

        print(f"H={H}:")
        print(f"  Cycles: {dict(by_length)}")
        print(f"  Total cycles: {len(cycles)}")
        print(f"  disj_33 = {disj_33}")
        print(f"  alpha (indep set counts): {dict(alpha)}")
        print(f"  I(Omega, 2) = {I_at_2}")
        print(f"  H == I(Omega,2): {H == I_at_2}")

        # Decompose H via OCF
        print(f"  H = 1", end="")
        for k in range(1, max(alpha.keys())+1):
            if alpha.get(k, 0) > 0:
                print(f" + {2**k}*{alpha[k]}", end="")
        print(f" = {I_at_2}")
        print()

# Now analyze the structure
print("\n" + "="*60)
print("ANALYSIS: Why H = disj^2 - 23*disj + 301?")
print("="*60)

for H, data in sorted(results_by_H.items()):
    alpha = data['alpha']
    disj = data['disj_33']

    # From OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
    # alpha_1 = c3 + c5 + c7 (all odd cycles)
    # But alpha_2 includes pairs of disjoint odd cycles of any length

    # Key question: what is alpha_2?
    # At n=7: disjoint cycle pairs need vertex-disjoint odd cycles
    # Possible: (3,3) = disj_33
    # (3,5): needs 8 vertices, but n=7. IMPOSSIBLE? Wait...
    # Actually at n=7, a 3-cycle uses 3 vertices and a 5-cycle uses 5 vertices.
    # If disjoint, need 8 > 7. So (3,5) disjoint is IMPOSSIBLE. ✓
    # Similarly all other pairs impossible.
    # So alpha_2 = disj_33. ✓

    # But then H = 1 + 2*alpha_1 + 4*disj_33
    # This is LINEAR in disj, not quadratic!
    # Unless alpha_1 also depends on disj...

    print(f"H={H}: alpha_1={alpha.get(1,0)}, alpha_2={alpha.get(2,0)}={disj}, "
          f"alpha_3={alpha.get(3,0)}")

    # Check: alpha_2 == disj_33?
    print(f"  alpha_2 == disj_33: {alpha.get(2,0) == disj}")

    # Check: does alpha_1 vary?
    c3 = data['cycles_by_length'].get(3, 0)
    c5 = data['cycles_by_length'].get(5, 0)
    c7 = data['cycles_by_length'].get(7, 0)
    print(f"  c3={c3}, c5={c5}, c7={c7}, alpha_1={c3+c5+c7}")

    # H = 1 + 2*(c3+c5+c7) + 4*alpha_2 + 8*alpha_3 + ...
    h_check = 1 + 2*(c3+c5+c7) + 4*alpha.get(2,0) + 8*alpha.get(3,0) + 16*alpha.get(4,0)
    print(f"  H from OCF = {h_check}")

    # The quadratic formula: H = disj^2 - 23*disj + 301
    h_quad = disj**2 - 23*disj + 301
    print(f"  H from quadratic = {h_quad}")
    print()

print("\n" + "="*60)
print("KEY INSIGHT: alpha_1 varies! c5 depends on the 3-cycle overlap structure.")
print("="*60)
print("""
For regular n=7:
  c3 = 14 (constant)
  c5 = varies with overlap structure
  c7 = varies
  alpha_1 = c3 + c5 + c7 = varies

The OCF says H = 1 + 2*alpha_1 + 4*alpha_2 + higher terms.
Since alpha_2 = disj_33 (only disjoint 3-cycle pairs contribute),
and alpha_k = 0 for k >= 3 at n=7,
we need c5 + c7 as a function of disj_33.

If H = disj^2 - 23*disj + 301 AND H = 1 + 2*(14+c5+c7) + 4*disj_33:
  disj^2 - 23*disj + 301 = 29 + 2*c5 + 2*c7 + 4*disj
  c5 + c7 = (disj^2 - 27*disj + 272) / 2

This gives a FORMULA for c5+c7 in terms of disj_33!

CHECK: disj=7:  c5+c7 = (49-189+272)/2 = 132/2 = 66
       disj=10: c5+c7 = (100-270+272)/2 = 102/2 = 51
       disj=14: c5+c7 = (196-378+272)/2 = 90/2 = 45
""")
