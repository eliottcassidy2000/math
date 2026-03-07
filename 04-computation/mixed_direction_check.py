#!/usr/bin/env python3
"""
Check whether mixed-direction cycle collections exist.

Finding from ocf_factor2_investigation.py: T-only = T^op-only = mixed.
This could mean:
(a) No mixed-direction collections exist, OR
(b) Mixed collections exist but cancel out in the weighted sum.

Let me check (a) directly: enumerate all valid permutations and flag
any that have BOTH a T-cycle and a T^op-cycle.

opus-2026-03-07-S37
"""
from itertools import permutations
from collections import defaultdict

def check_mixed_collections(n, edges):
    """Find permutations with both T-cycles and T^op-cycles."""
    edge_set = set(edges)
    opp_edges = set((j, i) for i in range(n) for j in range(n)
                    if i != j and (i, j) not in edge_set)

    mixed_count = 0
    mixed_examples = []

    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for start in range(n):
            if visited[start]:
                continue
            cycle = []
            curr = start
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = sigma[curr]
            cycles.append(tuple(cycle))

        # Check: all nontrivial cycles odd, and each is T or T^op
        has_T = False
        has_Top = False
        ok = True
        all_odd = True

        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                all_odd = False
                ok = False
                break

            is_T = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set
                       for i in range(len(cyc)))
            is_Top = all((cyc[i], cyc[(i+1) % len(cyc)]) in opp_edges
                        for i in range(len(cyc)))

            if not is_T and not is_Top:
                ok = False
                break

            if is_T and not is_Top:
                has_T = True
            elif is_Top and not is_T:
                has_Top = True
            else:
                # Both T and T^op — this means it's a cycle in both directions
                # Possible for self-converse cycles
                pass

        if ok and all_odd and has_T and has_Top:
            mixed_count += 1
            if len(mixed_examples) < 3:
                nontrivial = [(cyc, 'T' if all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set for i in range(len(cyc))) else 'T^op')
                              for cyc in cycles if len(cyc) > 1]
                mixed_examples.append(nontrivial)

    return mixed_count, mixed_examples

# Test C_5
n = 5
cyc_edges = [(i, (i+1) % n) for i in range(n)] + [(i, (i+2) % n) for i in range(n)]
mc, examples = check_mixed_collections(n, cyc_edges)
print(f"C_5: {mc} mixed-direction permutations")
for ex in examples:
    print(f"  Example: {ex}")

# Test Paley T_7
n = 7
QR = {1, 2, 4}
paley_edges = [(i, j) for i in range(n) for j in range(n)
               if i != j and (j - i) % n in QR]
mc, examples = check_mixed_collections(n, paley_edges)
print(f"\nPaley T_7: {mc} mixed-direction permutations")
for ex in examples:
    print(f"  Example: {ex}")

# Try a random tournament at n=7
import random
random.seed(42)
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1
rand_edges = [(i, j) for i in range(n) for j in range(n) if i != j and A[i][j]]
mc, examples = check_mixed_collections(n, rand_edges)
print(f"\nRandom T_7: {mc} mixed-direction permutations")
for ex in examples:
    print(f"  Example: {ex}")

# Now check: are there SELF-CONVERSE odd cycles (same cycle is both T and T^op)?
print(f"\n=== Self-converse odd cycles ===")
edge_set = set(paley_edges)
opp_set = set((j, i) for (i, j) in paley_edges)
# For Paley: T is self-converse, so T = T^op up to labeling
# Check: is Paley exactly self-converse?
print(f"Paley self-converse: {edge_set == opp_set}")
# No, T != T^op in general. But Paley has an anti-automorphism.

# Actually for Paley T_7: the map i -> -i mod 7 is an anti-automorphism
# (it reverses all edges). So T^op is isomorphic to T, but not equal.

# Can a directed 3-cycle (a->b->c->a) also be (a->c->b->a)?
# That would mean: a->b, b->c, c->a AND a->c, c->b, b->a.
# But a->b and b->a can't both hold in a tournament. So NO.
# Self-converse cycles DON'T EXIST.

print(f"\nConclusion: a directed cycle C cannot be both a T-cycle and a T^op-cycle")
print(f"(since that would require both a->b and b->a for some edge).")
print(f"Therefore: 'mixed' means genuinely having some T-cycles and some T^op-cycles.")
print(f"And the finding that mixed = T-only = T^op-only means:")
print(f"  EITHER no mixed collections exist (impossible for some tournaments)")
print(f"  OR mixed collections exist but their contributions exactly cancel.")
