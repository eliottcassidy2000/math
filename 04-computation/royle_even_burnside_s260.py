#!/usr/bin/env python3
"""
royle_even_burnside_s260.py — Verify the Royle equinumerosity
opus-2026-03-23-S260

The Royle theorem (arXiv:2204.01947): #even_graphs = #tournaments = A000568.

An "even graph" (Royle sense) has no automorphism reversing an odd number of
oriented edges. The Burnside formula for even graphs uses EVEN-ORDER permutations:
    #EvenGraphs = (1/n!) Σ_{|g| even} 2^{c(g_E)}
while tournaments use ODD-ORDER permutations:
    #Tournaments = (1/n!) Σ_{|g| odd} 2^{c(g_E)}
And #Graphs = #Tournaments + #EvenGraphs.

Wait — that's the definition. Let me verify:
- #Graphs = Σ_all 2^{c(g_E)} / n!
- #Tournaments = Σ_{|g| odd} 2^{c(g_E)} / n!  (odd-order g only fix tournaments)
- #RoyleEven = #Graphs - #Tournaments = Σ_{|g| even} 2^{c(g_E)} / n!

So Royle's theorem says: the number of "even graphs" (their definition) equals
the tournament count. The partition is:
    #Graphs = 2 * #Tournaments  ← ONLY IF THE COUNTS ARE EXACTLY EQUAL

Let me verify numerically.
"""

from math import comb, factorial, gcd, lcm
from collections import Counter

def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0: yield []; return
    for f in range(min(n, mx), 0, -1):
        for r in partitions(n - f, f): yield [f] + r

def conjugacy_class_size(n, ct):
    c = Counter(ct); r = factorial(n)
    for l, k in c.items(): r //= (l ** k) * factorial(k)
    return r

def edge_orbits(n, ct):
    """Count edge orbits under a permutation with cycle type ct."""
    # Standard formula: c(g_E) = sum_{i<=j} (if i<j: gcd(c_i,c_j), if i=j: floor(c_i/2))
    # where the sum is over pairs of cycles
    total = 0
    for i in range(len(ct)):
        total += ct[i] // 2  # self-pairing within cycle
        for j in range(i+1, len(ct)):
            total += gcd(ct[i], ct[j])  # pairing between cycles
    return total

def order_of_permutation(ct):
    """Order of a permutation with cycle type ct = lcm of cycle lengths."""
    result = 1
    for c in ct:
        result = lcm(result, c)
    return result

def compute_all(n):
    """Compute #Graphs, #Tournaments, #RoyleEven."""
    nf = factorial(n)
    graphs_sum = 0
    tourn_sum = 0
    royle_even_sum = 0

    for ct in partitions(n):
        ct_list = list(ct)
        ccs = conjugacy_class_size(n, ct_list)
        eo = edge_orbits(n, ct_list)
        fix = 2 ** eo
        ord_g = order_of_permutation(ct_list)

        graphs_sum += ccs * fix

        if ord_g % 2 == 1:  # odd order
            tourn_sum += ccs * fix
        else:  # even order
            royle_even_sum += ccs * fix

    return graphs_sum // nf, tourn_sum // nf, royle_even_sum // nf

if __name__ == '__main__':
    print("=" * 70)
    print("  ROYLE EQUINUMEROSITY VERIFICATION")
    print("  arXiv:2204.01947")
    print("  opus-2026-03-23-S260")
    print("=" * 70)

    print(f"\n{'n':>3} {'#Graphs':>12} {'#Tourn':>12} {'#RoyleEven':>12} "
          f"{'T==RE?':>7} {'G==T+RE?':>9} {'Ratio':>8}")

    for n in range(3, 20):
        G, T, RE = compute_all(n)
        match_te = "✓" if T == RE else "✗"
        match_sum = "✓" if G == T + RE else "✗"
        ratio = T / RE if RE > 0 else float('inf')
        print(f"{n:>3} {G:>12} {T:>12} {RE:>12} {match_te:>7} {match_sum:>9} {ratio:>8.4f}")

    print(f"\n{'='*70}")
    print("  CONCLUSION:")
    print("  If T == RE for all n: Royle equinumerosity holds.")
    print("  If T == RE: then #Graphs = 2 * #Tournaments (exact doubling!)")
    print(f"{'='*70}")

    # Also show the three key sequences
    print(f"\n  A000088 (graphs):     ", end="")
    for n in range(3, 16):
        G, _, _ = compute_all(n)
        print(f"{G}, ", end="")
    print()

    print(f"  A000568 (tournaments): ", end="")
    for n in range(3, 16):
        _, T, _ = compute_all(n)
        print(f"{T}, ", end="")
    print()

    print(f"  Royle-even graphs:     ", end="")
    for n in range(3, 16):
        _, _, RE = compute_all(n)
        print(f"{RE}, ", end="")
    print()

    # The degree-even (Euler graph) counts for comparison
    print(f"\n  A002854 (Euler graphs): checked separately")

    print("\nDONE.")
