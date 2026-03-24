#!/usr/bin/env python3
"""
trivial_vertex_s20gf.py — The trivial vertex bridge: Euler(n) ~ Tournament(n-1)
kind-pasteur-2026-03-25-S20gf

THE DIMENSIONAL CHAIN:
  Graph(n):         C(n,2) bits   [full edge set]
  Tournament(n):    C(n,2) arcs BUT C(n-1,2) free bits [tiling]
  Euler graph(n):   C(n-1,2) free bits [cycle space]
  Tournament(n-1):  C(n-1,2) arcs

All four have C(n-1,2) as their "effective dimension."
The trivial vertex absorbs the difference C(n,2) - C(n-1,2) = n-1.

HYPOTHESIS: Euler graphs on n vertices are equinumerous with
tournaments on n-1 vertices: A002854(n) = A000568(n-1)?

Let me check: A002854 = 1, 1, 2, 3, 7, 16, 54, 243, 2038
              A000568 = 1, 1, 1, 2, 4, 12, 56, 456, 6880

A002854(n) vs A000568(n-1):
  n=1: 1 vs A000568(0)=1?
  n=2: 1 vs 1 ✓
  n=3: 2 vs 1 ✗

Not a simple shift. But the DIMENSION match C(n-1,2) = C(n-1,2) is exact.

The Royle-Praeger count A000568(n) is BETWEEN A002854(n) and A000088(n).
It's the "tournament-sized" subset of graphs, bigger than Euler but smaller than all.

WHAT MAKES IT WORK: The sign representation.
  Euler: all degrees even (topological constraint)
  RP-even: sign(automorphism) = +1 (algebraic constraint)
  Tournament: complete + antisymmetric (structural constraint)

The user's insight: the tiling model bridges all three by making the
"trivial vertex" explicit. The base path pins one vertex, reducing
the tournament from n to effectively n-1 parameters. The Euler graph
also has n-1 free parameters (cycle space). The RP-even graph count
A000568(n) comes from a DIFFERENT identification of these n-1 parameters
that respects the sign representation.

TEST: Compute how the cycle space basis vectors relate to the
sign representation of Aut(G) for each even graph.
"""

import sys
from math import factorial, comb
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE TRIVIAL VERTEX BRIDGE")
print("  kind-pasteur-2026-03-25-S20gf")
print("=" * 80)

# First: verify the dimension match
print("\n  DIMENSION MATCH:")
print(f"  {'n':>3} {'C(n,2)':>8} {'C(n-1,2)':>10} {'n-1':>5} {'Euler(n)':>9} {'Tour(n)':>8} {'Tour(n-1)':>10} {'RP-even(n)':>11}")
A000568 = [1, 1, 1, 2, 4, 12, 56, 456, 6880]
A002854 = [1, 1, 1, 2, 3, 7, 16, 54, 243, 2038]
for n in range(1, 9):
    cn2 = comb(n, 2)
    cn12 = comb(n-1, 2)
    euler = A002854[n] if n < len(A002854) else '?'
    tour_n = A000568[n] if n < len(A000568) else '?'
    tour_n1 = A000568[n-1] if n-1 < len(A000568) else '?'
    print(f"  {n:3d} {cn2:8d} {cn12:10d} {n-1:5d} {euler:>9} {tour_n:>8} {tour_n1:>10} {tour_n:>11}")

# The key: A000568(n) = A002854(n) only at n=1,2.
# For n>=3: A000568(n) > A002854(n). Tournaments outnumber Euler graphs.
# But A000568(n) = |RP-even graphs| (different from Euler).

print(f"\n  RATIO COMPARISON:")
print(f"  {'n':>3} {'A000568/A002854':>16} {'A000568/A000568(n-1)':>22}")
for n in range(2, 9):
    if n < len(A000568) and n < len(A002854):
        r1 = A000568[n] / A002854[n]
        r2 = A000568[n] / A000568[n-1] if A000568[n-1] > 0 else 0
        print(f"  {n:3d} {r1:16.4f} {r2:22.4f}")

# Now: compute the sign representation for each graph automorphism
print(f"\n{'='*60}")
print(f"  SIGN REPRESENTATION AT n=4,5")
print(f"{'='*60}")

for n in [4, 5]:
    N = n
    edges = [(i,j) for i in range(N) for j in range(i+1,N)]
    m = len(edges)
    all_perms = list(permutations(range(N)))
    edge_idx = {e: k for k, e in enumerate(edges)}

    def graph_canon(bits):
        best = None
        for p in all_perms:
            nb = 0
            for k, (i,j) in enumerate(edges):
                pi, pj = min(p[i],p[j]), max(p[i],p[j])
                nk = edge_idx[(pi, pj)]
                if bits & (1 << k): nb |= (1 << nk)
            if best is None or nb < best: best = nb
        return best

    def graph_sign(bits, sigma):
        """Compute sgn_G(sigma) = product over edges {u,v} of (sigma(u)-sigma(v))/(u-v)."""
        sign = 1
        for k, (u, v) in enumerate(edges):
            if not (bits & (1 << k)): continue  # edge not present
            su, sv = sigma[u], sigma[v]
            # (su - sv) / (u - v)
            num = su - sv
            den = u - v
            if (num > 0) != (den > 0):  # different signs
                sign *= -1
        return sign

    def is_rp_even(bits):
        """RP-even: sgn_G(sigma) = +1 for ALL automorphisms."""
        # Find automorphisms
        for sigma in all_perms:
            # Check if sigma is an automorphism
            is_aut = True
            for k, (i,j) in enumerate(edges):
                pi, pj = min(sigma[i], sigma[j]), max(sigma[i], sigma[j])
                nk = edge_idx[(pi, pj)]
                if ((bits >> k) & 1) != ((bits >> nk) & 1):
                    is_aut = False
                    break
            if not is_aut: continue
            # Check sign
            s = graph_sign(bits, sigma)
            if s == -1:
                return False  # odd automorphism exists -> RP-odd
        return True  # all automorphisms have sign +1 -> RP-even

    def is_euler(bits):
        deg = [0] * N
        for k, (i,j) in enumerate(edges):
            if bits & (1 << k):
                deg[i] += 1; deg[j] += 1
        return all(d % 2 == 0 for d in deg)

    print(f"\n  n = {n}: classifying ALL {2**m} graphs")

    rp_even_classes = set()
    euler_classes = set()
    both_classes = set()
    rp_only_classes = set()
    euler_only_classes = set()
    neither_classes = set()

    for bits in range(1 << m):
        cn = graph_canon(bits)
        rp = is_rp_even(bits)
        eu = is_euler(bits)

        if rp and eu: both_classes.add(cn)
        elif rp and not eu: rp_only_classes.add(cn)
        elif eu and not rp: euler_only_classes.add(cn)
        else: neither_classes.add(cn)

        if rp: rp_even_classes.add(cn)
        if eu: euler_classes.add(cn)

    print(f"    RP-even classes:  {len(rp_even_classes)} (= A000568({n}) = {A000568[n]})")
    print(f"    Euler classes:    {len(euler_classes)} (= A002854({n}) = {A002854[n]})")
    print(f"    Both:             {len(both_classes)}")
    print(f"    RP-even only:     {len(rp_only_classes)}")
    print(f"    Euler only:       {len(euler_only_classes)}")
    print(f"    Neither:          {len(neither_classes)}")
    print(f"    Total classes:    {len(rp_even_classes | euler_classes | rp_only_classes | euler_only_classes | neither_classes)}")

    # THE KEY: Euler ⊂ RP-even? Or RP-even ⊂ Euler? Or neither?
    if euler_only_classes:
        print(f"    Euler NOT subset of RP-even (Euler-only exists)")
    if rp_only_classes:
        print(f"    RP-even NOT subset of Euler (RP-only exists)")
    if not euler_only_classes and not rp_only_classes:
        print(f"    RP-even = Euler (identical sets)")
    if euler_classes <= rp_even_classes:
        print(f"    Euler IS A SUBSET of RP-even")
    if rp_even_classes <= euler_classes:
        print(f"    RP-even IS A SUBSET of Euler")

print(f"\n{'='*60}")
print("SYNTHESIS")
print(f"{'='*60}")
print("""
The trivial vertex connects three counting families:
  Euler(n): all-even-degree graphs, counted by cycle space (dim C(n-1,2))
  RP-even(n): sign-even graphs, counted by Burnside (= A000568(n))
  Tournament(n): complete directed, counted by Burnside (= A000568(n))

The relationship Euler vs RP-even determines which graphs are
"reachable" from the tiling model:
  - If Euler subset RP-even: every Euler graph is also RP-even,
    and the extra RP-even graphs are non-Euler (have odd degrees)
  - This would mean: the cycle space bijection covers the Euler part,
    and the remaining A000568 - A002854 RP-even graphs need a
    DIFFERENT construction (using the trivial vertex's freedom)
""")

print("DONE.")
print("=" * 80)
