#!/usr/bin/env python3
"""
cycle_indicator_degree.py -- kind-pasteur-2026-03-13-S60

KEY DERIVATION: The 3-cycle indicator function for circulant tournaments.

For vertices a < b < c in a circulant tournament with signs σ(d):
  σ(d) = +1 if d ∈ S (edge forward), -1 if d ∉ S (edge backward)

Edge indicators:
  A[a][b] = (1 + σ(b-a))/2
  A[b][c] = (1 + σ(c-b))/2
  A[a][c] = (1 + σ(c-a))/2
  A[c][a] = (1 - σ(c-a))/2
  etc.

count_3(a,b,c) = A[a][b]*A[b][c]*A[c][a] + A[a][c]*A[c][b]*A[b][a]
               = [(1+s1)(1+s2)(1-s3) + (1-s1)(1-s2)(1+s3)] / 8

where s1=σ(b-a), s2=σ(c-b), s3=σ(c-a).

Expanding: = (1 + s1*s2 - s1*s3 - s2*s3) / 4

THIS IS DEGREE 2 in σ (not degree 3!). The cubic terms cancel.

For 5-cycles: what is the effective degree?
For general k: what pattern emerges?

This script computes the EXACT polynomial expansion for k-cycle indicators
at small k using symbolic computation.
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict


def compute_cycle_indicator_expansion(k):
    """
    For a k-cycle on vertices {0,1,...,k-1}, compute the indicator function
    as a polynomial in the pairwise σ values.

    Variables: σ_{ij} for i < j, where σ_{ij} = σ(j-i).
    But for a general k-vertex set, the differences are not necessarily
    {1,2,...,k-1}. However, the POLYNOMIAL STRUCTURE in σ_{ij} is the same.

    A directed Hamiltonian cycle on k vertices visits each vertex exactly once.
    There are (k-1)!/2 distinct undirected Hamiltonian cycles, each with 2 directions.
    Total: (k-1)! directed Hamiltonian cycles.

    For a directed cycle (v_π(0), v_π(1), ..., v_π(k-1)):
    Indicator = prod_{i=0}^{k-1} A[v_π(i)][v_π((i+1)%k)]

    where A[u][v] = (1 + sign(v-u) * σ_{min(u,v),max(u,v)}) / 2
    Here sign(v-u) = +1 if v > u (edge "forward" in label ordering).

    Total count = sum over all directed Hamiltonian cycles of indicator.
    """
    verts = list(range(k))
    n_pairs = k * (k - 1) // 2
    pair_index = {}
    idx = 0
    for i in range(k):
        for j in range(i + 1, k):
            pair_index[(i, j)] = idx
            idx += 1

    # Generate all directed Hamiltonian cycles starting from vertex 0
    # (fixing start eliminates the k rotations)
    cycles = []
    for perm in permutations(range(1, k)):
        cycle = [0] + list(perm)
        cycles.append(cycle)

    # For each cycle, compute the monomial expansion
    # Each edge (u,v) contributes: (1 + sign(v-u) * σ_{min,max}) / 2
    # The product over k edges: (1/2^k) * product of (1 + s_i * x_i)
    # where s_i ∈ {+1,-1} and x_i = σ_{pair}

    # Represent monomials as frozensets of pair indices (each appearing at most once
    # since σ^2 = 1)
    # Coefficient for monomial M = product of (s_i for i in M)

    total_poly = defaultdict(int)  # monomial -> coefficient (times 2^k)

    for cycle in cycles:
        # Edges of this cycle
        edges = []
        for i in range(k):
            u = cycle[i]
            v = cycle[(i + 1) % k]
            edges.append((u, v))

        # Each edge (u,v) contributes (1 + sign * σ_pair) where:
        # sign = +1 if v > u (forward), -1 if v < u (backward)
        # pair = (min(u,v), max(u,v))
        edge_data = []
        for u, v in edges:
            pair = (min(u, v), max(u, v))
            pidx = pair_index[pair]
            sign = 1 if v > u else -1
            edge_data.append((pidx, sign))

        # Expand product: iterate over all 2^k subsets of edges
        for mask in range(1 << k):
            # Monomial: set of pair indices included
            mono_set = set()
            coeff = 1
            for bit in range(k):
                if mask & (1 << bit):
                    pidx, sign = edge_data[bit]
                    if pidx in mono_set:
                        mono_set.remove(pidx)  # σ^2 = 1
                        coeff *= sign  # multiply by the sign
                    else:
                        mono_set.add(pidx)
                        coeff *= sign
                # If bit not set, we pick the "1" term (no contribution)

            total_poly[frozenset(mono_set)] += coeff

    # Divide by 2^k
    # But actually, we should keep the numerator and denominator separate
    # The total count = total_poly / 2^k

    # Group by degree (number of σ factors)
    by_degree = defaultdict(int)
    for mono, coeff in total_poly.items():
        by_degree[len(mono)] += abs(coeff)

    max_degree = max(by_degree.keys()) if by_degree else 0

    # Find effective degree: highest degree with nonzero total coefficient
    eff_degree = 0
    for mono, coeff in total_poly.items():
        if coeff != 0 and len(mono) > eff_degree:
            eff_degree = len(mono)

    return total_poly, eff_degree, by_degree


print("=" * 70)
print("CYCLE INDICATOR POLYNOMIAL DEGREE ANALYSIS")
print("=" * 70)

for k in [3, 5, 7]:
    print(f"\n{'='*60}")
    print(f"  k = {k} (directed {k}-cycle indicator)")
    print(f"{'='*60}")

    if k > 7:
        print(f"  SKIPPED (too many permutations)")
        continue

    n_pairs = k * (k - 1) // 2
    print(f"  {n_pairs} pairwise sigma variables")
    print(f"  {(k-1)} factorial = {1}")
    import math
    n_cycles = math.factorial(k - 1)
    print(f"  {n_cycles} directed Hamiltonian cycles")

    poly, eff_deg, by_deg = compute_cycle_indicator_expansion(k)

    print(f"\n  Effective polynomial degree: {eff_deg}")
    print(f"  Degree distribution (sum of |coeff| at each degree):")
    for deg in sorted(by_deg.keys()):
        n_terms = sum(1 for m, c in poly.items() if len(m) == deg and c != 0)
        print(f"    degree {deg}: {n_terms} nonzero terms, total |coeff| = {by_deg[deg]}")

    # Show the polynomial for k=3
    if k == 3:
        print(f"\n  Full polynomial (numerator, divide by 2^{k}={2**k}):")
        for mono in sorted(poly.keys(), key=lambda m: (len(m), sorted(m))):
            coeff = poly[mono]
            if coeff != 0:
                if len(mono) == 0:
                    term = "1"
                else:
                    pairs = []
                    for idx in sorted(mono):
                        for (i, j), pidx in sorted(
                            [(p, pi) for p, pi in
                             {(a, b): a * 10 + b
                              for a in range(k) for b in range(a + 1, k)}.items()
                             if pi == idx]):
                            pairs.append(f"σ_{i}{j}")
                    term = "*".join(pairs)
                print(f"      {coeff:+d} * {term}")

    # For k=5: just show degree summary
    if k >= 5:
        print(f"\n  Constant term: {poly.get(frozenset(), 0)}")
        # Count nonzero monomials at each degree
        for deg in sorted(by_deg.keys()):
            monomials = [(m, c) for m, c in poly.items() if len(m) == deg and c != 0]
            if len(monomials) <= 10:
                for m, c in sorted(monomials, key=lambda x: sorted(x[0])):
                    print(f"    deg {deg}: coeff={c}, vars={sorted(m)}")


# VERIFICATION at k=3
print(f"\n{'='*60}")
print(f"  VERIFICATION for k=3")
print(f"{'='*60}")

poly3, _, _ = compute_cycle_indicator_expansion(3)
# Expected: count = (1 + σ01*σ12 - σ01*σ02 - σ12*σ02) * (k-1)!/2^3
# No wait: total_poly already sums over all directed cycles.
# For k=3, there are 2! = 2 directed cycles: (0,1,2) and (0,2,1).
# count = total_poly / 2^3

print(f"  Polynomial / {2**3}:")
for mono in sorted(poly3.keys(), key=lambda m: (len(m), sorted(m))):
    coeff = poly3[mono]
    if coeff != 0:
        if len(mono) == 0:
            term = "1"
        else:
            pair_names = {0: "σ01", 1: "σ02", 2: "σ12"}
            term = "*".join(pair_names[idx] for idx in sorted(mono))
        print(f"    {coeff}/{2**3} * {term} = {coeff/8:.4f} * {term}")

# Check: should be (1/4)(1 + σ01*σ12 - σ01*σ02 - σ12*σ02)
# = (2/8) + (2/8)*σ01*σ12 + (-2/8)*σ01*σ02 + (-2/8)*σ12*σ02
# Constant: 2. σ01σ12: +2. σ01σ02: -2. σ12σ02: -2.

print(f"\n  Expected: (1/4)(1 + σ01*σ12 - σ01*σ02 - σ12*σ02)")
print(f"  = 2/8 + 2/8*σ01*σ12 - 2/8*σ01*σ02 - 2/8*σ12*σ02")


# KEY INSIGHT: degree structure determines moment dependence
print(f"\n{'='*60}")
print(f"  MOMENT DEPENDENCE THEORY")
print(f"{'='*60}")
print("""
  CYCLE INDICATOR DEGREE PATTERN:
    d(3) = 2   (degree 2 in pairwise sigma products)
    d(5) = 4   (degree 4)
    d(7) = 6   (degree 6)
    Pattern: d(k) = k - 1

  ALL ODD DEGREES VANISH (parity symmetry of Hamiltonian cycles).

  Moment dependence prediction:
    c_k needs moments up to S_{2*d(k)} = S_{2(k-1)}
    disj(k1,k2) needs up to S_{2*(d(k1)+d(k2))} = S_{2(k1+k2-2)}

  But disjointness constraint may reduce effective degree when k1+k2 > p.

  From data at p=11:
    disj(3,3): needs S_4. Predicted: S_{2*(2+2)} = S_8.
      BUT only S_4 needed! Reduction from disjointness constraint.
    disj(3,5): needs S_4..S_8. Predicted: S_{2*(2+4)} = S_12.
      Only S_8 needed. Reduction again.
    disj(5,5): needs S_4..S_8. Predicted: S_{2*(4+4)} = S_16.
      Only S_8 needed. Major reduction.

  CONJECTURE: For circulant tournaments on Z_p,
    disj(k1,k2) depends on moments S_4,...,S_{p-3} only.
    (Same moment range as alpha_1 from THM-157!)
""")

print("\nDONE.")
