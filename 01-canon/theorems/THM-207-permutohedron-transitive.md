# THM-207: Transitive Tournament Flip Graph = Permutohedron

**Status:** VERIFIED (n=3..6), PROOF COMPLETE
**Proved by:** opus-2026-03-14-S89b
**Verified in:** `04-computation/H_6t3_regression_89b.py`

## Statement

The subgraph of the tournament hypercube Q_m (m = n(n-1)/2) induced on the n! transitive tournaments is isomorphic to the **permutohedron** Π_n — the Cayley graph of S_n with generators = adjacent transpositions {s_1, ..., s_{n-1}}.

## Consequences

1. **Regularity**: The transitive flip graph is **(n-1)-regular**.
2. **Diameter**: The diameter is C(n,2) = m (the full hypercube dimension).
3. **Vertex count**: n! vertices, n!(n-1)/2 edges.
4. **H is constant**: H(T) = 1 for all transitive tournaments, so H is constant on the permutohedron.
5. **Adjacent flips create 3-cycles**: Flipping a non-adjacent pair in the score ordering always creates a 3-cycle.

## Proof

1. Each transitive tournament T_σ corresponds to a unique permutation σ ∈ S_n (the score ordering: σ(1) beats σ(2) beats ... beats σ(n)).

2. Flipping arc (i,j) in T_σ changes the relative order of i and j. If i and j are **adjacent in σ** (i.e., σ^{-1}(i) and σ^{-1}(j) differ by 1), then the result is another transitive tournament T_τ where τ = σ ∘ s_k for the appropriate adjacent transposition s_k.

3. If i and j are **not adjacent** in σ (there exists a vertex v between them in the score ordering), then flipping (i,j) creates a directed 3-cycle involving i, j, and v. So the result is NOT transitive.

4. Therefore T_σ and T_τ are at Hamming distance 1 in Q_m and both transitive ⟺ σ and τ differ by an adjacent transposition. This is exactly the permutohedron. □

## Computational Verification

| n | n! | Degree | Diameter | Expected C(n,2) | Match |
|---|-----|--------|----------|-----------------|-------|
| 3 | 6   | 2      | 3        | 3               | ✓     |
| 4 | 24  | 3      | 6        | 6               | ✓     |
| 5 | 120 | 4      | 10       | 10              | ✓     |
| 6 | 720 | 5      | 15       | 15              | ✓     |

## Connection to Other Results

- **THM-206**: Local minima of H = transitive tournaments = the permutohedron.
- **THM-205**: Cone preserves H, and the cone of a permutohedron vertex stays on the permutohedron of the larger tournament.
- The **Hamming distance distribution** between transitive tournaments is symmetric around m/2, matching the distribution of inversions in S_n.
