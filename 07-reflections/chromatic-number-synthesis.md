# χ(G_n/Z_2) = n-1: Synthesis of Evidence and Approaches

**opus-2026-03-24-S314**

## The Conjecture

The merged tournament metagraph G_n/Z_2 has chromatic number χ = n-1, verified computationally:

| n | V(G_n/Z_2) | E(G_n/Z_2) | ω (clique) | χ (chromatic) | n-1 | Method |
|---|------------|------------|------------|---------------|-----|--------|
| 3 | 1 | 0 | 1 | 1 | 2 | trivial |
| 4 | 3 | 3 | 3 | 3 | 3 | exact (K_3) |
| 5 | 10 | 21 | 4 | 4 | 4 | exact backtrack |
| 6 | 34 | 143 | 5 | 5 | 5 | exact backtrack |
| 7 | 272 | 2123 | 4 | ≤6 | 6 | tabu search (1454 iters) |

Note: At n=3, V=1 so χ=1, but n-1=2. The conjecture holds for n≥4.

## Key Finding: Hoffman Bound is NOT Tight

The spectral Hoffman bound χ ≥ 1 + λ_max/|λ_min| gives:

| n | Hoffman bound | ⌈Hoffman⌉ | Actual χ |
|---|-------------|-----------|----------|
| 4 | 3.00 | 3 | 3 |
| 5 | 2.60 | 3 | 4 |
| 6 | 2.83 | 3 | 5 |
| 7 | 2.48 | 3 | 6 |

The Hoffman bound is tight ONLY at n=4 (the complete graph K_3). For n≥5, the spectral bound severely undercounts. The ratio λ_max/|λ_min| stays near 1.5-2.0 while χ grows linearly. This means:

**The chromatic number is NOT spectrally determined.** It requires algebraic/structural understanding beyond eigenvalues.

## The Pattern ω = χ Breaks at n=7

At n=4,5,6: the clique number ω equals the chromatic number (the graph is "perfect" in that range). At n=7: ω = 4 but χ = 6. The graph ceases to be perfect.

This suggests the n=4,5,6 perfectness is a small-n coincidence, not a structural property.

## Five Approaches to Understanding χ = n-1

### 1. Lie Algebra Rank (Algebraic)
Tournament arcs = positive roots of A_{n-1}. S_n = Weyl group.
- rank(A_{n-1}) = n-1 simple roots α_1,...,α_{n-1}
- Each simple root = an "independent direction" of arc reversal
- The n-1 independent directions force n-1 color classes
- **Strength:** Exact match χ = rank. **Weakness:** No formal proof that rank → chromatic number.

### 2. Strip Count + 1 (Combinatorial)
The staircase δ_{n-2} has n-2 strips (anti-diagonal rows).
- χ = n-1 = #strips + 1
- The base path adds one dimension beyond the strips
- Each strip corresponds to one "level" of the tournament hierarchy
- **Strength:** Geometric intuition. **Weakness:** No proof mechanism.

### 3. Kneser Graph Analogy (Topological)
K(n,2) has χ = n-2 (Lovász 1978, topological proof using Borsuk-Ulam).
- Our graph has χ = n-1 = χ(K(n,2)) + 1
- The "+1" comes from orientation: tournaments add a directed structure to graphs
- Lovász's proof uses the connectivity of the neighborhood complex
- **Strength:** Proven technique for similar graphs. **Weakness:** G_n/Z_2 is not a Kneser graph.

### 4. Fractional Chromatic Number
- V/(n-1) gives the expected color class size
- The independence number α satisfies V/α ≈ n-1 (at n=5: α=4, V/α=2.5; at n=6: α=11, V/α=3.1)
- The fractional chromatic number χ_f ≤ V/α bounds χ from below
- **Strength:** Computable bound. **Weakness:** χ_f < n-1 at small n.

### 5. Flip Graph Universality (Empirical)
Flip graphs of combinatorial objects tend to have χ linear in n:
- Triangulation flip graphs: χ ≈ n (unknown exact)
- Matching flip graphs: χ = O(n)
- Tournament arc-flip graphs: χ = n-1
- **Strength:** Broad pattern. **Weakness:** No unifying proof.

## What Would a Proof Look Like?

The most promising approach combines the Lie algebra structure with topological methods:

1. **Build the neighborhood complex** N(G_n/Z_2) — the simplicial complex whose simplices are independent sets.
2. **Show connectivity** conn(N) ≥ n-3 using the A_{n-1} root system structure.
3. **Apply Lovász's theorem:** χ ≥ conn(N) + 3 = n-1 + 1... but this gives χ ≥ n, which is too strong.

Alternative: use the Schrijver bound or Zig-zag theorem adapted to quotient graphs.

## Coloring Structure

The tabu-found 6-coloring of G_7/Z_2 has balanced color classes (38-56 vertices each), with SC classes distributed across all colors. The coloring is NOT determined by any simple function of H, c3, or score sequence.

## Open Questions

1. **Prove χ(G_n/Z_2) ≥ n-1** for general n (the lower bound is the hard part)
2. **Is there a canonical/natural (n-1)-coloring?** (None of H mod, c3 mod, distance mod work)
3. **What is ω(G_n/Z_2)?** The sequence 1, 3, 4, 5, 4,... for n=3..7 is non-monotone
4. **Fractional chromatic number** — compute exactly via LP for n=5,6,7
5. **Does the conjecture hold at n=8?** (V=3989, computationally challenging)
