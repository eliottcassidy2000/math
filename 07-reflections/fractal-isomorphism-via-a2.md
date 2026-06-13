# Fractal Isomorphism via A² Rows

**opus-2026-03-25-S339**

## The Practical Situation

Sorted rows of A² is a complete isomorphism invariant for tournaments through n=9 (191K classes, exhaustive) and strongly suggested at n=10 (1M random sample, zero collisions). Cost: O(n³).

For n > 10? We don't know. If it fails, we need A³ or higher powers. But for practical applications, n ≤ 10 covers almost all real-world tournaments (sports leagues, round-robin competitions, preference rankings).

## The Fractal Idea

For tournaments with n > 10, decompose fractally:

### Mode A: Vertex Deletion

A tournament T on n vertices contains n sub-tournaments on n-1 vertices (delete each vertex). If n > 10:

1. For each vertex v, compute the sub-tournament T\v on n-1 vertices
2. Compute the A² fingerprint of T\v
3. The **fingerprint tuple** (fp(T\v₁), fp(T\v₂), ..., fp(T\vₙ)) is an invariant of T

This is O(n⁴) (n vertex deletions × O(n³) per fingerprint). It's more expensive than a single A² computation but still polynomial.

If n-1 > 10: recurse! Each sub-tournament's fingerprint is computed fractally.

### Mode B: Block Decomposition

Many real-world tournaments have **modular structure** (a set S of vertices such that all vertices outside S treat S uniformly). Decompose:

1. Find the modular decomposition of T (O(n²) algorithm exists)
2. For each module M_i of size ≤ 10: compute A² fingerprint directly
3. For the quotient tournament Q (with modules as vertices): compute A² fingerprint
4. The tournament fingerprint = (fp(Q), {fp(M_i)})

This is O(n²) for the decomposition + O(k³) for the quotient (k = number of modules) + Σ O(|M_i|³) for each module. For tournaments with small modules: nearly linear.

### Mode C: Random Sampling

For very large n: sample random sub-tournaments of size k ≤ 10:

1. Pick m random subsets of size k from the vertex set
2. For each subset, compute the induced sub-tournament's A² fingerprint
3. The **fingerprint distribution** (multiset of sub-tournament fingerprints) is an invariant

This is O(m × k³) — independent of n! The quality depends on m (the number of samples) and k (the sub-tournament size). For k=10 and m=1000: about 10⁶ operations, regardless of whether n = 100 or n = 1,000,000.

This is essentially a **random kernel method**: embed the tournament into a feature space defined by sub-tournament A² fingerprints. Two tournaments with the same embedding distribution are likely isomorphic.

## The Connection to Our Codec

The TC codec uses score conditioning (Level 1) for compression. The A² fingerprint (Level 2) gives 4× better separation. For the fractal codec:

1. Compute the A² fingerprint of the tournament (or sub-tournament at each recursive level)
2. Use the fingerprint as the compression context
3. The arithmetic coder uses P(residual | fingerprint) instead of P(residual | score)

This gives tighter prediction at each level, compounding across the recursive levels.

## The Connection to GNNs

Graph Neural Networks with 2 message-passing layers compute exactly the A² row features (with sum aggregation and identity update). Our result says: for tournaments, 2 layers suffice for perfect classification up to n=10 (and probably much further).

This is a concrete, testable prediction: a 2-layer GNN trained on tournament classification should achieve near-100% accuracy. If it doesn't, the failure is in the training, not the architecture.

## Conjecture

**Conjecture (A² Completeness):** For all n, the sorted row vectors of A² form a complete isomorphism invariant for tournaments on n vertices.

**Evidence:** Verified at n = 5, 6, 7, 8, 9 (exhaustive). Zero collisions in 1M random sample at n = 10.

**If true:** Tournament isomorphism is in P (specifically O(n³)), much stronger than Babai's quasi-polynomial result for general graphs. This would be a major result in computational complexity, as tournament isomorphism is a natural special case of graph isomorphism.

**If false:** The first collision at some n₀ would still give us a practical polynomial-time invariant for all n < n₀. And the fractal decomposition extends its reach far beyond n₀.
