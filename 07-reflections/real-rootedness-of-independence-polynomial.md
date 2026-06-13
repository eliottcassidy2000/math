# Real-Rootedness of the Odd-Cycle Independence Polynomial

**Session:** opus-2026-04-05-S24
**Status:** THEOREM (proved for n ≤ 8, conjectured for all n)

## The Discovery

For any tournament T, define Ω(T) as the conflict graph on all directed odd cycles (vertices = odd cycles, edges = pairs sharing a vertex). The independence polynomial is:

I(Ω(T), x) = Σ_k α_k x^k

where α_k = # independent sets of size k (= # collections of k pairwise vertex-disjoint odd cycles).

**Key fact:** H(T) = I(Ω(T), 2) (the OCF, proved by Grinberg-Stanley).

## Results

### Theorem (Real-Rootedness for n ≤ 8)

**For any tournament T on n ≤ 8 vertices, all zeros of I(Ω(T), x) are real and negative.**

**Proof:** By the Chudnovsky-Seymour theorem (2004), the independence polynomial of any claw-free graph has only real zeros. We show Ω(T) is claw-free for n ≤ 8.

A claw K_{1,3} in Ω(T) requires 3 pairwise non-adjacent vertices — i.e., 3 pairwise vertex-disjoint odd cycles C₁, C₂, C₃. Each odd cycle uses ≥ 3 vertices, so C₁ ∪ C₂ ∪ C₃ uses ≥ 9 vertices. But |V(T)| = n ≤ 8. Contradiction.

Therefore Ω(T) is claw-free for n ≤ 8, and by Chudnovsky-Seymour, I(Ω(T), x) has all real zeros. Since all coefficients α_k ≥ 0 and α_0 = 1 > 0, the zeros must be negative. □

### Computational Verification

Exhaustively verified for ALL iso classes at n = 3, 4, 5, 6 (74 classes total, 0 exceptions). Sampled 3000 tournaments at n = 7 — all have real negative zeros.

### Conjecture (All n)

**Conjecture:** I(Ω(T), x) has all real negative zeros for ALL tournaments T, at all n.

At n = 9, the conflict graph Ω(T) CAN have claws (three disjoint 3-cycles on 9 vertices), so the Chudnovsky-Seymour argument breaks. But real-rootedness may still hold — claw-free is sufficient but not necessary.

## Consequences

If I(Ω(T), x) has all real negative zeros, we can write:

**I(Ω(T), x) = Π_i (1 + x/r_i)** where r_i > 0

Then **H(T) = Π_i (1 + 2/r_i)**, decomposing the Hamiltonian path count as a product over "spectral modes" of the conflict graph.

### Log-concavity of α_k

Real-rootedness implies the coefficient sequence (α_0, α_1, ..., α_d) is **log-concave**: α_k² ≥ α_{k-1} α_{k+1}. This constrains which combinations (α_0, α_1, α_2, ...) can appear.

### The polynomial as a finer invariant

I(Ω(T), x) is strictly finer than H(T) = I(Ω(T), 2). At n = 6, tournaments with the same H can have different polynomials:
- H = 9: both I = 1+4x and I = 1+2x+x²
- H = 17: both I = 1+8x and I = 1+6x+x²
- H = 45: both I = 1+20x+x² and I = 1+14x+4x²

## Connection to Other Results

### Forbidden H Values
The α₁=3 gap is related to real-rootedness: when α₁=3 and the polynomial is real-rooted, the constraint on α�� becomes tighter. Specifically, log-concavity gives α₁² ≥ α_0 · α_2, i.e., 9 ≥ α_2. But the actual constraint is stronger: at n≤6, α₁=3 is impossible; at n=7, α₁=3 forces α₂=2.

### Hard-Core Lattice Gas
I(Ω(T), x) is the partition function of the hard-core lattice gas on Ω(T) at fugacity x. Real-rootedness means the partition function has no phase transition on the positive real axis — the system is always in the "gas phase" for x > 0. The evaluation at x = 2 (giving H(T)) is safely within this region.

### Heilmann-Lieb and the Matching Polynomial
The real-rootedness of independence polynomials of claw-free graphs generalizes the Heilmann-Lieb theorem for matching polynomials. This places tournament cycle-packing in the same algebraic framework as dimer models.
