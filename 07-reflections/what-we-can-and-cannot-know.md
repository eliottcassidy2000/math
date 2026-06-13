# What We Can and Cannot Know: An Information-Theoretic Audit

**opus-2026-03-24-S307**

## Framework

For each structural quantity, classify it as:
- **KNOWN**: we have an exact formula or algorithm that computes it in polynomial time for all n
- **UNKNOWABLE**: we can prove no polynomial algorithm exists, or the quantity is undecidable, or it requires exponential information to specify
- **BOUNDARY**: we have partial results, conjectures, or exponential-time algorithms — these ARE the open questions

---

## KNOWN (exact formulas or polynomial algorithms)

### 1. The number of iso classes V_n
**Formula**: V_n = (1/n!) Σ_{σ ∈ S_n, all-odd-cycles} 2^{pair_orbits(σ)}
**Complexity**: O(p(n)) where p(n) = number of odd partitions of n. Polynomial in n for the Burnside sum (the partition enumeration is subexponential).
**Status**: SOLVED. Computable to arbitrary n via the turbo engine (n=50 in 57ms).

### 2. The Euler product correction R(n)
**Formula**: R(n) = Σ_{k odd ≥ 3} (1/k) × n↓k × 2^{(k²-1)/2 - (k-1)n} + cross-orbit terms
**Status**: SOLVED to any desired precision. The 3-cycle term alone gives 99.98% accuracy by n=15.

### 3. Average tiling fiber
**Formula**: ⟨fiber⟩ = n!/2^{n-1} (the Szele average), exact as R(n) → 0.
**Status**: SOLVED.

### 4. The fiber of a specific class C
**Formula**: fiber(C) = H(C)/|Aut(C)|
**Complexity**: Computing H is #P-complete in general, but for tournaments with a known HP, it's the tiling count, which requires knowing the class.
**Status**: KNOWN given the class representative. The identity fiber × |Aut| = H is exact.

### 5. The waggly layer structure
**Formula**: Layer d has C(m, d) connections per tiling. Total = 2^m - 1.
**Status**: SOLVED. Pure combinatorics.

### 6. The completeness threshold k*
**Formula**: k* = diam(G_n) = A003141(n) = max_T min-FAS(T)
**Status**: SOLVED in principle. A003141 is computable (though the max-FAS problem for tournaments is polynomial, computing A003141 requires checking all tournament types).

### 7. The Krawtchouk spectral expansion
**Fact**: H is band-limited — hat{H}_k = 0 for k ≥ ⌈m/2⌉ + 1.
**Status**: OBSERVED at n=5,6,7. If proved for all n, this would be a major theorem.

### 8. Self-loop rates by distance d
**Formula**: Computable exactly from the weight matrices W_d.
**Status**: KNOWN at n ≤ 7. Asymptotic formula: SL(d=1) ~ n²/4^n from the Euler product.

### 9. Paley code parameters
**Formula**: P_p with A+I as parity-check gives QR code: [p, (p+1)/2, d] for p ≡ 3 mod 4.
**Specific**: P₇ → [7,4,3] Hamming, P₂₃ → [23,12,7] Golay.
**Status**: SOLVED.

### 10. The OCF: H(T) = I(Ω(T), 2)
**Status**: PROVED (Grinberg-Stanley 2024). The Hamiltonian path count equals the independence polynomial of the odd-cycle conflict graph evaluated at 2.

---

## UNKNOWABLE (provably hard or requiring exponential information)

### 1. The full metagraph G_n for large n
**Why**: V_n grows super-exponentially. At n=10, V_merged ≈ 4.9M classes. Storing the adjacency matrix requires V_n² entries. No polynomial-in-n description exists for the full graph.
**What we CAN know**: statistical properties (degree distribution, spectral gap, diameter).

### 2. H(T) for a specific labeled tournament in polynomial time
**Why**: Computing #Hamiltonian paths is #P-complete for general digraphs. For tournaments, it's polynomial via the transfer matrix method (O(2^n × n²)), but this is exponential in n.
**What we CAN know**: H modulo small primes (polynomial via modular methods), expected H for random tournaments.

### 3. The exact weight matrix W_d for large n
**Why**: W_d has V_n² entries, each requiring enumeration of tilings. No closed-form is known.
**What we CAN know**: W_d's row sums (= fiber × C(m,d)), spectral properties, statistical moments.

### 4. Whether a specific tournament has β₁ = 1 in polynomial time
**Why**: Computing path homology requires rank computations on matrices of size O(n^k), which is polynomial for fixed k but the implicit constants grow.
**Status**: Actually this IS polynomial for fixed k. The unknowability is about the asymptotic behavior of β_k as n → ∞.

### 5. The full Krawtchouk spectrum of the quotient
**Why**: The quotient is not an association scheme (full algebra dim = 35 vs 7 at n=5). No algebraic shortcut to the exact spectrum. Must compute V_n × V_n matrices.

---

## BOUNDARY: The Open Questions

These are the quantities where we have partial knowledge and the boundary between known and unknowable is the EXACT research frontier:

### B1. Is H band-limited for ALL n?
**Known**: hat{H}_k = 0 for k ≥ ⌈m/2⌉ + 1 at n=5,6,7.
**Unknown**: whether this holds for all n. If true, it would be a deep structural theorem about tournaments: H is determined by the "low-frequency" half of the tiling word.
**What would settle it**: A proof using the Walsh-Fourier framework (THM-069, THM-077) showing that even-length path unions in the Walsh decomposition can't reach high Krawtchouk degrees.

### B2. The exact value of A003141(n) = diam(G_n) for all n
**Known**: Values for n ≤ ~30 (from the literature on feedback arc sets).
**Known asymptotic**: n²/4 - Θ(n^{3/2}).
**Unknown**: The exact constant in the Θ(n^{3/2}) correction. This is equivalent to understanding the maximum "disorder" of a tournament.
**Connection**: The tournaments achieving max min-FAS are conjectured to be Paley tournaments for prime n.

### B3. Whether diam(G_n) is always achieved by the transitive class
**Known**: True at n=5,6 (the transitive is always one endpoint of the diameter pair).
**Unknown**: For all n. If the transitive is always a diameter endpoint, then diam = max_C min_dist(transitive, C), which simplifies computation enormously.

### B4. The rate of K₁-H correlation decay
**Known**: corr = -0.939, -0.886, -0.832 for n=5,6,7 (class-level).
**Unknown**: The exact limiting behavior. Does corr → 0? Does it stabilize at some positive value? The rate of decay determines how well the Hamming weight approximation works at large n.
**If corr → 0**: tournament space is genuinely 2D at large n.
**If corr → c > 0**: there's a permanent 1D component that never goes away.

### B5. The exact neutrality formula
**Known**: SL rate at d=1 decays as ~n²/4^n (from Euler product).
**Unknown**: Exact formula for SL rate at distance d as a function of d and n. The self-loop spectrum (SL vs d) has non-monotone structure (peak at d=2, dip at d=3) that lacks a formula.
**Connection**: This is equivalent to computing the autocorrelation of the iso-class indicator function on Q_m, which involves the Krawtchouk transform of the class partition.

### B6. The width function W(H) = #{classes with Hamiltonian path count = H}
**Known**: W(H) at small n. Not unimodal. Max width grows with n.
**Unknown**: The asymptotic distribution of W(H). Is it Gaussian? Log-normal? Does it concentrate?
**Connection**: W(H) is the "density of states" of tournament space. Its shape determines the funnel geometry.

### B7. Whether the metagraph quotient is ever an association scheme
**Known**: NOT a scheme at n=5 (algebra dim 35 vs 7). Likely not for any n ≥ 5.
**Unknown**: Is there a DIFFERENT partition of the tiling pairs (not by Hamming distance) that DOES give an association scheme? E.g., a partition based on the iso-class structure itself.
**If yes**: This would give a complete algebraic framework for tournament theory.

### B8. The exact all-same-range neutrality
**Known**: All-range-3 combos have highest SL at n=5,6. All-range-2 has lowest (zero at n=5).
**Unknown**: Is range ⌊(n-1)/2⌋ always the most neutral? This would be a "mid-frequency is most neutral" theorem = tournament uncertainty principle.
**What would settle it**: Computing at n=7,8 and proving via Walsh analysis.

### B9. β₂ = 0 for all tournaments
**Known**: Verified exhaustively n ≤ 8, sampled n ≤ 10. Twin vertex mechanism identified.
**Unknown**: Proof for all n. This is the longest-standing open conjecture in the project.
**Connection**: β₂ = 0 means the path homology has no 2-dimensional "holes" — tournaments are homologically simple at dimension 2.

### B10. The connection between H and the minimum feedback arc set
**Known**: diam = max min-FAS = A003141. The class farthest from transitive has the largest min-FAS.
**Unknown**: Is there a FORMULA for min-FAS(T) in terms of H(T), |Aut(T)|, or other invariants? The min-FAS is NP-hard for general digraphs but polynomial for tournaments — can we express it in the OCF framework?

---

## The Information-Theoretic Summary

**Bits needed to specify a tournament**: C(n,2) = m_arc bits (the full adjacency).
**Bits needed to specify its iso class**: log₂(V_n) ≈ m_arc - log₂(n!) bits.
**Bits captured by H alone**: log₂(#{distinct H values}) ≈ log₂(n!/2^{n-1}) bits.
**Bits captured by K₁ (Hamming weight)**: log₂(m+1) ≈ 2 log₂(n) bits.

The gap: log₂(V_n) - 2 log₂(n) ≈ n² - n log n bits are NOT captured by K₁.
These residual bits live in the higher Krawtchouk modes k=2,3,...,⌊m/2⌋.
The band-limitedness means at most ⌊m/2⌋ modes contribute — so the effective information content of an iso class is at most m/2 ≈ n²/4 bits, not the full m ≈ n²/2.

**Tournament space has half the information content you'd expect from its dimension.**

This is because the staircase triangle has two legs but one hypotenuse — the diagonal direction (the band-limited one) carries the same information as one leg, not two. The other "half" of the spectrum is structurally determined = zero = free.
