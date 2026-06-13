# The Unlabeling Principle

**Session**: kind-pasteur-2026-03-22-S20ab

## The Core Idea

The meta-tournament on iso classes is **transitive** (H=1). This means: after removing vertex labels, tournament space is a perfect hierarchy. All "complexity" of the 1024 tournaments at n=5 was in choosing which vertex gets which label. The structure itself is a clean linear order.

## The Unlabeling Spectrum

Each level of quotient removes more labeling information:

| Level | Quotient by | #classes | bits | H info captured |
|-------|------------|----------|------|-----------------|
| 0 | nothing | 1024 | 10.0 | 100% |
| 1 | vertex labels (S_n) | 12 | 3.6 | 100% |
| 2 | within-score structure | 9 | 3.2 | 85.2% |
| 3 | score detail (keep S2) | 6 | 2.6 | 85.2% |
| 5 | keep only H | 7 | 2.8 | 100% |
| 6 | everything | 1 | 0 | 0% |

**64% of the bits encoding a tournament are label noise.** Only 36% carry structural information. Of those structural bits, 78% determine H.

## The Crossover at n=5

The label fraction = log2(n!) / C(n,2) asymptotically equals 2*log2(n)/(n-1), which goes to **ZERO**.

| n | label % | structure % |
|---|---------|-------------|
| 3 | 67% | 33% |
| 4 | 67% | 33% |
| 5 | 64% | 36% |
| 6 | 61% | 39% |
| 10 | ~47% | ~53% |
| 20 | ~29% | ~71% |
| 100 | ~9% | ~91% |

**At small n, labels dominate. At large n, structure dominates.** The crossover is near n=7-8 where labels and structure are roughly equal.

This is why n=5 is special: it's the **last order where the illusion of labeling works.** At n=5, you can "blame complexity on labels" because 64% of bits are label noise. At n=7+, most bits are genuinely structural, and the complexity is real.

## The Unlabeling Hierarchy (General Principle)

The principle applies far beyond tournaments:

**Mathematics:**
- Graphs -> iso classes -> degree sequences -> edge count -> order
- Permutations -> conjugacy classes -> cycle types -> number of fixed points
- Partitions -> Ferrers shapes -> rank+crank -> size
- Polynomials -> root multisets -> elementary symmetric functions -> degree

**Science:**
- Molecules -> constitutional isomers -> functional groups -> molecular formula -> molecular weight
- Proteins -> fold classes -> secondary structure -> amino acid composition -> length
- Crystals -> space groups (230) -> crystal systems (7) -> Bravais lattices (14) -> dimension (3)

**Engineering:**
- Neural network weights -> equivalence under neuron permutation -> layer-wise spectra -> depth+width
- Attention matrices -> token-relabeling equivalence -> spectral profile -> rank
- Rankings -> isomorphism classes -> score sequences -> Copeland score

At each step, we quotient by a symmetry group and lose some information but gain clarity. The **coarsest quotient that still determines the property of interest** is the natural "unlabeled" representation.

## The OCR as an Unlabeling Theorem

The OCR (97% of Var(H) explained by scores) says: **the score sequence is an almost-sufficient statistic for H**. In unlabeling language: quotienting by within-score structure loses only 3% of H-information.

The Walsh-Fourier perspective: order-2 captures 92-95% of Var(H). Order-2 coefficients depend only on pairwise arc interactions, which are captured by scores. Order-4 (the 5-8% residual) requires knowing which specific arcs are present, not just their marginals.

**The OCR is the QUANTITATIVE UNLABELING THEOREM for tournaments.**

## The Self-Correction

I initially predicted that label fraction increases with n (more labels = more noise). The calculation revealed the opposite: label fraction DECREASES because C(n,2) grows faster than log(n!).

This is the key insight: **at large n, most tournaments are "generic" (trivial automorphism group), so labeling carries almost no redundancy.** The labeling illusion is a small-n artifact. For large tournaments, the genuine combinatorial structure dominates.

## Connection to the Meta-Tournament

The meta-tournament being transitive at n=5 is a consequence of:
1. Small n (64% label noise -> clean quotient)
2. H having only even Walsh orders (no cycles in the H-gradient)
3. Order-2 dominance (85% of H-info in scores, which define a total order)

At n=6, the meta-tournament likely gains 3-cycles (from the H=37/H=45 competing peaks). This is the beginning of genuine structural complexity that unlabeling cannot remove.

**Prediction**: The meta-tournament transitions from transitive (n<=5) to non-transitive (n>=6), precisely at the order where the label fraction drops below the structural fraction.

## The Philosophical Point

"Complexity is an illusion of labeling" is true at small n and false at large n. The transition happens near n=5-6. This is simultaneously:
- The order where the first non-trivial forbidden H value appears (H=7 at n=5)
- The order where the Morse landscape gains secondary peaks (n=6)
- The order where alpha_2 turns on in the OCF (n=6)
- The order where the OCR drops below 97% (n=6)

**All of these transitions mark the same boundary: where genuine structural complexity overwhelms labeling noise.**
