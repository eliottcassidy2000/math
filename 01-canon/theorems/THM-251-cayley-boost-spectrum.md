# THM-251: The Cayley Boost Spectrum and Functional Equation

**Status:** PROVED
**Session:** kind-pasteur-2026-03-17-S116n33

## Setup

Let T be a tournament on n vertices with a fixed canonical path 0 → 1 → ... → n-1. The non-path arcs form the tiling space {0,1}^m where m = C(n-1, 2).

The **random flip Markov chain** on the tiling space flips one uniformly random bit per step. Its eigenvalues are:

λ_k = (m - 2k)/m for k = 0, 1, ..., m

with multiplicity C(m, k) at each level.

## Theorem A (Cayley Boost Spectrum)

The Cayley transform Q(x) = (1+x)/(1-x) maps each eigenvalue to:

**Q(λ_k) = (m - k)/k** for k = 1, ..., m-1

with Q(λ_0) = ∞ and Q(λ_m) = 0.

### Proof

Q((m - 2k)/m) = (1 + (m-2k)/m) / (1 - (m-2k)/m) = (2m - 2k)/(2k) = (m-k)/k. ∎

## Theorem B (Functional Equation)

The Cayley boost spectrum satisfies:

**Q(λ_k) · Q(λ_{m-k}) = 1** for all k = 1, ..., m-1

### Proof

Q(λ_k) · Q(λ_{m-k}) = ((m-k)/k) · ((m-(m-k))/(m-k)) = ((m-k)/k) · (k/(m-k)) = 1. ∎

### Interpretation

This is a **duality pairing**: the Cayley boost of mode k and its "mirror" mode m-k are multiplicative inverses. The product Q(λ_k) · Q(λ_{m-k}) = 1 is the tournament analog of a functional equation relating a zeta function at s and 1-s.

## Theorem C (Hurwitz Factorization at n=6)

At n = 6 (m = 10), the Cayley boost values Q(λ_k) = (10-k)/k for k = 1,...,9 are:

| k | λ_k | Q(λ_k) | Factored |
|---|-----|---------|----------|
| 1 | 4/5 | 9 | 3² |
| 2 | 3/5 | 4 | 2² |
| 3 | 2/5 | 7/3 | 7/3 |
| 4 | 1/5 | 3/2 | 3/2 |
| 5 | 0 | 1 | 1 |
| 6 | -1/5 | 2/3 | 2/3 |
| 7 | -2/5 | 3/7 | 3/7 |
| 8 | -3/5 | 1/4 | 1/2² |
| 9 | -4/5 | 1/9 | 1/3² |

**Every Cayley boost is a ratio of powers of Hurwitz primes {2, 3, 7}.**

### Proof

The ratios (10-k)/k for k = 1,...,9 are 9/1, 8/2, 7/3, 6/4, 5/5, 4/6, 3/7, 2/8, 1/9. After reduction: 9, 4, 7/3, 3/2, 1, 2/3, 3/7, 1/4, 1/9. The primes appearing in numerators and denominators are: {2, 3, 7} (from 9=3², 4=2², 7, 3, 2). The prime 5 appears only in the eigenvalue denominators (λ_k has denominator 5) but NOT in the Cayley boosts, because 5 divides m/2 = 5 and cancels: (10-k)/k = (10-k)/k with 5 | 10 contributing to λ_k but not to Q(λ_k). ∎

## Remark

The appearance of {2, 3, 7} in the Cayley boost spectrum is specific to n=6 (m=10). At general n, the Cayley boost spectrum Q(λ_k) = (m-k)/k involves primes dividing values in {1, ..., m-1}, which depend on n.
