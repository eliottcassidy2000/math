# THM-252: The Rapidity Lattice of the Flip Chain

**Status:** PROVED
**Session:** kind-pasteur-2026-03-17-S116n33

## Statement

For the random flip Markov chain on the n-vertex canonical-path tiling space, the set of eigenvalue rapidities

R = {arctanh(λ_k) : 1 ≤ k ≤ m-1, λ_k ≠ 0}

is contained in the Q-vector space:

**R ⊆ Q · ln(p_1) + Q · ln(p_2) + ... + Q · ln(p_r)**

where {p_1, ..., p_r} are the prime factors of (m-1)! (equivalently, all primes ≤ m-1).

More precisely, the rapidity lattice is spanned by the logarithms of primes dividing the numbers {1, 2, ..., m-1}.

## Proof

arctanh(λ_k) = (1/2) · ln(Q(λ_k)) = (1/2) · ln((m-k)/k) by THM-251.

Since 1 ≤ k ≤ m-1 and λ_k ≠ 0 (i.e., k ≠ m/2 when m is even), the ratio (m-k)/k is a positive rational number whose numerator and denominator are both in {1, ..., m-1}.

Therefore:

ln((m-k)/k) = ln(m-k) - ln(k)

Each of ln(m-k) and ln(k) decomposes as a sum of ln(p) over the prime factorization of m-k and k respectively. All prime factors are ≤ m-1. ∎

## Specialization to n=6 (m=10)

At n=6, the nonzero-eigenvalue rapidities are:

| k | λ_k | (m-k)/k | arctanh(λ_k) | Expression in {ln 2, ln 3, ln 7} |
|---|-----|---------|---------------|-----------------------------------|
| 1 | 4/5 | 9 = 3² | ln(3) | **1 · ln(3)** |
| 2 | 3/5 | 4 = 2² | ln(2) | **1 · ln(2)** |
| 3 | 2/5 | 7/3 | (ln 7 - ln 3)/2 | **(-1/2) · ln(3) + (1/2) · ln(7)** |
| 4 | 1/5 | 3/2 | (ln 3 - ln 2)/2 | **(-1/2) · ln(2) + (1/2) · ln(3)** |

The rapidity lattice at n=6 is: **Q · ln(2) ⊕ Q · ln(3) ⊕ Q · ln(7)**

## Theorem (Rank and Independence)

By Baker's theorem (1966), the numbers ln(2), ln(3), ln(7) are **Q-linearly independent** (since 2, 3, 7 are multiplicatively independent integers). Therefore the rapidity lattice at n=6 is a free Q-module of **rank 3**.

## Corollary (Adelic Structure)

The rank-3 rapidity lattice Q · ln(2) ⊕ Q · ln(3) ⊕ Q · ln(7) is the **archimedean shadow** of the cuboid Z/2Z × Z/3Z × Z/7Z = Z/42Z.

Both are rank-3 structures built from {2, 3, 7}:
- **Cuboid** (p-adic completion): finite, periodic, lives in ∏ Z_p
- **Rapidity lattice** (archimedean completion): infinite, aperiodic, lives in R

This parallel is an instance of the **adelic product formula**: the local (p-adic) and global (archimedean) invariants of a number field are related by the product formula ∏_v |x|_v = 1. Here, the "number field" is Q(2,3,7) and the two completions encode complementary aspects of tournament structure.

## Remark (Conjugate Variables)

The eigenvalue λ_k = (5-k)/5 is **rational** (lives in Q).
Its arctanh = (1/2)ln((10-k)/k) is **transcendental** (lives in R \ Q̄, by Hermite-Lindemann).

The map λ ↦ arctanh(λ) is the **formal group logarithm** of the Cayley formal group F(x,y) = (x+y)/(1+xy). It crosses the boundary between the algebraic (eigenvalues) and transcendental (rapidities) worlds.

This crossing is IRREVERSIBLE in the following sense: the inverse map tanh: R → (-1,1) sends transcendental rapidities back to rational eigenvalues, but the information about WHICH transcendental value was the input is lost (since tanh(transcendental) = rational, there are uncountably many transcendental inputs mapping to each rational output).
