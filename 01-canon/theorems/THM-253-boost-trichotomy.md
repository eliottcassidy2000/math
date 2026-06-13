# THM-253: The Boost Trichotomy — Three Groups of Three

**Status:** PROVED
**Session:** kind-pasteur-2026-03-17-S116n33

## Observation

At n=6 (m=10), the 9 nontrivial Cayley boost values Q(λ_k) = (10-k)/k split into three groups of three:

| Group | k values | Boosts | Prime content | Role |
|-------|----------|--------|---------------|------|
| **I** | 1, 2, 3 | 9, 4, 7/3 | Hurwitz primes in NUMERATOR | PERSISTENT |
| **II** | 4, 5, 6 | 3/2, 1, 2/3 | Only {2,3} (no 7) | TRANSIENT |
| **III** | 7, 8, 9 | 3/7, 1/4, 1/9 | Hurwitz primes in DENOMINATOR | OSCILLATING |

## Theorem (Boost Trichotomy)

The Cayley boost spectrum Q(λ_k) for k = 1, ..., 2s-1 (where m = 2s = 10, s = 5) partitions into three groups:

**Group I** (k = 1, ..., s-2): Q(λ_k) > 1. The Hurwitz primes {2, 3, 7} appear in the numerators.

**Group II** (k = s-1, s, s+1): Q(λ_k) ≈ 1. Only the geometric base primes {2, 3} appear. The forbidden prime 7 is ABSENT.

**Group III** (k = s+2, ..., 2s-1): Q(λ_k) < 1. The Hurwitz primes appear in the denominators. Group III = 1 / Group I (by the functional equation).

### Proof

Q(λ_k) = (m-k)/k = (10-k)/k.

For k ≤ s-2 = 3: (10-k)/k ≥ 7/3 > 1. The numerator 10-k ∈ {9, 8, 7} contains factors of {2, 3, 7}.

For k ∈ {s-1, s, s+1} = {4, 5, 6}: (10-k)/k ∈ {6/4, 5/5, 4/6} = {3/2, 1, 2/3}. Only primes 2 and 3 appear. The prime 7 divides neither numerator nor denominator.

For k ≥ s+2 = 7: (10-k)/k ≤ 3/7 < 1. The denominator k ∈ {7, 8, 9} contains factors of {2, 3, 7}. Each value is the reciprocal of a Group I value: Q(λ_{10-k}) = 1/Q(λ_k). ∎

## The Three-by-Three Matrix

The 9 boosts form a 3×3 matrix indexed by (Hurwitz prime content, group):

```
              Group I (>1)    Group II (≈1)    Group III (<1)
              ----------      -----------      ------------
Prime 3:        9 = 3²         3/2              1/9 = 1/3²
Prime 2:        4 = 2²         2/3              1/4 = 1/2²
Prime 7:       7/3             1                3/7
```

Reading the columns:
- **Column I**: {3², 2², 7/3} — each Hurwitz prime appears SQUARED or as numerator
- **Column II**: {3/2, 2/3, 1} — the three ratios of {2, 3}, centered at 1
- **Column III**: {1/3², 1/2², 3/7} — reciprocals of Column I

Reading the rows:
- **Row "3"**: {9, 3/2, 1/9} — the prime 3 descends from 3² through 3/2 to 1/3²
- **Row "2"**: {4, 2/3, 1/4} — the prime 2 descends from 2² through 2/3 to 1/2²
- **Row "7"**: {7/3, 1, 3/7} — the prime 7 appears ONLY at the extremes, with 1 at center

## Interpretation: The Grand Trichotomy in Spectral Form

The three groups ARE the three types:

| Group | Eigenvalues | Decay | Boost | Trichotomy |
|-------|-------------|-------|-------|------------|
| **I** | λ near 1 | Slow (persistent) | > 1 (amplified) | **INERT** |
| **II** | λ near 0 | Fast (transient) | ≈ 1 (balanced) | **RAMIFIED** |
| **III** | λ near -1 | Slow (oscillating) | < 1 (damped) | **SPLIT** |

- **INERT modes** (Group I): The tournament REMEMBERS these longest. They correspond to individual arc orientations (degree 1-2). Hurwitz primes dominate the boost. These are the "bass notes" that ring longest.

- **RAMIFIED modes** (Group II): The tournament FORGETS these fastest. They correspond to the fine structure (degree 4-5). Only the geometric base {2,3} appears — the forbidden prime 7 is absent. These are the "high frequencies" that die quickly. The center k=5 (λ=0, boost=1) is the **fixed point** where the eigenvalue vanishes.

- **SPLIT modes** (Group III): The tournament OSCILLATES at these frequencies. They have the same decay rates as Group I (same |λ|) but alternate sign. Their boosts are reciprocals of Group I. These are the "out-of-phase" modes that create the overshoot phenomenon.

## The Absent Prime

The prime 7 is ABSENT from Group II. It appears only in Groups I and III — at the extremes of the spectrum. This mirrors the tournament structure:

- 7 is the **forbidden** prime (H=7 impossible)
- 7 appears in the **persistent** and **oscillating** modes but NOT in the **transient** modes
- The transient modes (middle of the spectrum) are governed by {2, 3} alone — the geometric base
- The forbidden prime controls the LONG-TERM behavior (persistence and oscillation) but NOT the short-term (transience)

## Corollary: Product Structure

The 9 boosts satisfy:

**Column I × Column III = 1** (entry-wise): 9 · (1/9) = 4 · (1/4) = (7/3) · (3/7) = 1.

**Row product** (across all three columns):
- Row "3": 9 · (3/2) · (1/9) = 3/2
- Row "2": 4 · (2/3) · (1/4) = 2/3
- Row "7": (7/3) · 1 · (3/7) = 1

The row products ARE Column II! The middle column is the GEOMETRIC MEAN of the extremes.

**Column product** (down all three rows):
- Column I: 9 · 4 · (7/3) = 84 = 2 · 42 = Hurwitz bound coefficient!
- Column II: (3/2) · (2/3) · 1 = 1
- Column III: (1/9) · (1/4) · (3/7) = 1/84

**Column I product = 84 = 2 × 42.** The product of the persistent boosts IS the Hurwitz automorphism bound.

## The Matrix Determinant

det(M) where M = [[9, 3/2, 1/9], [4, 2/3, 1/4], [7/3, 1, 3/7]]

= 9(2/3 · 3/7 - 1/4 · 1) - 3/2(4 · 3/7 - 1/4 · 7/3) + 1/9(4 · 1 - 2/3 · 7/3)

= 9(2/7 - 1/4) - 3/2(12/7 - 7/12) + 1/9(4 - 14/9)

= 9(8/28 - 7/28) - 3/2(144/84 - 49/84) + 1/9(36/9 - 14/9)

= 9(1/28) - 3/2(95/84) + 1/9(22/9)

= 9/28 - 285/168 + 22/81

= 9/28 - 95/56 + 22/81

This is a nonzero rational number — the matrix is INVERTIBLE. The three groups are genuinely independent.
