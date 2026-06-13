---
id: HYP-1582
name: Metallic Mean M_n - 1 as Constraint Propagation Threshold
status: EXPLORATORY
source: opus-2026-03-15-S72d
related: [THM-223, THM-217, THM-095, HYP-1581]
---

# HYP-1582: Metallic Mean and Tournament Constraint Propagation

## The Expression

x_n = (n-2 + √(n² + 4)) / 2 = M_n - 1

where M_n = (n + √(n² + 4)) / 2 is the n-th metallic mean (M₁ = φ golden, M₂ = 1+√2 silver, etc.).

## Key Algebraic Properties

- Satisfies x² = (n-2)x + n, equivalently x(x - (n-2)) = n
- Dominant eigenvalue of [[n-2, n], [1, 0]]
- Governs the recurrence a_{k+1} = (n-2)·a_k + n·a_{k-1}
- For large n: x_n ≈ n - 1 + 1/(n-1)
- Product of roots: x₊ · x₋ = -n; sum: x₊ + x₋ = n-2

## Established Connections to Tournaments

### 1. THM-217 Factorization (UNIQUE at n=6)

The THM-217 transfer matrix char poly λ³ - λ² - xλ - x factors as
(λ - (3-n))(λ² - (n-2)λ - n) at x = n(n-3).

This factorization is consistent ONLY at n=6, where x = 18 and the
eigenvalues are {-3, 2+√10, 2-√10}, with 2+√10 = x₆ = M₆ - 1.

### 2. Constraint Independence Fraction: 3/n

For the transitive tournament on n vertices:
- Constraint matrix C has C(n,3) rows (one per transitive triple)
- rank(C) = C(n,2) - n + 1 = (n-1)(n-2)/2
- max_rank / C(n,3) = 3/n exactly
- Redundancy ratio = (n-3)/n

### 3. C^T C Spectrum

For the transitive tournament: C^T C = n · P where P is the orthogonal
projection onto the column space of C^T. ALL nonzero eigenvalues equal n.

For ANY tournament (β₁=0 or β₁=1): the top eigenvalue of C^T C is ALWAYS n.
The metallic mean x_n does NOT appear as a direct eigenvalue of C^T C.

### 4. Simplicial Structure

The transitive tournament constraint system IS the simplicial chain complex:
- Constraints = 2-faces of the complete simplex
- Dependencies = 3-faces (4-vertex subsets)
- Each 4-vertex subset gives one dependency among 4 constraint rows
- Dependency dimension = C(n-1, 3) (matching C(n,4) - something)

### 5. β₁ and Rank Deficiency

β₁ = 0 ⟺ rank(C) = (n-1)(n-2)/2
β₁ = 1 ⟺ rank(C) = (n-1)(n-2)/2 - 1

For rank to drop: a constraint row must escape ALL dependency relations.
Each constraint participates in n-3 dependency relations (one per 4-subset).
For the relation to be "broken," at least one of 3 partner constraints must be
absent (cyclic triple). For ALL n-3 relations broken simultaneously, specific
global topology is required.

## Hypothesized Interpretation

The equation x(x - (n-2)) = n describes a BALANCE CONDITION:

- x = effective constraint propagation rate
- x - (n-2) = excess propagation beyond immediate neighborhood
  (n-2 = non-endpoint vertices = "free dimension")
- n = total system size

At the critical rate x_n, propagation exactly balances the system:
- Above x_n: over-determined → rank maximal → β₁ = 0
- Below x_n: under-determined → rank deficient → β₁ = 1

### 6. β₁=1 Count Sequence

Exact counts of β₁=1 tournaments:

| n | count | total | fraction | count/(n-1)! |
|---|-------|-------|----------|-------------|
| 3 | 2 | 8 | 0.2500 | 1 |
| 4 | 24 | 64 | 0.3750 | 4 |
| 5 | 304 | 1024 | 0.2969 | 38/3 ≈ 12.67 |
| 6 | 4800 | 32768 | 0.1465 | 40 |
| 7 | ~97000 | 2097152 | ~0.046 | ~135 (sampled) |

Successive ratios of count/(n-1)!: 4.000, 3.167, 3.158, ~3.36±0.3

The ratios **suggestively approach π ≈ 3.1416** but this is NOT confirmed.
If count(n)/(n-1)! ~ C·π^n, then the fraction f(n) ~ C·π^n·(n-1)!/2^{C(n,2)}
decays super-exponentially to 0.

**NOTE:** The prior session's claim of f(n=6) ≈ 0.297 was WRONG.
The correct exact value is 4800/32768 = 0.1465.

### 7. Simplicial C^T C Identity

For the transitive tournament on n vertices, the constraint Gram matrix
C^T C = n·P where P is the projection onto the constraint column space.
This means ALL nonzero eigenvalues are exactly n (multiplicity = rank(C)).

Proof sketch: Each column of C has exactly n-2 nonzero entries (edge (i,j)
appears in n-2 transitive triples). The constraint rows have 3 entries {+1,-1,-1}.
The identity C^T C = n·P follows from the incidence structure of the complete simplex.

## Status: EXPLORATORY

The metallic mean M_n - 1 was not found to appear directly in any
tournament spectral quantity. The n=6 THM-217 factorization is genuine but
isolated. The constraint propagation interpretation is conceptually
appealing but lacks a rigorous derivation.

The β₁=1 count/(n-1)! ratio sequence (4, 3.17, 3.16, ~3.4) is intriguing
and could converge to π, but more data (exact n=7) is needed.

## Open Questions

1. Does count(n)/(n-1)! truly grow as C·π^n? Need exact n=7 count.
2. Is there a Markov chain on tournament constraints with spectral gap 1/x_n?
3. Does the recurrence a_{k+1} = (n-2)a_k + n·a_{k-1} count any tournament object?
4. Can the 2×2 matrix [[n-2, n], [1, 0]] be identified with a specific
   endomorphism of the constraint chain complex?
5. Prove C^T C = n·P algebraically (not just computationally).

## Files

- `05-knowledge/results/metallic_beta1_decay_S72d.out`
- `05-knowledge/results/metallic_deeper_S72d.out`
- `05-knowledge/results/metallic_final_S72d.out`
- `05-knowledge/results/metallic_propagation_S72d.out`
