---
id: HYP-1593
name: The 3/n Bridge Between THM-217 and THM-224
status: ESTABLISHED (framework verified, some parts conjectural)
source: opus-2026-03-15-S72d
related: [THM-224, THM-217, THM-225, HYP-1582]
---

# HYP-1593: The 3/n Bridge — Unified Framework for Constraint Independence

## The Fraction 3/n

The fraction 3/n = rank(∂₂)/C(n,3) is the **independence ratio** of the
simplicial boundary operator on K_n. It connects:

### THM-224 (Spectral Side)
- C^TC = n·P, where rank(P) = (n-1)(n-2)/2 = C(n,3)·(3/n)
- The uniform eigenvalue n arises from:
  - trace = 3·C(n,3) [each triple has 3 edges]
  - rank = C(n,3)·(3/n) [H₁(K_n) = 0 forces max rank]
  - eigenvalue = trace/rank = n

### THM-217 (Algebraic Side)
- Transfer matrix char poly λ³-λ²-xλ-x factors at x = n(n-3) = n²·(1-3/n)
- Third eigenvalue α = 3-n = -n·(1-3/n) directly encodes the redundancy
- The metallic quadratic λ²-(n-2)λ-n gives the remaining roots

### The Bridge
Both encode the combinatorial identity: **"3 edges per triple, n-2 triples per edge."**

- THM-224 captures this spectrally (uniform eigenvalue n)
- THM-217 captures this algebraically (factorization at x = n²·redundancy)

## Spectral T/R Duality (PROVED)

For any tournament T on n vertices:

**C_T^TC_T + C_R^TC_R = n·P**

Eigenvalues pair as (λᵢ, n-λᵢ) between C_T^TC_T and C_R^TC_R on im(P).

**β₁(T) = 1 ⟺ C_R^TC_R has eigenvalue n on im(P)**
⟺ there exists v ∈ im(P) ∩ ker(C_T) (an edge function in the constraint
column space with zero flux through all transitive triples).

## The Two Quadratic Equations

| Property | Metallic (THM-217) | Self-dual (THM-224) |
|----------|-------------------|-------------------|
| Equation | μ²-(n-2)μ-n=0 | λ²-nλ+n=0 |
| Sum | n-2 | n |
| Product | -n | +n |
| Discriminant | n²+4 | n(n-4) |
| Center | (n-2)/2 | n/2 |
| Difference of discriminants | 4(n+1) | — |

The self-dual equation λ²-nλ+n=0 appears at the regular tournament
(t₃ = c₃), where T/R symmetry makes C_T^TC_T and C_R^TC_R isospectral.

## n=6: The Unique Consistency Point

n=6 is special because:

1. **Redundancy = Independence** (both = C(n,3)/2 = 10)
2. **THM-217 self-consistency** (gap = n-6 = 0)
   - The two constraints on x (from constant and linear terms) agree only at n=6
   - α·c = -(3-n)·n = n(n-3) = x (constant term)
   - α·b - c = (3-n)(n-2) - n = n²-4n+6 = x iff n=6 (linear term)
3. **x = n²/2** — the redundancy capacity equals half the total capacity

## Eigenvalue Patterns at n=5 (Exhaustive)

The spectra of C_T^TC_T are completely classified:

- **Top eigenvalue = n = 5** for ALL tournaments (THM-225)
- **Eigenvalues collapse toward n** as t₃ increases (more transitive triples)
- **Non-top eigenvalue sum** = 3·t₃ - n·(#eigs = n)
- **At t₃ = c₃**: self-dual eigenvalues (5±√5)/2, product = n
- **Integer eigenvalues** {2, 3, 4} appear for β₁=1 tournaments
- **Irrational eigenvalues** involve √5 (t₃=5,6) or √2 (t₃=7)

## C_R^TC_R Spectrum and β₁

β₁=1 ⟺ C_R^TC_R has eigenvalue n on im(P). Verified:

| t₃ | β₁ | C_R^TC_R eigs on im(P) |
|----|----|------------------------|
| 5  | 1  | [5, 3.618, 3.618, 1.382, 1.382, 0] |
| 6  | 0  | [4.618, 3.618, 2.382, 1.382, 0, 0] |
| 6  | 1  | [5, 3, 3, 1, 0, 0] or [5, 3, 2, 2, 0, 0] |
| 7  | 0  | [4.414, 3, 1.586, 0, 0, 0] |
| 7  | 1  | [5, 2, 2, 0, 0, 0] |

## Open Questions

1. **Prove THM-225 for general n** (top eigenvalue = n always)
2. **Does the golden ratio pair generalize?** At n=5 regular: (5±√5)/2.
   At n=7 regular: what are the C_T^TC_T eigenvalues?
3. **Physical meaning of the metallic mean** in constraint language?
4. **Does the count(n)/(n-1)! ratio for β₁=1 tournaments converge to π?**

## Files

- `05-knowledge/results/deep_3n_synthesis_S72d.out`
- `/tmp/deep_3n_synthesis.py`, `/tmp/deep_3n_proof.py`, `/tmp/n6_consistency.py`
