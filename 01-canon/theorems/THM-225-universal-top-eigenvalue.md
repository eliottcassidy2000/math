---
id: THM-225
name: Top Eigenvalue of Tournament Constraint Gram Matrix
status: PARTIALLY-TRUE
verified_computationally: n=5 exhaustive, n=6-8 sampled (20000 each, 0 violations)
counterexample: n=9 (~0.1% of tournaments have top eigenvalue < 9)
proved_by: opus-2026-03-15-S72d
related: [THM-224, THM-222]
---

# THM-225: Top Eigenvalue of C_T^TC_T

## Statement (CORRECTED)

**For n ≤ 8:** For any tournament T with t₃ ≥ 1, max eigenvalue of C_T^TC_T = n.

**For n ≥ 9:** NOT universal. Some tournaments have max eigenvalue < n.
The condition for max eigenvalue = n is: rank(C_R) < (n-1)(n-2)/2.
When rank(C_R) = (n-1)(n-2)/2 (the cyclic boundaries span im(C^T)),
the top eigenvalue drops below n.

At n=9: ~0.1% of tournaments violate top eig = n. Top eigenvalue ≈ 8.84-8.94.

**Original (FALSE) claim:** "For any tournament T on n ≥ 3 vertices with t₃ ≥ 1,
the largest eigenvalue of C_T^TC_T is exactly n." This fails at n=9.

Here C_T is the constraint matrix restricted to the transitive triples of T.

## Spectral Duality

From THM-224 (C^TC = n·P):

**C_T^TC_T + C_R^TC_R = n·P**

where C_R is the constraint matrix restricted to cyclic triples.

This gives the **spectral pairing**: if the eigenvalues of C_T^TC_T on im(P)
are λ₁ ≥ λ₂ ≥ ... ≥ λ_r (where r = (n-1)(n-2)/2), then the eigenvalues
of C_R^TC_R on im(P) are n-λ_r ≥ n-λ_{r-1} ≥ ... ≥ n-λ₁.

**Consequence:** β₁(T) = 1 iff C_R^TC_R has eigenvalue n on im(P), i.e.,
there exists an edge function v in the constraint column space with zero
flux through every transitive triple.

## Proof

**Upper bound (PROVED):** Since C_R^TC_R ≥ 0 and C_T^TC_T = n·P - C_R^TC_R,
all eigenvalues of C_T^TC_T are ≤ n.

**Equality criterion (PROVED):** max eigenvalue = n iff ∃ v ∈ im(P) ∩ ker(C_R).
Since im(C_R^T) ⊆ im(C^T) = im(P), this holds iff:

  dim(im(P) ∩ ker(C_R)) = r - rank(C_R) ≥ 1, i.e., **rank(C_R) < r**.

**For n ≤ 8 (PROVED):** rank(C_R) ≤ c₃ (since C_R has c₃ rows in ℝ^E).
By Kendall-Babington Smith, the maximum c₃ over all tournaments on n vertices is:
- Odd n: c₃ = C(n,3) - n·C((n-1)/2, 2) (achieved by regular tournament)
- Even n: c₃ = (n³-4n)/24 (achieved by near-regular tournament)

The key inequality is max(c₃) < r = (n-1)(n-2)/2 for 4 ≤ n ≤ 8:

| n | r | max c₃ | max c₃ < r |
|---|---|--------|------------|
| 3 | 1 | 1*     | (t₃≥1 ⟹ c₃=0) |
| 4 | 3 | 2      | ✓ |
| 5 | 6 | 5      | ✓ |
| 6 | 10| 8      | ✓ |
| 7 | 15| 14     | ✓ |
| 8 | 21| 20     | ✓ |

(*At n=3, the 3-cycle has c₃=1=r, but t₃≥1 forces c₃=0.)

Since rank(C_R) ≤ c₃ < r for all tournaments with t₃ ≥ 1 and n ≤ 8,
we have dim(ker(C_R) ∩ im(C^T)) = r - rank(C_R) ≥ 1.

**Note:** t₃ > c₃ is FALSE in general (regular tournament at n=5 has t₃=c₃=5).
The proof uses c₃ < r, NOT t₃ > c₃.

**For n ≥ 9 (FAILS):** At n=9, max c₃ = 30 > r = 28, and some tournaments
achieve rank(C_R) = r = 28. In these cases, ker(C_R) ∩ im(C^T) = {0},
and the top eigenvalue drops below n.

## Verified Data

### n=5 (exhaustive, all 1024 tournaments):

| t₃ | β₁ | Spectrum on im(P) | Count |
|----|----|-------------------|-------|
| 5  | 1  | [5, φ, φ, φ̄, φ̄, 0] | 24 |
| 6  | 0  | [5, 5, φ, φ-1, φ̄, φ̄-1] | 120 |
| 6  | 1  | [5, 5, 4, 2, 2, 0] | 120 |
| 6  | 1  | [5, 5, 3, 3, 2, 0] | 40 |
| 7  | 0  | [5, 5, 5, 2+√2, 2, 2-√2] | 120 |
| 7  | 1  | [5, 5, 5, 3, 3, 0] | 120 |
| 8  | 0  | [5, 5, 5, 5, 3, 1] | 240 |
| 9  | 0  | [5, 5, 5, 5, 5, 2] | 120 |
| 10 | 0  | [5, 5, 5, 5, 5, 5] | 120 |

where φ = (5+√5)/2 ≈ 3.618, φ̄ = (5-√5)/2 ≈ 1.382.

### n=6 (sampled, 10000 tournaments): ALL have top eigenvalue = 6.

## Key Identities

- **Trace:** trace(C_T^TC_T) = 3·t₃ (each transitive triple contributes 3 to diagonal)
- **Diagonal:** (C_T^TC_T)_{ee} = t₃(e) = number of transitive triples containing edge e
- **t₃(e) + c₃(e) = n-2** for every edge e (from THM-224 diagonal = n-2)

## Self-Dual Eigenvalues

At the regular tournament (n=5, t₃ = c₃ = 5), the T/R spectral duality
becomes self-duality: C_T^TC_T and C_R^TC_R have the SAME spectrum.
The non-top eigenvalues satisfy λ² - nλ + n = 0, giving the golden ratio
pair (5±√5)/2 with product = n and sum = n.

## Connection to THM-217

The self-dual equation λ² - nλ + n = 0 (from THM-224) and the metallic
equation μ² - (n-2)μ - n = 0 (from THM-217 factorization at x=n(n-3))
are related by:
- Discriminant difference: Δ_metallic - Δ_self-dual = (n²+4) - n(n-4) = 4(n+1)
- Both converge as n → ∞ (centered at n/2 with half-width n/2)
- Both encode the same 3/n independence fraction in different algebraic forms

## Open Questions

1. ~~Prove top eigenvalue = n for ALL n~~ RESOLVED: false at n ≥ 9 (MISTAKE-018)
2. Characterize which tournaments at n ≥ 9 have top eigenvalue < n (rank(C_R) = r)
3. Characterize which tournaments achieve specific spectral patterns
4. Does the number of distinct spectra at each t₃ follow a pattern?
5. What is the exact fraction of n=9 tournaments with top eigenvalue < n?

## Files

- `05-knowledge/results/deep_3n_synthesis_S72d.out`
- `/tmp/deep_3n_synthesis.py`
- `/tmp/deep_3n_proof.py`
