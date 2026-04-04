# THM-286: All Multilinear Coefficients of H Are Even

**Status:** PROVED (algebraic, for all n)
**Filed by:** opus-2026-04-04-S1
**Depends on:** Rédei's theorem (H(T) is always odd)

## Statement

Let H(t₁,...,tₘ) be the multilinear polynomial expressing the Hamiltonian path count as a function of tile bits t_i ∈ {0,1}. Write:

  H(t) = Σ_{S ⊆ [m]} c_S · Π_{i∈S} tᵢ

Then c_∅ = 1 (odd) and c_S is even for every non-empty S ⊆ [m].

Equivalently: **H(t) ≡ 1 (mod 2) for all t ∈ {0,1}^m.** This is Rédei's theorem restated in the tiling model.

## Proof

By induction on |S|.

**Base case |S|=1:** Let S = {i}. Setting all tiles to 0 except tile i:
  H(eᵢ) = c_∅ + c_{i} = 1 + c_{i}

Since H(eᵢ) is always odd (Rédei), c_{i} = H(eᵢ) - 1 is even. ∎

**Inductive step:** Assume c_T is even for all non-empty T with |T| < |S|.
Setting t = indicator vector of S:
  H(1_S) = Σ_{T ⊆ S} c_T = c_∅ + Σ_{∅ ≠ T ⊆ S} c_T = 1 + Σ_{∅ ≠ T ⊊ S} c_T + c_S

By induction, each c_T with ∅ ≠ T ⊊ S is even. So Σ_{∅ ≠ T ⊊ S} c_T is even.
Since H(1_S) is odd (Rédei), we need 1 + (even) + c_S to be odd, hence c_S is even. ∎

## Additional Computational Findings

The **minimum 2-adic valuation** of c_S (for |S| ≥ 1) is exactly 1 at all tested n (3-7). That is, while every coefficient is divisible by 2, NOT every coefficient is divisible by 4.

The **maximum-degree coefficients** have a different structure at odd vs even n:
- At odd n = 3, 5, 7: max-degree coefficients are ALL ±2
- At even n = 6: max-degree coefficients include ±4

**Conjecture (THM-286b):** At odd n, all maximum-degree multilinear coefficients are ±2.

## Computed Multilinear Spectra

| n | m | Degree | Total nonzero | Per degree (d0:d1:d2:d3:d4:d5:d6) |
|---|---|--------|---------------|--------------------------------------|
| 3 | 1 | 1 | 2 | 1:1 |
| 4 | 3 | 2 | 6 | 1:3:2 |
| 5 | 6 | 4 | 35 | 1:6:13:8:7 |
| 6 | 10 | 4 | 200 | 1:10:42:66:81 |
| 7 | 15 | 6 | 1782 | 1:15:101:309:626:407:323 |

## See Also
- THM-284 (linear coefficient = 2^(skip-1))
- THM-260 (Walsh degree = 2⌊(n-1)/2⌋)
- Rédei's theorem (H(T) always odd)
- Scripts: h_multilinear_complete.py, h_multilinear_analysis.py, max_degree_coeff_proof.py
