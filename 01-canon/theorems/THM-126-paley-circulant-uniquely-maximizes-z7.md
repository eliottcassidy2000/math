---
theorem_id: THM-126
title: Paley tournament uniquely maximizes H among circulants on Z_7
status: PROVED (computational — exhaustive over all 8 circulants on Z_7)
proved_by: opus-2026-03-12
date: 2026-03-12
related_hypotheses: [HYP-400, HYP-437]
related_theorems: [THM-125]
tags: [paley, circulant, maximizer, spectral, hamiltonian-paths]
---

## Statement

Among all 8 circulant tournaments on Z_7, the Paley tournament T_7 (connection set
S = QR_7 = {1,2,4} and its complement {3,5,6}) is the **unique** maximizer of H(T).

**Exact values:**
| Connection set S | H(T) | Spectral spread | Paley? |
|---|---|---|---|
| {1,2,4} | **189** | 0.0000 | YES |
| {3,5,6} | **189** | 0.0000 | YES |
| {1,2,3} | 175 | 1.6920 | no |
| {1,3,5} | 175 | 1.6920 | no |
| {1,4,5} | 175 | 1.6920 | no |
| {2,3,6} | 175 | 1.6920 | no |
| {2,4,6} | 175 | 1.6920 | no |
| {4,5,6} | 175 | 1.6920 | no |

(Note: {3,5,6} = complement of {1,2,4} = T_7^{op}; since T_7 ≅ T_7^{op} for p≡3 mod 4,
both give H=189.)

## Proof

Exhaustive computation via Held-Karp DP over all 8 = 2^{(p-1)/2} = 2^3 connection sets.
Each H value computed exactly. Script: `04-computation/paley_maximizer_circulant_test.py`.
Output: `05-knowledge/results/paley_maximizer_circulant_test.out`.

## Spectral Flatness Theorem (Z_7 case)

**Theorem:** Among circulant tournaments on Z_7, flat eigenvalue spectrum ↔ maximum H.

- Paley: all non-trivial eigenvalues |λ_k| = √2 for k=1,...,6 (spread = 0)
- All non-Paley: |λ_k| ∈ {2.2470, 0.8019, 0.5550} (spread = 1.6920)
- The flat-spectrum condition is EQUIVALENT to being Paley (or its complement) at p=7.

**Gauss sum computation:** For T_7 with S = QR_7 = {1,2,4}:
```
λ_k = sum_{s in QR_7} ω^(ks)   where ω = exp(2πi/7)
```
By the theory of Gauss sums for p ≡ 3 mod 4: |λ_k| = √((p+1)/4) = √2 for all k≠0.
(Correction to earlier conjecture: eigenvalue magnitude is √2, not √p/2 = √7/2 ≈ 1.3229.)

Proof: For p≡3 mod 4, the non-trivial Gauss sum satisfies |g|² = p.
The eigenvalue at k≠0 is λ_k = (η(k)·g - 1)/2 where η(k) ∈ {±1} is the Legendre symbol.
Then |λ_k|² = |η(k)·g - 1|²/4 = (p + 1)/4 = (7+1)/4 = 2. So |λ_k| = √2. ✓

## Significance

1. **Circulant reduction**: If the Paley maximizer conjecture can be reduced to circulants,
   this establishes the n=7 case.
2. **Spectral principle**: Flat eigenvalue spectrum = Ramanujan-optimal = H-optimal (for circulants).
3. **n=5 note**: At n=5, ALL 4 circulant tournaments achieve H=15 = OEIS max. The maximizer
   at n=5 is NOT unique. Paley does not apply (5≡1 mod 4, QR_5 not a tournament).

## Open Questions

- Does the spectral flatness ↔ H-maximization equivalence hold for ALL prime p≡3 mod 4?
- Is the max-H tournament among ALL (non-circulant) n-vertex tournaments always Paley for p prime?
- What is H for the Satake NDRT at q=13? (See INV-137.)
