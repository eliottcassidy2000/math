# THM-262: Seesaw Identity for Tournament Path Homology

**Status:** PROVED (algebraic identity). The seesaw conjecture β₁·β₃=0 is REDUCED to a rank bound but NOT fully proved.
**Filed by:** kind-pasteur-2026-03-25-S28
**Dependencies:** THM-108 (β₂=0 for all tournaments)
**Devil's advocate review (S29):** Identity derivation is CORRECT. Claims about "completing the proof" were premature — corrected below.

## Statement

For any tournament T on n ≥ 3 vertices:

**β₁(T) + β₃(T) = S(T) - dim(im ∂₄)**

where:
- S(T) = dim(ker ∂₁) + dim(Ω₃) - dim(Ω₂)
- dim(ker ∂₁) = C(n,2) - n + 1 (constant for all tournaments on n vertices)
- dim(Ω₂), dim(Ω₃) depend on the specific tournament T
- im ∂₄ = image of the boundary map ∂₄: Ω₄ → Ω₃

## Proof

From the GLMY chain complex:
```
Ω₄ →^{∂₄} Ω₃ →^{∂₃} Ω₂ →^{∂₂} Ω₁ →^{∂₁} Ω₀
```

**Step 1.** β₂ = 0 (THM-108) gives exactness at Ω₂: ker(∂₂) = im(∂₃).
Therefore: dim(Ω₂) - rank(∂₂) = rank(∂₃).

**Step 2.** β₁ = dim(ker ∂₁) - rank(∂₂), so rank(∂₂) = dim(ker ∂₁) - β₁.

**Step 3.** From Step 1: rank(∂₃) = dim(Ω₂) - rank(∂₂) = dim(Ω₂) - dim(ker ∂₁) + β₁.

**Step 4.** ker(∂₃) = dim(Ω₃) - rank(∂₃) = dim(Ω₃) - dim(Ω₂) + dim(ker ∂₁) - β₁.

**Step 5.** β₃ = ker(∂₃) - dim(im ∂₄) = dim(Ω₃) - dim(Ω₂) + dim(ker ∂₁) - β₁ - dim(im ∂₄).

**Step 6.** Rearranging: β₁ + β₃ = dim(ker ∂₁) + dim(Ω₃) - dim(Ω₂) - dim(im ∂₄) = S - dim(im ∂₄). QED.

**Note:** rank(∂₂) cancels in Steps 2-5 BECAUSE β₂ = 0. This is the essential role of THM-108.

## What This Proves and What It Doesn't

**PROVED:** The identity β₁ + β₃ = S - dim(im ∂₄) for all tournaments.

**NOT PROVED (the seesaw conjecture):** β₁·β₃ = 0.

The seesaw would follow from: **when β₁ = 1, dim(im ∂₄) ≥ S - 1.**
This is equivalent to: rank(∂₄: Ω₄ → Ω₃) ≥ dim(Ω₃) - dim(Ω₂) + C(n,2) - n.

This rank bound is verified computationally (300 samples at n=7, 3000 at n=8,
zero violations) but NOT proved algebraically. The bound says the boundary from
4-chains fills "almost all" of Ω₃, which is a RIGIDITY property of tournament
path complexes forced by tournament completeness.

## Computational Verification

| n | Samples | β₁+β₃ ∈ {0,1} | S range | Identity holds |
|---|---------|----------------|---------|----------------|
| 7 | 300 | ALL | [14, 35] | ALL |
| 8 | 3000 | ALL | — | ALL (p=0.0017 for seesaw) |

At n=7: ker(∂₃) - im(∂₄) = 0 for 273 tournaments (β₃=0), = 1 for 27 (β₃=1). Exact.
At n=8: β₃=2 occurs (1/3000). β₁+β₃ = 2 is allowed when β₁=0 (seesaw holds).

## Proof Strategy for the Full Seesaw (unfinished)

Via LES of (T, T\v) + strong induction:
1. Base: β₁·β₃ = 0 at n ≤ 6 (exhaustive verification).
2. If β₁(T) = 1, find v with β₁(T\v) = 1. By induction β₃(T\v) = 0.
3. LES gives: 0 → H₃(T) → H₃(T, T\v) → 0 (since H₃(T\v) = 0, H₂(T\v) = 0).
4. So β₃(T) = dim H₃(T, T\v).
5. **Missing piece:** Prove H₃(T, T\v) = 0 when β₁(T) = 1.

Step 5 is verified at n ≤ 7 (all β₁=1 tournaments have H₃(T,T\v) = 0 for all v).

## Related

- THM-108: β₂ = 0 for all tournaments (the prerequisite)
- THM-095: β₁·β₃ = 0 conditional on β₂ = 0 (the seesaw conjecture)
- OPEN-Q-024: Even Betti vanishing and seesaw questions
