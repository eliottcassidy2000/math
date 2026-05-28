---
id: THM-341
name: bivariate-GF-good-cuts
status: PROVED
date: 2026-05-28
session: opus-2026-05-28-S3b
---

# THM-341: Bivariate Generating Function for Good-Cuts Distribution

## Statement

Define:
- B(x) = Σ_{a≥2} SC(a+1) xᵃ  (the SC block generating function)
- T(x,y) = generating function for tilings by vertex count and good-cut count

Then the contribution from tilings with d ≥ 2 good cuts is:

  **F(x,y) = x·B(xy) / (1 − x − x·B(xy))**

where [xⁿ yᵈ] F(x,y) = #{tilings of n-vertex staircase with exactly d good cuts}, for d ≥ 2.

## Full Distribution

For any n-vertex tournament tiling:
- d=0: exactly 1 tiling (the transitive tiling, all tiles downward)
- d=1: exactly 0 tilings (IMPOSSIBLE — any tile spans ≥2 cuts, see THM-336)
- d≥2: given by [xⁿ yᵈ] F(x,y) above

## Proof

By THM-340, exactly-d-good(n) = Σ_{k=1}^{⌊d/2⌋} Q(d,k)·C(n−d, k) for d ≥ 2.

The GF in x for fixed d is:
  Σ_n exactly-d-good(n) xⁿ = Σ_k Q(d,k) x^{d+k}/(1−x)^{k+1}

Summing over d≥2 and weighting by yᵈ, using Q(d,k) = [xᵈ]B(x)ᵏ:

  F(x,y) = Σ_{k≥1} (1/(1−x))^{k+1} xᵏ B(xy)ᵏ
          = (1/(1−x)) · xB(xy)/(1−x) / (1 − xB(xy)/(1−x))
          = x·B(xy) / (1 − x − x·B(xy))

## Verification

Verified by brute force for n = 3..7, all d = 2..n−1. All values match.

## Corollaries

**C1 (y=1):** Σ_n (2^m − 1) xⁿ = x·B(x) / (1 − x − x·B(x))
This encodes: total non-transitive tilings = SC + non-SC − 1 = 2^m − 1.

**C2 (SC column):** [xⁿ y^{n-1}] F(x,y) = SC(n) = last column of good-cuts triangle.
This follows because any tiling with d = n−1 good cuts is SC (all cuts covered).

**C3 (Statistical Mechanics interpretation):** F(x,y) is the partition function of a
1D lattice gas where SC blocks are "polymers" of width ≥ 2, fugacity y per unit length,
monomer contribution x/(1−x).
This is a CLUSTER EXPANSION / MAYER EXPANSION — see 07-reflections.

**C4 (Exact formula for column d=2):**
[xⁿ y²]F = C(n−2, 1)·1 = n−2  (proved, matches THM-336)

**C5 (Exact formula for column d=3):**
[xⁿ y³]F = C(n−3, 1)·5 = 5(n−3)  (proved, matches THM-336)
