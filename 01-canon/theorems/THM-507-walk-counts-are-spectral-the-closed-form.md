---
id: THM-507
title: The walk counts of a tournament are spectral — the exact closed form 1ᵀadj(xI−A)1 = (−1)ⁿcharA(−x−1) − charA(x), and hence the entire A-affine pencil determinant is an explicit function of charA
status: PROVED (elementary; matrix-determinant lemma + the tournament identity A−J=−(Aᵀ+I)). VERIFIED with exact integer arithmetic: per-tournament polynomial identity exhaustive n≤6 (all labeled), random n=7,8,9,10,12 (0 failures); walk-count constancy across all 28 cospectral classes at n=6 (0 splits).
source: monad-explorer-2026-06-15-S8
depends_on:
  - THM-506   # permanental companion; the det/per face-by-face program
  - THM-505   # OCF non-spectral defect; spectral skeleton + carriers
related:
  - HYP-2517  # the whole A-affine pencil is spectral (this theorem PROVES its linchpin)
  - HYP-2518  # (this session) the walk-count closed form
  - OPEN-Q-096
  - reflection: the-skew-determinant-is-the-signed-even-face-and-it-is-spectral-monad-s7
  - reflection: the-walk-function-is-the-complement-shift-monad-s8
---

# THM-507 — the walk counts of a tournament are spectral, with an exact closed form

**One line.** For any tournament `T` with 0/1 adjacency `A` (so `A + Aᵀ = J − I`), the
numerator of the walk-generating function is an explicit function of the characteristic
polynomial:

> **`N(x) := 1ᵀ adj(xI − A) 1 = (−1)ⁿ charA(−x−1) − charA(x) = char_{A−J}(x) − charA(x)`**,

equivalently the walk-generating function is the ratio of two characteristic polynomials,

> **`F(x) := Σ_{k≥0} w_k x^{−k−1} = 1ᵀ(xI−A)⁻¹1 = ∏ᵢ(x+1+λᵢ)/∏ᵢ(x−λᵢ) − 1`**,   `w_k = 1ᵀAᵏ1`,

with `{λᵢ}` the eigenvalues of `A`. Hence **all walk counts `w_k` are spectral** (functions of
`charA` alone), for every `k` and every `n`. This PROVES the linchpin left open by S7
(HYP-2517 / the skew-determinant reflection): since every A-affine signing's determinant
reduces to charA times this walk function, **the whole A-affine pencil `P(α,β,γ)=αA+β(J−I)+γI`
has a spectral characteristic polynomial** — upgraded from VERIFIED to PROVED.

---

## Statement

Let `A` be the adjacency matrix of a tournament on `n` vertices, `J = 11ᵀ`, and write
`charA(x) = det(xI − A) = ∏ᵢ(x − λᵢ)`. Define the *walk counts* `w_k = 1ᵀAᵏ1` (= number of
directed walks of length `k`) and the walk-generating function `F(x) = Σ_{k≥0} w_k x^{−k−1}`.

**(a)** `N(x) := 1ᵀ adj(xI − A) 1 = det(xI − A + J) − charA(x)`.

**(b)** `det(xI − A + J) = det((x+1)I + A) = ∏ᵢ(x+1+λᵢ) = (−1)ⁿ charA(−x−1) = char_{A−J}(x)`.

**(c)** Therefore `N(x) = (−1)ⁿ charA(−x−1) − charA(x)` and
`F(x) = N(x)/charA(x) = ∏ᵢ(x+1+λᵢ)/∏ᵢ(x−λᵢ) − 1`.
In particular each `w_k` is a fixed polynomial in the coefficients of `charA`.

**(d)** *(reciprocity)* `(1 + F(x))(1 + F(−x−1)) = 1`.

## Proof

**(a)** is the matrix-determinant lemma in its adjugate (inverse-free) form:
`det(M + ab ᵀ) = det M + bᵀ adj(M) a`. Take `M = xI − A`, `a = b = 1`. ∎

**(b)** Here is the only place the tournament structure enters. From `A + Aᵀ = J − I` we get
`J − A = Aᵀ + I`, so
`xI − A + J = xI + (J − A) = (x+1)I + Aᵀ`.
Determinant is transpose-invariant, so `det(xI − A + J) = det((x+1)I + A)`. The eigenvalues of
`(x+1)I + A` are `x+1+λᵢ`, giving `∏ᵢ(x+1+λᵢ)`. Finally
`∏ᵢ(x+1+λᵢ) = (−1)ⁿ∏ᵢ(−x−1−λᵢ) = (−1)ⁿ∏ᵢ((−x−1)−λᵢ) = (−1)ⁿ charA(−x−1)`, and the same matrix
`(x+1)I + Aᵀ = xI − (A − J)` shows it equals `char_{A−J}(x)`. ∎

**(c)** Combine (a),(b): `charA(x) + N(x) = (−1)ⁿ charA(−x−1)`, i.e.
`N(x) = (−1)ⁿ charA(−x−1) − charA(x)`. The Neumann series
`(xI−A)⁻¹ = Σ_{k≥0} Aᵏ x^{−k−1}` gives `F(x) = 1ᵀ(xI−A)⁻¹1 = Σ_k w_k x^{−k−1}`, while
`(xI−A)⁻¹ = adj(xI−A)/charA(x)` gives `F(x) = N(x)/charA(x)`. Both numerator and denominator
are functions of `charA`, so every Taylor coefficient `w_k` is a function of `charA`. ∎

**(d)** Put `G(x) = 1 + F(x) = (charA(x)+N(x))/charA(x) = (−1)ⁿ charA(−x−1)/charA(x)`. Then
`G(x)G(−x−1) = [(−1)ⁿ charA(−x−1)/charA(x)]·[(−1)ⁿ charA(x)/charA(−x−1)] = 1`. ∎

## The mechanism — why tournaments escape the cospectral walk obstruction

For a *general* graph the walk counts are **not** spectral: the smallest cospectral pair
`C₄ ⊔ K₁` and `K_{1,4}` share the spectrum `{2,0,0,0,−2}` but have `w₂ = Σ degᵢ² = 16` vs `20`.
The walk counts of a general matrix depend on the *main eigenvalues* and the angles of `1`
against the eigenspaces — data not in the spectrum. The all-ones rank-1 perturbation
`A ↦ A − J` generically scrambles the spectrum in an angle-dependent way (this is the entire
"main-eigenvalue / walk-matrix" theory).

The tournament miracle is step (b): `A − J = −(Aᵀ + I)`. So for a tournament the all-ones
perturbation is *not* angle-dependent at all — it is a **transpose-and-shift**, whose
eigenvalues `{−1−λᵢ}` are forced by the spectrum alone. The cospectral obstruction never gets
to act. The single identity `A + Aᵀ = J − I` (complement = converse) is exactly what is
needed and used.

## Corollaries

1. **The whole A-affine pencil is spectral (PROVED).** For `P(α,β,γ) = αA + β(J−I) + γI`, with
   `y = x + β − γ`, the matrix-determinant lemma gives
   `det(xI − P) = det(yI − αA − βJ) = det(yI − αA)·(1 − β·1ᵀ(yI−αA)⁻¹1)`.
   `det(yI − αA) = αⁿ charA(y/α)` is spectral, and `1ᵀ(yI−αA)⁻¹1 = α⁻¹ F_α(y)` where
   `F_α` is the walk function with `A→αA` (eigenvalues `αλᵢ`), spectral by (c). Hence
   `det(xI−P)` is a function of `charA`. This contains: `A` itself `(1,0,0)`; the **skew**
   `S = 2A−(J−I)` `(2,−1,0)` (THM-507 ⇒ the S7 skew char poly `det(xI−S)=∏(x²+μⱼ²)` is
   spectral, now PROVED, not just verified); and every Hermitian length-mod-`r` signing
   `M_ω = ωA + ω̄Aᵀ`. **No determinant of an A-affine matrix can see non-spectral content.**

2. **An exact spectral formula for the walk counts.** Writing `charA(x)=Σ_k(−1)ᵏeₖ x^{n−k}`
   (`eₖ` = elementary symmetric in the `λᵢ`; `e₁ = tr A = 0`), `N(x)=Σ_k eₖ(x+1)^{n−k} −
   Σ_k(−1)ᵏeₖx^{n−k}`. The leading term gives `w₀ = n`; the next, using `e₁=0`,
   `w₁ = C(n,2)`; and `w₂ = (n−1)C(n,2) − Σsᵢ²` collapses to a function of `tr A³` by Moon's
   identity `Σsᵢ² = 2C(n,3) − 2c₃ + C(n,2)`, `c₃ = tr(A³)/3` — recovering the S7 low-order
   computations as special cases of one formula.

3. **The striking refinement is explained.** The score sequence and the power-moments `Σsᵢᵖ`
   (`p ≥ 3`) are *not* spectral (they split 14 cospectral classes at n=6, 85 at n=7), yet the
   walk counts `w_k` are. Reason: `w_k` is the *symmetric* contraction `1ᵀAᵏ1`, which (c)
   expresses through `charA` via the complement shift; the score moments are the *diagonal*
   data `Σ(Aᵏ1)ᵢᵖ`, which retain the angle information the contraction discards.

## What this resolves / does not

- **Resolves:** HYP-2517's open mechanism ("a clean general proof that `w_k` is spectral")
  and the S7 handoff item #1. The entire signed/determinantal world (the A-affine pencil) is
  now PROVABLY spectral, not merely verified to n=7.
- **Does not touch the non-spectral side.** `H = I(Ω,2)`, the unsigned/permanental faces
  (THM-506), and `I(Ω_even,·)` remain non-spectral. THM-507 sharpens *which* objects are on the
  easy (spectral, determinant, P) shore; the Valiant det/per wall is unchanged. The
  odd-length face still has no determinantal home.

## Verification

`04-computation/walk_counts_spectral_proof_monad_s8.py` (+ `.out`), exact integer arithmetic:
the per-tournament identity `N(x) = (−1)ⁿ charA(−x−1) − charA(x)`, the spectral collapse
`char_{A−J}(x) = (−1)ⁿ charA(−x−1)`, and the series match `series(N/charA)_k = w_k = 1ᵀAᵏ1`
(all `k ≤ 2n`) hold with **0 failures** over: exhaustive n=3,4,5,6 (8/64/1024/32768 labeled
tournaments); random n=7,8,9,10,12 (400 each). At n=6 all 28 cospectral classes have a single
walk-vector (0 splits) — the constancy S7 observed, now proved by the per-tournament identity.

## Provenance

This is the tournament case of the *generalized-spectral / walk-matrix* circle of ideas
(complement = converse for tournaments ⟹ ordinary spectrum determines the generalized
spectrum; Wang–Xu walk-matrix theory; *Spectral characterizations of tournaments*,
ScienceDirect 2022). The contribution here is the **explicit closed form**
`F = ∏(x+1+λ)/∏(x−λ) − 1` via the matrix-determinant lemma, the **`A−J = −(Aᵀ+I)`
no-angle-dependence mechanism**, and the placement inside the Φ / det-per / spectral-defect
program — which turns a spectral-DS fact into the PROOF that closes the A-affine-pencil
question (THM-506 / HYP-2517).
