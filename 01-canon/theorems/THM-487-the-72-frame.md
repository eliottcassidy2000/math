# THM-487: the [72,36,16] frame — its extremal enumerator is positive, its lattice exists, so the obstruction is purely code-combinatorial (and the Paley gauge route stalls at d=12)

**Status:** VERIFIED (the enumerator facts, exact) + LITERATURE-PINNED (SOTA, with
DOIs from the recon round) + framing (the obstruction localization). No new
theorem about existence is claimed — this is the precise statement of WHERE the
problem lives, sharpening OPEN-Q-061/HYP-2415.
**Source:** kind-pasteur-2026-06-11-S3 (HYP-2421). Uses THM-486 (the η-bridge and
n = 3696), THM-481 (the Paley/eQR gauge), THM-484/claudebox (the 72-as-3·Δ gate).

## A. The enumerator is not the obstruction (VERIFIED)

The extremal Type II enumerator W₇₂ (the unique Gleason-invariant degree-72
enumerator with A₄ = … = A₁₂ = 0, minimum weight 16) has
  A₀ = 1, **A₁₆ = 249849** (Sloane's value), A₂₀ = 18106704, …, all coefficients
  POSITIVE (built exactly from W_{ê₈}, P₂₄; full vector in the .out).
By THM-486 C the first negative coefficient in the extremal Type II family appears
only at **n = 3696** (Zhang's exact threshold, reproduced). Length 72 = m·24 with
m = 3 sits FAR below — so:

> The obstruction to a [72,36,16] code is NOT weight-enumerator positivity, NOT
> the Mallows–Sloane bound (it is the bound, and attainable as a polynomial), and
> NOT Rains' shadow bound (72 ≡ 0 mod 24, no shadow penalty). The extremal W₇₂ is
> a bona-fide non-negative integer enumerator. The obstruction, if the code does
> not exist, is purely **code-combinatorial existence** — a question about ℤ₂
> realizability, not modular forms.

This is the asymmetry with the LATTICE side, which is SOLVED: Nebe's extremal
even unimodular Γ₇₂ (min 8) exists (arXiv:1008.2862; Crelle 673 (2012) 237–247),
built as Barnes ⊗ Leech over ℤ[(1+√−7)/2] with (PSL₂(7)×SL₂(25)):2 ⊆ Aut. The
modular-form-side enumerator/theta is non-negative and realized by a lattice;
only the binary-code realization is open.

## B. The Paley/eQR gauge route reaches exactly d = 12 (VERIFIED, THM-481)

72 = q+1 at the prime q = 71 ≡ 7 (mod 8), so THM-481's tournament-gauge identity
applies: C(border(Paley₇₁)) = eQR(72), a Type II self-dual [72,36] code — but with
**d = 12, not 16** (eQR(72) is the classical d=12 code, Aut ⊇ PSL₂(71)). This is
the FIRST failure of the eQR extremal ladder: extremal at lengths 8, 24, 32, 48
(d = 4,8,8,12 = 4⌊n/24⌋+4), the Paley gauge stops being extremal at 72. The
symmetry that makes the eQR codes natural (a transitive PSL₂ action) is exactly
what the SOTA says the extremal code CANNOT have (Borello: |Aut| ∈ {1,2,3,4,5},
abelian, "almost a rigid object"). The gauge construction and extremality pull in
opposite directions at 72: **arithmetic symmetry caps the gauge distance at 12;
extremality (if achievable) requires near-rigidity.** This is the precise
content of HYP-2415's "asymmetry conjecture."

## C. The 72 = 3·24 = Δ³ gate (the modular reading)

Extremality at 24m kills the first m coefficients = the Δ^{≤m} corrections
(THM-486 B). At m = 3 (length 72) extremality kills exactly the Δ, Δ², Δ³ strata.
The lattice side realizes the Δ³ gate (Γ₇₂); the code side is the open question of
whether the ℤ₂-Gleason ring's Δ³-extremal point is a genuine code. The
deterministic secular machinery (THM-486 C) is silent at m = 3 (no negativity
until m = 154), so — unlike the lattice nonexistence at dimension 163264 — the 72
problem gets NO help from modular-form positivity. It is the rare case where the
analytic obstruction theory is vacuous and the answer must be combinatorial.

## D. State of the art (recon, DOIs)

- Open since Sloane 1973 (IEEE-IT 19, 251); prizes Sloane \$10 + Dougherty \$100
  (existence) + Harada \$200 (nonexistence) = \$310; no Conway \$1000 prize (that
  list is the 99-graph, thrackle, etc.).
- Aut(C) ∈ {1, C₂, C₃, C₂×C₂, C₅}, all elements prime order, involutions and
  order-3 fixed-point-free, order-5 with 2 fixed points (Borello survey
  arXiv:1311.3868 Thm 6.3, unchanged since the Yorgov–Yorgov 2014 order-4
  elimination, arXiv:1310.2570). CAUTION: arXiv:1303.4920 ("|Aut| ∈ {1,3,5}") was
  WITHDRAWN — orders 2 and 4 (Klein) remain open.
- Best known [72,36] self-dual codes have d = 12 (thousands, incl. eQR(72)); even
  d = 14 has never been reached at length 72.
- Lattice analogue SOLVED (Γ₇₂, Nebe 2010/2012).

## Honesty

- No existence/nonexistence progress is claimed. The value of this file is the
  localization: the modular-form/enumerator/lattice machinery that controls other
  lengths is non-obstructing at 72, isolating the question as combinatorial, and
  the gauge-vs-rigidity tension (B) as the structural reason the natural
  constructions stall.
- A₁₆ = 249849 and the all-positive W₇₂ are exact; the Aut and SOTA facts are
  cited, not re-derived.

**Cross-refs:** THM-486 (n = 3696, the η-bridge), THM-481 (eQR gauge),
THM-484/claudebox (involution modulus, the 72 gate), OPEN-Q-061, HYP-2415,
HYP-2421 (sharpened here).
