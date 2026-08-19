---
id: THM-1330
title: THE KELLER MONOID NECESSARY-STRUCTURE ATLAS — units = degree 1, the ideal {deg ≥ 3}, finite factorization, universal non-properness/ramification-at-infinity, the inverse Jelonek problem, and the Zariski cusp selection rule
status: PROVED-with-citations FOR THE INDIVIDUAL NECESSARY-STRUCTURE CLAIMS; NOT A CLASSIFICATION OF ALL COUNTEREXAMPLES. The monoid/ideal/factorization laws are elementary given Keller 1939 [birational case], Białynicki-Birula–Rosenlicht [injective ⟹ automorphism], and Campbell [Galois ⟹ invertible]. Non-properness is a covering-space argument; the selection rule uses Zariski–Lefschetz and Deligne–Fulton; the displayed example's cusp computations are exact/numeric-certified. The inverse Jelonek realization problem, classification beyond the verified weighted irreducible family, and JC(2) remain open.
source: mac-mini-2026-07-20-S135 (owner: can we get an exact picture of ALL counterexamples?)
depends_on: [THM-1300, THM-1305, THM-1310, THM-1315, HYP-8115]
related: [THM-1320, THM-1325 (Drużkowski strata), THM-3438, klein-S325 (Smith selection rule), HYP-8125]
output: 05-knowledge/results/jc_keller_monoid_cusp_rule_macmini_S135.out
---

# THM-1330 — necessary structure of the Keller nonunit locus

> **Scope correction (MISTAKE-236).** Earlier versions called the results below
> an "exact picture" and called the covering data a classification. They are a
> strong list of universal invariants and necessary conditions. They do not
> enumerate irreducible Keller maps, prove that every abstract covering triple
> is realizable, or give a complete invariant up to polynomial equivalence.

Fix n and let 𝒦ₙ = {polynomial F: ℂⁿ → ℂⁿ with det JF ∈ ℂ*} (Keller maps),
𝒳ₙ = 𝒦ₙ ∖ Aut(ℂⁿ) (the counterexamples; nonempty for n ≥ 3).

## 1. The algebraic shape (exact)

- **(monoid)** 𝒦ₙ is closed under composition (chain rule: the determinant of a composite is
  a product of constants). Units of 𝒦ₙ = Aut(ℂⁿ).
- **(degree homomorphism)** deg F := [ℂ(x₁..xₙ) : F*ℂ(x₁..xₙ)] is multiplicative
  (function-field tower), giving deg: 𝒦ₙ → (ℕ≥1, ×).
- **(units = degree 1)** deg F = 1 ⟺ F ∈ Aut. (⟸ trivial; ⟹ is Keller's own 1939 theorem:
  birational Keller maps are invertible.)
- **(no degree 2)** deg F = 2 forces a Galois cover; Campbell: Galois Keller ⟹ invertible.
  Hence **𝒳ₙ = deg⁻¹({3, 4, 5, …})** — the exact boundary of the counterexample set.
- **(two-sided ideal)** if G∘F ∈ Aut then F is injective, hence (Białynicki-Birula–
  Rosenlicht) F ∈ Aut, hence G ∈ Aut. So 𝒳ₙ is a two-sided ideal: counterexamples absorb
  composition.
- **(finite factorization)** deg multiplicativity + units = deg⁻¹(1) ⟹ every F ∈ 𝒳ₙ is a
  FINITE composition of irreducible Keller maps. The realizable-degree image `R_n=deg(𝒦ₙ)`
  is a multiplicative submonoid of `{1,3,4,5,…}`.  For every `n>=3`, THM-3438's
  weighted lifts prove the exact value spectrum `R_n={1,3,4,5,…}`; for `n=1` only degree
  one occurs, while the `n=2` image remains tied to `JC(2)`.  Since every nonunit factor
  has degree at least three, every map of degree `3` through `8` is automatically a
  composition atom; more generally this holds whenever the degree is not a product
  `ab` with `a,b>=3`.  THM-3438 further gives an `S_d` atom in every grade `d>=3`,
  while compositions realize every product grade `ab`; hence every such grade is
  mixed and degree nine is the first.  Degree one is the unit, not an atom.
- **(degree-3 monodromy)** deg 3 non-Galois ⟹ monodromy exactly S₃ (THM-1310/1315 for F;
  forced generally).

## 2. The geometric shape (universal anatomy — every counterexample, provably)

- **(étale)** det JF ∈ ℂ* ⟹ F is everywhere a local isomorphism: NO finite ramification,
  ever. All branching of every counterexample is at infinity.
- **(non-proper, three lines)** if F were proper, étale + proper over simply-connected ℂⁿ
  makes F a finite covering of ℂⁿ, hence trivial, hence deg 1, hence a unit. So every
  counterexample is non-proper: **the Jelonek set A(F) ≠ ∅ always** — escape at infinity is
  not a feature of one example (THM-1315) but a law of the whole set.
- **(covering constraint)** over ℂⁿ ∖ A(F), F is proper + étale = an honest deg-F-sheeted
  covering. Hence every counterexample is exactly a triple
  **(A ⊂ ℂⁿ closed, a connected non-Galois d-cover of ℂⁿ∖A given by a monodromy surjection
  π₁(ℂⁿ∖A) ↠ (transitive, non-regular) G ≤ S_d, degeneration data over A)**
  subject to polynomial-with-constant-Jacobian gluing. This records necessary
  data, but no theorem here says that it is a complete invariant or that an
  abstract allowed triple glues. Those realization and equivalence questions form the
  **inverse Jelonek problem**. (For F: A = the quartic {K = 0}, d = 3, G = S₃, degeneration
  = one sheet escaping — all explicit.)

## 3. The selection rule (which Jelonek hypersurfaces are even allowed)

π₁(ℂⁿ∖A) must SURJECT onto a nonabelian group already at d = 3 (S₃). By Zariski–Lefschetz,
π₁(ℂ³∖A) ≅ π₁(plane ∖ (generic plane section of A)); by Deligne–Fulton, a NODAL section
curve has abelian complement group. Therefore:

> **the generic plane section of the Jelonek hypersurface of any degree-3 counterexample
> must carry worse-than-nodal singularities.**

**Verified on F (frozen out):** the section of {K = 0} by a generic plane is a quartic with
EXACTLY THREE singular points, ALL of Hessian rank 1 (cusps), ZERO nodes — a three-cuspidal
quartic: **Zariski's 1929 curve**, the first curve known to have nonabelian complement group
(metacyclic of order 12, which surjects onto S₃ — precisely the required monodromy supply).
The singular locus of {K = 0} itself is the Gröbner-exact cuspidal edge
{16a = b³c, (3bc−4)² = 0}, i.e. the rational curve **(a,b,c) = (b²/12, b, 4/(3b))** — the
A₂-edge the fleet's ADE reading (kps/klein) predicted. Also verified: |F⁻¹| = 3 and
|(F∘F)⁻¹| = 9 at a generic target (degree multiplicativity in the flesh).

## 4. What is exactly known vs open about 𝒳ₙ (the honest ledger)

KNOWN EXACTLY: the boundary (units = deg 1); for `n>=3`, the degree-value spectrum
`{1,3,4,5,…}` (THM-3438); the degree stratification (no 2; ideal = deg ≥ 3);
finite factorization; universal étaleness/non-properness/ramification-at-infinity; the
covering constraint; the cusp selection rule; near-surjectivity (the image misses at
most codimension 2 — van den Essen circle; the fixed F has
`F(A^3)=A^3 minus E` and fibre spectrum `3/1/0`, never `2`, by the
MISTAKE-282 correction to THM-1315); at
(n, d) = (3, 3): rigidity of F (THM-1305: deformation tangent = orbit tangent) and uniqueness
in every tested design box (THM-1310, boxeph-S142).

CONJECTURAL: Mount-Everest (HYP-8115(d)): at (3,3) there is exactly ONE irreducible
counterexample up to two-sided Aut-equivalence — the strongest "exact picture" statement;
current evidence is rigidity + multi-box uniqueness, not a proof.

OPEN: classification within each realized degree, decomposition in the numerically
factorable grades (beginning at degree 9), `JC₂` (hence whether `𝒳₂=∅`) and the Dixmier
floor `DC₂ ⟹ JC₂ ⟹ DC₁`; unique factorization (up to units) in 𝒦ₙ; the full
inverse-Jelonek realization theory (which allowed triples actually glue polynomially).
