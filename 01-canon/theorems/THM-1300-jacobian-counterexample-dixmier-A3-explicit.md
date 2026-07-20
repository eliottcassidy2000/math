---
id: THM-1300
title: The Jacobian Conjecture counterexample (owner-supplied, dim 3, det = −2) — independently verified exactly; the EXPLICIT Dixmier counterexample φ: A₃ → A₃ (all 18 Weyl/flatness identities verified symbolically; proper self-embedding via the module one-liner); the ℂ*-equivariant structure (weights (1,−1,−2) → (−2,−1,1)) exposing the triple collision as 1 fixed-branch + 2 orbit-branch points of the torus doubling λ ↦ λ²; the formal inverse's unbounded dyadic coefficient ladder.
status: >
  VERIFIED-EXACT — det JF ≡ −2 and the triple collision F(0,0,−1/4) =
  F(1,−3/2,13/2) = F(−1,3/2,13/2) = (−1/4,0,0) (exact rational arithmetic,
  dependency-free script). Hence JC is FALSE for n = 3, and for all n ≥ 3 by
  identity padding. PROVED (modulo classical cited facts: A_n simple; the Weyl
  presentation's universal property) — φ(X_i) = F_i, φ(D_j) = Σ_k B_jk D_k with
  B = (JF^T)^{-1} is a well-defined injective NON-surjective endomorphism of A₃:
  the Dixmier conjecture is FALSE for A_n, n ≥ 3, constructively. VERIFIED-EXACT —
  the ℂ*-equivariance and the orbit-branch collision law. PROVENANCE — the map is
  owner-supplied (2026-07-19); two independent in-repo verifications (this file;
  kind-pasteur S128c97); literature/web search finds no public source yet.
  OPEN at the bottom of the tower: JC₂, DC₁, DC₂ are NOT decided by this.
source: death-star-2026-07-19-S59m (HYP-8075; owner prompt). Concurrent: kind-pasteur S128c97 (HYP-8070 — verification, σ-equivariance, Groebner degree, mod-p statistics).
depends_on: []
related:
  - HYP-8075 (this session), HYP-8070 (kind-pasteur's complementary streams)
  - 07-reflections/jacobian-dixmier-through-the-repos-eyes-deathstar-S59m.md
scripts:
  - 04-computation/jacobian_counterexample_verify_deathstar_S59m.py -> 05-knowledge/results/jacobian_counterexample_verify_deathstar_S59m.out
  - 04-computation/dixmier_explicit_endomorphism_A3_deathstar_S59m.py -> 05-knowledge/results/dixmier_explicit_endomorphism_A3_deathstar_S59m.out (+ .verdict.txt)
  - 04-computation/jacobian_torus_equivariance_deathstar_S59m.py -> 05-knowledge/results/jacobian_torus_equivariance_deathstar_S59m.out
  - 04-computation/jacobian_formal_inverse_2adic_deathstar_S59m.py -> 05-knowledge/results/jacobian_formal_inverse_2adic_deathstar_S59m.out
---

# THM-1300 — the JC counterexample, its explicit Dixmier transfer, and its torus anatomy

## 0. The map

F = (F₁, F₂, F₃): ℂ³ → ℂ³, with u := 1 + xy:

- F₁ = u³z + y²u(4 + 3xy)   (degree 7)
- F₂ = y + 3xu²z + 3xy²(4 + 3xy)   (degree 6)
- F₃ = 2x − 3x²y − x³z   (degree 4)

**Verified exactly:** det JF ≡ −2 (constant), and the three points (0, 0, −1/4),
(1, −3/2, 13/2), (−1, 3/2, 13/2) all map to (−1/4, 0, 0). A unit-Jacobian
polynomial self-map of ℂ³ that is not injective: **the Jacobian Conjecture is
false at n = 3**, hence (append identity coordinates — det and non-injectivity
persist) **for every n ≥ 3**. n = 2 is untouched.

## 1. The explicit Dixmier counterexample (the session's construction)

Let B := (JF^T)^{-1} = adj(JF)^T/(−2) — a polynomial matrix over ℤ[1/2] because
det JF is the unit −2. Entries: degrees 6–11, only 5–14 terms each (stored in
full in the results file). Define on generators of the Weyl algebra A₃:

  φ(X_i) = F_i(X₁, X₂, X₃),  φ(D_j) = Σ_k B_jk(X) D_k.

**Symbolically verified (exact, 18 identities):**
- (R1) B·(JF)^T = I — equivalently [φD_j, φX_i] = δ_ij (9 identities);
- (R2) Σ_k (B_ik ∂_k B_jl − B_jk ∂_k B_il) = 0 for all i < j, l — equivalently
  [φD_i, φD_j] = 0 (9 identities; the flatness/commuting-fields condition).
- [φX_i, φX_j] = 0 is automatic (both lie in ℂ[X]).

By the universal property of the Weyl presentation, φ extends to an
endomorphism of A₃; since A₃ is simple and φ ≠ 0, **φ is injective**.

**φ is not surjective — the module one-liner.** Every element of im(φ) is a
ℂ-combination of φ(X^β D^α) = F^β V^α, where V_j := Σ_k B_jk D_k acts on the
standard module ℂ[x,y,z] as a derivation (so V_j(1) = 0). If X₁ ∈ im(φ), apply
both sides of X₁ = Σ c_{β,α} F^β V^α to 1 ∈ ℂ[x,y,z]: every α ≠ 0 term dies,
leaving x₁ = Σ_β c_{β,0} F^β ∈ ℂ[F]. The same for x₂, x₃ gives ℂ[F] = ℂ[X],
i.e. a polynomial left inverse G with G ∘ F = id — forcing F injective,
contradicting the verified triple collision. Hence X₁ ∉ im(φ):

**φ is an explicit proper self-embedding of A₃. The Dixmier conjecture is false
for A₃ — hence for A_n, n ≥ 3 (pad φ with identity on the extra generators).**
This is the classical DC_n ⟹ JC_n transfer run constructively on a concrete map;
nothing in it depends on the map's provenance, only on the verified arithmetic.

## 2. The ℂ*-anatomy of the collision (verified)

Under the weighted action λ·(x, y, z) = (λx, λ⁻¹y, λ⁻²z):

- Each component is **weighted-homogeneous**: w(F₁) = −2, w(F₂) = −1, w(F₃) = +1
  (the building blocks u and 4+3xy have weight 0). So F is ℂ*-equivariant:
  F(λ·p) = (λ⁻²F₁, λ⁻¹F₂, λF₃)(p). kind-pasteur's involution σ = (−x,−y,z) is
  the λ = −1 element of this torus.
- The target a-axis {(a,0,0)} has (at least) two preimage branches:
  the **fixed branch** x = y = 0, on which F(0,0,z) = (z,0,0) — bijective onto
  the axis; and the **orbit branch** ℂ*·(1, −3/2, 13/2), on which (verified as
  an exact Laurent identity) F(λ, −3/(2λ), 13/(2λ²)) = (−1/(4λ²), 0, 0) —
  the orbit maps **2:1** onto the punctured axis via λ ↦ λ².
- The verified triple collision is exactly λ = ±1 of the orbit branch plus the
  fixed point: **3 = 1 + 2, a fixed point plus a doubled orbit** — the
  "Rédei-shaped fiber" of kind-pasteur's reading, now with the doubling map
  λ ↦ λ² identified as the engine. Geometric degree of F is ≥ 3 (Groebner
  confirmation of "= 3" is kind-pasteur's stream).

## 3. The formal inverse's dyadic ladder (measured)

F(0) = 0 and JF(0) = (x,y,z) ↦ (z, y, 2x), so the compositional inverse G
exists as a formal power series — and can never be polynomial (a polynomial
G with G∘F = id would force injectivity). Computed to total degree 12:

- G is nonzero at **every** degree, but **sparse**: 3,2,2,4,3,4,5,5,5,7,6,7
  terms per degree (vs 91 possible at degree 12) — the sparsity is the
  weighted grading of §2 acting on the inverse.
- Coefficients live in ℤ[1/2] and their minimal 2-adic valuation **descends
  without bound, in pairs**: −1,−1,−3,−3,−4,−4,−7,−7,−8,−8,−10,−10 at degrees
  1..12 (the pairing is the λ ↦ −λ parity). Every mod-2^k truncation of the
  inverse dies at finite degree; the obstruction to inverting F is a
  **dyadic staircase** with no bottom. (Consistently: det ≡ 0 mod 2 — the
  counterexample degenerates at the prime 2 itself.)

## 4. What stands after the fall

- FALSE: JC_n (n ≥ 3), DC_n (n ≥ 3) — the latter constructively (§1).
- UNTOUCHED: JC₂; DC₁ (Dixmier's original question), DC₂. The classical
  bridges JC_{2n} ⟹ DC_n (Tsuchimoto; Belov-Kanel–Kontsevich) with 2n ≥ 4 are
  now vacuous at the refuted end; DC₂ ⟹ JC₂ remains a live route to JC₂.
  The open frontier has moved to the BOTTOM of the doubling tower.
- The repo-facing reading (Mode A/B, the doubling ladder, the observer
  principle, the dyadic layer) is in the companion reflection.

## 5. Fleet integration (added at close — six sessions, one prompt)

Concurrent independent streams: kind-pasteur S128c97 (σ-equivariance, mod-p plan),
klein S323 (Groebner fiber-complete, degree-3 étale, the SYMPLECTIC ℂ⁶ cotangent
lift Φ = (F, J^{-T}ξ) with Φ*ω = ω — the explicit n-vs-2n object; S₃ monodromy),
boxeph S140 (full S₃ Chebotarev (1/2, 1/6, 1/3); T(p) collision counts), mac-mini
S129 (Groebner; odd-fiber parity; tournament-Dixmier testable), opus S413 (parity
lemma; second fiber F⁻¹(1,0,0) = {(0,0,1)} ∪ {(±i/2, ±3i, −26)}; JC₂ ⟹ DC₁ the
last bridge). This file's unique adds interlock with two of their filed items:

- **klein's named task "a concrete element outside im(φ_F)" is CLOSED by §1**:
  X₁ ∉ im(φ), with the self-contained module one-liner (no citation needed).
  boxeph's "unowned, high value" lead (the explicit endomorphism) is this file.
- **opus's second fiber is a corollary of §2's orbit law**: the Laurent identity
  F(λ, −3/(2λ), 13/(2λ²)) = (−1/(4λ²), 0, 0) is a polynomial identity over ℚ,
  hence holds in every extension; the fiber over (a, 0, 0), a ≠ 0, is
  {(0, 0, a)} ∪ {the two λ with λ² = −1/(4a)} — at a = −1/4: λ = ±1 (the
  original collision); at a = 1: λ = ±i/2, giving exactly (±i/2, ±3i, −26).
  Over ℝ the orbit pair is real iff a < 0 — real non-injectivity lives only on
  the negative half-axis. The fleet's three measured fibers (a = −1/4, 1, and
  the mod-p 1-vs-3 dichotomy on the axis, which is the QR split of −1/(4a))
  are one formula.
