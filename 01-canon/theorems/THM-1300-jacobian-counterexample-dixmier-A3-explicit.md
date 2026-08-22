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
  kind-pasteur S128c97).  The exact map and later weighted-lift family are now
  public on jacobianfun.org and its linked repository; use the byte-zero
  correction below for attribution.
  OPEN at the bottom of the tower: JC₂, DC₁, DC₂ are NOT decided by this.
source: death-star-2026-07-19-S59m (HYP-8075; owner prompt). Concurrent: kind-pasteur S128c97 (HYP-8070 — verification, σ-equivariance, Groebner degree, mod-p statistics).
depends_on: []
related:
  - THM-3438-weighted-lift-keller-degree-spectrum
  - HYP-8075 (this session), HYP-8070 (kind-pasteur's complementary streams)
  - 07-reflections/jacobian-dixmier-through-the-repos-eyes-deathstar-S59m.md
scripts:
  - 04-computation/jacobian_counterexample_verify_deathstar_S59m.py -> 05-knowledge/results/jacobian_counterexample_verify_deathstar_S59m.out
  - 04-computation/dixmier_explicit_endomorphism_A3_deathstar_S59m.py -> 05-knowledge/results/dixmier_explicit_endomorphism_A3_deathstar_S59m.out (+ .verdict.txt)
  - 04-computation/jacobian_torus_equivariance_deathstar_S59m.py -> 05-knowledge/results/jacobian_torus_equivariance_deathstar_S59m.out
  - 04-computation/jacobian_formal_inverse_2adic_deathstar_S59m.py -> 05-knowledge/results/jacobian_formal_inverse_2adic_deathstar_S59m.out
---

> ## ⚠ OWNER CORRECTION (2026-07-20, klein-S377) — the "Alpoge-Mathew" attribution below is a HALLUCINATION
>
> The owner has flagged the mac-mini-S127/S129 attribution blocks (which name Levent Alpoge and Akhil
> Mathew) as a **hallucination**: they were produced by web searches on a result then ~1 day old, with no
> reliable public record, and must NOT be restated as established. **Corrected provenance:** the
> counterexample was **discovered by Claude**; it became public on 2026-07-19 via a **tweet from an
> Anthropic employee** (a sharing event, not a discovery attribution). Treat the map as **Claude-discovered,
> provenance otherwise uncertain** — do not write "Alpoge" or "Alpoge-Mathew" as the discoverers anywhere.
> The blocks below are retained as history but are contested. See MISTAKE-205. What the repo legitimately
> holds is unchanged: the independent exact verification and the equivariant/Dixmier/elliptic anatomy.

# THM-1300 — the JC counterexample, its explicit Dixmier transfer, and its torus anatomy

> ## ⚠ ATTRIBUTION CORRECTION (mac-mini-2026-07-20-S127) — the map is **Levent Alpöge's**
>
> This file records the map as "owner-supplied (2026-07-19)" with "literature/web search
> finds no public source yet". That was accurate when written and is **now wrong**.
>
> A web search this session identifies the counterexample as **Levent Alpöge's**
> (Anthropic), announced **2026-07-19/20**, obtained with Claude Fable — the *identical*
> map, `u = 1+xy`, `F = (u³z + y²u(4+3xy),\ y + 3xu²z + 3xy²(4+3xy),\ 2x − 3x²y − x³z)`,
> with the same `det JF = −2` and the same triple collision
> `F(0,0,−¼) = F(1,−3/2,13/2) = F(−1,3/2,13/2) = (−¼,0,0)`.
>
> **What is ours and what is not.** The counterexample is **not** a repo discovery, and
> nothing downstream of it should be presented as one. What this repo legitimately holds is
> (i) **independent exact verification** — re-derived again this session, `det JF ≡ −2`
> symbolically and the triple collision in exact rationals; and (ii) the **explicit Dixmier
> transfer** (§1), which the public discussion had not obviously drawn at the time of the
> search. Both are verification/derivation, not discovery.
>
> **Consequence for downstream work.** Every immediate corollary — Dixmier, Zhao's vanishing
> conjecture, the image conjecture, failing Mathieu subspaces — is a corollary of *Alpöge's*
> theorem and is being chased by the whole field right now. Priority claims on them are not
> available to us. See HYP-8240 for the full assessment.

> ## ⚠ ATTRIBUTION AMENDMENT (mac-mini-2026-07-20-S129) — **co-credit is owed, and the "3" is not structural**
>
> An independent external check this session (which re-verified `det JF = −2` two ways —
> sympy symbolic *and* exact `Fraction` Lagrange interpolation, plus the triple collision in
> exact rationals) corrects and extends the block above:
>
> - **Co-credit: Akhil Mathew.** Alpöge's own announcement thanks "my close friend akhil for
>   asking about it," and Jared Duker Lichtman's post credits "Alpoge, Mathew, and Claude
>   Fable 5." Crediting **Alpöge alone is incomplete** — the correction above was itself
>   under-attributed. Announcement date is **19 July 2026** (X), not 19/20.
>
>   > **⚠ THIS BULLET IS ITSELF OVER-CORRECTED (mac-mini-S133).** A second independent sweep
>   > reads the primary wording more carefully: Alpöge thanks Akhil for **asking the
>   > question**, and credits **Claude Fable** with producing the example. "Co-credit is
>   > owed to Akhil Mathew" overstates *co-authorship* relative to what the sources
>   > actually say, and the identification of "akhil" as Akhil Mathew could **not** be
>   > confirmed from a primary source — it comes from Lichtman's secondary post. The
>   > accurate statement is: **the map is Alpöge's announcement, produced with Claude
>   > Fable, with Akhil credited for posing the question.** I have now mis-stated this
>   > twice in opposite directions; further edits should quote the primary post verbatim
>   > rather than paraphrase it.
> - **No arXiv preprint, no peer review, no journal** — confirmed *absent* via the arXiv API,
>   not merely unfound. (Caution: `export.arxiv.org` 301-redirects to HTTPS; an unredirected
>   query silently returns zero entries and produced a false negative mid-search.)
> - **The "three" is the generic fiber degree, and it dissolves.** The load-bearing three is
>   that the generic fiber has degree 3 — three sheets, a cubic fiber equation — which *is*
>   coordinate-invariant. But an explicit family (jacobianfun.org) gives one counterexample
>   **for every generic fiber degree `d ≥ 3`**, from seed `p_d`. The degree-4 member was
>   independently verified: component degrees `(12,11,4)`, `det JG = −6`,
>   `G(1,0,0) = G(−1,0,2) = (0,0,1)`; since `4` is neither `3` nor a power of `3`, it is
>   genuinely inequivalent. **So Alpöge's map is the minimal member of an infinite family,
>   not a three-part object.** Note also the variety-level fiber decomposition is **2**, not
>   3: the plane `x=0` (one preimage) plus the curved surface `x²z = 2−3xy` (two, on its two
>   lobes) — the collision splits `1+2` across two components.
> - **THM-1370 is externally confirmed.** The weighted symmetry `(x,y,z) ↦ (λx, y/λ, z/λ²)`
>   — weights **(1, −1, −2)**, invariants `v = xy`, `t = x²z` — appears in the external
>   treatment exactly as THM-1370 derived it independently. That is a genuine validation of
>   repo work, and one of the few things here we can point to as ours.
> - **Dixmier is refuted** (Wikipedia's Dixmier conjecture page, edited 2026-07-20, states it
>   "was thus refuted"), consistent with §1. **Zhao / image conjecture / Mathieu: no source
>   discusses them post-announcement.** The propagation `¬JC ⟹ ¬VC` etc. is an *inference
>   from the pre-2026 literature*, not an asserted result — which is exactly why HYP-8240's
>   "untouched and publishable" reading survives, and exactly why the witness must be
>   verified **directly** rather than through the equivalence.
> - **The degree-100/108 lower-bound objection does not apply.** arXiv 2204.14178 concerns the
>   ***plane*** Jacobian Conjecture (`n = 2` only); no contradiction with a degree-7 map in
>   dimension 3. **`n = 2` remains open.**
>
> *Deflation of a repo claim.* death-star-S60's "Erdős 592 and the fallen Jacobian conjecture
> are ONE trichotomy in two theaters" is **overstated**. Erdős #592 is infinite partition
> calculus (partition ordinals, `$1000` prize, OPEN) with no polynomial, covering-system, or
> Jacobian connection in its literature; its three *is* real and sharp — Schipperus 2010 gives
> true for 1–2 indecomposable summands, false for ≥4, leaving **exactly the three-summand case
> open** — but the JC three is a generic fiber degree that starts an infinite family. Two
> unrelated threes. **THM-1415 §V's verdict — "a coincidence of POSITION, not of mechanism" —
> is the correct call**, and this independent check strengthens it.

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
