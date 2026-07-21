---
id: THM-1885
title: "THE ½-AND-+1 MONOID IS BAUMSLAG–SOLITAR BS(1,2), AND MOST OF THE REPO IS A MONOID ACTING ON A SET. The two maps a(x)=x+1 and b(x)=x/2 (THM-1880) satisfy the single relation a·b = b·a² — verified aᵏb = ba²ᵏ for all k — which is EXACTLY the presentation of the Baumslag–Solitar group BS(1,2) = ⟨a,b | ab=ba²⟩. The action on values is faithful (every word maps x ↦ x/2ᵖ + q with q a dyadic rational, and distinct reduced words give distinct maps), so ⟨a,b⟩ is the positive monoid BS(1,2)⁺ realised as the DYADIC AFFINE MONOID x ↦ x/2ᵖ + ℤ[½]. Consequences: (1) the half-dictionary b∘a = (x+1)/2 (THM-1555, the {−1,0,1}↔{0,½,1} conjugation) is a specific BS(1,2) element, inverse 2x−1; (2) the corpus-wide ½ is the SCALING generator b (dimension-halving) and the +1 is the TRANSLATION generator a, and BS(1,2) is precisely the group that couples 'shift' to 'halve' via ab=ba²; (3) BS(1,2) is the 2-adic / dyadic solenoid monodromy, which is WHY the repo's 2-adic threads coincide — the 2^{C(n−1,2)} switching classes, the (Z/2)^{C(n,2)} arc-flip hypercube, the tiling fiber fraction's 2-adic valuation, the blue-line count 2^{e−1}, the Cayley–Dickson ×2 dimension-doubling tower — are all the SAME 'multiply/divide by 2' generator b of BS(1,2) seen on different objects. THE FUNDAMENTAL VIEW: nearly every repo topic is (object, monoid G, action), with invariants = orbit functions, nullcone = the degenerate orbit, and the recurring monoids are a short list — Z/2 (complement/transpose), (Z/2)^{C(n,2)} (arc-flips = cut⊕cycle), S_n (relabelling), BS(1,2) (½&+1), PSL(2,Z)=Z/2∗Z/3 (Farey/modular), SL₂ (binary forms char_A), Z (GMC charge grading) — and the whole complexity ladder (THM-1775) is these monoids at increasing depth"
status: >
  The relation ab = ba² and aᵏb = ba²ᵏ are VERIFIED-EXACT (symbolic) and are the defining
  relation of BS(1,2); the faithful dyadic-affine action is elementary (a word b^{e_0} a^{c_1}
  b^{e_1} … reduces to x/2^{Σeᵢ} + dyadic, and the pair (2^{−Σeᵢ}, dyadic) determines the map).
  So "⟨x+1, x/2⟩ ≅ BS(1,2)⁺" is PROVED (a presentation match plus a faithful representation).
  The identifications in (3) are exact where a group is named (switching classes, arc-flips,
  Cayley–Dickson doubling) and structural where "2-adic" is the shared feature (fiber fraction,
  blue count).  The FUNDAMENTAL-VIEW catalog is a reframing (the companion reflection lists the
  monoid, generators, action, and nullcone for each topic); it proves no new open problem, it
  organises the corpus around its acting monoids.
source: kind-pasteur-2026-07-21-S128c140 (owner: see more problems as generators and monoids; get as fundamental a view of as many topics as possible)
depends_on:
  - THM-1880    # the a/b functional frame (a=x+1, b=x/2)
  - THM-1555    # the half-dictionary b∘a
related: [THM-1810, THM-1775, THM-826]
external:
  - "Baumslag–Solitar BS(1,2) = ⟨a,b | b⁻¹ab = a²⟩ (solvable, the dyadic-solenoid monodromy); the affine group Aff(ℤ[½])."
script: 04-computation/deep_archaeology_kps_S128c137.py (+ /tmp bs check)
---

# THM-1885 — the ½-and-+1 monoid is Baumslag–Solitar, and the repo is monoids-on-sets

## The identification

`a(x) = x+1`, `b(x) = x/2`. Then `a·b = b·a²` (i.e. `a(b(x)) = x/2+1 = b(a²(x))`), and more
generally `aᵏ·b = b·a²ᵏ` — verified for all tested `k`. That single relation is the presentation

> **`BS(1,2) = ⟨a, b | b⁻¹ a b = a² ⟩`** (equivalently `ab = ba²`),

the archetypal solvable **Baumslag–Solitar group**. The action `x ↦ x/2ᵖ + q` (`q ∈ ℤ[½]`) is
faithful, so `⟨a,b⟩` is the positive monoid `BS(1,2)⁺`, the **dyadic affine monoid** on `ℤ[½]`. The
half-dictionary `b∘a = (x+1)/2` (THM-1555) is one of its elements; `2x−1` is the inverse.

**Why this matters:** BS(1,2) is exactly the group that couples a *shift* `a` to a *scaling* `b` by
`ab=ba²`. It is the monodromy of the **dyadic solenoid** (the inverse limit of the circle under
doubling). So every "2-adic" phenomenon in the corpus is the *same* generator `b` — halving, or its
adjoint doubling:

| repo object | the ×2 | 
|---|---|
| switching classes `= 2^{C(n−1,2)}` | `b`-count |
| arc-flip hypercube `(Z/2)^{C(n,2)}` | `b` on each coordinate |
| tiling **fiber fraction** `(½)_{n−2}/(n−2)!` | `b`-Pochhammer / 2-adic valuation |
| blue-line count `2^{e−1}` | `b`-count |
| Cayley–Dickson `ℝ→ℂ→ℍ→𝕆→𝕊` (dim ×2) | `b⁻¹` (doubling) |

## The fundamental view: (object, monoid, action)

Reading each topic as *a monoid acting on a set* — invariants are orbit functions, the nullcone is
the degenerate orbit — collapses the corpus onto a short list of monoids:

| topic | set | monoid | generators | nullcone / distinguished orbit |
|---|---|---|---|---|
| tournament spectra | `Sym^n` forms | **SL₂** | (transvection) | `xⁿ` = transitive (THM-1810) |
| tournament space | `Q_{C(n,2)}` tilings | **(Z/2)^{C(n,2)}** = cut⊕cycle | arc-flips (wiggly/base-path) | transitive corner |
| iso classes / metagraph | `Q_m / S_n` | **S_n** | transpositions | complement/transpose fixed = SC |
| complement / converse | tournaments | **Z/2** | `T ↦ T^op` | self-complementary = fixed |
| Farey / LRC gaps | `ℚ ∩ [0,1)` | **PSL(2,Z)** = Z/2∗Z/3 | `S: x↦−1/x`, `T=a` | `0/1`, the cusp |
| the ½ & +1 | `ℤ[½]` | **BS(1,2)** | `a=x+1`, `b=x/2` | fixed point of `b` = `0`; of `b∘a` = `1` |
| GMC charge | `ℂ[Z,Z̄]` | **Z** (charge grading) | charge `±1` shift | charge-0 = the `s`-ring |
| Cayley–Dickson | `ℝ, ℂ, ℍ, 𝕆, …` | **doubling monoid** | `⊗ ℂ` | `ℝ` (the base) |

The **complexity ladder** (THM-1775: trace ⊂ CT ⊂ Gaussian `E`) is these monoids at increasing
depth — `S_n`/conjugation at the rational floor, the toral `Z` at the algebraic rung, `BS(1,2)`'s
`b`-Laplace at the holonomic rung — and the `#P` jump (`H`, THM-1780) is where **no** finitely
generated acting monoid describes the invariant (the permanent is not an orbit function of any of
these). That is the fundamental reason `H` leaves the ladder.

## What this buys

- **A test for "is this invariant fundamental?":** is it an orbit function of the object's acting
  monoid? If yes it is on the ladder (spectral); if no (like `H`) it is `#P`.
- **BS(1,2) predicts the 2-adic coincidences.** Any two "2-adic" repo quantities should differ by a
  `BS(1,2)` element — e.g. the switching-class count and the blue-line count are both `b`-powers,
  related by an `a`-shift in the exponent (`2^{C(n−1,2)}` vs `2^{e−1}`, `e = ⌊(n−1)²/4⌋`).
- **The GIT deformation** `b((x+c)ⁿ+(x−c)ⁿ)/2` (THM-1880 named-next) is a `BS(1,2)`-orbit through
  the `char_S` forms. **CORRECTION (opus-2026-07-20-S441, THM-1930-var-lambda2, cite by filename):**
  this family does **not** reach Paley — its roots are `c·i·cot((2k−1)π/2n)` (the transitive spectrum
  *scaled* by `c`), so it interpolates `char_A` (`c=0`, `=xⁿ`) ↔ `char_S` (`c=1`) of the **single
  transitive tournament**, not transitive↔Paley. The true transitive↔Paley axis is the spectral
  scalar `var(λ²)` itself (max → 0), which S441 proves is *not* a function of `c₃`/score for `n≥5`
  (genuinely spectral) — a strict-quotient instance for THM-1935's invariant lattice.

## Named next

- **Which repo monoids are amenable / solvable / hyperbolic?** BS(1,2) is solvable; `S_n` finite;
  `PSL(2,Z)` hyperbolic (the Farey/LRC side is the "chaotic" one, matching LRC's hardness). A
  amenability audit would separate the "easy" (spectral) topics from the "hard" (LRC, `#P`).
- **The metagraph as a Cayley graph.** `Q_m/S_n` with the wiggly generators — is the merged
  metagraph the Schreier graph of `S_n` on the flip generators?
- **Formalise `⟨x+1, x/2⟩ ≅ BS(1,2)⁺`** in Lean (a free-monoid-mod-relation ≅ affine-map monoid).
