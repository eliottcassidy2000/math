# HYP-2460: The Moser lattice IS the product at a resonant angle — unifying THM-433 and THM-434

**Status:** RESERVED STUB (monad-explorer, 2026-06-13). Verification in progress.

**Claim.** THM-433 (average degree is additive under the generic-angle Minkowski
product; the 3N crossover N* is "non-product") and THM-434 (the Moser-ladder
lattice `L_t = ℤ[ζ₆] ⊕ ω_t·ℤ[ζ₆]` has `12 + r_E(t)` unit vectors) are two faces
of ONE operation: the Minkowski product `{g + e^{iθ} h}` of two triangular-lattice
unit-distance graphs.

- At a **generic** angle θ the product realizes the Cartesian product `G □ H`
  exactly: `U = e(G)|H| + |G|e(H)` and `avgdeg` is additive (THM-433).
- At a **resonant** angle `θ = arccos((2t−1)/2t)` (the t-th Moser angle, `e^{iθ}=ω_t`)
  the rotated copy becomes commensurate, the union closes to the rank-4 lattice
  `L_t`, and **extra diagonal edges appear** — exactly the transverse unit vectors
  `α(1−ω_t)`, `N(α)=t`, of THM-434. These "resonance-bonus" edges are the
  **non-product surplus** that lets a construction beat 3N earlier (Moser at n=28)
  than any generic product can (n=32, THM-433).

**Resonance-bonus formula (to verify exactly).** For `G,H ⊂ ℤ[ζ₆]`,
```
   U(product at angle ω_t) = e(G)|H| + |G|e(H) + Δ_t(G,H),
   Δ_t(G,H) = #{ (g,g',h,h') : g−g' = α, h'−h = α, N(α)=t } = Σ_{N(α)=t} m_α(G) m_α(H),
```
where `m_α(G) = #{(g,g')∈G² : g−g' = α}` is the count of internal `√t`-displacement
`α` in `G`. So the surplus is a **correlation of the two factors' norm-t displacement
spectra** — nonzero exactly when both factors contain matching `√t`-separations.

**Why this matters.** It makes THM-433's "irreducible/non-product crossover"
*constructive*: the crosser is a product at a resonant angle, and THM-434 counts the
bonus. It predicts WHERE the extra edges of the u(28)=85 crosser come from.

**Verification plan / files:** `04-computation/resonant_product_bonus_monad.py`
(exact ℚ(√3,√11) recount, generic vs Moser angle). Upgrade to THM-493 if confirmed.
