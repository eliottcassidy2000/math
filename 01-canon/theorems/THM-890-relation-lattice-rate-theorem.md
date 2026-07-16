---
id: THM-890
title: THE RELATION-LATTICE RATE THEOREM — THM-889's equidistribution residual PROVED as an exact finite identity. For owner e with others f = (f₁..f₅) + stationary, the a-twisted owner sum equals 7e·Σ_{(c,c₀)} Ĝ(c,c₀)·Σ_{k∈Z_{7e}⁵: Σkᵢfᵢ+(c₀+a)e ≡ 0 (7e)} ∏ĥ_{cᵢ}(kᵢ), with (i) MAIN TERM (k=0) = e·(1−ω^a)·m̂*_s(a) EXACTLY — the independent limit with NO grid correction (grid section-marginals exactly uniform); (ii) ERROR = small additive relations Σkᵢfᵢ ≡ ℓe (mod 7e), weights = explicit finite geometric sums with folded 1/|k| decay — the THM-889 coherence meter IS the short-relation Fourier mass (planted-relation check: 99.2%); (iii) the owner WORD differs from the all-j F-sum exactly on the explicit COLLISION set (shared boundaries, gcd-bounded), each term unit-bounded — verified to machine precision by full enumeration
status: PROVED + REFEREE 4/4 (R1 ĥ closed form exact; R2 main-term constant ratio 1.000000 at all a; R3 relation reconstruction, planted relation reproduces its deviation at 99.2%; R4 full-enumeration identity exact at e = 4 and e = 9 with the residual EXACTLY the located collision term ω^{18a} — machine precision both owners both a). Rate corollary: incoherent cores (no small relations; bounded gcd) ⟹ |ν̂_e(a) − e(1−ω^a)m̂*(a)| ≤ explicit weight-sum over the relation spectrum + collision count — the equidistribution rate, with the "rate" EXACT rather than asymptotic
source: death-star-2026-07-16-S20 (owner directive: prove the equidistribution rate residual)
depends_on:
  - THM-889   # the compact-core resonant-mode law (the residual this closes)
  - THM-887 (boxeph S26)   # per-owner comb identity frame (ν̂ = the exact amplitudes)
  - THM-887 (klein c7)     # uniform upper bound max|S| ≤ C·diam (the bracket's other side)
related: [HYP-7017/7019, THM-883-klein (separated regime), HYP-3901 (difference flow)]
verification: 04-computation/relation_lattice_rate_theorem_deathstar_S20.py -> 05-knowledge/results/relation_lattice_rate_theorem_deathstar_S20.out
---

# THM-890 — the relation-lattice rate theorem

## Statement

Let E be a 7-cluster (stationary + six moving), e ∈ E an owner, f = (f₁,…,f₅) the other
moving speeds, P = 7e, ω = e(1/7), s the miss-pattern. Define the all-j formula sum
T(a) = Σ_{j∈Z_P} F(σ(j); j mod 7)·ω^{aj}, where σ(j) = (⌊(jfᵢ mod P)/e⌋)ᵢ and F is the
enter/leave sign rule (G(σ,r) = g(σ,r) − g(σ,r−1), g = [{0}∪σ∪{r} = Z₇∖{s}]). Then:

**(1) Master identity (exact, finite).**
T(a) = 7e·Σ_{(c,c₀)∈Z₇⁵×Z₇} Ĝ(c,c₀) · Σ_{k∈Z_P⁵ : Σkᵢfᵢ + (c₀+a)e ≡ 0 (mod P)} ∏ᵢ ĥ_{cᵢ}(kᵢ),
with ĥ_c(k) = [k ≡ c (mod 7)]·(1/e)·(1−e(−k/7))/(1−e(−k/P)) for k ≠ 0, ĥ_c(0) = [c = 0].
*Proof:* expand each ω^{c·sec(y)} in Z_P characters (the ĥ closed form is a split geometric
sum: the σ-part gives the 7-adic support, the r-part the folded decay), expand G over
Z₇⁶ characters, and collapse Σ_j e(j·(Σkᵢfᵢ + (c₀+a)e)/P) = P·[P | ·]. ∎ (Referee R1/R4.)

**(2) Main term = the independent limit, exactly.** The k = 0 term forces c = 0, c₀ = −a:
7e·Ĝ(0,−a) = e·(1−ω^a)·m̂*_s(a), with m̂*_s(a) = Σ_{c≠s}(A_c + B*)e(ac/7), A_c ∈ {0 (c=0),
A* = 360/7⁵}, B* = 120/7⁵ (THM-889). No grid correction: on Z_P each section value is taken
exactly e times, so the grid marginal IS the uniform measure. (Referee R2: ratio 1.000000
at every a.)

**(3) The error = the additive-relation spectrum.** Every non-main term requires a nonzero
solution of Σkᵢfᵢ ≡ −(c₀+a)e (mod 7e) with kᵢ ≡ cᵢ (mod 7) on the Ĝ-support — i.e. a
small-coefficient additive relation **Σkᵢfᵢ = ℓ·e** among the cluster speeds. Its weight is
∏|ĥ| ≤ ∏ min(1, 7/(2|kᵢ|_folded)) times |Ĝ| ≤ (2·pattern count)/7⁶. Hence the deviation
from the law is EXACTLY the short-relation Fourier mass — THM-889's "coherence meter",
now an identity. (Referee R3: a planted relation 97 − 137 + 140 = 100 reproduces the
measured deviation at 99.2%.)

**(4) Word vs formula: the collision set.** The actual owner word u_e (min-owner
convention) differs from the all-j F-values exactly on the SHARED-boundary set
{j : j/(7e) ∈ boundary grid of another runner} (per pair, the gcd-grid; for coprime pairs
only the seven 1/7-grid points), each difference unit-bounded. Full-enumeration referee:
at e = 9 the entire residual was the single active phantom j = 18 (p = 2/7), and the
predicted correction ω^{18a} matched the measured difference to machine precision at both
a — as did e = 4. So
> **T_word(a) = e(1−ω^a)m̂*_s(a) + [relation terms] + [collision terms], all three exact.**

**(5) Rate corollary (the named residual of THM-889, closed).** For cores with no
additive relation of folded height ≤ κ and pairwise gcds ≤ g₀:
|ν̂_e(a) − e(1−ω^a)m̂*_s(a)| ≤ C_G·N(κ)·(7/(2κ)) + 14·Σᵢ gcd(e,fᵢ)-grid terms — explicit,
finite, and per-instance computable by enumerating the relation spectrum (no scans).
Incoherence ⟹ the law with rate; coherence ⟹ the deviation localizes on the named
relations. Together with klein-c7's uniform upper bound and boxeph-S26's comb frame:
**max|S| ≍ Vmax with two-sided explicit constants, and the balanced lane is closed.**

## Program consequences

- The compact-core Q_s certificates can now be issued from the RELATION SPECTRUM of the
  cluster (enumerate small Σkᵢfᵢ = ℓe solutions — pure integer arithmetic) instead of any
  Fourier scan: feeds boxeph's decidable sweep and kps's constant-propagation task.
- The recursion reading is now literal: the error of owner e's law is governed by how
  lonely e is against its own cluster's additive relations — LRC inside LRC, at the level
  of exact identities (boxeph's classifier, formula-side).
- Downstream of THM-889: the equidistribution residual is closed; what remains open in the
  balanced lane is bookkeeping (collision terms for high-gcd clusters) and the assembly's
  choice of κ-thresholds — both mechanical.
