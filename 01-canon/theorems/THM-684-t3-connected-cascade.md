---
id: THM-684
title: The t≥3 bootstrap — the orthogonality identity (character layer sums = t-fold multiplicative box counts M_t(U) = #{y ∈ Bᵗ : Πy ≡ Πu}/(q−1), EXACT) + the CS cascade per support (|M_t − bᵗ/q| ≤ √b·√(orbit pair-variance), verified DOMINATING 3–4.6×) + THE HONEST CORRECTION: raw M_t deviations scale ~q because they CONTAIN the lower layers (b × pair-deviations — measured exactly), so the vanishing cascade must run on the CONNECTED (cumulant) counts; the absolute per-support assembly over-counts the signed layer total 10–40× (the 11th confirmation of the standing law)
status: MIXED, honestly scoped. PROVED: (I) the orthogonality identity — for any support U (|U| = t), Σ_{χ_l on U, Πχ = χ₀} Π ĉ(χ_l)χ_l(u_l) = M_t(U)/(q−1) with M_t(U) = #{y ∈ Bᵗ : Π y ≡ Π u (mod q)} (trivial characters included; the pure layer follows by inclusion–exclusion over sub-supports); (II) the CS cascade — M_t(T) − bᵗ/q = Σ_{y₁∈B}[M_{t−1}(T·y₁⁻¹) − bᵗ⁻¹/q], hence |dev_t| ≤ √b·(Σ_{y₁}dev_{t−1}²)^{1/2}, recursing to the THM-683/ACZ pair objects. VERIFIED: the t = 3 CS bound dominates every sampled support (322 vs 108 @q=139; 745 vs 230 @q=239; 1522 vs 329 @q=383 — generic AND quarantined). MEASURED (the correction): raw M₃ deviations ≈ 0.8–1.0·q — exactly b × (typical pair-deviation): the t=2-within-t=3 mass dominates; the PURE (connected) layer-3 deviations are the S226-measured small remainders (0.1–5). The vanishing assembly therefore requires the CONNECTED-count cascade (inclusion–exclusion before CS) — stated below, its full verification the named next step. The absolute per-support assembly of raw bounds gives ~40 vs the signed truth ~1: the 11th documented confirmation that absolute assemblies over signed structure fail — expected, recorded.
source: klein-2026-07-09-S232 (HYP-5875; the S226-spec bootstrap, executed with the correction found)
depends_on:
  - THM-683   # the pair object (Ostrowski + ACZ) the cascade recurses to
related:
  - HYP-5845/5865/5870 (the character program), THM-680/681
---

# THM-684 — the t≥3 bootstrap: identity, cascade, and the connected-form correction

## I. The orthogonality identity (proved)

For a support U = {u₁,…,u_t} of speeds and characters χ_l on U with Πχ_l = χ₀
(trivial included), orthogonality of characters gives EXACTLY

> Σ_{Πχ=χ₀} Π_l ĉ(χ_l)χ_l(u_l) = (1/(q−1))·M_t(U),
> M_t(U) := #{(y₁,…,y_t) ∈ Bᵗ : y₁⋯y_t ≡ u₁⋯u_t (mod q)},

so every character layer is a t-fold multiplicative box count — the natural
generalization of THM-683's N_w (t = 2). The PURE layer (all χ_l nontrivial)
follows by inclusion–exclusion over sub-supports: the connected count
C_t(U) := Σ_{V⊆U}(−1)^{|U∖V|}·b^{|U∖V|}-weighted M_{|V|}(V)-terms. ∎

## II. The CS cascade (proved; verified dominating)

Peeling one coordinate: M_t(T) − bᵗ/q = Σ_{y₁∈B}[M_{t−1}(T y₁⁻¹) − bᵗ⁻¹/q],
so by Cauchy–Schwarz |M_t(T) − bᵗ/q| ≤ √b·(Σ_{y₁∈B} dev_{t−1}(T y₁⁻¹)²)^{1/2},
recursing to t = 2 where the orbit variance is the ACZ object (THM-683 V).
Verified at t = 3 on generic and quarantined instances: the bound dominates
every sampled support with 3–4.6× slack.

## III. The connected-form correction (measured; the honest finding)

Raw M₃ deviations measure ≈ 0.8–1.0·q — exactly b × (typical pair deviation):
the mass is the t = 2 layers LIVING INSIDE M₃ (one trivial character), not the
pure layer. The S226 exact decomposition already showed pure-t≥3 remainders of
only 0.1–5. Hence the VANISHING assembly must run the cascade on the CONNECTED
counts C_t (inclusion–exclusion first, then peel-and-CS): the connected pair
object is THM-683's centered N_w, and the connected cascade's verification is
the named next step. Corollary of the attempt: the absolute per-support
assembly of raw bounds gives ~40 against a signed truth of ~1 — the 11th
documented instance of the standing law (absolute assemblies over signed
resonance structure over-count by orders); recorded, expected, and the reason
the connected form is not optional.

## Files

`05-knowledge/results/lrc14_t3_cs_bootstrap_klein_S232.out` (the cascade
domination + the raw-vs-connected scale measurement).
