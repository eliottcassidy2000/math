---
id: THM-684
title: The t≥3 bootstrap — the orthogonality identity [CORRECTED S233: the box object is the COMMON-MULTIPLIER (partial-live) count A_t(U) = #{c : c·u ∈ B ∀u∈U}, with A₂ = THM-683's ratio object and A₁₃ = LM(q) itself — NOT the product count M_t] + the CS cascade (a true statement about the product counts) + the connected (cumulant) cascade EXECUTED (S233): exact Möbius forms, the RELATION-TRIPLE LAW (dev₃/q → exact Kronecker-line constants; Schur = −17/1372), the exact layer closure of LM/Q, and the collapse of the absolute over-count from 10–40× to 1.5×
status: PROVED with a CORRECTION (see MISTAKE-136). CORRECTED (I): the layer-sum identity is Σ_{χ_l on U, Πχ = χ₀} Π ĉ(χ_l)χ_l(u_l) = A_t(U)/(q−1) with **A_t(U) = #{c ∈ (ℤ/q)ˣ : c·u mod q ∈ B for every u ∈ U}** — locked by direct character sums at q=61 (t=2,3; 1e-9) and the Möbius pure forms; the originally stated product count M_t = #{y ∈ Bᵗ : Πy ≡ Πu} is a DIFFERENT object (M₂ ≠ A₂ at 60/78 pairs, q=139; the near-agreement of the two at small q masked the error at S232). (II) stands as proved — but about the product counts M_t; the layer-side cascade is the character-side Hölder chain (S226). (III)'s qualitative finding (raw counts contain the lower layers; the vanishing assembly needs connected counts) SURVIVES and is now executed: see the S233 addendum — exact connected forms, the relation-triple law with closed-form torus constants (verified to 4–5 significant figures), the exact layer closure, and the measured q₀ split (noise decays; the relation part is a finite exact list).
source: klein-2026-07-09-S232 (HYP-5875; the S226-spec bootstrap, executed with the correction found); corrected and completed klein-2026-07-09-S233 (HYP-5885)
depends_on:
  - THM-683   # the pair object (Ostrowski + ACZ) the cascade recurses to
related:
  - HYP-5845/5865/5870 (the character program), THM-680/681
---

# THM-684 — the t≥3 bootstrap: identity, cascade, and the connected-form correction

## I. The orthogonality identity (proved; CORRECTED S233 — MISTAKE-136)

For a support U = {u₁,…,u_t} of speeds and characters χ_l on U with Πχ_l = χ₀
(trivial included), with ĉ(χ) = (1/(q−1))Σ_{y∈B}χ̄(y), orthogonality gives EXACTLY

> Σ_{Πχ=χ₀} Π_l ĉ(χ_l)χ_l(u_l) = (1/(q−1))·**A_t(U)**,
> **A_t(U) := #{c ∈ (ℤ/q)ˣ : c·u_l mod q ∈ B for every l}**

— the COMMON-MULTIPLIER (partial-live) count. Proof: expand each ĉ, sum the
t−1 free characters; each forces u_l·y_l⁻¹ equal to a common value d, i.e.
y_l = c·u_l with c = d⁻¹. Locked numerically at q=61 (t = 2, 3, two supports,
1e-9). Consequences: **A₂({u₁,u₂}) = N_{u₂/u₁}** (THM-683's ratio object,
verbatim), **A₁ = b exactly**, **A₀ = q−1**, and **A₁₃(full support) = LM(q)**
— the box counts ARE the partial live counts, interpolating the character
program into the certificate program at t = 13.

**The original statement's box object M_t(U) = #{y ∈ Bᵗ : Πy ≡ Πu} is a
different object** (its expansion carries a SINGLE character on all
coordinates, the diagonal tuple family): M₂ = A₂ fails at 60/78 pairs
(q = 139). The two agree in main term b²/q, which masked the error at S232.

The PURE layer (all χ_l nontrivial) follows by inclusion–exclusion over
sub-supports (Pure(∅) = 1):

> Pure(U) = Σ_{V⊆U} (−b/(q−1))^{|U∖V|} · A_{|V|}(V)/(q−1),

integer forms (Q = q−1): Pure₂ = Q·A₂ − b², Pure₃ = Q²·A₃ − b·Q·Σ_{3 pairs}A₂
+ 2b³ (verified = the direct all-nontrivial character sums, 1e-9). ∎

## II. The CS cascade (proved; verified dominating — a statement about M_t)

Peeling one coordinate: M_t(T) − bᵗ/q = Σ_{y₁∈B}[M_{t−1}(T y₁⁻¹) − bᵗ⁻¹/q],
so by Cauchy–Schwarz |M_t(T) − bᵗ/q| ≤ √b·(Σ_{y₁∈B} dev_{t−1}(T y₁⁻¹)²)^{1/2},
recursing to t = 2 where the orbit variance is the ACZ object (THM-683 V).
Verified at t = 3 on generic and quarantined instances: the bound dominates
every sampled support with 3–4.6× slack. [S233 scope note: this peel is
correct algebra for the PRODUCT counts M_t; the layer objects are the A_t of
the corrected (I), whose cascade is the character-side Hölder chain (S226) —
the connected execution below supersedes this as the program's route.]

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

## Addendum (S233) — the connected cascade EXECUTED: the relation-triple law,
## exact torus constants, and the exact layer closure

All 78 pairs and 286 triples exact (bitmask popcounts), q = 139…5003, on the
generic adversary GEN and the quarantined near-dilation DIL. Notation:
dev_t := Pure_t/Qᵗ⁻¹ (so dev₂ = A₂ − b²/Q, the centered THM-683 object);
layer_t = (b/Q)^{13−t}·(1/Q)·Σ_{|U|=t} dev_t(U); LM/Q = (b/Q)¹³ + Σ_t layer_t
EXACTLY (closure verified against LM by 13-way popcount).

**(A) THE RELATION-TRIPLE LAW.** The top-|dev₃| triples are EXACTLY the
instance's relation triples: GEN's 9 additive triples (Schur v_a+v_b=v_c and
AP 2v_b=v_a+v_c) occupy ranks 1–9 of 286, all NEGATIVE; DIL: 18 of the top 20
are relation triples (it has 67). And dev₃/q converges to EXACT Kronecker-line
constants — A₃(U)/q is the fraction of the line α·(v_i,v_j,v_k) ⊂ T³ inside
B̃³ (B̃ = [1/14, 13/14]); a low-height relation confines the line to a rational
subtorus, and the connected combination of the section integrals is the limit:

- Schur, generic direction: c₃ = T₃ − (6/7)³ = 121/196 − 216/343 =
  **−17/1372 = −0.012391** (T₃ = ∫∫_{B̃²}1_{B̃}(s+t), the band's additive
  triple self-correlation). Measured: −0.012679 (3% gap = pair corrections).
- AP, generic direction: c₃ = −0.007288 (numeric integral).
- Stacked relations (DIL): line (1,2,3) [80+160=240 AND ratios 1:2:3]:
  predicted −0.029640, measured −0.029654; line (1,3,4): predicted −0.024538,
  measured −0.024452 — **4–5 significant figures**. Relation stacking deepens
  the well (the covering/deep-well phenomenology, quantified).

The additive coordinate (Schur/E3) thus SURFACES inside the multiplicative
cumulant layer with computable rational constants: the two-coordinate
quarantine (additive E3/R(v) + multiplicative small-ratio) is ONE object at
the connected level, and THM-681's W₀ exact-relation ledger is the
non-vanishing part of this cascade seen from the Fourier side.

**(B) The exact layer closure.** GEN: LM/Q = 0.097–0.103 at EVERY q (139 →
5003); the deficit ≈ −0.038 vs (6/7)¹³ MIGRATES from layer 2 (−0.049 at
q=139) to layer 3 (−0.065 at q=5003), with R_{≥4} back-compensating (+0.029
at q=5003 — relation QUADRUPLES enter the inclusion–exclusion sign-flipped).
DIL: LM/Q ≈ 0.026 (5× suppressed); layer₃^sig → −0.44, R_{≥4} → +0.25, no
truncation convergence at any order — the quarantined family's coherence
persists in every cumulant layer (the ladder owns it, as designed).

**(C) Scales and the q₀ split.** max|dev₂| ~ q^{0.55} (ACZ-consistent);
non-relation dev₃ stays noise (median 2.87 at q=5003). The absolute-assembly
over-count COLLAPSES from 10–40× (raw counts, S232-III) to **1.5×**
(connected counts): most of the historical over-counting was self-inflicted
layer mixing, not the standing law's irreducible signed content. Measured
q₀(t≤3, absolute ≤ budget) = 1009 for GEN; absolute ≤ budget/2 is NOT reached
by q = 5003 — because the relation part does not decay and should not be
pushed to: it is a FINITE EXACT LIST with closed-form constants. The honest
architecture: **handle the relation lattice exactly (it is THM-681's W₀
object, computable per instance), let only the noise decay.** This also
explains death-star-S11's concurrent findings (per-triple layer constant in
q; closure requires cross-triple signed cancellation): the constant part is
the deterministic relation term; the cross-triple cancellation is the
relation-lattice resummation (R_{≥4} back-compensation), now with named
constants.

## Files

`05-knowledge/results/lrc14_t3_cs_bootstrap_klein_S232.out` (the cascade
domination + the raw-vs-connected scale measurement). S233:
`04-computation/lrc14_connected_cascade_klein_S233.py` (convention lock +
sweep), `04-computation/lrc14_cascade_closure_klein_S233.py` (M₂≠A₂ forensic,
exact closure, relation-triple classification), and
`05-knowledge/results/lrc14_torus_constants_klein_S233.out` (the subtorus /
Kronecker-line constants vs measured dev₃/q), with `.out` companions.
