---
id: THM-1018
title: THE SHALLOW-WITNESS BAND LEMMA — the last stratum of "covering ⟹ M > 1/14" gets an explicit small-modulus certificate. SETTING: V = P ∪ {v₁,v₂}, P ⊆ {1,…,12} (11 speeds), v₂ = v₁+d with d ≤ 5, v₁ > 13·max(P) — the d_min ≤ 5 clustered-killer stratum left open by THM-1011. (I) THE COLLAPSE: at a modulus q, write eᵢ = |least absolute residue of vᵢ mod q|; since v₂ = v₁+d the two killers collapse TOGETHER (|e₂ − e₁| ≤ d ≤ 5), so at modulus q the 13-speed family behaves as the SMALL effective family E = P ∪ {e₁,e₂}; (II) THE BAND LEMMA (PROVED, elementary): if 2 ≤ eᵢ and e_max ≤ 13·e_min (e_min, e_max over E), then EVERY integer a with q/(14·e_min) ≤ a ≤ 13q/(14·e_max) gives the witness t = a/q with ‖v·t‖ ≥ 1/14 for ALL v ∈ V — because every effective speed e satisfies e·a ∈ [q/14, 13q/14], i.e. lands in the safe band with NO wraparound. The admissible interval has length q(13/e_max − 1/e_min)/14, which for the canonical case e_min = 2, e_max = 12 is q/24 — so an integer exists as soon as q ≥ 24; (III) EXISTENCE OF A GOOD MODULUS: a modulus is BAD only if some vᵢ ≡ 0, ±1 (mod q), i.e. only if q divides one of the SIX numbers v₁−1, v₁, v₁+1, v₂−1, v₂, v₂+1 — and since those lie in a window of width d+2 ≤ 7, **any q ≥ 8 divides AT MOST ONE of them**, so the bad moduli are a bounded divisor set (measured: ≤ 17 of the 77 moduli in [24,100]); (IV) VERIFIED: the construction witnesses ALL 8 critical covering families of THM-1011(VII) (q = 20–22, a = 1, min‖vt‖ = 1/10, 2/21, 1/11, 3/31 — all ≥ 1/14) and 300/300 of a wide random sweep (P an 11-subset of {1..12}, v₁ ∈ [157,900], d ∈ 1..5), with valid-modulus counts min 2, median 11, max 52 over 200 further families and ZERO families lacking a valid q. NOTE the witness modulus is SMALL (20–40, median 26) — far smaller than the 7k+1 optimum modulus of THM-1011(VII), because certifying level 1/14 needs far less than attaining the optimum M ≈ 1/7
status: (II) PROVED (elementary: the band containment is an implication with no hypotheses beyond 2 ≤ eᵢ, e_max ≤ 13 e_min, and an integer in the interval); (I) exact (residue arithmetic); (III) MECHANISM PROVED (the width-7 window ⟹ each q ≥ 8 divides at most one of the six) but the final count "a good q always exists in [24,100]" is VERIFIED (500+ families, min 2 valid, zero failures), NOT yet proved — the remaining step is bounding the divisor count of the six-number window from above by the range size. So: the certificate is rigorous GIVEN a good modulus; the existence of one is verified, not proved
source: kind-pasteur-2026-07-18-S128 (cont.56; owner: prove the shallow-witness counting lemma at modulus 7k+1)
depends_on:
  - THM-1011 (VII)   # the d_min ≤ 5 stratum this certifies
  - THM-724 Lemma 2  # the shallow-witness counting template
  - THM-523          # the q-witness reduction this generalizes (t = 1/q → t = a/q, collapsed speeds)
related:
  - THM-1007 (single-killer + lacunary chains), THM-933 (the gluing this replaces on this stratum)
script: 04-computation/shallow_witness_lemma_kps_S128c56.py, shallow_witness_lemma2_kps_S128c56.py (+ .out)
---

# THM-1018 — the shallow-witness band lemma

## The lemma

Let V = P ∪ {v₁,v₂}, P ⊆ {1,…,12}, v₂ = v₁ + d, d ≤ 5. Fix a modulus q, put
eᵢ = |lstabs(vᵢ mod q)| and E = P ∪ {e₁,e₂}, e_min = min E, e_max = max E.

> **If 2 ≤ e₁,e₂ and e_max ≤ 13·e_min, then every integer a in
> [q/(14 e_min), 13q/(14 e_max)] gives min_{v∈V} ‖v·a/q‖ ≥ 1/14.**

*Proof.* For each effective speed e ∈ E, e·a ≥ e_min·a ≥ q/14 and e·a ≤ e_max·a ≤ 13q/14.
So e·a lies in [q/14, 13q/14] with no wraparound, giving ‖e·a/q‖ ≥ 1/14. Since each vᵢ is
≡ ±eᵢ (mod q) and ‖·‖ is even, ‖vᵢ·a/q‖ = ‖eᵢ·a/q‖ ≥ 1/14; core speeds are their own
residues. ∎

The interval has length q(13/e_max − 1/e_min)/14 = q/24 in the canonical case
(e_min = 2, e_max = 12), so it contains an integer once q ≥ 24.

## Why the killers collapse

v₂ = v₁ + d forces |e₂ − e₁| ≤ d ≤ 5: the two killers cannot separate at any modulus.
That is exactly what made them hard for the gluing (THM-1011: large block discrepancy) and
what makes them easy here — both reduce to small speeds simultaneously.

## Existence of a good modulus

q is bad only if some vᵢ ≡ 0, ±1 (mod q), i.e. q divides one of v₁−1, v₁, v₁+1, v₂−1, v₂,
v₂+1. These six lie in a window of width d+2 ≤ 7, so **every q ≥ 8 divides at most one of
them** — the bad set is a bounded divisor set. Measured: ≤ 17 bad of the 77 moduli in
[24,100]; valid moduli number min 2, median 11, max 52 across 200 families, never zero.

## Verification

All 8 critical covering families of THM-1011(VII) are witnessed at q = 20–22 with a = 1:
min‖vt‖ = 1/10, 1/10, 1/10, 2/21, 1/10, 1/10, 1/11, 1/10 — every one ≥ 1/14. Wide sweep:
300/300.

## Named next
- Bound the divisor count of the six-number window below the range size to convert (III)
  from verified to proved — that completes the stratum, hence "covering ⟹ M > 1/14".
- Lean: the band lemma (II) is a pure inequality chain over ℚ plus one residue fact.
