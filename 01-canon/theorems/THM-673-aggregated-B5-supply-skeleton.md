---
id: THM-673
title: THM-671 part 6, the a-priori resolved-modulus supply — corrected to the AGGREGATED form. (A) the dispersal lemma — the m≠0 low-height envelope mass over all moduli q ∈ (V, 2V] totals W₁(H), V-INDEPENDENT (each (j,m) pins exactly one q), so per-modulus mass vanishes on average as V grows; (B) the per-modulus absolute-envelope prediction is REFUTED (corr(deficit, R_H) ≈ 0 — signs matter; the Mertens lesson at the modular level): per-q a-priori resolution via envelopes is dead, the target must be the AGGREGATE avg_q B5 > 0; (C) VERIFIED: avg_q B5/(q−1) ∈ [0.05, 0.12] > 0 on generic covering instances (V = 91..280, adversarial included), negative exactly on the (near-)dilated-interval family, which the exact-relation budget E_H separates by 7× (1.81 vs ≤ 0.25) — and which primitivity thins to perturbations (c·{1..13} has gcd c: exact dilations are imprimitive)
status: MIXED, honestly scoped. (A) PROVED (the counting is one line: q = (j·v)/m is determined by (j, m); the weighted total W₁(H) = Σ_j (#divisors of |j·v| in (V,2V]) · w(j) ≤ Σ_j ‖j‖₁ w(j) is V-free; VERIFIED: Σ_q R_H = 100–143 across V = 91..280 and instances). (B) REFUTED as a per-modulus predictor (measured corr −0.28..+0.40 ≈ 0 across 9 instances; the deficit is governed by SIGNED exact K̂ products, not absolute envelopes — second independent instance of the "cancellation is essential / don't bound absolutely" lesson: klein-S198 Mertens warning, kps-S89 box). (C) VERIFIED (the table below); the two remaining a-priori items are NAMED: (C1) the signed box of the m=0 low-relation contribution to avg_q S_d (closed-form per relation, kps-S89-box pattern); (C2) quantitative E_H-rigidity for the branch (E_H > E₀ ⟹ near-dilated ⟹ imprimitive core + perturbation census). LRC(14) NOT closed by this file; the certificate (THM-671) remains empirically universal.
source: klein-2026-07-09-S211 (HYP-5761; owner-directed "run THM-671 part 6")
depends_on:
  - THM-671   # the B5 certificate whose a-priori supply this addresses
  - THM-604   # the iid values; the truncation
  - LEM-015 / kps-S114 LRCSchurRigidity  # the E3-rigidity top of the branch
related:
  - HYP-5732, HYP-5761, HYP-5758
  - kps-S89 box bound (the signed-box pattern for (C1))
  - klein-S198 Mertens warning (absolute bounds vs cancellation)
---

# THM-673 — the aggregated-B5 supply skeleton (part 6, corrected)

## (A) The dispersal lemma (PROVED)

For support-≤5 vectors j (‖j‖∞ ≤ H, weights w(j) = Π min(1/7, 1/(2|j_l|)) — the
K̂ envelope) and m ≠ 0 with Σ j_l v_l = m·q: the pair (j, m) determines q. Hence

> Σ_{q ∈ (V,2V]} R_H(S, q) ≤ W₁(H) := Σ_j #{m : (j·v)/m ∈ (V,2V] ∩ Z} · w(j),

and W₁(H) ≤ Σ_j ‖j‖₁ w(j) — independent of V. VERIFIED: Σ_q R_H = 100–143
across all instances and scales tested. Consequence (Markov): the fraction of
moduli with R_H > δ is ≤ W₁/(δV) → 0. With W₁ ≈ 10², per-modulus envelope
control needs V ~ 10³⁺ — large; this motivated testing (B) at realistic V.

## (B) The refutation (measured; load-bearing negative)

Prediction tested: deficit(q) := B₅-iid − B5(S,q)/(q−1) ≈ c·(E_H(S) + R_H(S,q)).
Measured per-instance corr(deficit, R_H) over q ∈ (V, 2V]: −0.28..+0.40 ≈ 0
(9 instances). **The absolute-envelope mass does not govern the per-modulus
deficit** — the low-height resonances enter with signed exact K̂ products whose
cancellation is essential (|j| ≤ 3 coordinates all carry envelope 1/7, erasing
the height information that the true K̂ values and signs carry). Per-modulus
a-priori resolution via absolute envelopes is DEAD; certification per modulus
stays computational (THM-671's B5 itself). This is the modular-frame instance
of the project's standing lesson (Mertens warning; the D3 box).

## (C) The corrected target and its verification

**Target:** Σ_{q∈(V,2V]} B5(S,q) > 0 (⟹ some q has B5 > 0 ⟹ M(S) ≥ 1/14).

| instance | V | avg B5/q | %q with B5>0 | max B5/q | min B5/q |
|---|---|---|---|---|---|
| adv-worst-120 | 120 | +0.1001 | 100.0% | 0.189 | +0.012 |
| adv-120b | 120 | +0.0639 | 79.2% | 0.183 | −0.081 |
| adv-200 | 200 | +0.0506 | 61.5% | 0.200 | −0.256 |
| @91 7-struct | 91 | +0.1206 | 98.9% | 0.328 | 0.000 |
| adv-280 | 280 | +0.0964 | 99.3% | 0.183 | −0.010 |
| dilAP 2·{1..13} | 26 | **−2.194** | 3.8% | 0.222 | −15.57 |
| near-dilAP | 28 | **−2.836** | 10.7% | 0.316 | −14.94 |

Positive with margin on every generic covering instance; negative exactly on
the (near-)dilated-interval family. The exact-relation budget E_H(S) (support
≤ 5, height ≤ 2) separates the branch by 7×: E_H = 1.81 (dilated) vs
0.07–0.25 (everything else). Even branch instances keep max B5/q > 0
(certified moduli exist — THM-671's empirical universality is untouched).

**Structural thinning of the branch:** an exact dilated interval c·{1..13} has
gcd = c > 1 — IMPRIMITIVE, excluded by normalization. The Lemma-C branch is
only the PERTURBED family (near-dilations breaking the gcd), where E_H is
lowered by every broken relation. Quantitative version = (C2) below.

## What remains (named, elementary-shaped)

- **(C1) the signed low-relation box:** the m = 0 exact relations' contribution
  to avg_q S_d in closed form per relation (exact K̂ products, signs kept),
  boxed as kps-S89 boxed D3 — replaces the refuted envelope bound.
- **(C2) E_H-rigidity:** E_H(S) > E₀ ⟹ S within explicit distance of c·{1..13}
  (stability version of LEM-015/LRCSchurRigidity) ⟹ imprimitive core +
  perturbation census. E₀ ≈ 0.3–0.5 empirically.
- **(C3)** the finite band below the (C1)+(C2) threshold: B5-witness-search
  enumeration (kps-S115 pattern; per-set cost trivial, enumeration is the cost).

## Files

`04-computation/lrc14_part6_supply_klein_S211.py` (+out; the binary-resolution
misfire — kept as the record of why the weighted form is needed),
`lrc14_part6_weighted_mass_klein_S211.py` (+out; (A) verification, (B)
refutation), `05-knowledge/results/lrc14_avg_B5_klein_S211.out` ((C) table).
