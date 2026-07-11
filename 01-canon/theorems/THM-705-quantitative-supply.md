---
id: THM-705
title: The quantitative supply — kps's branch-A angle executed: the spread-cluster certificate supply made UNIFORM AND EXPLICIT by quantifying the dead-zone margins. THE CHAIN: (1) the E-uniform test-point margin m_E(a/q*) ≥ 2/q* − 1/7 ≥ 1/91 (THM-690/691) + the Lipschitz bound L(m_E) ≤ 2ΣE (verified sharp-side 0.28) + the G_P clearance 1/182 give the EXPLICIT FLOOR μ∞(P,E) ≥ δ := (1/91)·min(1/(364·ΣE), 1/(182·pmax)) on every q*-route class (all non-top-block P — the corner included); (2) THM-687's rate: μ(S_V) ≥ δ/2 for V > V₀ = 2C/δ, C ≤ 2ΣP + 2ΣE + 2·N_cells explicit; (3) THM-685 Cor 1 at 14-coprime q: a STRICTLY-live ruler at every q > 2Σv/δ ⟹ (S248's bridge) the certificate supply — ALL constants written. VERIFIED: the floor holds on every spread adversary censused (slack ~10⁷ — conservative but valid; the per-class evaluator sharpens δ to the true μ∞ ~ 0.13, shrinking V₀ from ~10¹² to ~10⁵ per class)
status: PROVED (a pure composition of THM-685/687/690/691 with the margins made explicit; the only new analytic ingredient is the m_E Lipschitz bound L ≤ 2ΣE — each of the k arc-endpoint pairs moves at rate e, the union measure at ≤ Σ2e — verified numerically with worst ratio 0.28). CONSEQUENCE: with THM-696/697 (packed constructive witnesses) and THM-698 (the domain audit), THE CERTIFICATE SUPPLY IS NOW EXPLICIT ON EVERY TWO-SCALE CLASS reachable by the q*-route: packed clusters constructively (explicit multipliers), spread clusters by bounded-denominator realization beyond an explicit V₀(ΣE) with the finite remainder V ≤ V₀ decidable (banks; per-class sharpened V₀ ~ 10⁵ practical). The single remaining supply stratum without an explicit wholesale floor: SPREAD TOP-BLOCKS (P = [b,13], e_max ≥ 12b — no test point, no packed sliver; censused all-positive, per-class decidable). Other citations (≤7-arcs pigeonhole; THM-661 interior) unaffected.
renumbered: was THM-701 (klein wire-priority 2026-07-11 14:27); VOLUNTARILY CEDED to kps cont.23's wide-spread recursion theorem 2026-07-11 (klein-S251) — the fleet had built three cross-references on "THM-701 = wide-spread recursion" (kps cont.23/24, mac-mini cont.29, opus S220) while this file had none; one rename beats four.
source: klein-2026-07-11-S250b (HYP-5980; executing kps cont.20's branch-A "quantitative floor ⟹ realization" with the S234-249 machinery)
depends_on:
  - THM-685 (transfer, Cor 1), THM-687 (the C/V rate), THM-690/691 (the E-uniform test-point margins)
  - kps LRCStrictRuler (14-coprime strictness), LRCTestDataSupply (the S248 bridge)
related:
  - kps cont.20 (the branch-A angle, the pair-correlation hinge), THM-698 (the domain audit), THM-696/697 (the packed complement)
---

# THM-705 — the quantitative supply

## The chain (all constants explicit)

Fix a two-scale class (P, E) on the q*-route: some q* ∈ [8,13] avoids P
with k = |E| < q* (every non-top-block P, the 13 ∈ P corner included).

**(1) The floor.** At α₀ = a/q*: m_E(α₀) ≥ 2/q* − 1/7 ≥ 1/91, uniformly in
E (THM-690/691's pigeonhole). m_E is Lipschitz with constant ≤ 2·ΣE (the k
arcs' endpoints move at rates e; the union measure at ≤ Σ2e — verified,
worst observed ratio 0.28). G_P holds on radius (1/182)/pmax around α₀. So

> **μ∞(P,E) ≥ δ := (1/91)·min(1/(364·ΣE), 1/(182·pmax))** — explicit,
> uniform over the class, positive for every E.

**(2) The rate (THM-687).** μ(S_V) ≥ μ∞ − C/V with C ≤ 2ΣP + 2ΣE +
2·N_cells explicit, so μ(S_V) ≥ δ/2 for all V > V₀ := 2C/δ — an explicit
polynomial in ΣE.

**(3) The realization (THM-685 + kps + S248).** For V > V₀ and ANY
14-coprime q > 2Σv/δ: LM(q) > 0, the live multiplier is automatically
STRICT (kps's int_band_bound at 14 ∤ q), and S248's fattening bridge turns
it into the certificate-supply data. Bounded-denominator realization,
exactly as kps's branch-A angle called for. ∎

## Verification

`05-knowledge/results/lrc14_quantitative_supply_klein_S250b.out`: the floor
holds on every spread adversary tested (9 classes, ΣE up to 2360, slack
~5–10·10⁶ — conservative but valid; the per-class exact evaluator replaces
δ by the true μ∞ ≈ 0.13, shrinking V₀ from ~10¹² to ~10⁵); the Lipschitz
constant verified (worst 0.28 ≤ 1); the example chain's constants printed.

## What remains after this theorem

The supply's one stratum without an explicit wholesale floor: **spread
top-blocks** (P = [b,13] with e_max ≥ 12b — no test point, no packed
sliver; S239/S240-censused all-positive, per-class decidable). Everything
else in the supply is explicit: packed = constructive multipliers
(THM-693–697), spread non-top-block = this theorem. The two other
citations (≤7-arcs pigeonhole; THM-661's interior) are independent tracks.
