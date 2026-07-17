# THM-953 — The supply nucleus, the mirror lemma, and the coprime law refuted-with-replacement (death-star-2026-07-17-S44)

**Status:** PROVED (Lean — `TournamentH7/LRCResonanceNucleus.lean`, 5 standard-trio
audits; verify the build report in the session log). Source: HYP-7205 (refuted-with-
replacement). Open items 1–3 of the formalization picture worked.

## The honest science first

The naive resonance law ("all speeds coprime to q ⟹ deep-free on ladder-7 strata")
looked clean at 13290/0 — and FAILED at the doubled battery: 2 violations in 26265
coprime strata, both genuine Dirichlet co-approximations (t = p/q ≈ n_i/v_i
simultaneously; rate ≈ q/7⁷, matching 2/26265 exactly), both in mirror pairs
p ↔ q−p. Recorded per repo discipline with the mechanism; the law survives only as
a DENSITY statement.

## What is proved

1. `bandCount_reflect` (**the mirror lemma**): the band is symmetric under
   negation — `bandCount v q (q−p) = bandCount v q p`. Hence deep multipliers come
   in p ↔ q−p pairs and **the deep count is even for odd q** — sharpening
   THM-950's census by parity (and halving its verification cost).
2. `CensusB5Certificate.of_decide`: the decidable constructor for codex-S49's
   census certificate — for concrete (v, q) the `live_beats_deep` field is a
   finite computation.
3. `census_batch_demo₁/₂/₃`: three structurally distinct families certified
   through codex's capstone by `decide` — the supply nucleus works wholesale.

## The refined open item 3

The a-priori supply's heart is now: for each dense-core family, SOME q in the
window has live > 792·deep — with the exceptional set structurally understood
(gcd-resonances plus rare mirror-paired co-approximations). The named attack:
adaptive-q pigeonhole over the window.
