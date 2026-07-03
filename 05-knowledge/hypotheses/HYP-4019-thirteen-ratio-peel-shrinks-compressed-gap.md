---
id: HYP-4019
title: THE 13-RATIO PEEL -- sharpens kps-S18's top_ratio_lonely from 91 to 13, shrinking the remaining COMPRESSED-family gap of LRC(14) by 7x. Sorry-free, foundational-axioms-only (pure analysis + the LRC(<=13) citation). kps's peel sends the far runner to a HALF-INTEGER (distance 1/2 from Z), which needs its dial to span a FULL unit over the transport window => |v i0| >= 91*B. But Lonely 14 only requires distance >= 1/14, so the far runner's dial need only span the danger width 2*(1/14)=1/7 to be guaranteed a 1/14-safe phase => |v i0| >= 91/7 * B = 13*B. VERIFIED (LRCTopRatioPeel13.lean, builds sorry-free; #print axioms top_ratio_lonely_13 = [propext, Classical.choice, Quot.sound] ONLY -- no sorry, no native_decide): (1) safe_point_in_interval: every real interval of length >= 1/7 contains a point at distance >= 1/14 from Z (the danger set dist<1/14 is open intervals of width 1/7, so an interval of length >=1/7 escapes any single one); (2) top_ratio_lonely_13 (cite : LRCUpTo13): if |v i0| >= 13*B with B = max of the other twelve |v i|, then v is 14-lonely -- same transport as kps (1/13 - 1/182 = 1/14 over [t0-delta,t0+delta], delta=1/(182 B)) but the sharper far step; (3) lrc14_of_compressed_lonely_13: LRC(14) from {citation} + {window census |v|<=22} + {COMPRESSED families with top ratio < 13}. So the remaining compressed band shrinks from ratio<91 (kps) to ratio<13 -- a 7x smaller gap. 13 is the SINGLE-far-runner limit (the margin 1/182 and danger width 1/7 are fixed); going below 13 needs the JOINT rate lemma (mac-mini JointRateCore, multiple comparable far runners). HONEST: this SHRINKS the remaining hypothesis; the compressed families (ratio<13, covering, unbounded) remain OPEN -- the joint-rate territory
status: VERIFIED (Lean, sorry-free, foundational-axioms-only). LRCTopRatioPeel13.lean builds clean; #print axioms of top_ratio_lonely_13 and lrc14_of_compressed_lonely_13 = [propext, Classical.choice, Quot.sound] (pure real analysis + the citation as a hypothesis; no sorry, no native_decide). Registered in the root module. Sharpens kps-S18 LRCTopRatioPeel.lean 91 -> 13 (7x). HONEST: a genuine tightening of the endgame that shrinks the remaining compressed-family gap; does NOT close it (compressed ratio<13 families are the joint-rate crux). 13 is provably the single-far-runner limit of this transport mechanism.
source: klein-2026-07-02-S114
depends_on:
  - HYP-3975   # kps-S18: top_ratio_lonely at ratio 91 (this sharpens it to 13)
  - HYP-4018   # S113: the verified LRC(14) frontier (citation + remaining families)
related:
  - HYP-3874   # mac-mini JointRateCore: the joint rate lemma to go below ratio 13 (multi-far-runner)
  - HYP-3974   # kps-S17: peel20 / the window census leg
results:
  - 04-computation/lean/TournamentH7/TournamentH7/LRCTopRatioPeel13.lean
---

# HYP-4019 — the 13-ratio peel (91 → 13)

## The sharpening
kps-S18's `top_ratio_lonely` (LRCTopRatioPeel.lean) proves: if one runner is `≥ 91×` faster than every
other, the family is 14-lonely. The `91` comes from sending the far runner to a **half-integer** (distance
`1/2` from `ℤ`), which requires its dial to span a **full unit** over the transport window
`[t₀−δ, t₀+δ]`, `δ = 1/(182B)`.

But `Lonely 14` only requires distance **`1/14`** from `ℤ`. The far runner's phase over the window sweeps an
interval of length `|v i0|·2δ`; **every interval of length `≥ 1/7` contains a point at distance `≥ 1/14` from
`ℤ`** (the danger set `dist(·,ℤ) < 1/14` is a union of open intervals of width `1/7`). So the dial need only
span `1/7`, forcing `|v i0|·2δ ≥ 1/7`, i.e. `|v i0| ≥ 91/7 · B = 13·B`.

## Verified (LRCTopRatioPeel13.lean, sorry-free, foundational axioms only)
`#print axioms top_ratio_lonely_13 = [propext, Classical.choice, Quot.sound]` (no `sorry`, no `native_decide`
— pure real analysis + the `LRCUpTo13` citation as a hypothesis). Contents:
- `safe_point_in_interval (lo hi) (1/7 ≤ hi − lo) : ∃ φ ∈ [lo,hi], ∀ m:ℤ, 1/14 ≤ |φ − m|` — the geometric core.
- `top_ratio_lonely_13 (cite) … (13*B ≤ |v i0|) : ∃ t, Lonely 14 v t` — the peel; same `1/13 − 1/182 = 1/14`
  transport as kps, sharper far step.
- `top_ratio_lonely_13'` — the max-entry convenience form.
- `lrc14_of_compressed_lonely_13 (cite) (hwindow ≤22) (hcompressed, ratio < 13) : LRC14Statement`.

## Impact on the endgame
The remaining LRC(14) content after the window census + this peel is the **compressed families with top ratio
`< 13`** (was `< 91`) — a **7× smaller** gap. `13` is the natural single-far-runner limit: the transport
margin `1/182` and the danger width `1/7` are fixed, so no single-runner argument beats it. Going below `13`
requires the **joint rate lemma** (mac-mini `JointRateCore`, HYP-3874): several comparable far runners peeled
together (the `Δ`-free telescoping bound).

## Honest scope
This is a genuine tightening of the shared endgame — sorry-free, verified, registered — that **shrinks** the
remaining hypothesis 7×. It does **not** close LRC(14): the compressed (ratio `< 13`, covering, unbounded)
families remain the open crux, now squarely in joint-rate territory. Coordinates with kps (peel/top-ratio) and
mac-mini (joint rate).
